/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#include "scene/image_oiio.h"
#include "scene/image.h"

#include "util/image.h"
#include "util/image_maketx.h"
#include "util/path.h"
#include "util/progress.h"
#include "util/string.h"
#include "util/texture.h"
#include "util/thread.h"
#include "util/types_base.h"
#include "util/unique_ptr.h"

#include <OpenImageIO/filesystem.h>

CCL_NAMESPACE_BEGIN

OIIOImageLoader::OIIOImageLoader(const string &filepath) : original_filepath(filepath) {}

OIIOImageLoader::~OIIOImageLoader() = default;

bool OIIOImageLoader::resolve_texture_cache(const bool auto_generate,
                                            const string &texture_cache_path,
                                            const ustring &colorspace,
                                            const ImageAlphaType alpha_type,
                                            Progress &progress)
{
  const std::string &filepath = get_filepath();
  const bool found = resolve_tx(filepath,
                                texture_cache_path,
                                colorspace,
                                alpha_type,
                                IMAGE_FORMAT_PLAIN,
                                texture_cache_filepath);

  if (found) {
    return true;
  }

  if (!auto_generate) {
    texture_cache_filepath.clear();
    return false;
  }

  progress.set_status("Generating tx cache", path_filename(texture_cache_filepath));

  if (!make_tx(filepath, texture_cache_filepath, colorspace, alpha_type, IMAGE_FORMAT_PLAIN)) {
    texture_cache_filepath.clear();
  }

  return true;
}

bool OIIOImageLoader::load_metadata(ImageMetaData &metadata)
{
  return metadata.load_metadata(get_filepath());
}

bool OIIOImageLoader::load_pixels_full(const ImageMetaData &metadata, uint8_t *pixels)
{
  if (!metadata.load_pixels(get_filepath(), pixels)) {
    return false;
  }

  /* TODO: skip when loading tx file without texture cache? */
  metadata.conform_pixels(pixels);
  return true;
}

template<typename StorageType>
static bool oiio_load_pixels_tile(const unique_ptr<ImageInput> &in,
                                  const ImageMetaData &metadata,
                                  const int miplevel,
                                  const int64_t x,
                                  const int64_t y,
                                  const int64_t w,
                                  const int64_t h,
                                  const OIIO::TypeDesc typedesc,
                                  const int64_t x_stride,
                                  const int64_t y_stride,
                                  StorageType *pixels)

{
  const int channels = metadata.channels;
  const int64_t num_pixels = w * h;

  /* Read pixels through OpenImageIO. */
  StorageType *readpixels = pixels;
  vector<StorageType> tmppixels;
  if (channels > 4) {
    tmppixels.resize(num_pixels * channels);
    readpixels = &tmppixels[0];
  }

  if (!in->read_tiles(0,
                      miplevel,
                      x,
                      x + w,
                      y,
                      y + h,
                      0,
                      1,
                      0,
                      channels,
                      typedesc,
                      readpixels,
                      x_stride,
                      y_stride))
  {
    return false;
  }

  if (channels > 4) {
    for (int64_t i = num_pixels - 1, pixel = 0; pixel < num_pixels; pixel++, i--) {
      pixels[i * 4 + 3] = tmppixels[i * channels + 3];
      pixels[i * 4 + 2] = tmppixels[i * channels + 2];
      pixels[i * 4 + 1] = tmppixels[i * channels + 1];
      pixels[i * 4 + 0] = tmppixels[i * channels + 0];
    }
    tmppixels.clear();
  }

  return true;
}

static bool oiio_load_pixels_tile(const unique_ptr<ImageInput> &in,
                                  const ImageMetaData &metadata,
                                  const int64_t /*height*/,
                                  const int miplevel,
                                  const int64_t x,
                                  const int64_t y,
                                  const int64_t w,
                                  const int64_t h,
                                  const int64_t x_stride,
                                  const int64_t y_stride,
                                  uint8_t *pixels)
{
  /* Note: Here we don't flip vertical pixel order, because this was already done in the tx file.
   * TODO: Blender reading should also unflip them, or we should not flip them in the tx file.
   * Another possibility would be to change tile indexing or the convention in Cycles. */
  switch (metadata.type) {
    case IMAGE_DATA_TYPE_BYTE:
    case IMAGE_DATA_TYPE_BYTE4:
      return oiio_load_pixels_tile<uint8_t>(
          in, metadata, miplevel, x, y, w, h, TypeDesc::UINT8, x_stride, y_stride, pixels);
    case IMAGE_DATA_TYPE_USHORT:
    case IMAGE_DATA_TYPE_USHORT4:
      return oiio_load_pixels_tile<uint16_t>(in,
                                             metadata,
                                             miplevel,
                                             x,
                                             y,
                                             w,
                                             h,
                                             TypeDesc::USHORT,
                                             x_stride,
                                             y_stride,
                                             reinterpret_cast<uint16_t *>(pixels));
      break;
    case IMAGE_DATA_TYPE_HALF:
    case IMAGE_DATA_TYPE_HALF4:
      return oiio_load_pixels_tile<half>(in,
                                         metadata,
                                         miplevel,
                                         x,
                                         y,
                                         w,
                                         h,
                                         TypeDesc::HALF,
                                         x_stride,
                                         y_stride,
                                         reinterpret_cast<half *>(pixels));
    case IMAGE_DATA_TYPE_FLOAT:
    case IMAGE_DATA_TYPE_FLOAT4:
      return oiio_load_pixels_tile<float>(in,
                                          metadata,
                                          miplevel,
                                          x,
                                          y,
                                          w,
                                          h,
                                          TypeDesc::FLOAT,
                                          x_stride,
                                          y_stride,
                                          reinterpret_cast<float *>(pixels));
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT3:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT4:
    case IMAGE_DATA_TYPE_NANOVDB_FPN:
    case IMAGE_DATA_TYPE_NANOVDB_FP16:
    case IMAGE_DATA_TYPE_NANOVDB_EMPTY:
    case IMAGE_DATA_NUM_TYPES:
      return false;
  }

  return false;
}

static bool oiio_load_pixels_tile_adjacent(const unique_ptr<ImageInput> &in,
                                           const ImageMetaData &metadata,
                                           const int64_t width,
                                           const int64_t height,
                                           const int miplevel,
                                           const int64_t x,
                                           const int64_t y,
                                           const int64_t w,
                                           const int64_t h,
                                           const int64_t x_stride,
                                           const int64_t y_stride,
                                           const int x_adjacent,
                                           const int y_adjacent,
                                           const int64_t padding,
                                           const ExtensionType extension,
                                           uint8_t *pixels)
{
  const int64_t tile_size = metadata.tile_size;

  int64_t x_new = x + x_adjacent * tile_size;
  int64_t y_new = y + y_adjacent * tile_size;

  const bool in_range = x_new >= 0 && x_new < width && y_new >= 0 && y_new < height;

  const int64_t pad_x = (x_adjacent < 0) ? 0 : (x_adjacent == 0) ? padding : padding + w;
  const int64_t pad_y = (y_adjacent < 0) ? 0 : (y_adjacent == 0) ? padding : padding + h;
  const int64_t pad_w = (x_adjacent == 0) ? w : padding;
  const int64_t pad_h = (y_adjacent == 0) ? h : padding;

  if (!in_range) {
    /* Adjacent tile does not exist, fill in padding depending on extension mode. */
    if (extension == EXTENSION_EXTEND) {
      /* Duplicate pixels from border of tile. */
      for (int64_t j = 0; j < pad_h; j++) {
        for (int64_t i = 0; i < pad_w; i++) {
          const int64_t source_x = (x_adjacent < 0) ? 0 : (x_adjacent == 0) ? i : w - 1;
          const int64_t source_y = (y_adjacent < 0) ? 0 : (y_adjacent == 0) ? j : h - 1;

          std::copy_n(pixels + (padding + source_x) * x_stride + (padding + source_y) * y_stride,
                      x_stride,
                      pixels + (pad_x + i) * x_stride + (pad_y + j) * y_stride);
        }
      }
      return true;
    }
    if (extension == EXTENSION_MIRROR) {
      /* Mirror pixels from border of tile. */
      for (int64_t j = 0; j < pad_h; j++) {
        for (int64_t i = 0; i < pad_w; i++) {
          const int64_t source_x = (x_adjacent < 0)  ? padding - 1 - i :
                                   (x_adjacent == 0) ? i :
                                                       w - 1 - i;
          const int64_t source_y = (y_adjacent < 0)  ? padding - 1 - j :
                                   (y_adjacent == 0) ? j :
                                                       h - 1 - j;

          std::copy_n(pixels + (padding + source_x) * x_stride + (padding + source_y) * y_stride,
                      x_stride,
                      pixels + (pad_x + i) * x_stride + (pad_y + j) * y_stride);
        }
      }
      return true;
    }
    if (extension == EXTENSION_CLIP) {
      /* Fill with zeros. */
      for (int64_t j = 0; j < pad_h; j++) {
        std::fill_n(pixels + pad_x * x_stride + (pad_y + j) * y_stride, x_stride * pad_w, 0);
      }
      return true;
    }
    if (extension == EXTENSION_REPEAT) {
      /* Wrap around for repeat mode. */
      // TODO: fails if not multiple of tile size
      if (x_new < 0) {
        x_new = (divide_up(width, tile_size) - 1) * tile_size;
      }
      else if (x_new >= width) {
        x_new = 0;
      }

      if (y_new < 0) {
        y_new = (divide_up(height, tile_size) - 1) * tile_size;
      }
      else if (y_new >= height) {
        y_new = 0;
      }
    }
  }

  /* Load adjacent tiles. */
  vector<uint8_t> tile_pixels(tile_size * tile_size * x_stride, 0);
  if (!oiio_load_pixels_tile(in,
                             metadata,
                             height,
                             miplevel,
                             x_new,
                             y_new,
                             tile_size,
                             tile_size,
                             x_stride,
                             tile_size * x_stride,
                             tile_pixels.data()))
  {
    return false;
  }

  /* Copy pixels from adjacent tiles. */
  // TODO: verify this works for very small and non-pow2 images
  const int64_t tile_x = (x_adjacent < 0) ? std::min(width, tile_size) - padding : 0;
  const int64_t tile_y = (y_adjacent < 0) ? std::min(height, tile_size) - padding : 0;

  for (int64_t j = 0; j < pad_h; j++) {
    std::copy_n(tile_pixels.data() + tile_x * x_stride + (tile_y + j) * (tile_size * x_stride),
                x_stride * pad_w,
                pixels + pad_x * x_stride + (pad_y + j) * y_stride);
  }

  return true;
}

bool OIIOImageLoader::load_pixels_tile(const ImageMetaData &metadata,
                                       const int miplevel,
                                       const int64_t x,
                                       const int64_t y,
                                       const int64_t w,
                                       const int64_t h,
                                       const int64_t x_stride,
                                       const int64_t y_stride,
                                       const int64_t padding,
                                       const ExtensionType extension,
                                       uint8_t *pixels)
{
  assert(metadata.tile_size != 0);

  if (filehandle_failed) {
    return false;
  }
  if (!filehandle) {
    thread_scoped_lock lock(mutex);
    if (filehandle_failed) {
      return false;
    }
    if (!filehandle) {
      const string &filepath = get_filepath();
      unique_ptr<ImageInput> in = unique_ptr<ImageInput>(ImageInput::create(filepath));
      if (!in) {
        filehandle_failed = true;
        return false;
      }

      ImageSpec spec = ImageSpec();
      ImageSpec config;
      config.attribute("oiio:UnassociatedAlpha", 1);

      if (!in->open(filepath, spec, config)) {
        filehandle_failed = true;
        return false;
      }

      filehandle = std::move(in);
    }
  }

  const int64_t width = divide_up(metadata.width, 1 << miplevel);
  const int64_t height = divide_up(metadata.height, 1 << miplevel);

  /* Load center pixels. */
  bool ok = oiio_load_pixels_tile(filehandle,
                                  metadata,
                                  height,
                                  miplevel,
                                  x,
                                  y,
                                  w,
                                  h,
                                  x_stride,
                                  y_stride,
                                  pixels + padding * x_stride + padding * y_stride);

  /* Pad tile borders from adjacent tiles. */
  if (padding > 0) {
    for (int j = -1; j <= 1; j++) {
      for (int i = -1; i <= 1; i++) {
        if (i == 0 && j == 0) {
          continue;
        }
        ok &= oiio_load_pixels_tile_adjacent(filehandle,
                                             metadata,
                                             width,
                                             height,
                                             miplevel,
                                             x,
                                             y,
                                             w,
                                             h,
                                             x_stride,
                                             y_stride,
                                             i,
                                             j,
                                             padding,
                                             extension,
                                             pixels);
      }
    }
  }

  return ok;
}

void OIIOImageLoader::drop_file_handle()
{
  filehandle.reset();
}

string OIIOImageLoader::name() const
{
  return path_filename(get_filepath());
}

const string &OIIOImageLoader::get_filepath() const
{
  return (texture_cache_filepath.empty()) ? original_filepath : texture_cache_filepath;
}

bool OIIOImageLoader::equals(const ImageLoader &other) const
{
  const OIIOImageLoader &other_loader = (const OIIOImageLoader &)other;
  return original_filepath == other_loader.original_filepath;
}

CCL_NAMESPACE_END
