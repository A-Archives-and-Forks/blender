/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#include "scene/image_oiio.h"
#include "scene/image.h"

#include "util/image.h"
#include "util/image_maketx.h"
#include "util/log.h"
#include "util/path.h"
#include "util/string.h"
#include "util/texture.h"
#include "util/thread.h"
#include "util/types_base.h"
#include "util/unique_ptr.h"

#include <OpenImageIO/filesystem.h>

CCL_NAMESPACE_BEGIN

OIIOImageLoader::OIIOImageLoader(const string &filepath) : original_filepath(filepath) {}

OIIOImageLoader::~OIIOImageLoader() = default;

static bool texture_cache_file_outdated(const string &filepath, const string &tx_filepath)
{
  if (!path_is_file(tx_filepath)) {
    return true;
  }

  std::time_t in_time = OIIO::Filesystem::last_write_time(filepath);
  std::time_t out_time = OIIO::Filesystem::last_write_time(tx_filepath);

  /* TODO: Compare metadata? maketx:full_command_line? */

  if (in_time == out_time) {
    LOG_INFO << "Using texture cache file: " << tx_filepath;
    return false;
  }

  LOG_INFO << "Texture cache file is outdated: " << tx_filepath;
  return true;
}

bool OIIOImageLoader::resolve_texture_cache(const bool auto_generate,
                                            const string &texture_cache_path,
                                            const ImageAlphaType alpha_type)
{
  /* Nothing to do if file doesn't even exist. */
  const string &filepath = get_filepath();

  if (!path_exists(filepath)) {
    return false;
  }

  /* TODO: progress display for users. */
  /* TODO: delay auto generating in case image is not used. */
  /* TODO: check if it's is a texture cache file we can actually use? */
  /* TODO: different filenames for different wrap modes, colorspace, etc? */
  /* TODO: avoid overwriting other file types? */
  const char *ext = ".tx";
  if (string_endswith(filepath, ext)) {
    return true;
  }

  /* TODO: check if path_is_relative function properly handles things like network drives. */
  /* TODO: add hash of full path to filename when using an absolute path, to avoid conflicts?
   * Though this would not be portable? */
  const string tx_filename = path_filename(filepath) + ext;
  const string tx_filepath = path_join(path_is_relative(texture_cache_path) ?
                                           path_join(path_dirname(filepath), texture_cache_path) :
                                           texture_cache_path,
                                       tx_filename);
  if (!texture_cache_file_outdated(filepath, tx_filepath)) {
    texture_cache_filepath = tx_filepath;
    return true;
  }

  /* Check in the same directory. */
  if (!texture_cache_path.empty()) {
    const string tx_local_filepath = filepath + ext;
    if (!texture_cache_file_outdated(filepath, tx_local_filepath)) {
      texture_cache_filepath = tx_local_filepath;
      return true;
    }
  }

  /* Check in default subdirectory. */
  /* TODO: not sure if we should do this. */
  const char *default_texture_cache_dir = "texture_cache";
  if (texture_cache_path != default_texture_cache_dir) {
    const string tx_default_filepath = path_join(
        path_join(path_dirname(filepath), default_texture_cache_dir), tx_filename);
    if (!texture_cache_file_outdated(filepath, tx_default_filepath)) {
      texture_cache_filepath = tx_default_filepath;
      return true;
    }
  }

  if (!auto_generate) {
    return false;
  }

  /* Auto generate. */
  LOG_INFO << "Auto generating texture cache file: " << tx_filepath;

  /* TODO: explicitly check for write permission? And even write somewhere else? */

  if (!path_create_directories(tx_filepath)) {
    LOG_WARNING << "Failed to create directory for texture cache: " << path_dirname(tx_filepath);
    return false;
  }

  /* TODO: Create ImageCache to limit memory usage, enable forcefloat? */
  /* TODO: MakeTxEnvLatl support. */
  ImageSpec configspec;
  configspec.attribute("maketx:constant_color_detect", true);
  configspec.attribute("maketx:monochrome_detect", true);
  configspec.attribute("maketx:compute_average", true);
  configspec.attribute("maketx:fixnan", true);
  /* TODO: temporarily resize to power of two until we can load other resolutions. */
  configspec.attribute("maketx:resize", true);
  /* TODO: configspec.attribute("maketx:filtername", filtername); */
  /* TODO: configspec.attribute("maketx:prman_options", true); */
  /* TODO: configspec.attribute("maketx:unpremult", associate_alpha); */
  /* TODO: configspec.attribute("maketx:incolorspace", colorspace); */
  /* TODO: configspec.attribute("maketx:wrapmodes", "black,black"); */
  /* TODO: configspec.attribute("maketx:full_command_line", full_command_line); */
  /* TODO: configspec.attribtue("maketx:set_full_to_pixels", true); */

  /* TODO: make associate alpha match Blender logic exactly. */
  ImageMetaData metadata;
  load_metadata(metadata);
  metadata.finalize(alpha_type);
  configspec.attribute("maketx:ignore_unassoc", !metadata.associate_alpha);

  OIIO::ImageBufAlgo::MakeTextureMode mode = OIIO::ImageBufAlgo::MakeTxTexture;
  std::stringstream outstream;

  if (!make_tx(mode, filepath, tx_filepath, configspec, &outstream)) {
    /* TODO: this will contain non-errors as well. OIIO::geterror() gets just the errors but is not
     * thread safe. */
    LOG_WARNING << "Failed to generate tx file: " << outstream.str();
    return false;
  }

  /* Stamp with same time as input image file to detect updates. */
  OIIO::Filesystem::last_write_time(tx_filepath, OIIO::Filesystem::last_write_time(filepath));
  assert(path_is_file(tx_filepath));
  texture_cache_filepath = tx_filepath;
  return true;
}

bool OIIOImageLoader::load_metadata(ImageMetaData &metadata)
{
  return metadata.load_metadata(get_filepath());
}

template<TypeDesc::BASETYPE FileFormat, typename StorageType>
static bool oiio_load_pixels_full(const ImageMetaData &metadata,
                                  const unique_ptr<ImageInput> &in,
                                  StorageType *pixels)
{
  const int64_t width = metadata.width;
  const int64_t height = metadata.height;
  const int channels = metadata.channels;
  const int64_t num_pixels = width * height;

  /* Read pixels through OpenImageIO. */
  StorageType *readpixels = pixels;
  vector<StorageType> tmppixels;
  if (channels > 4) {
    tmppixels.resize(num_pixels * channels);
    readpixels = &tmppixels[0];
  }

  const int64_t scanlinesize = width * channels * sizeof(StorageType);
  if (!in->read_image(0,
                      0,
                      0,
                      channels,
                      FileFormat,
                      (uchar *)readpixels + (height - 1) * scanlinesize,
                      AutoStride,
                      -scanlinesize,
                      AutoStride))
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
                                  const int64_t height,
                                  const int miplevel,
                                  const int64_t x,
                                  const int64_t y,
                                  const int64_t w,
                                  const int64_t h,
                                  const int64_t x_stride,
                                  const int64_t y_stride,
                                  uint8_t *pixels)
{
  /* Flip vertical pixel order from OIIO to Cycles convention. */
  // TODO: this fails if image is not multiple of tile size
  const int64_t flip_y = height - h - y;
  const int64_t flip_y_stride = -y_stride;
  uint8_t *flip_pixels = pixels + (h - 1) * y_stride;

  switch (metadata.type) {
    case IMAGE_DATA_TYPE_BYTE:
    case IMAGE_DATA_TYPE_BYTE4:
      return oiio_load_pixels_tile<uint8_t>(in,
                                            metadata,
                                            miplevel,
                                            x,
                                            flip_y,
                                            w,
                                            h,
                                            TypeDesc::UINT8,
                                            x_stride,
                                            flip_y_stride,
                                            flip_pixels);
    case IMAGE_DATA_TYPE_USHORT:
    case IMAGE_DATA_TYPE_USHORT4:
      return oiio_load_pixels_tile<uint16_t>(in,
                                             metadata,
                                             miplevel,
                                             x,
                                             flip_y,
                                             w,
                                             h,
                                             TypeDesc::USHORT,
                                             x_stride,
                                             flip_y_stride,
                                             reinterpret_cast<uint16_t *>(flip_pixels));
      break;
    case IMAGE_DATA_TYPE_HALF:
    case IMAGE_DATA_TYPE_HALF4:
      return oiio_load_pixels_tile<half>(in,
                                         metadata,
                                         miplevel,
                                         x,
                                         flip_y,
                                         w,
                                         h,
                                         TypeDesc::HALF,
                                         x_stride,
                                         flip_y_stride,
                                         reinterpret_cast<half *>(flip_pixels));
    case IMAGE_DATA_TYPE_FLOAT:
    case IMAGE_DATA_TYPE_FLOAT4:
      return oiio_load_pixels_tile<float>(in,
                                          metadata,
                                          miplevel,
                                          x,
                                          flip_y,
                                          w,
                                          h,
                                          TypeDesc::FLOAT,
                                          x_stride,
                                          flip_y_stride,
                                          reinterpret_cast<float *>(flip_pixels));
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
      filehandle = unique_ptr<ImageInput>(ImageInput::create(filepath));
      if (!filehandle) {
        filehandle_failed = true;
        return false;
      }

      ImageSpec spec = ImageSpec();
      if (!filehandle->open(filepath, spec)) {
        filehandle_failed = true;
        filehandle.reset();
        return false;
      }
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
