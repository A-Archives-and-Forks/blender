/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#include "scene/image_oiio.h"
#include "scene/image.h"

#include "util/colorspace.h"
#include "util/image.h"
#include "util/log.h"
#include "util/maketx.h"
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

static bool load_metadata_color(const ImageSpec &spec, const char *name, float4 &r_color)
{
  string_view metadata_color = spec.get_string_attribute(name);
  if (metadata_color.size() == 0) {
    return false;
  }

  vector<float> color;
  while (metadata_color.size()) {
    float val;
    if (!OIIO::Strutil::parse_float(metadata_color, val)) {
      break;
    }
    color.push_back(val);
    if (!OIIO::Strutil::parse_char(metadata_color, ',')) {
      break;
    }
  }

  if (color.size() != size_t(spec.nchannels)) {
    return false;
  }

  switch (spec.nchannels) {
    case 1:
      r_color = make_float4(color[0], color[0], color[0], 1.0f);
      return true;
    case 2:
      r_color = make_float4(color[0], color[0], color[0], color[1]);
      return true;
    case 3:
      r_color = make_float4(color[0], color[1], color[2], 1.0f);
      return true;
    case 4:
      r_color = make_float4(color[0], color[1], color[2], color[3]);
      return true;
    default:
      return false;
  }
}

bool OIIOImageLoader::load_metadata(ImageMetaData &metadata)
{
  /* Perform preliminary checks, with meaningful logging. */
  const string &filepath = get_filepath();
  if (!path_exists(filepath)) {
    LOG_WARNING << "File " << filepath << " does not exist.";
    return false;
  }
  if (path_is_directory(filepath)) {
    LOG_WARNING << "File " << filepath << " is a directory, can't use as image.";
    return false;
  }

  unique_ptr<ImageInput> in(ImageInput::create(filepath));
  if (!in) {
    LOG_WARNING << "File " << filepath << " failed to open.";
    return false;
  }

  ImageSpec spec;
  ImageSpec config = ImageSpec();

  /* Load without automatic OIIO alpha conversion. */
  config.attribute("oiio:UnassociatedAlpha", 1);

  if (!in->open(filepath, spec, config)) {
    LOG_WARNING << "File " << filepath << " failed to open.";
    return false;
  }

  metadata.width = spec.width;
  metadata.height = spec.height;
  metadata.compress_as_srgb = false;

  /* Check the main format, and channel formats. */
  size_t channel_size = spec.format.basesize();

  bool is_float = false;
  bool is_half = false;

  if (spec.format.is_floating_point()) {
    is_float = true;
  }

  for (size_t channel = 0; channel < spec.channelformats.size(); channel++) {
    channel_size = max(channel_size, spec.channelformats[channel].basesize());
    if (spec.channelformats[channel].is_floating_point()) {
      is_float = true;
    }
  }

  /* check if it's half float */
  if (spec.format == TypeDesc::HALF) {
    is_half = true;
  }

  /* set type and channels */
  metadata.channels = spec.nchannels;

  if (is_half) {
    metadata.type = (metadata.channels > 1) ? IMAGE_DATA_TYPE_HALF4 : IMAGE_DATA_TYPE_HALF;
  }
  else if (is_float) {
    metadata.type = (metadata.channels > 1) ? IMAGE_DATA_TYPE_FLOAT4 : IMAGE_DATA_TYPE_FLOAT;
  }
  else if (spec.format == TypeDesc::USHORT) {
    metadata.type = (metadata.channels > 1) ? IMAGE_DATA_TYPE_USHORT4 : IMAGE_DATA_TYPE_USHORT;
  }
  else {
    metadata.type = (metadata.channels > 1) ? IMAGE_DATA_TYPE_BYTE4 : IMAGE_DATA_TYPE_BYTE;
  }

  metadata.colorspace_file_format = in->format_name();
  metadata.colorspace_file_hint = spec.get_string_attribute("oiio:ColorSpace");

  metadata.associate_alpha = spec.get_int_attribute("oiio:UnassociatedAlpha", 0);

  if (!metadata.associate_alpha && spec.alpha_channel != -1) {
    /* Workaround OIIO not detecting TGA file alpha the same as Blender (since #3019).
     * We want anything not marked as premultiplied alpha to get associated. */
    if (strcmp(in->format_name(), "targa") == 0) {
      metadata.associate_alpha = spec.get_int_attribute("targa:alpha_type", -1) != 4;
    }
    /* OIIO DDS reader never sets UnassociatedAlpha attribute. */
    if (strcmp(in->format_name(), "dds") == 0) {
      metadata.associate_alpha = true;
    }
    /* Workaround OIIO bug that sets oiio:UnassociatedAlpha on the last layer
     * but not composite image that we read. */
    if (strcmp(in->format_name(), "psd") == 0) {
      metadata.associate_alpha = true;
    }
  }

  /* Load constant or average color. */
  if (load_metadata_color(spec, "oiio:ConstantColor", metadata.average_color)) {
    // TODO: avoid loading tiles entirely
  }
  else {
    load_metadata_color(spec, "oiio:AverageColor", metadata.average_color);
  }

  if (spec.tile_width) {
    // TODO: only do for particular file formats?
    if (!is_power_of_two(spec.tile_width)) {
      LOG_DEBUG << "Image " << name() << "has tiles, but tile size not power of two ("
                << spec.tile_width << ")";
    }
    else if (spec.tile_width != spec.tile_height) {
      LOG_DEBUG << "Image " << name() << "has tiles, but tile size is not square ("
                << spec.tile_width << "x" << spec.tile_height << ")";
    }
    else if (spec.tile_depth != 1) {
      LOG_DEBUG << "Image " << name() << "has tiles, but depth is not 1";
    }
    else if (spec.tile_width < KERNEL_IMAGE_TEX_PADDING * 4) {
      LOG_DEBUG << "Image " << name() << "has tiles, but tile size too small (found "
                << spec.tile_width << ", minimum " << KERNEL_IMAGE_TEX_PADDING * 4 << ")";
    }
    else if (metadata.width < spec.tile_width && metadata.height < spec.tile_width) {
      // TODO: there are artifacts loading images smaller than tile size, because
      // the repeat mode is not respected for padding. Fix and and remove this exception.
      LOG_DEBUG << "Image " << name()
                << "has tiles, but image resolution is smaller than tile size";
    }
    else {
      metadata.tile_size = spec.tile_width;
    }
  }

  LOG_DEBUG << "Image " << name() << ", " << metadata.width << "x" << metadata.height << ", "
            << (metadata.tile_size ? "tiled" : "untiled");

  return true;
}

template<typename StorageType>
static void oiio_conform_alpha_cmyk(const ImageMetaData &metadata,
                                    const unique_ptr<ImageInput> &in,
                                    StorageType *pixels,
                                    const size_t num_pixels)
{
  /* CMYK to RGBA. */
  const bool cmyk = strcmp(in->format_name(), "jpeg") == 0 && metadata.channels == 4;
  if (cmyk) {
    const StorageType one = util_image_cast_from_float<StorageType>(1.0f);

    for (int64_t i = num_pixels - 1, pixel = 0; pixel < num_pixels; pixel++, i--) {
      const float c = util_image_cast_to_float(pixels[i * 4 + 0]);
      const float m = util_image_cast_to_float(pixels[i * 4 + 1]);
      const float y = util_image_cast_to_float(pixels[i * 4 + 2]);
      const float k = util_image_cast_to_float(pixels[i * 4 + 3]);
      pixels[i * 4 + 0] = util_image_cast_from_float<StorageType>((1.0f - c) * (1.0f - k));
      pixels[i * 4 + 1] = util_image_cast_from_float<StorageType>((1.0f - m) * (1.0f - k));
      pixels[i * 4 + 2] = util_image_cast_from_float<StorageType>((1.0f - y) * (1.0f - k));
      pixels[i * 4 + 3] = one;
    }
  }

  if (metadata.channels == 4 && metadata.associate_alpha) {
    for (int64_t i = num_pixels - 1, pixel = 0; pixel < num_pixels; pixel++, i--) {
      const StorageType alpha = pixels[i * 4 + 3];
      pixels[i * 4 + 0] = util_image_multiply_native(pixels[i * 4 + 0], alpha);
      pixels[i * 4 + 1] = util_image_multiply_native(pixels[i * 4 + 1], alpha);
      pixels[i * 4 + 2] = util_image_multiply_native(pixels[i * 4 + 2], alpha);
    }
  }
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

  oiio_conform_alpha_cmyk(metadata, in, pixels, width * height);

  return true;
}

bool OIIOImageLoader::load_pixels_full(const ImageMetaData &metadata, uint8_t *pixels)
{
  /* load image from file through OIIO */
  const string &filepath = get_filepath();
  unique_ptr<ImageInput> in = unique_ptr<ImageInput>(ImageInput::create(filepath));
  if (!in) {
    return false;
  }

  ImageSpec spec = ImageSpec();
  ImageSpec config = ImageSpec();

  /* Load without automatic OIIO alpha conversion, we do it ourselves. OIIO
   * will associate alpha in the 8bit buffer for PNGs, which leads to too
   * much precision loss when we load it as half float to do a color-space transform. */
  config.attribute("oiio:UnassociatedAlpha", 1);

  if (!in->open(filepath, spec, config)) {
    return false;
  }

  switch (metadata.type) {
    case IMAGE_DATA_TYPE_BYTE:
    case IMAGE_DATA_TYPE_BYTE4:
      return oiio_load_pixels_full<TypeDesc::UINT8, uchar>(metadata, in, (uchar *)pixels);
    case IMAGE_DATA_TYPE_USHORT:
    case IMAGE_DATA_TYPE_USHORT4:
      return oiio_load_pixels_full<TypeDesc::USHORT, uint16_t>(metadata, in, (uint16_t *)pixels);
    case IMAGE_DATA_TYPE_HALF:
    case IMAGE_DATA_TYPE_HALF4:
      return oiio_load_pixels_full<TypeDesc::HALF, half>(metadata, in, (half *)pixels);
    case IMAGE_DATA_TYPE_FLOAT:
    case IMAGE_DATA_TYPE_FLOAT4:
      return oiio_load_pixels_full<TypeDesc::FLOAT, float>(metadata, in, (float *)pixels);
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT3:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT4:
    case IMAGE_DATA_TYPE_NANOVDB_FPN:
    case IMAGE_DATA_TYPE_NANOVDB_FP16:
    case IMAGE_DATA_TYPE_NANOVDB_EMPTY:
    case IMAGE_DATA_NUM_TYPES:
      break;
  }

  return false;
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

  oiio_conform_alpha_cmyk<StorageType>(metadata, in, pixels, w * h);
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
