/* SPDX-FileCopyrightText: 2011-2025 Blender Authors
 * SPDX-License-Identifier: Apache-2.0 */

#include <cstdio>

#include <OpenImageIO/filesystem.h>
#include <OpenImageIO/typedesc.h>

#include "util/color.h"
#include "util/colorspace.h"
#include "util/image.h"
#include "util/image_maketx.h"
#include "util/image_metadata.h"
#include "util/log.h"
#include "util/param.h"
#include "util/texture.h"

CCL_NAMESPACE_BEGIN

ImageMetaData::ImageMetaData() = default;

bool ImageMetaData::operator==(const ImageMetaData &other) const
{
  return channels == other.channels && width == other.width && height == other.height &&
         use_transform_3d == other.use_transform_3d &&
         (!use_transform_3d || transform_3d == other.transform_3d) && type == other.type &&
         colorspace == other.colorspace && compress_as_srgb == other.compress_as_srgb;
}

bool ImageMetaData::is_float() const
{
  return (type == IMAGE_DATA_TYPE_FLOAT || type == IMAGE_DATA_TYPE_FLOAT4 ||
          type == IMAGE_DATA_TYPE_HALF || type == IMAGE_DATA_TYPE_HALF4);
}

bool ImageMetaData::is_half() const
{
  return (type == IMAGE_DATA_TYPE_HALF || type == IMAGE_DATA_TYPE_HALF4);
}

bool ImageMetaData::is_rgba() const
{
  return (type == IMAGE_DATA_TYPE_BYTE4 || type == IMAGE_DATA_TYPE_USHORT4 ||
          type == IMAGE_DATA_TYPE_HALF4 || type == IMAGE_DATA_TYPE_FLOAT4);
}

TypeDesc ImageMetaData::typedesc() const
{
  switch (type) {
    case IMAGE_DATA_TYPE_BYTE4:
    case IMAGE_DATA_TYPE_BYTE:
      return OIIO::TypeUInt8;
    case IMAGE_DATA_TYPE_USHORT:
    case IMAGE_DATA_TYPE_USHORT4:
      return OIIO::TypeUInt16;
    case IMAGE_DATA_TYPE_HALF4:
    case IMAGE_DATA_TYPE_HALF:
      return OIIO::TypeHalf;
    case IMAGE_DATA_TYPE_FLOAT4:
    case IMAGE_DATA_TYPE_FLOAT:
      return OIIO::TypeFloat;
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT3:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT4:
    case IMAGE_DATA_TYPE_NANOVDB_FPN:
    case IMAGE_DATA_TYPE_NANOVDB_FP16:
    case IMAGE_DATA_TYPE_NANOVDB_EMPTY:
    case IMAGE_DATA_NUM_TYPES:
      break;
  }
  return TypeUnknown;
}

void ImageMetaData::finalize(const ImageAlphaType alpha_type)
{
  /* Convert used specified color spaces to one we know how to handle. */
  colorspace = ColorSpaceManager::detect_known_colorspace(
      colorspace, colorspace_file_hint.c_str(), colorspace_file_format, is_float());

  if (colorspace == u_colorspace_scene_linear || colorspace == u_colorspace_data) {
    /* Nothing to do. */
  }
  else if (colorspace == u_colorspace_scene_linear_srgb) {
    /* Keep sRGB colorspace stored as sRGB, to save memory and/or loading time
     * for the common case of 8bit sRGB images like PNG. */
    compress_as_srgb = true;
  }
  else {
    /* If colorspace conversion needed, use half instead of short so we can
     * represent HDR values that might result from conversion. */
    if (type == IMAGE_DATA_TYPE_BYTE || type == IMAGE_DATA_TYPE_USHORT) {
      type = IMAGE_DATA_TYPE_HALF;
    }
    else if (type == IMAGE_DATA_TYPE_BYTE4 || type == IMAGE_DATA_TYPE_USHORT4) {
      type = IMAGE_DATA_TYPE_HALF4;
    }
  }

  ignore_alpha = alpha_type == IMAGE_ALPHA_IGNORE;

  /* For typical RGBA images we let OIIO convert to associated alpha,
   * but some types we want to leave the RGB channels untouched. */
  associate_alpha = associate_alpha && !(ColorSpaceManager::colorspace_is_data(colorspace) ||
                                         alpha_type == IMAGE_ALPHA_IGNORE ||
                                         alpha_type == IMAGE_ALPHA_CHANNEL_PACKED);

  /* Convert average color to scene linear colorspace. */
  if (!is_zero(average_color) && colorspace != u_colorspace_data &&
      colorspace != u_colorspace_scene_linear)
  {
    if (colorspace == u_colorspace_scene_linear_srgb) {
      average_color = color_srgb_to_linear_v4(average_color);
    }
    else {
      ColorSpaceManager::to_scene_linear(colorspace, &average_color.x, 1, 1, 1, true, false);
    }
  }
}

void ImageMetaData::make_float()
{
  switch (type) {
    case IMAGE_DATA_TYPE_BYTE:
    case IMAGE_DATA_TYPE_USHORT:
    case IMAGE_DATA_TYPE_HALF:
      type = IMAGE_DATA_TYPE_FLOAT;
      break;
    case IMAGE_DATA_TYPE_BYTE4:
    case IMAGE_DATA_TYPE_USHORT4:
    case IMAGE_DATA_TYPE_HALF4:
      type = IMAGE_DATA_TYPE_FLOAT4;
      break;
    case IMAGE_DATA_TYPE_FLOAT:
    case IMAGE_DATA_TYPE_FLOAT4:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT3:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT4:
    case IMAGE_DATA_TYPE_NANOVDB_FPN:
    case IMAGE_DATA_TYPE_NANOVDB_FP16:
    case IMAGE_DATA_TYPE_NANOVDB_EMPTY:
    case IMAGE_DATA_NUM_TYPES:
      break;
  }
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

bool ImageMetaData::load_metadata(OIIO::string_view filepath, OIIO::ImageSpec *r_spec)
{
  /* Perform preliminary checks, with meaningful logging. */
  if (!OIIO::Filesystem::exists(filepath)) {
    LOG_WARNING << "File " << filepath << " does not exist.";
    return false;
  }
  if (OIIO::Filesystem::is_directory(filepath)) {
    LOG_WARNING << "File " << filepath << " is a directory, cannot use as image.";
    return false;
  }

  std::unique_ptr<ImageInput> in(ImageInput::create(filepath));
  if (!in) {
    LOG_WARNING << "File " << filepath << " failed to open.";
    return false;
  }

  ImageSpec spec;
  ImageSpec config;

  /* Load without automatic OIIO alpha conversion. */
  config.attribute("oiio:UnassociatedAlpha", 1);

  if (!in->open(filepath, spec, config)) {
    LOG_WARNING << "File " << filepath << " failed to open.";
    return false;
  }

  width = spec.width;
  height = spec.height;
  compress_as_srgb = false;

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
  channels = spec.nchannels;

  if (is_half) {
    type = (channels > 1) ? IMAGE_DATA_TYPE_HALF4 : IMAGE_DATA_TYPE_HALF;
  }
  else if (is_float) {
    type = (channels > 1) ? IMAGE_DATA_TYPE_FLOAT4 : IMAGE_DATA_TYPE_FLOAT;
  }
  else if (spec.format == TypeDesc::USHORT) {
    type = (channels > 1) ? IMAGE_DATA_TYPE_USHORT4 : IMAGE_DATA_TYPE_USHORT;
  }
  else {
    type = (channels > 1) ? IMAGE_DATA_TYPE_BYTE4 : IMAGE_DATA_TYPE_BYTE;
  }

  colorspace_file_format = in->format_name();
  colorspace_file_hint = spec.get_string_attribute("oiio:ColorSpace");

  associate_alpha = spec.get_int_attribute("oiio:UnassociatedAlpha", 0);

  if (!associate_alpha && spec.alpha_channel != -1) {
    /* Workaround OIIO not detecting TGA file alpha the same as Blender (since #3019).
     * We want anything not marked as premultiplied alpha to get associated. */
    if (strcmp(in->format_name(), "targa") == 0) {
      associate_alpha = spec.get_int_attribute("targa:alpha_type", -1) != 4;
    }
    /* OIIO DDS reader never sets UnassociatedAlpha attribute. */
    if (strcmp(in->format_name(), "dds") == 0) {
      associate_alpha = true;
    }
    /* Workaround OIIO bug that sets oiio:UnassociatedAlpha on the last layer
     * but not composite image that we read. */
    if (strcmp(in->format_name(), "psd") == 0) {
      associate_alpha = true;
    }
  }

  /* Workaround OIIO bug that sets oiio:UnassociatedAlpha on the last layer
   * but not composite image that we read. */
  is_cmyk = strcmp(in->format_name(), "jpeg") == 0 && channels == 4;

  /* Load constant or average color. */
  if (load_metadata_color(spec, "oiio:ConstantColor", average_color)) {
    // TODO: avoid loading tiles entirely
  }
  else {
    load_metadata_color(spec, "oiio:AverageColor", average_color);
  }

  if (spec.tile_width) {
    const std::string software = spec.get_string_attribute("Software");
    int tx_file_format_version = INT_MAX;
    sscanf(software.c_str(), "Blender maketx v%d", &tx_file_format_version);

    if (tx_file_format_version == INT_MAX) {
      LOG_DEBUG << "Image " << OIIO::Filesystem::filename(filepath)
                << " has tiles, but is missing blender:TxFileFormatVersion";
    }
    else if (tx_file_format_version < 0 || tx_file_format_version > TX_FILE_FORMAT_VERSION) {
      LOG_DEBUG << "Image " << OIIO::Filesystem::filename(filepath)
                << " has tiles, but file format version " << tx_file_format_version
                << " is not supported by this version of Cycles";
    }
    else if (!(channels == 1 || channels == 4)) {
      LOG_DEBUG << "Image " << OIIO::Filesystem::filename(filepath)
                << " has tiles, but expected 1 or 4 channels, found " << channels;
    }
    else if (!is_power_of_two(spec.tile_width)) {
      LOG_DEBUG << "Image " << OIIO::Filesystem::filename(filepath)
                << " has tiles, but tile size not power of two (" << spec.tile_width << ")";
    }
    else if (spec.tile_width != spec.tile_height) {
      LOG_DEBUG << "Image " << OIIO::Filesystem::filename(filepath)
                << " has tiles, but tile size is not square (" << spec.tile_width << "x"
                << spec.tile_height << ")";
    }
    else if (spec.tile_depth != 1) {
      LOG_DEBUG << "Image " << OIIO::Filesystem::filename(filepath)
                << " has tiles, but depth is not 1";
    }
    else if (spec.tile_width < KERNEL_IMAGE_TEX_PADDING * 4) {
      LOG_DEBUG << "Image " << OIIO::Filesystem::filename(filepath)
                << " has tiles, but tile size too small (found " << spec.tile_width << ", minimum "
                << KERNEL_IMAGE_TEX_PADDING * 4 << ")";
    }
    else if (width < spec.tile_width && height < spec.tile_width) {
      // TODO: there are artifacts loading images smaller than tile size, because
      // the repeat mode is not respected for padding. Fix and and remove this exception.
      LOG_DEBUG << "Image " << OIIO::Filesystem::filename(filepath)
                << " has tiles, but image resolution is smaller than tile size";
    }
    else {
      tile_size = spec.tile_width;
    }
  }

  LOG_DEBUG << "Image " << OIIO::Filesystem::filename(filepath) << ", " << width << "x" << height
            << ", " << (tile_size ? "tiled" : "untiled");

  if (r_spec) {
    *r_spec = spec;
  }

  return true;
}

template<typename StorageType>
static void conform_pixels_to_metadata_type(const ImageMetaData &metadata,
                                            StorageType *pixels,
                                            const int64_t width,
                                            const int64_t height,
                                            const int64_t x_stride,
                                            const int64_t y_stride)
{
  /* The kernel can handle 1 and 4 channel images. Anything that is not a single
   * channel image is converted to RGBA format. */
  const int channels = metadata.channels;
  const bool is_rgba = metadata.is_rgba();

  /* CMYK to RGBA. */
  if (metadata.is_cmyk && is_rgba) {
    const StorageType one = util_image_cast_from_float<StorageType>(1.0f);

    for (int64_t j = 0; j < height; j++) {
      StorageType *pixel = pixels + j * y_stride * 4;
      for (int64_t i = 0; i < width; i++, pixel += 4) {
        const float c = util_image_cast_to_float(pixel[0]);
        const float m = util_image_cast_to_float(pixel[1]);
        const float y = util_image_cast_to_float(pixel[2]);
        const float k = util_image_cast_to_float(pixel[3]);
        pixel[0] = util_image_cast_from_float<StorageType>((1.0f - c) * (1.0f - k));
        pixel[1] = util_image_cast_from_float<StorageType>((1.0f - m) * (1.0f - k));
        pixel[2] = util_image_cast_from_float<StorageType>((1.0f - y) * (1.0f - k));
        pixel[3] = one;
      }
    }
  }

  /* Associate alpha. */
  if (channels == 4 && metadata.associate_alpha) {
    for (int64_t j = 0; j < height; j++) {
      StorageType *pixel = pixels + j * y_stride * 4;
      for (int64_t i = 0; i < width; i++, pixel += 4) {
        const StorageType alpha = pixel[3];
        pixel[0] = util_image_multiply_native(pixel[0], alpha);
        pixel[1] = util_image_multiply_native(pixel[1], alpha);
        pixel[2] = util_image_multiply_native(pixel[2], alpha);
      }
    }
  }

  if (is_rgba) {
    const StorageType one = util_image_cast_from_float<StorageType>(1.0f);

    if (channels == 2) {
      /* Grayscale + alpha to RGBA. */
      for (int64_t j = height - 1; j >= 0; j--) {
        StorageType *out_pixels = pixels + j * y_stride * 4;
        StorageType *in_pixels = pixels + j * y_stride * x_stride;
        for (int64_t i = width - 1; i >= 0; i--) {
          out_pixels[i * 4 + 3] = in_pixels[i * x_stride + 1];
          out_pixels[i * 4 + 2] = in_pixels[i * x_stride + 0];
          out_pixels[i * 4 + 1] = in_pixels[i * x_stride + 0];
          out_pixels[i * 4 + 0] = in_pixels[i * x_stride + 0];
        }
      }
    }
    else if (channels == 3) {
      /* RGB to RGBA. */
      for (int64_t j = height - 1; j >= 0; j--) {
        StorageType *out_pixels = pixels + j * y_stride * 4;
        StorageType *in_pixels = pixels + j * y_stride * x_stride;
        for (int64_t i = width - 1; i >= 0; i--) {
          out_pixels[i * 4 + 3] = one;
          out_pixels[i * 4 + 2] = in_pixels[i * x_stride + 2];
          out_pixels[i * 4 + 1] = in_pixels[i * x_stride + 1];
          out_pixels[i * 4 + 0] = in_pixels[i * x_stride + 0];
        }
      }
    }
    else if (channels == 1) {
      /* Grayscale to RGBA. */
      for (int64_t j = height - 1; j >= 0; j--) {
        StorageType *out_pixels = pixels + j * y_stride * 4;
        StorageType *in_pixels = pixels + j * y_stride * x_stride;
        for (int64_t i = width - 1; i >= 0; i--) {
          out_pixels[i * 4 + 3] = one;
          out_pixels[i * 4 + 2] = in_pixels[i * x_stride];
          out_pixels[i * 4 + 1] = in_pixels[i * x_stride];
          out_pixels[i * 4 + 0] = in_pixels[i * x_stride];
        }
      }
    }

    /* Disable alpha if requested by the user. */
    if (metadata.ignore_alpha) {
      for (int64_t j = 0; j < height; j++) {
        StorageType *out_pixels = pixels + j * y_stride * 4;
        for (int64_t i = 0; i < width; i++) {
          out_pixels[i * 4 + 3] = one;
        }
      }
    }
  }

  if (metadata.colorspace != u_colorspace_scene_linear &&
      metadata.colorspace != u_colorspace_scene_linear_srgb &&
      metadata.colorspace != u_colorspace_data)
  {
    /* Convert to scene linear. */
    ColorSpaceManager::to_scene_linear(
        metadata.colorspace, pixels, width, height, y_stride, is_rgba, metadata.compress_as_srgb);
  }

  /* Make sure we don't have buggy values. */
  if constexpr (std::is_same_v<float, StorageType> || std::is_same_v<half, StorageType>) {
    /* For RGBA buffers we put all channels to 0 if either of them is not
     * finite. This way we avoid possible artifacts caused by fully changed
     * hue. */
    if (is_rgba) {
      for (int64_t j = 0; j < height; j++) {
        StorageType *pixel = pixels + j * y_stride * 4;
        for (int64_t i = 0; i < width; i++, pixel += 4) {
          if (!util_image_is_finite(pixel[0]) || !util_image_is_finite(pixel[1]) ||
              !util_image_is_finite(pixel[2]) || !util_image_is_finite(pixel[3]))
          {
            pixel[0] = 0;
            pixel[1] = 0;
            pixel[2] = 0;
            pixel[3] = 0;
          }
        }
      }
    }
    else {
      for (int64_t j = 0; j < height; j++) {
        StorageType *pixel = pixels + j * y_stride;
        for (int64_t i = 0; i < width; i++, pixel++) {
          if (!util_image_is_finite(pixel[0])) {
            pixel[0] = 0;
          }
        }
      }
    }
  }
}

void ImageMetaData::conform_pixels(void *pixels,
                                   const int64_t width,
                                   const int64_t height,
                                   const int64_t x_stride,
                                   const int64_t y_stride) const
{
  switch (type) {
    case IMAGE_DATA_TYPE_BYTE4:
    case IMAGE_DATA_TYPE_BYTE:
      conform_pixels_to_metadata_type<uchar>(
          *this, static_cast<uchar *>(pixels), width, height, x_stride, y_stride);
      break;
    case IMAGE_DATA_TYPE_USHORT:
    case IMAGE_DATA_TYPE_USHORT4:
      conform_pixels_to_metadata_type<uint16_t>(
          *this, static_cast<uint16_t *>(pixels), width, height, x_stride, y_stride);
      break;
    case IMAGE_DATA_TYPE_HALF4:
    case IMAGE_DATA_TYPE_HALF:
      conform_pixels_to_metadata_type<half>(
          *this, static_cast<half *>(pixels), width, height, x_stride, y_stride);
      break;
    case IMAGE_DATA_TYPE_FLOAT4:
    case IMAGE_DATA_TYPE_FLOAT:
      conform_pixels_to_metadata_type<float>(
          *this, static_cast<float *>(pixels), width, height, x_stride, y_stride);
      break;
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT3:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT4:
    case IMAGE_DATA_TYPE_NANOVDB_FPN:
    case IMAGE_DATA_TYPE_NANOVDB_FP16:
    case IMAGE_DATA_TYPE_NANOVDB_EMPTY:
    case IMAGE_DATA_NUM_TYPES:
      break;
  }
}

void ImageMetaData::conform_pixels(void *pixels) const
{
  conform_pixels(pixels, width, height, channels, width);
}

template<TypeDesc::BASETYPE FileFormat, typename StorageType>
static bool load_pixels_tx(const ImageMetaData &metadata,
                           const std::unique_ptr<ImageInput> &in,
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

bool ImageMetaData::load_pixels(OIIO::string_view filepath, void *pixels) const
{
  /* load image from file through OIIO */
  std::unique_ptr<ImageInput> in = ImageInput::create(filepath);
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

  switch (type) {
    case IMAGE_DATA_TYPE_BYTE:
    case IMAGE_DATA_TYPE_BYTE4:
      return load_pixels_tx<TypeDesc::UINT8, uchar>(*this, in, (uchar *)pixels);
    case IMAGE_DATA_TYPE_USHORT:
    case IMAGE_DATA_TYPE_USHORT4:
      return load_pixels_tx<TypeDesc::USHORT, uint16_t>(*this, in, (uint16_t *)pixels);
    case IMAGE_DATA_TYPE_HALF:
    case IMAGE_DATA_TYPE_HALF4:
      return load_pixels_tx<TypeDesc::HALF, half>(*this, in, (half *)pixels);
    case IMAGE_DATA_TYPE_FLOAT:
    case IMAGE_DATA_TYPE_FLOAT4:
      return load_pixels_tx<TypeDesc::FLOAT, float>(*this, in, (float *)pixels);
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

CCL_NAMESPACE_END
