/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#include "util/colorspace.h"
#include "util/color.h"
#include "util/image.h"
#include "util/log.h"
#include "util/map.h"
#include "util/math.h"
#include "util/thread.h"
#include "util/transform.h"
#include "util/vector.h"

#ifdef WITH_OCIO
#  include <OpenColorIO/OpenColorIO.h>
namespace OCIO = OCIO_NAMESPACE;
#endif

CCL_NAMESPACE_BEGIN

/* Builtin colorspaces. */
ustring u_colorspace_auto;
ustring u_colorspace_raw("__builtin_raw");
ustring u_colorspace_srgb("__builtin_srgb");

/* Cached data. */
#ifdef WITH_OCIO
static thread_mutex cache_colorspaces_mutex;
static thread_mutex cache_processors_mutex;
static unordered_map<ustring, ustring> cached_colorspaces;
static unordered_map<ustring, OCIO::ConstProcessorRcPtr> cached_processors;
#endif

ColorSpaceProcessor *ColorSpaceManager::get_processor(ustring colorspace)
{
#ifdef WITH_OCIO
  /* Only use this for OpenColorIO color spaces, not the builtin ones. */
  assert(colorspace != u_colorspace_srgb && colorspace != u_colorspace_auto);

  if (colorspace == u_colorspace_raw) {
    return nullptr;
  }

  OCIO::ConstConfigRcPtr config = nullptr;
  try {
    config = OCIO::GetCurrentConfig();
  }
  catch (const OCIO::Exception &exception) {
    LOG_WARNING << "OCIO config error: " << exception.what();
    return nullptr;
  }

  if (!config) {
    return nullptr;
  }

  /* Cache processor until free_memory(), memory overhead is expected to be
   * small and the processor is likely to be reused. */
  const thread_scoped_lock cache_processors_lock(cache_processors_mutex);
  if (cached_processors.find(colorspace) == cached_processors.end()) {
    try {
      cached_processors[colorspace] = config->getProcessor(colorspace.c_str(), "scene_linear");
    }
    catch (const OCIO::Exception &exception) {
      cached_processors[colorspace] = OCIO::ConstProcessorRcPtr();
      LOG_WARNING << "Colorspace " << colorspace.c_str()
                  << " can't be converted to scene_linear: " << exception.what();
    }
  }

  const OCIO::Processor *processor = cached_processors[colorspace].get();
  return (ColorSpaceProcessor *)processor;
#else
  /* No OpenColorIO. */
  (void)colorspace;
  return nullptr;
#endif
}

bool ColorSpaceManager::colorspace_is_data(ustring colorspace)
{
  if (colorspace == u_colorspace_auto || colorspace == u_colorspace_raw ||
      colorspace == u_colorspace_srgb)
  {
    return false;
  }

#ifdef WITH_OCIO
  OCIO::ConstConfigRcPtr config = nullptr;
  try {
    config = OCIO::GetCurrentConfig();
  }
  catch (const OCIO::Exception &exception) {
    LOG_WARNING << "OCIO config error: " << exception.what();
    return false;
  }

  if (!config) {
    return false;
  }

  try {
    const OCIO::ConstColorSpaceRcPtr space = config->getColorSpace(colorspace.c_str());
    return space && space->isData();
  }
  catch (const OCIO::Exception &) {
    return false;
  }
#else
  return false;
#endif
}

ustring ColorSpaceManager::detect_known_colorspace(ustring colorspace,
                                                   const char *file_colorspace,
                                                   const char *file_format,
                                                   bool is_float)
{
  if (colorspace == u_colorspace_auto) {
    /* Auto detect sRGB or raw if none specified. */
    if (is_float) {
      const bool srgb = (strcmp(file_colorspace, "sRGB") == 0 ||
                         strcmp(file_colorspace, "GammaCorrected") == 0 ||
                         (file_colorspace[0] == '\0' &&
                          (strcmp(file_format, "png") == 0 || strcmp(file_format, "jpeg") == 0 ||
                           strcmp(file_format, "tiff") == 0 || strcmp(file_format, "dpx") == 0 ||
                           strcmp(file_format, "jpeg2000") == 0)));
      return srgb ? u_colorspace_srgb : u_colorspace_raw;
    }
    return u_colorspace_srgb;
  }

  /* Builtin colorspaces. */
  if (colorspace == u_colorspace_srgb || colorspace == u_colorspace_raw) {
    return colorspace;
  }

  /* Use OpenColorIO. */
#ifdef WITH_OCIO
  {
    const thread_scoped_lock cache_lock(cache_colorspaces_mutex);
    /* Cached lookup. */
    if (cached_colorspaces.find(colorspace) != cached_colorspaces.end()) {
      return cached_colorspaces[colorspace];
    }
  }

  /* Detect if it matches a simple builtin colorspace. */
  bool is_scene_linear;
  bool is_srgb;
  is_builtin_colorspace(colorspace, is_scene_linear, is_srgb);

  const thread_scoped_lock cache_lock(cache_colorspaces_mutex);
  if (is_scene_linear) {
    LOG_INFO << "Colorspace " << colorspace.string() << " is no-op";
    cached_colorspaces[colorspace] = u_colorspace_raw;
    return u_colorspace_raw;
  }
  if (is_srgb) {
    LOG_INFO << "Colorspace " << colorspace.string() << " is sRGB";
    cached_colorspaces[colorspace] = u_colorspace_srgb;
    return u_colorspace_srgb;
  }

  /* Verify if we can convert from the requested color space. */
  if (!get_processor(colorspace)) {
    OCIO::ConstConfigRcPtr config = nullptr;
    try {
      config = OCIO::GetCurrentConfig();
    }
    catch (const OCIO::Exception &exception) {
      LOG_WARNING << "OCIO config error: " << exception.what();
      return u_colorspace_raw;
    }

    if (!config || !config->getColorSpace(colorspace.c_str())) {
      LOG_WARNING << "Colorspace " << colorspace.c_str() << " not found, using raw instead";
    }
    else {
      LOG_WARNING << "Colorspace " << colorspace.c_str()
                  << " can't be converted to scene_linear, using raw instead";
    }
    cached_colorspaces[colorspace] = u_colorspace_raw;
    return u_colorspace_raw;
  }

  /* Convert to/from colorspace with OpenColorIO. */
  LOG_INFO << "Colorspace " << colorspace.string() << " handled through OpenColorIO";
  cached_colorspaces[colorspace] = colorspace;
  return colorspace;
#else
  LOG_WARNING << "Colorspace " << colorspace.c_str()
              << " not available, built without OpenColorIO";
  return u_colorspace_raw;
#endif
}

void ColorSpaceManager::is_builtin_colorspace(ustring colorspace,
                                              bool &is_scene_linear,
                                              bool &is_srgb)
{
#ifdef WITH_OCIO
  const OCIO::Processor *processor = (const OCIO::Processor *)get_processor(colorspace);
  if (!processor) {
    is_scene_linear = false;
    is_srgb = false;
    return;
  }

  const OCIO::ConstCPUProcessorRcPtr device_processor = processor->getDefaultCPUProcessor();
  is_scene_linear = true;
  is_srgb = true;
  for (int i = 0; i < 256; i++) {
    const float v = i / 255.0f;

    float cR[3] = {v, 0, 0};
    float cG[3] = {0, v, 0};
    float cB[3] = {0, 0, v};
    float cW[3] = {v, v, v};
    device_processor->applyRGB(cR);
    device_processor->applyRGB(cG);
    device_processor->applyRGB(cB);
    device_processor->applyRGB(cW);

    /* Make sure that there is no channel crosstalk. */
    if (fabsf(cR[1]) > 1e-5f || fabsf(cR[2]) > 1e-5f || fabsf(cG[0]) > 1e-5f ||
        fabsf(cG[2]) > 1e-5f || fabsf(cB[0]) > 1e-5f || fabsf(cB[1]) > 1e-5f)
    {
      is_scene_linear = false;
      is_srgb = false;
      break;
    }
    /* Make sure that the three primaries combine linearly. */
    if (!compare_floats(cR[0], cW[0], 1e-6f, 64) || !compare_floats(cG[1], cW[1], 1e-6f, 64) ||
        !compare_floats(cB[2], cW[2], 1e-6f, 64))
    {
      is_scene_linear = false;
      is_srgb = false;
      break;
    }
    /* Make sure that the three channels behave identically. */
    if (!compare_floats(cW[0], cW[1], 1e-6f, 64) || !compare_floats(cW[1], cW[2], 1e-6f, 64)) {
      is_scene_linear = false;
      is_srgb = false;
      break;
    }

    const float out_v = average(make_float3(cW[0], cW[1], cW[2]));
    if (!compare_floats(v, out_v, 1e-6f, 64)) {
      is_scene_linear = false;
    }
    if (!compare_floats(color_srgb_to_linear(v), out_v, 1e-4f, 64)) {
      is_srgb = false;
    }
  }
#else
  (void)colorspace;
  is_scene_linear = false;
  is_srgb = false;
#endif
}

#ifdef WITH_OCIO

template<typename T> inline float4 cast_to_float4(T *data)
{
  return make_float4(util_image_cast_to_float(data[0]),
                     util_image_cast_to_float(data[1]),
                     util_image_cast_to_float(data[2]),
                     util_image_cast_to_float(data[3]));
}

template<typename T> inline void cast_from_float4(T *data, const float4 value)
{
  data[0] = util_image_cast_from_float<T>(value.x);
  data[1] = util_image_cast_from_float<T>(value.y);
  data[2] = util_image_cast_from_float<T>(value.z);
  data[3] = util_image_cast_from_float<T>(value.w);
}

/* Slower versions for other all data types, which needs to convert to float and back. */
template<typename T, bool compress_as_srgb = false>
inline void processor_apply_pixels_rgba(const OCIO::Processor *processor,
                                        T *pixels,
                                        const int64_t width,
                                        const int64_t height,
                                        const int64_t y_stride)
{
  /* TODO: implement faster version for when we know the conversion
   * is a simple matrix transform between linear spaces. In that case
   * un-premultiply is not needed. */
  const OCIO::ConstCPUProcessorRcPtr device_processor = processor->getDefaultCPUProcessor();

  /* Process large images in chunks to keep temporary memory requirement down. */
  const int64_t chunk_rows = divide_up(std::min((int64_t)(16 * 1024 * 1024), width * height),
                                       width);
  vector<float4> float_pixels(chunk_rows * width);

  for (int64_t row = 0; row < height; row += chunk_rows) {
    const int64_t num_rows = std::min(chunk_rows, height - row);

    for (int64_t j = 0; j < num_rows; j++) {
      T *pixel = pixels + (row + j) * y_stride * 4;
      float4 *float_pixel = float_pixels.data() + j * width;
      for (int64_t i = 0; i < width; i++, pixel += 4, float_pixel++) {
        float4 value = cast_to_float4(pixel);

        if (!(value.w <= 0.0f || value.w == 1.0f)) {
          const float inv_alpha = 1.0f / value.w;
          value.x *= inv_alpha;
          value.y *= inv_alpha;
          value.z *= inv_alpha;
        }

        *float_pixel = value;
      }
    }

    const OCIO::PackedImageDesc desc((float *)float_pixels.data(), num_rows * width, 1, 4);
    device_processor->apply(desc);

    for (int64_t j = 0; j < num_rows; j++) {
      T *pixel = pixels + (row + j) * y_stride * 4;
      float4 *float_pixel = float_pixels.data() + j * width;
      for (int64_t i = 0; i < width; i++, pixel += 4, float_pixel++) {
        float4 value = *float_pixel;

        if (compress_as_srgb) {
          value = color_linear_to_srgb_v4(value);
        }

        if (!(value.w <= 0.0f || value.w == 1.0f)) {
          value.x *= value.w;
          value.y *= value.w;
          value.z *= value.w;
        }

        cast_from_float4(pixel, value);
      }
    }
  }
}

template<typename T, bool compress_as_srgb = false>
inline void processor_apply_pixels_grayscale(const OCIO::Processor *processor,
                                             T *pixels,
                                             const int64_t width,
                                             const int64_t height,
                                             const int64_t y_stride)
{
  const OCIO::ConstCPUProcessorRcPtr device_processor = processor->getDefaultCPUProcessor();

  /* Process large images in chunks to keep temporary memory requirement down. */
  const int64_t chunk_rows = divide_up(std::min((int64_t)(16 * 1024 * 1024), width * height),
                                       width);
  vector<float> float_pixels(chunk_rows * width * 3);

  for (int64_t row = 0; row < height; row += chunk_rows) {
    const int64_t num_rows = std::min(chunk_rows, height - row);

    /* Convert to 3 channels, since that's the minimum required by OpenColorIO. */
    for (int64_t j = 0; j < num_rows; j++) {
      T *pixel = pixels + (row + j) * y_stride;
      float *float_pixel = float_pixels.data() + j * width * 3;
      for (int64_t i = 0; i < width; i++, pixel++, float_pixel += 3) {
        const float f = util_image_cast_to_float<T>(*pixel);
        float_pixel[0] = f;
        float_pixel[1] = f;
        float_pixel[2] = f;
      }
    }

    const OCIO::PackedImageDesc desc((float *)float_pixels.data(), num_rows * width, 1, 3);
    device_processor->apply(desc);

    for (int64_t j = 0; j < num_rows; j++) {
      T *pixel = pixels + (row + j) * y_stride;
      float *float_pixel = float_pixels.data() + j * width * 3;
      for (int64_t i = 0; i < width; i++, pixel++, float_pixel += 3) {
        float f = average(make_float3(float_pixel[0], float_pixel[1], float_pixel[2]));
        if (compress_as_srgb) {
          f = color_linear_to_srgb(f);
        }
        *pixel = util_image_cast_from_float<T>(f);
      }
    }
  }
}

#endif

template<typename T>
void ColorSpaceManager::to_scene_linear(ustring colorspace,
                                        T *pixels,
                                        const int64_t width,
                                        const int64_t height,
                                        const int64_t y_stride,
                                        bool is_rgba,
                                        bool compress_as_srgb)
{
#ifdef WITH_OCIO
  const OCIO::Processor *processor = (const OCIO::Processor *)get_processor(colorspace);

  if (processor) {
    if (is_rgba) {
      if (compress_as_srgb) {
        /* Compress output as sRGB. */
        processor_apply_pixels_rgba<T, true>(processor, pixels, width, height, y_stride);
      }
      else {
        /* Write output as scene linear directly. */
        processor_apply_pixels_rgba<T>(processor, pixels, width, height, y_stride);
      }
    }
    else {
      if (compress_as_srgb) {
        /* Compress output as sRGB. */
        processor_apply_pixels_grayscale<T, true>(processor, pixels, width, height, y_stride);
      }
      else {
        /* Write output as scene linear directly. */
        processor_apply_pixels_grayscale<T>(processor, pixels, width, height, y_stride);
      }
    }
  }
#else
  (void)colorspace;
  (void)pixels;
  (void)width;
  (void)height;
  (void)y_stride;
  (void)is_rgba;
  (void)compress_as_srgb;
#endif
}

void ColorSpaceManager::to_scene_linear(ColorSpaceProcessor *processor_,
                                        float *pixel,
                                        const int channels)
{
#ifdef WITH_OCIO
  const OCIO::Processor *processor = (const OCIO::Processor *)processor_;

  if (processor) {
    const OCIO::ConstCPUProcessorRcPtr device_processor = processor->getDefaultCPUProcessor();
    if (channels == 1) {
      float3 rgb = make_float3(pixel[0], pixel[0], pixel[0]);
      device_processor->applyRGB(&rgb.x);
      pixel[0] = average(rgb);
    }
    if (channels == 3) {
      device_processor->applyRGB(pixel);
    }
    else if (channels == 4) {
      if (pixel[3] == 1.0f || pixel[3] == 0.0f) {
        /* Fast path for RGBA. */
        device_processor->applyRGB(pixel);
      }
      else {
        /* Un-associate and associate alpha since color management should not
         * be affected by transparency. */
        const float alpha = pixel[3];
        const float inv_alpha = 1.0f / alpha;

        pixel[0] *= inv_alpha;
        pixel[1] *= inv_alpha;
        pixel[2] *= inv_alpha;

        device_processor->applyRGB(pixel);

        pixel[0] *= alpha;
        pixel[1] *= alpha;
        pixel[2] *= alpha;
      }
    }
  }
#else
  (void)processor_;
  (void)pixel;
  (void)channels;
#endif
}

void ColorSpaceManager::free_memory()
{
#ifdef WITH_OCIO
  map_free_memory(cached_colorspaces);
  map_free_memory(cached_processors);
#endif
}

void ColorSpaceManager::init_fallback_config()
{
#ifdef WITH_OCIO
  OCIO::SetCurrentConfig(OCIO::Config::CreateRaw());
#endif
}

#ifdef WITH_OCIO
static bool to_scene_linear_transform(OCIO::ConstConfigRcPtr &config,
                                      const char *colorspace,
                                      Transform &to_scene_linear)
{
  OCIO::ConstProcessorRcPtr processor;
  try {
    processor = config->getProcessor("scene_linear", colorspace);
  }
  catch (OCIO::Exception &) {
    return false;
  }

  if (!processor) {
    return false;
  }

  const OCIO::ConstCPUProcessorRcPtr device_processor = processor->getDefaultCPUProcessor();
  if (!device_processor) {
    return false;
  }

  to_scene_linear = transform_identity();
  device_processor->applyRGB(&to_scene_linear.x.x);
  device_processor->applyRGB(&to_scene_linear.y.x);
  device_processor->applyRGB(&to_scene_linear.z.x);
  to_scene_linear = transform_transposed_inverse(to_scene_linear);
  return true;
}
#endif

bool ColorSpaceManager::get_xyz_to_scene_linear_rgb(Transform &xyz_to_rgb)
{
#ifdef WITH_OCIO
  /* Get from OpenColorO config if it has the required roles. */
  OCIO::ConstConfigRcPtr config = nullptr;
  try {
    config = OCIO::GetCurrentConfig();
  }
  catch (OCIO::Exception &exception) {
    LOG_WARNING << "OCIO config error: " << exception.what();
    return false;
  }

  if (!(config && config->hasRole("scene_linear"))) {
    return false;
  }

  if (config->hasRole("aces_interchange")) {
    /* Standard OpenColorIO role, defined as ACES AP0 (ACES2065-1). */
    Transform aces_to_rgb;
    if (!to_scene_linear_transform(config, "aces_interchange", aces_to_rgb)) {
      return false;
    }

    /* This is the OpenColorIO builtin transform:
     * UTILITY - ACES-AP0_to_CIE-XYZ-D65_BFD. */
    const Transform ACES_AP0_to_xyz_D65 = make_transform(0.938280f,
                                                         -0.004451f,
                                                         0.016628f,
                                                         0.000000f,
                                                         0.337369f,
                                                         0.729522f,
                                                         -0.066890f,
                                                         0.000000f,
                                                         0.001174f,
                                                         -0.003711f,
                                                         1.091595f,
                                                         0.000000f);
    const Transform xyz_to_aces = transform_inverse(ACES_AP0_to_xyz_D65);
    xyz_to_rgb = aces_to_rgb * xyz_to_aces;
    return true;
  }

  if (config->hasRole("XYZ")) {
    /* Custom role used before the standard existed. */
    if (to_scene_linear_transform(config, "XYZ", xyz_to_rgb)) {
      return true;
    }
  }

  /* No reference role found to determine XYZ. */
#endif
  return false;
}

/* Template instantiations so we don't have to inline functions. */
template void ColorSpaceManager::to_scene_linear(
    ustring, uchar *, int64_t, int64_t, int64_t, bool, bool);
template void ColorSpaceManager::to_scene_linear(
    ustring, ushort *, int64_t, int64_t, int64_t, bool, bool);
template void ColorSpaceManager::to_scene_linear(
    ustring, half *, int64_t, int64_t, int64_t, bool, bool);
template void ColorSpaceManager::to_scene_linear(
    ustring, float *, int64_t, int64_t, int64_t, bool, bool);

CCL_NAMESPACE_END
