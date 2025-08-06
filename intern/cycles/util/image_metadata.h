/* SPDX-FileCopyrightText: 2011-2025 Blender Authors
 * SPDX-License-Identifier: Apache-2.0 */

#pragma once

#include <cstdlib>

#include <OpenImageIO/imageio.h>
#include <OpenImageIO/string_view.h>

#include "util/colorspace.h"
#include "util/string.h"
#include "util/texture.h"

CCL_NAMESPACE_BEGIN

/* Image MetaData
 *
 * Information about the image that is available before the image pixels are loaded. */
class ImageMetaData {
 public:
  /* Set by ImageLoader.load_metadata(). */
  int channels = 0;
  int64_t width = 0, height = 0;
  int64_t byte_size = 0;
  ImageDataType type = IMAGE_DATA_NUM_TYPES;

  /* Optional color space, defaults to scene linear. */
  ustring colorspace = u_colorspace_scene_linear;
  string colorspace_file_hint;
  const char *colorspace_file_format = "";

  /* Optional transform for 3D images. */
  bool use_transform_3d = false;
  Transform transform_3d = transform_identity();

  /* Automatically set. */
  bool compress_as_srgb = false;
  bool associate_alpha = false;
  bool ignore_alpha = false;
  bool is_cmyk = false;

  /* Tiling */
  uint32_t tile_size = 0;
  float4 average_color = zero_float4();

  ImageMetaData();
  bool operator==(const ImageMetaData &other) const;

  bool load_metadata(OIIO::string_view filepath, OIIO::ImageSpec *r_spec = nullptr);
  void finalize(const ImageAlphaType alpha_type);
  void make_float();

  bool is_float() const;
  bool is_half() const;
  bool is_rgba() const;
  TypeDesc typedesc() const;

  bool load_pixels(OIIO::string_view filepath, void *pixels) const;
  void conform_pixels(void *pixels) const;
  void conform_pixels(void *pixels,
                      const int64_t width,
                      const int64_t height,
                      const int64_t x_stride,
                      const int64_t y_stride) const;
};

CCL_NAMESPACE_END
