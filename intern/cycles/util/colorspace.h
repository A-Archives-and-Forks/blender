/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#pragma once

#include "util/param.h"
#include "util/string.h"

CCL_NAMESPACE_BEGIN

/* Automatically determine colorspace. */
extern ustring u_colorspace_auto;
/* Non-color data. */
extern ustring u_colorspace_data;
/* Scene linear colorspace used for rendering. */
extern ustring u_colorspace_scene_linear;
/* Scene linear + sRGB transfer function . */
extern ustring u_colorspace_scene_linear_srgb;
/* sRGB. */
extern ustring u_colorspace_srgb;

class ColorSpaceProcessor;
struct Transform;

class ColorSpaceManager {
 public:
  /* Convert used specified colorspace to a colorspace that we are able to
   * convert to and from. If the colorspace is u_colorspace_auto, we auto
   * detect a colospace. */
  static ustring detect_known_colorspace(ustring colorspace,
                                         const char *file_colorspace,
                                         const char *file_format,
                                         bool is_float);

  /* Test if colorspace is for non-color data. */
  static bool colorspace_is_data(ustring colorspace);

  /* Convert pixels in the specified colorspace to scene linear color for
   * rendering. Must be a colorspace returned from detect_known_colorspace. */
  template<typename T>
  static void to_scene_linear(ustring colorspace,
                              T *pixels,
                              const int64_t width,
                              const int64_t height,
                              const int64_t y_stride,
                              bool is_rgba,
                              bool compress_as_srgb);

  /* Clear memory when the application exits. Invalidates all processors. */
  static void free_memory();

  /* Create a fallback color space configuration.
   *
   * This may be useful to allow regression test to create a configuration which is considered
   * valid without knowing the actual configuration used by the final application. */
  static void init_fallback_config();

  /* Compute matrix to convert from XYZ to scene linear RGB, based on the config. */
  static Transform get_xyz_to_scene_linear_rgb();
  static Transform get_xyz_to_rec709();
  /* Compute unique string for texture cache hashing and metadata. */
  static const string &get_xyz_to_scene_linear_rgb_string();
  /* Determine if scene linear is linear Rec.709 */
  static bool get_scene_linear_is_rec709();

 private:
  static ColorSpaceProcessor *get_processor(ustring colorspace);

  static void is_builtin_colorspace(ustring colorspace,
                                    bool &is_scene_linear,
                                    bool &is_scene_linear_srgb);
};

CCL_NAMESPACE_END
