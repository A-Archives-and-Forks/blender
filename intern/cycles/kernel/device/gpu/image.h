/* SPDX-FileCopyrightText: 2017-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#pragma once

#include "kernel/globals.h"
#include "kernel/util/image.h"

CCL_NAMESPACE_BEGIN

ccl_device_inline float frac(const float x, ccl_private int *ix)
{
  int i = float_to_int(x) - ((x < 0.0f) ? 1 : 0);
  *ix = i;
  return x - (float)i;
}

/* w0, w1, w2, and w3 are the four cubic B-spline basis functions. */
ccl_device float cubic_w0(const float a)
{
  return (1.0f / 6.0f) * (a * (a * (-a + 3.0f) - 3.0f) + 1.0f);
}
ccl_device float cubic_w1(const float a)
{
  return (1.0f / 6.0f) * (a * a * (3.0f * a - 6.0f) + 4.0f);
}
ccl_device float cubic_w2(const float a)
{
  return (1.0f / 6.0f) * (a * (a * (-3.0f * a + 3.0f) + 3.0f) + 1.0f);
}
ccl_device float cubic_w3(const float a)
{
  return (1.0f / 6.0f) * (a * a * a);
}

/* g0 and g1 are the two amplitude functions. */
ccl_device float cubic_g0(const float a)
{
  return cubic_w0(a) + cubic_w1(a);
}
ccl_device float cubic_g1(const float a)
{
  return cubic_w2(a) + cubic_w3(a);
}

/* h0 and h1 are the two offset functions */
ccl_device float cubic_h0(const float a)
{
  return (cubic_w1(a) / cubic_g0(a)) - 1.0f;
}
ccl_device float cubic_h1(const float a)
{
  return (cubic_w3(a) / cubic_g1(a)) + 1.0f;
}

/* Fast bicubic texture lookup using 4 bilinear lookups, adapted from CUDA samples. */
template<typename T>
ccl_device_noinline T kernel_image_interp_bicubic(const ccl_global KernelImageInfo &info,
                                                  const float2 uv)
{
  ccl_gpu_image_object_2D tex = (ccl_gpu_image_object_2D)info.data;

  const float x = (uv.x * info.width) - 0.5f;
  const float y = (uv.y * info.height) - 0.5f;

  const float px = floorf(x);
  const float py = floorf(y);
  const float fx = x - px;
  const float fy = y - py;

  const float g0x = cubic_g0(fx);
  const float g1x = cubic_g1(fx);
  /* Note +0.5 offset to compensate for CUDA linear filtering convention. */
  // TODO: turn division into multiplication?
  const float x0 = (px + cubic_h0(fx) + 0.5f) / info.width;
  const float x1 = (px + cubic_h1(fx) + 0.5f) / info.width;
  const float y0 = (py + cubic_h0(fy) + 0.5f) / info.height;
  const float y1 = (py + cubic_h1(fy) + 0.5f) / info.height;

  return cubic_g0(fy) * (g0x * ccl_gpu_image_object_read_2D<T>(tex, x0, y0) +
                         g1x * ccl_gpu_image_object_read_2D<T>(tex, x1, y0)) +
         cubic_g1(fy) * (g0x * ccl_gpu_image_object_read_2D<T>(tex, x0, y1) +
                         g1x * ccl_gpu_image_object_read_2D<T>(tex, x1, y1));
}

ccl_device float4 kernel_image_interp(KernelGlobals kg,
                                      ccl_private ShaderData *sd,
                                      const int tex_id,
                                      float2 uv,
                                      const differential2 duv)
{
  if (tex_id == KERNEL_IMAGE_NONE) {
    return IMAGE_TEXTURE_MISSING_RGBA;
  }

  const ccl_global KernelImageTexture &tex = kernel_data_fetch(image_textures, tex_id);
  const ccl_global KernelImageInfo *info;

  if (tex.tile_descriptor_offset != UINT_MAX) {
    /* Wrapping. */
    if (!kernel_image_tile_wrap(ExtensionType(tex.extension), uv)) {
      return zero_float4();
    }

    /* Tile mapping */
    float2 xy = zero_float2();
    const KernelTileDescriptor tile_descriptor = kernel_image_tile_map(kg, sd, tex, uv, duv, xy);

    if (!kernel_tile_descriptor_loaded(tile_descriptor)) {
      if (tile_descriptor == KERNEL_TILE_LOAD_FAILED) {
        return IMAGE_TEXTURE_MISSING_RGBA;
      }
      // TODO: cancel shader execution
      return tex.average_color;
    }

    info = &kernel_data_fetch(image_info, kernel_tile_descriptor_slot(tile_descriptor));

    /* Convert to normalized space again. */
    // TODO: avoid this, or at least turn division into multiplication
    uv = make_float2(xy.x / info->width, xy.y / info->height);
  }
  else {
    /* Full image sampling. */
    if (tex.slot == KERNEL_IMAGE_NONE) {
      return IMAGE_TEXTURE_MISSING_RGBA;
    }

    info = &kernel_data_fetch(image_info, tex.slot);
  }

  /* float4, byte4, ushort4 and half4 */
  const int texture_type = info->data_type;
  if (texture_type == IMAGE_DATA_TYPE_FLOAT4 || texture_type == IMAGE_DATA_TYPE_BYTE4 ||
      texture_type == IMAGE_DATA_TYPE_HALF4 || texture_type == IMAGE_DATA_TYPE_USHORT4)
  {
    if (info->interpolation == INTERPOLATION_CUBIC || info->interpolation == INTERPOLATION_SMART) {
      return kernel_image_interp_bicubic<float4>(*info, uv);
    }
    else {
      ccl_gpu_image_object_2D tex = (ccl_gpu_image_object_2D)info->data;
      return ccl_gpu_image_object_read_2D<float4>(tex, uv.x, uv.y);
    }
  }
  /* float, byte and half */
  else {
    float f;

    if (info->interpolation == INTERPOLATION_CUBIC || info->interpolation == INTERPOLATION_SMART) {
      f = kernel_image_interp_bicubic<float>(*info, uv);
    }
    else {
      ccl_gpu_image_object_2D tex = (ccl_gpu_image_object_2D)info->data;
      f = ccl_gpu_image_object_read_2D<float>(tex, uv.x, uv.y);
    }

    return make_float4(f, f, f, 1.0f);
  }
}

ccl_device_forceinline float4 kernel_image_interp_with_udim(KernelGlobals kg,
                                                            ccl_private ShaderData *sd,
                                                            const int image_id,
                                                            float2 uv,
                                                            const differential2 duv)
{
  const int tex_id = kernel_image_udim_map(kg, image_id, uv);
  if (tex_id == KERNEL_IMAGE_NONE) {
    return IMAGE_TEXTURE_MISSING_RGBA;
  }

  return kernel_image_interp(kg, sd, tex_id, uv, duv);
}

CCL_NAMESPACE_END
