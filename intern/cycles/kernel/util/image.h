
/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#pragma once

#include "kernel/globals.h"
#include "kernel/sample/lcg.h"

#include "util/defines.h"
#include "util/math_fast.h"
#include "util/texture.h"

CCL_NAMESPACE_BEGIN

ccl_device_forceinline int kernel_image_udim_map(KernelGlobals kg,
                                                 const int id,
                                                 ccl_private float2 &uv)
{
  if (id >= 0) {
    return id;
  }

  const int tx = (int)uv.x;
  const int ty = (int)uv.y;
  if (tx < 0 || ty < 0 || tx >= 10) {
    return KERNEL_IMAGE_NONE;
  }

  int udim_id = -id - 1;
  const int num_udims = kernel_data_fetch(image_texture_udims, udim_id++).tile;
  const int tile = 1001 + 10 * ty + tx;

  for (int i = 0; i < num_udims; i++) {
    const KernelImageUDIM udim = kernel_data_fetch(image_texture_udims, udim_id++);
    if (udim.tile == tile) {
      /* If we found the tile, offset the UVs to be relative to it. */
      uv.x -= tx;
      uv.y -= ty;
      return udim.tex_id;
    }
  }

  return KERNEL_IMAGE_NONE;
}

ccl_device_forceinline bool kernel_image_tile_wrap(const ExtensionType extension,
                                                   ccl_private float2 &uv)
{
  /* Wrapping. */
  switch (extension) {
    case EXTENSION_REPEAT:
      uv = uv - floor(uv);
      return true;
    case EXTENSION_CLIP:
      // TODO: this is not exactly the same as full image wrapping
      return (uv.x >= 0.0f && uv.x <= 1.0f && uv.y >= 0.0f && uv.y <= 1.0f);
    case EXTENSION_EXTEND:
      uv = clamp(uv, zero_float2(), one_float2());
      return true;
    case EXTENSION_MIRROR: {
      // TODO: replace fmod with uv.x - floor(uv.x)?
      uv.x = fmodf(fabsf(uv.x), 2.0f);
      if (uv.x >= 1.0f) {
        uv.x = 2.0f - uv.x;
      }
      uv.y = fmodf(fabsf(uv.y), 2.0f);
      if (uv.y >= 1.0f) {
        uv.y = 2.0f - uv.y;
      }
      return true;
    }
    default:
      break;
  }

  return false;
}

/* From UV coordinates in 0..1 range, compute tile and pixel coordinates. */
ccl_device_forceinline KernelTileDescriptor
kernel_image_tile_map(KernelGlobals kg,
                      ccl_private ShaderData *sd,
                      const ccl_global KernelImageTexture &tex,
                      const float2 uv,
                      const differential2 duv,
                      ccl_private float2 &xy)
{
  /* Find mipmap level. */

  // TODO: make this faster
  // TODO: find good ratio, check how OIIO handles anisotropy
  const float dudxy = len(make_float2(duv.dx.x, duv.dy.x)) * float(tex.width);
  const float dvdxy = len(make_float2(duv.dx.y, duv.dy.y)) * float(tex.height);

  /* Limit max anisotropy ratio, to avoid loading too high mip resolutions
   * for stretched UV coordinates, which don't really benefit from it anyway. */
  const float maxdxy = max(dudxy, dvdxy);
  const float mindxy = min(dudxy, dvdxy);
  const float max_aniso_ratio = 16.0f;
  const float sampledxy = max(mindxy, maxdxy * (1.0f / max_aniso_ratio));

  /* Select mipmap level. */
  float flevel = fast_log2f(sampledxy);
  if (sd->lcg_state != 0) {
    /* Randomize mip level, except for some cases like displacement or background. */
    /* TODO: What is the appropriate value here? Smaller means fewer tiles loaded but also
     * potential artifacts. */
    const float transition = 0.5f;
    flevel += (lcg_step_float(&sd->lcg_state) - 0.5f) * transition;
  }
  const int level = clamp(int(flevel), 0, tex.tile_levels - 1);

  /* Compute width of this mipmap level. */
  const int width = divide_up_by_shift(tex.width, level);
  const int height = divide_up_by_shift(tex.height, level);

  /* Convert coordinates to pixel space. */
  xy = uv * make_float2(width, height);

  /* Tile mapping */
  const int ix = clamp((int)xy.x, 0, width - 1);
  const int iy = clamp((int)xy.y, 0, height - 1);
  const int tile_size_shift = tex.tile_size_shift;
  const int tile_size_padded = (1 << tile_size_shift) + KERNEL_IMAGE_TEX_PADDING * 2;

  const int tile_x = ix >> tile_size_shift;
  const int tile_y = iy >> tile_size_shift;
  const int tile_offset = kernel_data_fetch(image_texture_tile_descriptors,
                                            tex.tile_descriptor_offset + level) +
                          tile_x + tile_y * divide_up_by_shift(width, tile_size_shift);

  const KernelTileDescriptor tile_descriptor = kernel_data_fetch(
      image_texture_tile_descriptors, tex.tile_descriptor_offset + tile_offset);

  if (!kernel_tile_descriptor_loaded(tile_descriptor)) {
#ifdef __KERNEL_GPU__
    /* For GPU, mark load requested and cancel shader execution. */
    if (tile_descriptor == KERNEL_TILE_LOAD_NONE) {
      kernel_data_assign(image_texture_tile_descriptors,
                         tex.tile_descriptor_offset + tile_offset,
                         KERNEL_TILE_LOAD_REQUEST);
    }
    return tile_descriptor;
#else
    kg->image_load_tile(
        tex.slot,
        level,
        tile_x * (1 << tile_size_shift),
        tile_y * (1 << tile_size_shift),
        &kg->image_texture_tile_descriptors.data[tex.tile_descriptor_offset + tile_offset]);

    /* For CPU, load tile immediately. */
    if (!kernel_tile_descriptor_loaded(tile_descriptor)) {
      return tile_descriptor;
    }
#endif
  }

  /* Remap coordinates into tiled image space. */
  const int offset = kernel_tile_descriptor_offset(tile_descriptor);
  xy += make_float2(KERNEL_IMAGE_TEX_PADDING - (tile_x << tile_size_shift) +
                        offset * tile_size_padded,
                    KERNEL_IMAGE_TEX_PADDING - (tile_y << tile_size_shift));

  return tile_descriptor;
}

CCL_NAMESPACE_END
