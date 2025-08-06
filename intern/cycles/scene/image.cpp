/* SPDX-FileCopyrightText: 2011-2022 Blender Foundation
 *
 * SPDX-License-Identifier: Apache-2.0 */

#include "scene/image.h"
#include "scene/image_oiio.h"
#include "scene/image_vdb.h"
#include "scene/scene.h"
#include "scene/stats.h"

#include "util/atomic.h"
#include "util/colorspace.h"
#include "util/image.h"
#include "util/image_impl.h"
#include "util/log.h"
#include "util/progress.h"
#include "util/task.h"
#include "util/texture.h"
#include "util/types_base.h"

#include <OpenImageIO/thread.h>

CCL_NAMESPACE_BEGIN

/* Image Handle */

ImageHandle::ImageHandle() = default;

ImageHandle::ImageHandle(ImageSlot *image_slot, ImageManager *manager)
    : image_slot(image_slot), manager(manager)
{
  if (image_slot) {
    image_slot->users++;
  }
}

ImageHandle::ImageHandle(const ImageHandle &other)
    : image_slot(other.image_slot), manager(other.manager)
{
  if (image_slot) {
    image_slot->users++;
  }
}

ImageHandle::ImageHandle(ImageHandle &&other)
    : image_slot(other.image_slot), manager(other.manager)
{
  if (&other == this) {
    abort();
  }
  other.image_slot = nullptr;
  other.manager = nullptr;
}

ImageHandle &ImageHandle::operator=(const ImageHandle &other)
{
  clear();
  image_slot = other.image_slot;
  manager = other.manager;

  if (image_slot) {
    image_slot->users++;
  }

  return *this;
}

ImageHandle &ImageHandle::operator=(ImageHandle &&other)
{
  clear();
  image_slot = other.image_slot;
  manager = other.manager;
  other.image_slot = nullptr;
  other.manager = nullptr;

  return *this;
}

ImageHandle::~ImageHandle()
{
  clear();
}

void ImageHandle::clear()
{
  /* Don't remove immediately, rather do it all together later on. one of
   * the reasons for this is that on shader changes we add and remove nodes
   * that use them, but we do not want to reload the image all the time. */
  if (image_slot) {
    assert(image_slot->users >= 1);
    image_slot->users--;
    if (image_slot->users == 0) {
      manager->tag_update();
    }
    image_slot = nullptr;
  }

  manager = nullptr;
}

bool ImageHandle::empty() const
{
  return image_slot == nullptr;
}

int ImageHandle::num_tiles() const
{
  if (image_slot && image_slot->type == ImageSlot::UDIM) {
    ImageUDIM *udim = static_cast<ImageUDIM *>(image_slot);
    return udim->tiles.size();
  }

  return 0;
}

ImageMetaData ImageHandle::metadata(Progress &progress)
{
  if (image_slot) {
    if (image_slot->type == ImageSlot::SINGLE) {
      ImageSingle *img = static_cast<ImageSingle *>(image_slot);
      manager->load_image_metadata(img, progress);
      return img->metadata;
    }
    if (image_slot->type == ImageSlot::UDIM) {
      ImageUDIM *udim = static_cast<ImageUDIM *>(image_slot);
      return udim->tiles[0].second.metadata(progress);
    }
  }

  return ImageMetaData();
}

int ImageHandle::kernel_id() const
{
  return (image_slot) ? image_slot->id : KERNEL_IMAGE_NONE;
}

device_image *ImageHandle::vdb_image_memory() const
{
  if (image_slot == nullptr || image_slot->type != ImageSlot::SINGLE) {
    return nullptr;
  }

  ImageSingle *img = static_cast<ImageSingle *>(image_slot);
  return img->vdb_memory;
}

VDBImageLoader *ImageHandle::vdb_loader() const
{
  if (image_slot == nullptr || image_slot->type != ImageSlot::SINGLE) {
    return nullptr;
  }

  ImageSingle *img = static_cast<ImageSingle *>(image_slot);
  ImageLoader *loader = img->loader.get();
  if (loader == nullptr) {
    return nullptr;
  }

  if (loader->is_vdb_loader()) {
    return dynamic_cast<VDBImageLoader *>(loader);
  }

  return nullptr;
}

ImageManager *ImageHandle::get_manager() const
{
  return manager;
}

bool ImageHandle::operator==(const ImageHandle &other) const
{
  return image_slot == other.image_slot && manager == other.manager;
}

/* Image Loader */

ImageLoader::ImageLoader() = default;

int ImageLoader::get_tile_number() const
{
  return 0;
}

bool ImageLoader::equals(const ImageLoader *a, const ImageLoader *b)
{
  if (a == nullptr && b == nullptr) {
    return true;
  }
  return (a && b && typeid(*a) == typeid(*b) && a->equals(*b));
}

bool ImageLoader::is_vdb_loader() const
{
  return false;
}

/* Image Manager */

ImageManager::ImageManager(const DeviceInfo & /*info*/, const SceneParams &params)
{
  use_texture_cache = params.use_texture_cache;
  auto_texture_cache = params.auto_texture_cache;
  texture_cache_path = params.texture_cache_path;
}

ImageManager::~ImageManager()
{
  for (size_t slot = 0; slot < images.size(); slot++) {
    assert(!images[slot]);
  }
}

bool ImageManager::set_animation_frame_update(const int frame)
{
  if (frame != animation_frame) {
    const thread_scoped_lock device_lock(images_mutex);
    animation_frame = frame;

    for (size_t slot = 0; slot < images.size(); slot++) {
      if (images[slot] && images[slot]->params.animated) {
        return true;
      }
    }
  }

  return false;
}

void ImageManager::load_image_metadata(ImageSingle *img, Progress &progress)
{
  if (!img->need_metadata) {
    return;
  }

  const thread_scoped_lock image_lock(img->mutex);
  if (!img->need_metadata) {
    return;
  }

  /* Change image to use tx file if supported. */
  /* TODO: Can we further delay auto generating until we are certain the image is used.
   * On the other hand, shold we move this earlier so explicit use of the same tx file
   * and source iamge is deduplicated? */
  if (use_texture_cache) {
    img->loader->resolve_texture_cache(auto_texture_cache,
                                       texture_cache_path,
                                       img->params.colorspace,
                                       img->params.alpha_type,
                                       progress);
  }

  ImageMetaData &metadata = img->metadata;
  metadata = ImageMetaData();
  metadata.colorspace = img->params.colorspace;

  if (img->loader->load_metadata(metadata)) {
    assert(metadata.type != IMAGE_DATA_NUM_TYPES);
  }
  else {
    metadata.type = IMAGE_DATA_TYPE_BYTE4;
  }

  metadata.finalize(img->params.alpha_type);

  img->need_metadata = false;
}

ImageHandle ImageManager::add_image(const string &filename, const ImageParams &params)
{
  ImageSingle *image = add_image_slot(make_unique<OIIOImageLoader>(filename), params, false);
  return ImageHandle(image, this);
}

ImageHandle ImageManager::add_image(const string &filename,
                                    const ImageParams &params,
                                    const array<int> &tiles)
{
  if (tiles.empty()) {
    return add_image(filename, params);
  }

  vector<std::pair<int, ImageHandle>> udim_tiles;
  for (const int tile : tiles) {
    string tile_filename = filename;

    /* Since we don't have information about the exact tile format used in this code location,
     * just attempt all replacement patterns that Blender supports. */
    string_replace(tile_filename, "<UDIM>", string_printf("%04d", tile));

    const int u = ((tile - 1001) % 10);
    const int v = ((tile - 1001) / 10);
    string_replace(tile_filename, "<UVTILE>", string_printf("u%d_v%d", u + 1, v + 1));

    ImageSingle *image = add_image_slot(
        make_unique<OIIOImageLoader>(tile_filename), params, false);
    udim_tiles.emplace_back(tile, ImageHandle(image, this));
  }

  ImageUDIM *udim = add_image_slot(std::move(udim_tiles));
  return ImageHandle(udim, this);
}

ImageHandle ImageManager::add_image(unique_ptr<ImageLoader> &&loader,
                                    const ImageParams &params,
                                    const bool builtin)
{
  ImageSingle *image = add_image_slot(std::move(loader), params, builtin);
  return ImageHandle(image, this);
}

ImageHandle ImageManager::add_image(vector<unique_ptr<ImageLoader>> &&loaders,
                                    const ImageParams &params)
{
  vector<std::pair<int, ImageHandle>> udim_tiles;

  for (unique_ptr<ImageLoader> &loader : loaders) {
    unique_ptr<ImageLoader> local_loader;
    std::swap(loader, local_loader);
    ImageSingle *image = add_image_slot(std::move(local_loader), params, true);
    udim_tiles.emplace_back(image->loader->get_tile_number(), ImageHandle(image, this));
  }

  ImageUDIM *udim = add_image_slot(std::move(udim_tiles));
  return ImageHandle(udim, this);
}

/* ImageManager */

ImageSingle *ImageManager::add_image_slot(unique_ptr<ImageLoader> &&loader,
                                          const ImageParams &params,
                                          const bool builtin)
{
  const thread_scoped_lock device_lock(images_mutex);

  /* Find existing image. */
  size_t slot;
  for (slot = 0; slot < images.size(); slot++) {
    ImageSingle *img = images[slot].get();
    if (img && ImageLoader::equals(img->loader.get(), loader.get()) && img->params == params) {
      return img;
    }
  }

  /* Find free slot. */
  for (slot = 0; slot < images.size(); slot++) {
    if (!images[slot]) {
      break;
    }
  }

  if (slot == images.size()) {
    images.resize(images.size() + 1);
  }

  /* Add new image. */
  unique_ptr<ImageSingle> img = make_unique<ImageSingle>();
  img->type = ImageSlot::SINGLE;
  img->id = slot;
  img->params = params;
  img->loader = std::move(loader);
  img->builtin = builtin;

  images[slot] = std::move(img);

  tag_update();

  return images[slot].get();
}

ImageUDIM *ImageManager::add_image_slot(vector<std::pair<int, ImageHandle>> &&tiles)
{
  const thread_scoped_lock device_lock(images_mutex);

  /* Find existing UDIM. */
  size_t slot;
  for (slot = 0; slot < image_udims.size(); slot++) {
    ImageUDIM *udim = image_udims[slot].get();
    if (udim && udim->tiles == tiles) {
      return udim;
    }
  }

  /* Find free slot. */
  for (slot = 0; slot < image_udims.size(); slot++) {
    if (!image_udims[slot]) {
      break;
    }
  }

  if (slot == image_udims.size()) {
    image_udims.resize(image_udims.size() + 1);
  }

  /* Add new image. */
  unique_ptr<ImageUDIM> img = make_unique<ImageUDIM>();
  img->type = ImageSlot::UDIM;
  img->id = -num_udim_tiles - 1;
  img->tiles = std::move(tiles);

  num_udim_tiles += img->tiles.size() + 1;

  image_udims[slot] = std::move(img);

  tag_update();

  return image_udims[slot].get();
}

template<TypeDesc::BASETYPE FileFormat, typename StorageType>
bool ImageManager::file_load_image(Device *device, ImageSingle *img, const int texture_limit)
{
  /* Ignore empty images. */
  if (!(img->metadata.channels > 0)) {
    return false;
  }

  /* Get metadata. */
  const int width = img->metadata.width;
  const int height = img->metadata.height;

  /* Read pixels. */
  vector<StorageType> pixels_storage;
  StorageType *pixels;
  const int64_t max_size = max(width, height);
  if (max_size == 0) {
    /* Don't bother with empty images. */
    return false;
  }

  /* Allocate memory as needed, may be smaller to resize down. */
  if (texture_limit > 0 && max_size > texture_limit) {
    pixels_storage.resize(((int64_t)width) * height * 4);
    pixels = &pixels_storage[0];
  }
  else {
    pixels = image_cache
                 .alloc_full(device,
                             img->metadata.type,
                             img->params.interpolation,
                             img->params.extension,
                             width,
                             height,
                             img->texture_slot)
                 .data<StorageType>();
  }

  if (pixels == nullptr) {
    /* Could be that we've run out of memory. */
    return false;
  }

  if (!img->loader->load_pixels_full(img->metadata, (uint8_t *)pixels)) {
    return false;
  }

  /* Scale image down if needed. */
  if (!pixels_storage.empty()) {
    float scale_factor = 1.0f;
    while (max_size * scale_factor > texture_limit) {
      scale_factor *= 0.5f;
    }
    LOG_WORK << "Scaling image " << img->loader->name() << " by a factor of " << scale_factor
             << ".";
    vector<StorageType> scaled_pixels;
    int64_t scaled_width;
    int64_t scaled_height;

    util_image_resize_pixels(pixels_storage,
                             width,
                             height,
                             img->metadata.is_rgba() ? 4 : 1,
                             scale_factor,
                             &scaled_pixels,
                             &scaled_width,
                             &scaled_height);

    StorageType *texture_pixels = image_cache
                                      .alloc_full(device,
                                                  img->metadata.type,
                                                  img->params.interpolation,
                                                  img->params.extension,
                                                  scaled_width,
                                                  scaled_height,
                                                  img->texture_slot)
                                      .data<StorageType>();
    std::copy_n(scaled_pixels.data(), scaled_pixels.size(), texture_pixels);
  }

  return true;
}

void ImageManager::device_resize_image_textures(Scene *scene)
{
  const thread_scoped_lock device_lock(device_mutex);
  DeviceScene &dscene = scene->dscene;

  if (dscene.image_textures.size() < images.size()) {
    dscene.image_textures.resize(images.size());
  }
}

void ImageManager::device_copy_image_textures(Device *device, Scene *scene)
{
  image_cache.copy_to_device_if_modified();

  const thread_scoped_lock device_lock(device_mutex);
  DeviceScene &dscene = scene->dscene;

  dscene.image_textures.copy_to_device_if_modified();
  dscene.image_texture_tile_descriptors.copy_to_device_if_modified();
  dscene.image_texture_udims.copy_to_device_if_modified();

  device->set_image_cache_func(
      /* This is called by the CPU kernel to immediately load a tile. */
      /* TODO: resizing image_info is not thread safe for CPU rendering, there is
       * a workaround in the CPUDevice constructor. */
      /* TODO: drop file handle. */
      [this, device, scene](
          size_t slot, int miplevel, int x, int y, KernelTileDescriptor *tile_descriptor) {
        /* If we can atomically set KERNEL_TILE_LOAD_REQUEST, this thread is responsible
         * for loading the tile. */
        KernelTileDescriptor tile_descriptor_old = *tile_descriptor;
        if (tile_descriptor_old != KERNEL_TILE_LOAD_REQUEST &&
            tile_descriptor_old ==
                atomic_cas_uint32(tile_descriptor, tile_descriptor_old, KERNEL_TILE_LOAD_REQUEST))
        {
          KernelTileDescriptor tile_descriptor_new = device_update_tile_requested(
              device, scene, images[slot].get(), miplevel, x, y);
          /* TODO: this is inefficient, and with CPU + GPU rendering we only want to update CPU. */
          image_cache.copy_to_device_if_modified();
          *tile_descriptor = tile_descriptor_new;
          return;
        }

        /* Wait for other thread to load the tile. */
        OIIO::atomic_backoff backoff;
        while (*tile_descriptor == KERNEL_TILE_LOAD_REQUEST) {
          backoff();
        }
      },
      [this, device, scene] { this->device_update_requested(device, scene); });
}

void ImageManager::device_load_image_full(Device *device, Scene *scene, const size_t slot)
{
  ImageSingle *img = images[slot].get();

  const ImageDataType type = img->metadata.type;
  const int texture_limit = scene->params.texture_limit;

  /* Create new texture. */
  switch (type) {
    case IMAGE_DATA_TYPE_FLOAT4:
      file_load_image<TypeDesc::FLOAT, float>(device, img, texture_limit);
      break;
    case IMAGE_DATA_TYPE_FLOAT:
      file_load_image<TypeDesc::FLOAT, float>(device, img, texture_limit);
      break;
    case IMAGE_DATA_TYPE_BYTE4:
      file_load_image<TypeDesc::UINT8, uchar>(device, img, texture_limit);
      break;
    case IMAGE_DATA_TYPE_BYTE:
      file_load_image<TypeDesc::UINT8, uchar>(device, img, texture_limit);
      break;
    case IMAGE_DATA_TYPE_HALF4:
      file_load_image<TypeDesc::HALF, half>(device, img, texture_limit);
      break;
    case IMAGE_DATA_TYPE_HALF:
      file_load_image<TypeDesc::HALF, half>(device, img, texture_limit);
      break;
    case IMAGE_DATA_TYPE_USHORT:
      file_load_image<TypeDesc::USHORT, uint16_t>(device, img, texture_limit);
      break;
    case IMAGE_DATA_TYPE_USHORT4:
      file_load_image<TypeDesc::USHORT, uint16_t>(device, img, texture_limit);
      break;
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT3:
    case IMAGE_DATA_TYPE_NANOVDB_FLOAT4:
    case IMAGE_DATA_TYPE_NANOVDB_FPN:
    case IMAGE_DATA_TYPE_NANOVDB_FP16:
    case IMAGE_DATA_TYPE_NANOVDB_EMPTY: {
#ifdef WITH_NANOVDB
      img->vdb_memory = &image_cache.alloc_full(device,
                                                type,
                                                img->params.interpolation,
                                                img->params.extension,
                                                img->metadata.byte_size,
                                                0,
                                                img->texture_slot);

      uint8_t *pixels = img->vdb_memory->data<uint8_t>();
      if (pixels) {
        img->loader->load_pixels_full(img->metadata, pixels);
      }
#endif
      break;
    }
    case IMAGE_DATA_NUM_TYPES:
      break;
  }
}

void ImageManager::device_load_image_tiled(Scene *scene, const size_t slot)
{
  ImageSingle *img = images[slot].get();

  vector<KernelTileDescriptor> levels;
  const int max_miplevels = img->params.interpolation != INTERPOLATION_CLOSEST ? 1 : INT_MAX;
  const int tile_size = img->metadata.tile_size;

  int num_tiles = 0;

  for (int miplevel = 0; max_miplevels; miplevel++) {
    const int width = divide_up(img->metadata.width, 1 << miplevel);
    const int height = divide_up(img->metadata.height, 1 << miplevel);

    levels.push_back(num_tiles);

    num_tiles += divide_up(width, tile_size) * divide_up(height, tile_size);

    if (width <= tile_size && height <= tile_size) {
      break;
    }
  }

  {
    // TODO: make this more efficient
    const thread_scoped_lock device_lock(device_mutex);
    const int tile_descriptor_offset = scene->dscene.image_texture_tile_descriptors.size();
    scene->dscene.image_texture_tile_descriptors.resize(tile_descriptor_offset + levels.size() +
                                                        num_tiles);

    KernelTileDescriptor *descr_data = scene->dscene.image_texture_tile_descriptors.data() +
                                       tile_descriptor_offset;

    for (int i = 0; i < levels.size(); i++) {
      descr_data[i] = levels.size() + levels[i];
    }
    std::fill_n(descr_data + levels.size(), num_tiles, KERNEL_TILE_LOAD_NONE);

    img->tile_descriptor_offset = tile_descriptor_offset;
    img->tile_descriptor_levels = levels.size();
    img->tile_descriptor_num = num_tiles;
  }
}

KernelTileDescriptor ImageManager::device_update_tile_requested(Device *device,
                                                                Scene *scene,
                                                                ImageSingle *img,
                                                                const int miplevel,
                                                                const size_t x,
                                                                const size_t y)
{
  const int width = divide_up(img->metadata.width, 1 << miplevel);
  const int height = divide_up(img->metadata.height, 1 << miplevel);
  const size_t tile_size = img->metadata.tile_size;
  const size_t w = min(width - x, tile_size);
  const size_t h = min(height - y, tile_size);
  const size_t tile_size_padded = tile_size + KERNEL_IMAGE_TEX_PADDING * 2;

  KernelTileDescriptor tile_descriptor;

  device_image &mem = image_cache.alloc_tile(device,
                                             img->metadata.type,
                                             img->params.interpolation,
                                             img->params.extension,
                                             tile_size_padded,
                                             tile_descriptor);

  const size_t pixel_bytes = mem.data_elements * datatype_size(mem.data_type);
  /* TODO: Handle case where channels > 4. */
  const size_t x_stride = pixel_bytes;
  const size_t y_stride = mem.data_width * pixel_bytes;
  const size_t x_offset = kernel_tile_descriptor_offset(tile_descriptor) * tile_size_padded *
                          pixel_bytes;

  uint8_t *pixels = mem.data<uint8_t>() + x_offset;

  /* TODO: opening file handles is not thread safe! */
  const bool ok = img->loader->load_pixels_tile(img->metadata,
                                                miplevel,
                                                x,
                                                y,
                                                w,
                                                h,
                                                x_stride,
                                                y_stride,
                                                KERNEL_IMAGE_TEX_PADDING,
                                                img->params.extension,
                                                pixels);

  scene->dscene.image_texture_tile_descriptors.tag_modified();

  if (ok) {
    LOG_DEBUG << "Load image tile: " << img->loader->name() << ", mip level " << miplevel << " ("
              << x << " " << y << ")";
  }
  else {
    LOG_WARNING << "Failed to load image tile: " << img->loader->name() << ", mip level "
                << miplevel << " (" << x << " " << y << ")";
  }

  return (ok) ? tile_descriptor : KERNEL_TILE_LOAD_FAILED;
}

void ImageManager::device_update_image_requested(Device *device, Scene *scene, ImageSingle *img)
{
  const size_t tile_size = img->metadata.tile_size;

  KernelTileDescriptor *tile_descriptors = scene->dscene.image_texture_tile_descriptors.data() +
                                           img->tile_descriptor_offset +
                                           img->tile_descriptor_levels;

  size_t i = 0;
  for (int miplevel = 0; miplevel < img->tile_descriptor_levels; miplevel++) {
    const int width = divide_up(img->metadata.width, 1 << miplevel);
    const int height = divide_up(img->metadata.height, 1 << miplevel);

    for (size_t y = 0; y < height; y += tile_size) {
      for (size_t x = 0; x < width; x += tile_size, i++) {
        assert(i < img->tile_descriptor_num);

        if (tile_descriptors[i] != KERNEL_TILE_LOAD_REQUEST) {
          continue;
        }

        tile_descriptors[i] = device_update_tile_requested(device, scene, img, miplevel, x, y);
      }
    }
  }

  img->loader->drop_file_handle();
}

void ImageManager::device_load_image(Device *device,
                                     Scene *scene,
                                     const size_t slot,
                                     Progress &progress)
{
  if (progress.get_cancel()) {
    return;
  }

  ImageSingle *img = images[slot].get();

  progress.set_status("Updating Images", "Loading " + img->loader->name());

  load_image_metadata(img, progress);

  KernelImageTexture tex;
  tex.width = img->metadata.width;
  tex.height = img->metadata.height;
  tex.interpolation = img->params.interpolation;
  tex.extension = img->params.extension;
  tex.use_transform_3d = img->metadata.use_transform_3d;
  tex.transform_3d = img->metadata.transform_3d;
  tex.average_color = img->metadata.average_color;

  if (use_texture_cache && img->metadata.tile_size) {
    assert(is_power_of_two(img->metadata.tile_size));

    device_load_image_tiled(scene, slot);

    tex.tile_descriptor_offset = img->tile_descriptor_offset;
    tex.tile_size_shift = __bsr(img->metadata.tile_size);
    tex.tile_levels = img->tile_descriptor_levels;
    tex.slot = slot;
  }
  else {
    device_load_image_full(device, scene, slot);
    tex.slot = img->texture_slot;
  }

  /* Update image texture device data. */
  scene->dscene.image_textures[slot] = tex;
  scene->dscene.image_textures.tag_modified();

  /* Cleanup memory in image loader. */
  img->loader->cleanup();
  img->need_load = false;
}

void ImageManager::device_free_image(Scene *scene, size_t slot)
{
  ImageSingle *img = images[slot].get();
  if (img == nullptr) {
    return;
  }

  if (img->texture_slot != KERNEL_IMAGE_NONE) {
    image_cache.free_full(img->texture_slot);
  }
  if (img->tile_descriptor_offset != KERNEL_IMAGE_NONE) {
    // TODO: shrink image_texture_tile_descriptors
    KernelTileDescriptor *tile_descriptors = scene->dscene.image_texture_tile_descriptors.data() +
                                             img->tile_descriptor_offset +
                                             img->tile_descriptor_levels;

    for (int i = 0; i < img->tile_descriptor_num; i++) {
      if (kernel_tile_descriptor_loaded(tile_descriptors[i])) {
        image_cache.free_tile(tile_descriptors[i]);
      }
    }
  }

  images[slot].reset();
}

void ImageManager::device_update_requested(Device *device, Scene *scene)
{
  // TODO: not supported for MEM_GLOBAL
  // TODO: only do if modified
  // scene->dscene.image_texture_tile_descriptors.copy_from_device();

  parallel_for(blocked_range<size_t>(0, images.size(), 1), [&](const blocked_range<size_t> &r) {
    for (size_t i = r.begin(); i != r.end(); i++) {
      unique_ptr<ImageSingle> &img = images[i];
      if (img && img->tile_descriptor_offset != KERNEL_IMAGE_NONE) {
        device_update_image_requested(device, scene, img.get());
      }
    }
  });

  device_copy_image_textures(device, scene);
}

void ImageManager::device_update_udims(Device * /*device*/, Scene *scene)
{
  // TODO: Shrink image_texture_udims
  const thread_scoped_lock device_lock(device_mutex);
  device_vector<KernelImageUDIM> &device_udims = scene->dscene.image_texture_udims;
  if (device_udims.size() == num_udim_tiles) {
    return;
  }

  device_udims.resize(num_udim_tiles);

  for (size_t slot = 0; slot < image_udims.size(); slot++) {
    ImageUDIM *udim = image_udims[slot].get();
    if (udim) {
      if (udim->users == 0) {
        image_udims[slot].reset();
      }
      else if (udim->need_load) {
        const uint udim_offset = -udim->id - 1;
        KernelImageUDIM *udim_data = device_udims.data() + udim_offset;

        udim_data[0] = KernelImageUDIM{int(udim->tiles.size()), 0};
        for (int i = 0; i < udim->tiles.size(); i++) {
          const auto &tile = udim->tiles[i];
          udim_data[i + 1] = KernelImageUDIM{tile.first, tile.second.kernel_id()};
        }
        udim->need_load = false;
      }
    }
  }
}

void ImageManager::device_update(Device *device, Scene *scene, Progress &progress)
{
  if (!need_update()) {
    return;
  }

  const scoped_callback_timer timer([scene](double time) {
    if (scene->update_stats) {
      scene->update_stats->image.times.add_entry({"device_update", time});
    }
  });

  /* Update UDIM ids. */
  device_update_udims(device, scene);

  /* Resize devices arrays to match. */
  device_resize_image_textures(scene);

  /* Free and load images. */
  TaskPool pool;
  for (size_t slot = 0; slot < images.size(); slot++) {
    ImageSingle *img = images[slot].get();
    if (img && img->users == 0) {
      device_free_image(scene, slot);
    }
    else if (img && img->need_load) {
      pool.push([this, device, scene, slot, &progress] {
        device_load_image(device, scene, slot, progress);
      });
    }
  }

  pool.wait_work();

  /* Copy device arrays. */
  device_copy_image_textures(device, scene);

  need_update_ = false;
}

void ImageManager::device_load_handles(Device *device,
                                       Scene *scene,
                                       Progress &progress,
                                       const set<const ImageHandle *> &handles)
{
  /* Update UDIM ids. */
  device_update_udims(device, scene);

  /* Resize devices arrays to match number of images. */
  device_resize_image_textures(scene);

  /* Load handles. */
  TaskPool pool;
  auto load_image = [&](const ImageSingle *img) {
    pool.push([this, device, scene, img, &progress] {
      assert(img != nullptr);
      if (img->users == 0) {
        device_free_image(scene, img->id);
      }
      else if (img->need_load) {
        device_load_image(device, scene, img->id, progress);
      }
    });
  };

  for (const ImageHandle *handle : handles) {
    if (handle->empty()) {
      continue;
    }

    if (handle->image_slot->type == ImageSlot::SINGLE) {
      load_image(static_cast<const ImageSingle *>(handle->image_slot));
    }
    else {
      for (const auto &tile : static_cast<const ImageUDIM *>(handle->image_slot)->tiles) {
        load_image(static_cast<const ImageSingle *>(tile.second.image_slot));
      }
    }
  }
  pool.wait_work();

  /* Copy device arrays. */
  device_copy_image_textures(device, scene);
}

void ImageManager::device_load_builtin(Device *device, Scene *scene, Progress &progress)
{
  /* Load only builtin images, Blender needs this to load evaluated
   * scene data from depsgraph before it is freed. */
  if (!need_update()) {
    return;
  }

  device_resize_image_textures(scene);

  TaskPool pool;
  for (size_t slot = 0; slot < images.size(); slot++) {
    ImageSingle *img = images[slot].get();
    if (img && img->need_load && img->builtin) {
      pool.push([this, device, scene, slot, &progress] {
        device_load_image(device, scene, slot, progress);
      });
    }
  }

  pool.wait_work();

  device_copy_image_textures(device, scene);
}

void ImageManager::device_free_builtin(Scene *scene)
{
  image_udims.clear();
  for (size_t slot = 0; slot < images.size(); slot++) {
    ImageSingle *img = images[slot].get();
    if (img && img->builtin) {
      device_free_image(scene, slot);
    }
  }
}

void ImageManager::device_free(Scene *scene)
{
  image_udims.clear();
  for (size_t slot = 0; slot < images.size(); slot++) {
    device_free_image(scene, slot);
  }
  images.clear();
  image_cache.device_free();
  scene->dscene.image_textures.free();
  scene->dscene.image_texture_tile_descriptors.free();
  scene->dscene.image_texture_udims.free();
}

void ImageManager::collect_statistics(RenderStats *stats)
{
  for (const unique_ptr<ImageSingle> &image : images) {
    if (!image) {
      /* Image may have been freed due to lack of users. */
      continue;
    }

    // TODO: collection from image cache
    stats->image.textures.add_entry(NamedSizeEntry(image->loader->name(), 0));
  }
}

void ImageManager::tag_update()
{
  need_update_ = true;
}

bool ImageManager::need_update() const
{
  return need_update_;
}

CCL_NAMESPACE_END
