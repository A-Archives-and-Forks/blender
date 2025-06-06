/* SPDX-FileCopyrightText: 2024 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup edasset
 */

#pragma once

#include "ED_asset_indexer.hh"

struct AssetMetaData;
namespace blender {
class StringRefNull;
}  // namespace blender
namespace blender::io::serialize {
class DictionaryValue;
class Value;
}  // namespace blender::io::serialize

namespace blender::ed::asset::index {

struct RemoteListingAssetEntry;

std::unique_ptr<io::serialize::Value> read_contents(StringRefNull filepath);

AssetMetaData *asset_metadata_from_dictionary(const io::serialize::DictionaryValue &entry);

enum class ReadingResult {
  Success,
  Failure,
  Cancelled,
};

/**
 * Reading of API schema version 1. See #read_remote_listing() on \a process_fn.
 * \param version_root_dirpath: Absolute path to the remote listing root directory.
 * \param version_listing_filepath: Absolute path to the remote listing meta file of this version
 *                                  (e.g. `_v1/asset-listing.json`).
 */
ReadingResult read_remote_listing_v1(StringRefNull listing_root_dirpath,
                                     StringRefNull version_listing_filepath,
                                     RemoteListingEntryProcessFn process_fn);

}  // namespace blender::ed::asset::index
