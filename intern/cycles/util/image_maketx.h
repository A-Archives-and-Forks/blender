/* SPDX-FileCopyrightText: 2025 OpenImageIO project
 * SPDX-FileCopyrightText: 2025 Blender Authors
 * SPDX-License-Identifier: Apache-2.0
 *
 * This is a modified version of maketexture.cpp from OpenImageIO, to add a few
 * features missing in the native implementation. */

#include <cstdlib>
#include <iostream>

#include <OpenImageIO/imagebufalgo.h>

CCL_NAMESPACE_BEGIN

bool make_tx(OIIO::ImageBufAlgo::MakeTextureMode mode,
             OIIO::string_view filename,
             OIIO::string_view outputfilename,
             const OIIO::ImageSpec &configspec,
             std::ostream *outstream);

bool make_tx(OIIO::ImageBufAlgo::MakeTextureMode mode,
             const OIIO::ImageBuf &input,
             OIIO::string_view outputfilename,
             const OIIO::ImageSpec &configspec,
             std::ostream *outstream);

CCL_NAMESPACE_END
