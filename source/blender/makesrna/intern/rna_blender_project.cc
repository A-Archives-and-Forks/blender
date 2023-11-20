/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup RNA
 */

#include "RNA_define.hh"

#include "rna_internal.h"

#ifdef RNA_RUNTIME

#  include "BLI_string_ref.hh"

#  include "BKE_asset_library_custom.h"
#  include "BKE_blender_project.hh"

#  include "BLT_translation.h"

#  include "WM_api.hh"

using namespace blender;

static void rna_BlenderProject_update(Main * /*bmain*/, Scene * /*scene*/, PointerRNA * /*ptr*/)
{
  /* TODO evaluate which props should send which notifiers. */
  /* Force full redraw of all windows. */
  WM_main_add_notifier(NC_WINDOW, nullptr);
}

static void rna_BlenderProject_name_get(PointerRNA *ptr, char *value)
{
  const bke::BlenderProject *project = static_cast<bke::BlenderProject *>(ptr->data);
  if (!project) {
    value[0] = '\0';
    return;
  }

  strcpy(value, project->project_name().c_str());
}

static int rna_BlenderProject_name_length(PointerRNA *ptr)
{
  const bke::BlenderProject *project = static_cast<bke::BlenderProject *>(ptr->data);
  if (!project) {
    return 0;
  }

  return project->project_name().size();
}

static void rna_BlenderProject_name_set(PointerRNA *ptr, const char *value)
{
  bke::BlenderProject *project = static_cast<bke::BlenderProject *>(ptr->data);

  if (!project) {
    return;
  }

  project->set_project_name(value);
}

static void rna_BlenderProject_root_path_get(PointerRNA *ptr, char *value)
{
  const bke::BlenderProject *project = static_cast<bke::BlenderProject *>(ptr->data);
  if (!project) {
    value[0] = '\0';
    return;
  }

  strcpy(value, project->root_path().c_str());
}

static int rna_BlenderProject_root_path_length(PointerRNA *ptr)
{
  const bke::BlenderProject *project = static_cast<bke::BlenderProject *>(ptr->data);
  if (!project) {
    return 0;
  }

  return project->root_path().size();
}

static void rna_BlenderProject_root_path_set(PointerRNA * /*ptr*/, const char * /*value*/)
{
  /* Property is not editable, see #rna_BlenderProject_root_path_editable(). */
  BLI_assert_unreachable();
}

static int rna_BlenderProject_root_path_editable(PointerRNA * /*ptr*/, const char **r_info)
{
  /* Path is never editable (setting up a project is an operation), but return a nicer disabled
   * hint. */
  *r_info = N_("Project location cannot be changed, displayed for informal purposes only");
  return 0;
}

static void rna_BlenderProject_asset_libraries_begin(CollectionPropertyIterator *iter,
                                                     PointerRNA *ptr)
{
  bke::BlenderProject *project = static_cast<bke::BlenderProject *>(ptr->data);
  ListBase &asset_libraries = project->asset_library_definitions();
  rna_iterator_listbase_begin(iter, &asset_libraries, nullptr);
}

static bool rna_BlenderProject_is_dirty_get(PointerRNA *ptr)
{
  const bke::BlenderProject *project = static_cast<bke::BlenderProject *>(ptr->data);
  return project->has_unsaved_changes();
}

#else

void RNA_def_blender_project(BlenderRNA *brna)
{
  StructRNA *srna = RNA_def_struct(brna, "BlenderProject", nullptr);
  RNA_def_struct_ui_text(srna, "Blender Project", "");

  PropertyRNA *prop;

  prop = RNA_def_property(srna, "name", PROP_STRING, PROP_NONE);
  RNA_def_property_string_funcs(prop,
                                "rna_BlenderProject_name_get",
                                "rna_BlenderProject_name_length",
                                "rna_BlenderProject_name_set");
  RNA_def_property_ui_text(prop, "Name", "The identifier for the project");
  RNA_def_struct_name_property(srna, prop);
  RNA_def_property_update(prop, 0, "rna_BlenderProject_update");

  prop = RNA_def_property(srna, "root_path", PROP_STRING, PROP_NONE);
  RNA_def_property_string_funcs(prop,
                                "rna_BlenderProject_root_path_get",
                                "rna_BlenderProject_root_path_length",
                                "rna_BlenderProject_root_path_set");
  RNA_def_property_editable_func(prop, "rna_BlenderProject_root_path_editable");
  RNA_def_property_ui_text(prop, "Location", "The location of the project on disk");

  prop = RNA_def_property(srna, "asset_libraries", PROP_COLLECTION, PROP_NONE);
  RNA_def_property_struct_type(prop, "CustomAssetLibraryDefinition");
  RNA_def_property_collection_funcs(prop,
                                    "rna_BlenderProject_asset_libraries_begin",
                                    "rna_iterator_listbase_next",
                                    "rna_iterator_listbase_end",
                                    "rna_iterator_listbase_get",
                                    nullptr,
                                    nullptr,
                                    nullptr,
                                    nullptr);
  RNA_def_property_ui_text(prop, "Asset Libraries", "");

  prop = RNA_def_property(srna, "is_dirty", PROP_BOOLEAN, PROP_NONE);
  RNA_def_property_boolean_funcs(prop, "rna_BlenderProject_is_dirty_get", nullptr);
  RNA_def_property_clear_flag(prop, PROP_EDITABLE);
  RNA_def_property_ui_text(
      prop,
      "Dirty",
      "Project settings have changed since read from disk. Save the settings to keep them");
}

#endif
