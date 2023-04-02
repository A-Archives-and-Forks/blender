/* SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup bke
 */

/*

TODO:

Convergence improvements:
1. DONE: Limit number of edges processed per run.
2. DONE: Scale split steps by ratio of long to short edges to
   prevent runaway tesselation.
3. DONE: Detect and dissolve three and four valence vertices that are surrounded by
   all tris.
4. DONE: Use different (coarser) brush spacing for applying dyntopo

Drawing improvements:
4. PARTIAL DONE: Build and cache vertex index buffers, to reduce GPU bandwidth

Topology rake:
5. DONE: Enable new curvature topology rake code and add to UI.
6. DONE: Add code to cache curvature data per vertex in a CD layer.

*/

#include "MEM_guardedalloc.h"

#include "BLI_alloca.h"
#include "BLI_buffer.h"
#include "BLI_ghash.h"
#include "BLI_hash.h"
#include "BLI_heap_simple.h"
#include "BLI_math.h"
#include "BLI_memarena.h"
#include "BLI_rand.h"
#include "BLI_sort_utils.h"
#include "BLI_task.h"
#include "BLI_utildefines.h"

#include "BLI_index_range.hh"
#include "BLI_map.hh"
#include "BLI_math_vector_types.hh"
#include "BLI_set.hh"
#include "BLI_vector.hh"

#include "PIL_time.h"
#include "atomic_ops.h"

#include "DNA_material_types.h"
#include "DNA_mesh_types.h"

#include "BKE_DerivedMesh.h"
#include "BKE_ccg.h"
#include "BKE_context.h"
#include "BKE_global.h"
#include "BKE_paint.h"
#include "BKE_pbvh.h"

#include "DRW_pbvh.hh"

#include "atomic_ops.h"
#include "bmesh.h"
#include "bmesh_log.h"
#include "dyntopo_intern.hh"
#include "pbvh_intern.hh"

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <cstdarg>

using blender::float2;
using blender::float3;
using blender::IndexRange;
using blender::Map;
using blender::Set;
using blender::Vector;

template<typename T> T *c_array_from_vector(Vector<T> &array)
{
  T *ret = MEM_cnew_array<T>(array.size(), __func__);
  memcpy(static_cast<void *>(ret), static_cast<void *>(array.data()), sizeof(T) * array.size());
  return ret;
}

static void _debugprint(const char *fmt, ...)
{
  va_list args;
  va_start(args, fmt);
  vprintf(fmt, args);
  va_end(args);
}

#ifdef PBVH_BMESH_DEBUG
void pbvh_bmesh_check_nodes_simple(PBVH *pbvh)
{
  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;
    BMFace *f;

    if (!(node->flag & PBVH_Leaf)) {
      continue;
    }

    TGSET_ITER (f, node->bm_faces) {
      if (!f || f->head.htype != BM_FACE) {
        _debugprint("Corrupted (freed?) face in node->bm_faces\n");
        continue;
      }

      if (BM_ELEM_CD_GET_INT(f, pbvh->cd_face_node_offset) != i) {
        _debugprint("Face in more then one node\n");
      }
    }
    TGSET_ITER_END;
  }
}

ATTR_NO_OPT void pbvh_bmesh_check_nodes(PBVH *pbvh)
{
  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (node->flag & PBVH_Leaf) {
      pbvh_bmesh_check_other_verts(node);
    }
  }

  BMVert *v;
  BMIter iter;

  BM_ITER_MESH (v, &iter, pbvh->header.bm, BM_VERTS_OF_MESH) {
    int ni = BM_ELEM_CD_GET_INT(v, pbvh->cd_vert_node_offset);

    if (ni >= 0 && (!v->e || !v->e->l)) {
      _debugprint("wire vert had node reference: %p (type %d)\n", v, v->head.htype);
      // BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, DYNTOPO_NODE_NONE);
    }

    if (ni < -1 || ni >= pbvh->totnode) {
      _debugprint("vert node ref was invalid: %p (type %d)\n", v, v->head.htype);
      continue;
    }

    if (ni == -1) {
      continue;
    }

    PBVHNode *node = pbvh->nodes + ni;
    if (!(node->flag & PBVH_Leaf) || !node->bm_unique_verts) {
      _debugprint("vert node ref was in non leaf node");
      continue;
    }

    if (!BLI_table_gset_haskey(node->bm_unique_verts, v)) {
      _debugprint("vert not in node->bm_unique_verts\n");
    }

    if (BLI_table_gset_haskey(node->bm_other_verts, v)) {
      _debugprint("vert in node->bm_other_verts");
    }

    MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, v);
    BKE_pbvh_bmesh_check_valence(pbvh, (PBVHVertRef){.i = (intptr_t)v});

    if (BM_vert_edge_count(v) != mv->valence) {
      _debugprint("cached vertex valence mismatch; old: %d, should be: %d\n",
                  mv->valence,
                  BM_vert_edge_count(v));
    }
  }

  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;
    BMVert *v;
    BMFace *f;

    // delete nodes should
    if (node->flag & PBVH_Delete) {
      _debugprint("orphaned delete node\n");
    }

    if (!(node->flag & PBVH_Leaf)) {
      if (node->bm_unique_verts || node->bm_other_verts || node->bm_faces) {
        _debugprint("dangling leaf pointers in non-leaf node\n");
      }

      continue;
    }

    TGSET_ITER (v, node->bm_unique_verts) {
      int ni = BM_ELEM_CD_GET_INT(v, pbvh->cd_vert_node_offset);

      if (ni != i) {
        if (ni >= 0 && ni < pbvh->totnode) {
          PBVHNode *node2 = pbvh->nodes + ni;
          _debugprint("v node offset is wrong, %d\n",
                      !node2->bm_unique_verts ? 0 :
                                                BLI_table_gset_haskey(node2->bm_unique_verts, v));
        }
        else {
          _debugprint("v node offset is wrong\n");
        }
      }

      if (!v || v->head.htype != BM_VERT) {
        _debugprint("corruption in pbvh! bm_unique_verts\n");
      }
      else if (BLI_table_gset_haskey(node->bm_other_verts, v)) {
        _debugprint("v in both unique and other verts\n");
      }
    }
    TGSET_ITER_END;

    TGSET_ITER (f, node->bm_faces) {
      if (!f || f->head.htype != BM_FACE) {
        _debugprint("corruption in pbvh! bm_faces\n");
        continue;
      }

      int ni = BM_ELEM_CD_GET_INT(f, pbvh->cd_face_node_offset);
      if (pbvh->nodes + ni != node) {
        _debugprint("face in multiple nodes!\n");
      }
    }
    TGSET_ITER_END;

    TGSET_ITER (v, node->bm_other_verts) {
      if (!v || v->head.htype != BM_VERT) {
        _debugprint("corruption in pbvh! bm_other_verts\n");
      }
      else if (BLI_table_gset_haskey(node->bm_unique_verts, v)) {
        _debugprint("v in both unique and other verts\n");
      }
    }
    TGSET_ITER_END;
  }
}

ATTR_NO_OPT void pbvh_bmesh_pbvh_bmesh_check_nodes(PBVH *pbvh)
{
  pbvh_bmesh_check_nodes(pbvh);
}
#endif

/** \} */

/****************************** Vertex/Face APIs ******************************/
namespace blender::dyntopo {

void pbvh_kill_vert(PBVH *pbvh, BMVert *v, bool log_vert, bool log_edges)
{
  BMEdge *e = v->e;
  bm_logstack_push();

  if (log_edges) {
    if (e) {
      do {
        BM_log_edge_removed(pbvh->header.bm, pbvh->bm_log, e);
      } while ((e = BM_DISK_EDGE_NEXT(e, v)) != v->e);
    }
  }

  if (log_vert) {
    BM_log_vert_removed(pbvh->header.bm, pbvh->bm_log, v);
  }

  BM_idmap_release(pbvh->bm_idmap, (BMElem *)v, true);
  BM_vert_kill(pbvh->header.bm, v);
  bm_logstack_pop();
}

static BMVert *pbvh_bmesh_vert_create(PBVH *pbvh,
                                      int node_index,
                                      const float co[3],
                                      const float no[3],
                                      BMVert *v_example,
                                      const int cd_vert_mask_offset)
{
  PBVHNode *node = &pbvh->nodes[node_index];

  BLI_assert((pbvh->totnode == 1 || node_index) && node_index <= pbvh->totnode);

  /* avoid initializing customdata because its quite involved */
  BMVert *v = BM_vert_create(pbvh->header.bm, co, nullptr, BM_CREATE_NOP);
  MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, v);

  pbvh_boundary_update_bmesh(pbvh, v);
  MV_ADD_FLAG(mv, SCULPTVERT_NEED_DISK_SORT | SCULPTVERT_NEED_VALENCE);

  if (v_example) {
    v->head.hflag = v_example->head.hflag;

    CustomData_bmesh_copy_data(
        &pbvh->header.bm->vdata, &pbvh->header.bm->vdata, v_example->head.data, &v->head.data);

    /* This value is logged below */
    copy_v3_v3(v->no, no);

    // keep MSculptVert copied from v_example as-is
  }
  else {
    MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, v);

    copy_v3_v3(mv->origco, co);
    copy_v3_v3(mv->origno, no);
    mv->origmask = 0.0f;

    /* This value is logged below */
    copy_v3_v3(v->no, no);
  }

  BLI_table_gset_insert(node->bm_unique_verts, v);
  BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, node_index);

  node->flag |= PBVH_UpdateDrawBuffers | PBVH_UpdateBB | PBVH_UpdateTris | PBVH_UpdateOtherVerts;

  /* Log the new vertex */
  BM_log_vert_added(pbvh->header.bm, pbvh->bm_log, v);
  v->head.index = pbvh->header.bm->totvert;  // set provisional index

  return v;
}

static BMFace *bmesh_face_create_edge_log(PBVH *pbvh,
                                          BMVert *v_tri[3],
                                          BMEdge *e_tri[3],
                                          const BMFace *f_example)
{
  BMFace *f;

  if (!e_tri) {
    BMEdge *e_tri2[3];

    for (int i = 0; i < 3; i++) {
      BMVert *v1 = v_tri[i];
      BMVert *v2 = v_tri[(i + 1) % 3];

      BMEdge *e = BM_edge_exists(v1, v2);

      if (!e) {
        e = BM_edge_create(pbvh->header.bm, v1, v2, nullptr, BM_CREATE_NOP);
        BM_log_edge_added(pbvh->header.bm, pbvh->bm_log, e);
      }

      e_tri2[i] = e;
    }

    // f = BM_face_create_verts(pbvh->header.bm, v_tri, 3, f_example, BM_CREATE_NOP, true);
    f = BM_face_create(pbvh->header.bm, v_tri, e_tri2, 3, f_example, BM_CREATE_NOP);
  }
  else {
    f = BM_face_create(pbvh->header.bm, v_tri, e_tri, 3, f_example, BM_CREATE_NOP);
  }

  if (f_example) {
    f->head.hflag = f_example->head.hflag;
  }

  return f;
}

/**
 * \note Callers are responsible for checking if the face exists before adding.
 */
BMFace *pbvh_bmesh_face_create(PBVH *pbvh,
                               int node_index,
                               BMVert *v_tri[3],
                               BMEdge *e_tri[3],
                               const BMFace *f_example,
                               bool ensure_verts,
                               bool log_face)
{
  PBVHNode *node = &pbvh->nodes[node_index];

  /* ensure we never add existing face */
  BLI_assert(!BM_face_exists(v_tri, 3));

  BMFace *f = bmesh_face_create_edge_log(pbvh, v_tri, e_tri, f_example);

  BLI_table_gset_insert(node->bm_faces, f);
  BM_ELEM_CD_SET_INT(f, pbvh->cd_face_node_offset, node_index);

  /* mark node for update */
  node->flag |= PBVH_UpdateDrawBuffers | PBVH_UpdateNormals | PBVH_UpdateTris |
                PBVH_UpdateOtherVerts;
  node->flag &= ~PBVH_FullyHidden;

  /* Log the new face */
  if (log_face) {
    BM_log_face_added(pbvh->header.bm, pbvh->bm_log, f);
  }

  int cd_vert_node = pbvh->cd_vert_node_offset;

  if (ensure_verts) {
    BMLoop *l = f->l_first;
    do {
      int ni = BM_ELEM_CD_GET_INT(l->v, cd_vert_node);

      if (ni == DYNTOPO_NODE_NONE) {
        BLI_table_gset_add(node->bm_unique_verts, l->v);
        BM_ELEM_CD_SET_INT(l->v, cd_vert_node, node_index);

        node->flag |= PBVH_UpdateDrawBuffers | PBVH_UpdateBB | PBVH_UpdateTris |
                      PBVH_UpdateOtherVerts;
      }

      pbvh_boundary_update_bmesh(pbvh, l->v);

      MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, l->v);
      MV_ADD_FLAG(mv, SCULPTVERT_NEED_DISK_SORT | SCULPTVERT_NEED_VALENCE);

      l = l->next;
    } while (l != f->l_first);
  }
  else {
    BMLoop *l = f->l_first;
    do {
      pbvh_boundary_update_bmesh(pbvh, l->v);

      MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, l->v);
      MV_ADD_FLAG(mv, SCULPTVERT_NEED_DISK_SORT | SCULPTVERT_NEED_VALENCE);
    } while ((l = l->next) != f->l_first);
  }

  return f;
}

BMVert *BKE_pbvh_vert_create_bmesh(
    PBVH *pbvh, float co[3], float no[3], PBVHNode *node, BMVert *v_example)
{
  if (!node) {
    for (int i = 0; i < pbvh->totnode; i++) {
      PBVHNode *node2 = pbvh->nodes + i;

      if (!(node2->flag & PBVH_Leaf)) {
        continue;
      }

      /* Ensure we have at least some node somewhere picked. */
      node = node2;

      bool ok = true;

      for (int j = 0; j < 3; j++) {
        if (co[j] < node2->vb.bmin[j] || co[j] >= node2->vb.bmax[j]) {
          continue;
        }
      }

      if (ok) {
        break;
      }
    }
  }

  BMVert *v;

  if (!node) {
    printf("possible pbvh error\n");
    v = BM_vert_create(pbvh->header.bm, co, v_example, BM_CREATE_NOP);
    BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, DYNTOPO_NODE_NONE);

    pbvh_boundary_update_bmesh(pbvh, v);
    MSculptVert *mv = (MSculptVert *)BM_ELEM_CD_GET_VOID_P(v, pbvh->cd_sculpt_vert);
    MV_ADD_FLAG(mv, SCULPTVERT_NEED_VALENCE | SCULPTVERT_NEED_DISK_SORT);

    copy_v3_v3(mv->origco, co);

    return v;
  }

  return pbvh_bmesh_vert_create(
      pbvh, node - pbvh->nodes, co, no, v_example, pbvh->cd_vert_mask_offset);
}

PBVHNode *BKE_pbvh_node_from_face_bmesh(PBVH *pbvh, BMFace *f)
{
  return pbvh->nodes + BM_ELEM_CD_GET_INT(f, pbvh->cd_face_node_offset);
}

BMFace *BKE_pbvh_face_create_bmesh(PBVH *pbvh,
                                   BMVert *v_tri[3],
                                   BMEdge *e_tri[3],
                                   const BMFace *f_example)
{
  int ni = DYNTOPO_NODE_NONE;

  for (int i = 0; i < 3; i++) {
    BMVert *v = v_tri[i];
    BMLoop *l;
    BMIter iter;

    BM_ITER_ELEM (l, &iter, v, BM_LOOPS_OF_VERT) {
      int ni2 = BM_ELEM_CD_GET_INT(l->f, pbvh->cd_face_node_offset);
      if (ni2 != DYNTOPO_NODE_NONE) {
        ni = ni2;
        break;
      }
    }
  }

  if (ni == DYNTOPO_NODE_NONE) {
    BMFace *f;

    /* No existing nodes? Find one. */
    for (int i = 0; i < pbvh->totnode; i++) {
      PBVHNode *node = pbvh->nodes + i;

      if (!(node->flag & PBVH_Leaf)) {
        continue;
      }

      for (int j = 0; j < 3; j++) {
        BMVert *v = v_tri[j];

        bool ok = true;

        for (int k = 0; k < 3; k++) {
          if (v->co[k] < node->vb.bmin[k] || v->co[k] >= node->vb.bmax[k]) {
            ok = false;
          }
        }

        if (ok &&
            (ni == DYNTOPO_NODE_NONE || BLI_table_gset_len(node->bm_faces) < pbvh->leaf_limit)) {
          ni = i;
          break;
        }
      }

      if (ni != DYNTOPO_NODE_NONE) {
        break;
      }
    }

    if (ni == DYNTOPO_NODE_NONE) {
      /* Empty pbvh? */
      f = bmesh_face_create_edge_log(pbvh, v_tri, e_tri, f_example);

      BM_ELEM_CD_SET_INT(f, pbvh->cd_face_node_offset, DYNTOPO_NODE_NONE);

      return f;
    }
  }

  return pbvh_bmesh_face_create(pbvh, ni, v_tri, e_tri, f_example, true, true);
}

#define pbvh_bmesh_node_vert_use_count_is_equal(pbvh, node, v, n) \
  (pbvh_bmesh_node_vert_use_count_at_most(pbvh, node, v, (n) + 1) == n)

static int pbvh_bmesh_node_vert_use_count_at_most(PBVH *pbvh,
                                                  PBVHNode *node,
                                                  BMVert *v,
                                                  const int count_max)
{
  int count = 0;
  BMFace *f;

  BM_FACES_OF_VERT_ITER_BEGIN (f, v) {
    PBVHNode *f_node = pbvh_bmesh_node_from_face(pbvh, f);
    if (f_node == node) {
      count++;
      if (count == count_max) {
        return count;
      }
    }
  }
  BM_FACES_OF_VERT_ITER_END;

  return count;
}

/* Return a node that uses vertex 'v' other than its current owner */
static PBVHNode *pbvh_bmesh_vert_other_node_find(PBVH *pbvh, BMVert *v)
{
  PBVHNode *current_node = pbvh_bmesh_node_from_vert(pbvh, v);
  BMFace *f;

  BM_FACES_OF_VERT_ITER_BEGIN (f, v) {
    PBVHNode *f_node = pbvh_bmesh_node_from_face(pbvh, f);

    if (f_node != current_node) {
      return f_node;
    }
  }
  BM_FACES_OF_VERT_ITER_END;

  return nullptr;
}

static void pbvh_bmesh_vert_ownership_transfer(PBVH *pbvh, PBVHNode *new_owner, BMVert *v)
{
  PBVHNode *current_owner = pbvh_bmesh_node_from_vert(pbvh, v);
  /* Mark node for update. */

  if (current_owner) {
    current_owner->flag |= PBVH_UpdateDrawBuffers | PBVH_UpdateBB;

    BLI_assert(current_owner != new_owner);

    /* Remove current ownership. */
    BLI_table_gset_remove(current_owner->bm_unique_verts, v, nullptr);
  }

  /* Set new ownership. */
  BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, new_owner - pbvh->nodes);
  BLI_table_gset_insert(new_owner->bm_unique_verts, v);

  /* Mark node for update. */
  new_owner->flag |= PBVH_UpdateDrawBuffers | PBVH_UpdateBB | PBVH_UpdateOtherVerts;
}

void pbvh_bmesh_vert_remove(PBVH *pbvh, BMVert *v)
{
  /* never match for first time */
  int f_node_index_prev = DYNTOPO_NODE_NONE;
  const int updateflag = PBVH_UpdateDrawBuffers | PBVH_UpdateBB | PBVH_UpdateTris |
                         PBVH_UpdateNormals | PBVH_UpdateOtherVerts;

  PBVHNode *v_node = pbvh_bmesh_node_from_vert(pbvh, v);

  if (v_node) {
    BLI_table_gset_remove(v_node->bm_unique_verts, v, nullptr);
    v_node->flag |= (PBVHNodeFlags)updateflag;
  }

  BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, DYNTOPO_NODE_NONE);

  /* Have to check each neighboring face's node */
  BMFace *f;
  BM_FACES_OF_VERT_ITER_BEGIN (f, v) {
    const int f_node_index = pbvh_bmesh_node_index_from_face(pbvh, f);

    if (f_node_index == DYNTOPO_NODE_NONE) {
      continue;
    }

    /* faces often share the same node,
     * quick check to avoid redundant #BLI_table_gset_remove calls */
    if (f_node_index_prev != f_node_index) {
      f_node_index_prev = f_node_index;

      PBVHNode *f_node = &pbvh->nodes[f_node_index];
      f_node->flag |= (PBVHNodeFlags)updateflag;  // flag update of bm_other_verts

      BLI_assert(!BLI_table_gset_haskey(f_node->bm_unique_verts, v));
    }
  }
  BM_FACES_OF_VERT_ITER_END;
}

void pbvh_bmesh_face_remove(
    PBVH *pbvh, BMFace *f, bool log_face, bool check_verts, bool ensure_ownership_transfer)
{
  PBVHNode *f_node = pbvh_bmesh_node_from_face(pbvh, f);

  if (!f_node || !(f_node->flag & PBVH_Leaf)) {
    printf("pbvh corruption\n");
    fflush(stdout);
    return;
  }

  bm_logstack_push();

  /* Check if any of this face's vertices need to be removed
   * from the node */
  if (check_verts) {
    BMLoop *l_first = BM_FACE_FIRST_LOOP(f);
    BMLoop *l_iter = l_first;
    do {
      BMVert *v = l_iter->v;
      if (pbvh_bmesh_node_vert_use_count_is_equal(pbvh, f_node, v, 1)) {
        if (BM_ELEM_CD_GET_INT(v, pbvh->cd_vert_node_offset) == f_node - pbvh->nodes) {
          // if (BLI_table_gset_haskey(f_node->bm_unique_verts, v)) {
          /* Find a different node that uses 'v' */
          PBVHNode *new_node;

          new_node = pbvh_bmesh_vert_other_node_find(pbvh, v);
          // BLI_assert(new_node || BM_vert_face_count_is_equal(v, 1));

          if (new_node) {
            pbvh_bmesh_vert_ownership_transfer(pbvh, new_node, v);
          }
          else if (ensure_ownership_transfer && !BM_vert_face_count_is_equal(v, 1)) {
            pbvh_bmesh_vert_remove(pbvh, v);

            f_node->flag |= PBVH_RebuildNodeVerts | PBVH_UpdateOtherVerts;
            // printf("failed to find new_node\n");
          }
        }
      }
    } while ((l_iter = l_iter->next) != l_first);
  }

  /* Remove face from node and top level */
  BLI_table_gset_remove(f_node->bm_faces, f, nullptr);
  BM_ELEM_CD_SET_INT(f, pbvh->cd_face_node_offset, DYNTOPO_NODE_NONE);

  /* Log removed face */
  if (log_face) {
    BM_log_face_removed(pbvh->header.bm, pbvh->bm_log, f);
  }

  /* mark node for update */
  f_node->flag |= PBVH_UpdateDrawBuffers | PBVH_UpdateNormals | PBVH_UpdateTris |
                  PBVH_UpdateOtherVerts;

  bm_logstack_pop();
}
}  // namespace blender::dyntopo

/****************************** Building ******************************/

/* Update node data after splitting */
static void pbvh_bmesh_node_finalize(PBVH *pbvh,
                                     const int node_index,
                                     const int cd_vert_node_offset,
                                     const int cd_face_node_offset,
                                     bool add_orco)
{
  PBVHNode *n = &pbvh->nodes[node_index];
  bool has_visible = false;

  n->draw_batches = nullptr;

  /* Create vert hash sets */
  if (!n->bm_unique_verts) {
    n->bm_unique_verts = BLI_table_gset_new("bm_unique_verts");
  }
  n->bm_other_verts = BLI_table_gset_new("bm_other_verts");

  BB_reset(&n->vb);
  BB_reset(&n->orig_vb);
  BMFace *f;

  TGSET_ITER (f, n->bm_faces) {
    /* Update ownership of faces */
    BM_ELEM_CD_SET_INT(f, cd_face_node_offset, node_index);

    /* Update vertices */
    BMLoop *l_first = BM_FACE_FIRST_LOOP(f);
    BMLoop *l_iter = l_first;

    do {
      BMVert *v = l_iter->v;
      MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, v);

      int *flags = BM_ELEM_CD_PTR<int *>(v, pbvh->cd_boundary_flag);
      *flags |= SCULPT_BOUNDARY_NEEDS_UPDATE;

      if (!BLI_table_gset_haskey(n->bm_unique_verts, v)) {
        if (BM_ELEM_CD_GET_INT(v, cd_vert_node_offset) != DYNTOPO_NODE_NONE) {
          BLI_table_gset_add(n->bm_other_verts, v);
        }
        else {
          BLI_table_gset_insert(n->bm_unique_verts, v);
          BM_ELEM_CD_SET_INT(v, cd_vert_node_offset, node_index);
        }
      }
      /* Update node bounding box */
      BB_expand(&n->vb, v->co);
      BB_expand(&n->orig_vb, mv->origco);
    } while ((l_iter = l_iter->next) != l_first);

    if (!BM_elem_flag_test(f, BM_ELEM_HIDDEN)) {
      has_visible = true;
    }
  }
  TGSET_ITER_END

  BLI_assert(n->vb.bmin[0] <= n->vb.bmax[0] && n->vb.bmin[1] <= n->vb.bmax[1] &&
             n->vb.bmin[2] <= n->vb.bmax[2]);

  /* Build GPU buffers for new node and update vertex normals */
  BKE_pbvh_node_mark_rebuild_draw(n);

  BKE_pbvh_node_fully_hidden_set(n, !has_visible);
  n->flag |= PBVH_UpdateNormals | PBVH_UpdateCurvatureDir | PBVH_UpdateTris;
  n->flag |= PBVH_UpdateBB | PBVH_UpdateOriginalBB;

  if (add_orco) {
    BKE_pbvh_bmesh_check_tris(pbvh, n);
  }
}

static void pbvh_print_mem_size(PBVH *pbvh)
{
  BMesh *bm = pbvh->header.bm;
  CustomData *cdatas[4] = {&bm->vdata, &bm->edata, &bm->ldata, &bm->pdata};

  int tots[4] = {bm->totvert, bm->totedge, bm->totloop, bm->totface};
  int sizes[4] = {
      (int)sizeof(BMVert), (int)sizeof(BMEdge), (int)sizeof(BMLoop), (int)sizeof(BMFace)};

  float memsize1[4] = {0.0f, 0.0f, 0.0f, 0.0f};
  float memsize2[4] = {0.0f, 0.0f, 0.0f, 0.0f};
  float tot = 0.0f;

  for (int i = 0; i < 4; i++) {
    CustomData *cdata = cdatas[i];

    memsize1[i] = (float)(sizes[i] * tots[i]) / 1024.0f / 1024.0f;
    memsize2[i] = (float)(cdata->totsize * tots[i]) / 1024.0f / 1024.0f;

    tot += memsize1[i] + memsize2[i];
  }

  printf("base sizes:\n");
  printf("  v: %.2fmb e: %.2fmb l: %.2fmb f: %.2fmb\n",
         memsize1[0],
         memsize1[1],
         memsize1[2],
         memsize1[3]);

  printf("custom attribute sizes:\n");
  printf("  v: %.2fmb e: %.2fmb l: %.2fmb f: %.2fmb\n",
         memsize2[0],
         memsize2[1],
         memsize2[2],
         memsize2[3]);

  int ptrsize = (int)sizeof(void *);

  float memsize3[3] = {(float)(ptrsize * bm->idmap.map_size) / 1024.0f / 1024.0f,
                       (float)(ptrsize * bm->idmap.freelist_len) / 1024.0f / 1024.0f,
                       (float)(4 * bm->idmap.free_ids_size) / 1024.0f / 1024.0f};

  printf("idmap sizes:\n  map_size: %.2fmb freelist_len: %.2fmb free_ids_size: %.2fmb\n",
         memsize3[0],
         memsize3[1],
         memsize3[2]);

  tot += memsize3[0] + memsize3[1] + memsize3[2];

  printf("total: %.2f\n", tot);

#ifdef WITH_BM_ID_FREELIST
  if (bm->idmap.free_idx_map) {
    printf("freelist length: %d\n", bm->idmap.freelist_len);
    /* printf("free_idx_map: nentries %d, size %d: nfreecells: %d\n",
           bm->idmap.free_idx_map->nentries,
           bm->idmap.free_idx_map->nbuckets,
           bm->idmap.free_idx_map->nfreecells);*/
  }
#endif
}

/* Recursively split the node if it exceeds the leaf_limit */
static void pbvh_bmesh_node_split(
    PBVH *pbvh, const BBC *bbc_array, int node_index, bool add_orco, int depth)
{
  const int cd_vert_node_offset = pbvh->cd_vert_node_offset;
  const int cd_face_node_offset = pbvh->cd_face_node_offset;
  PBVHNode *n = &pbvh->nodes[node_index];

#ifdef PROXY_ADVANCED
  BKE_pbvh_free_proxyarray(pbvh, n);
#endif

  if (n->depth >= PBVH_STACK_FIXED_DEPTH || BLI_table_gset_len(n->bm_faces) <= pbvh->leaf_limit) {
    /* Node limit not exceeded */
    pbvh_bmesh_node_finalize(pbvh, node_index, cd_vert_node_offset, cd_face_node_offset, add_orco);
    return;
  }

  /* Calculate bounding box around primitive centroids */
  BB cb;
  BB_reset(&cb);
  BMFace *f;

  TGSET_ITER (f, n->bm_faces) {
    const BBC *bbc = &bbc_array[BM_elem_index_get(f)];

    BB_expand(&cb, bbc->bcentroid);
  }
  TGSET_ITER_END

  /* Find widest axis and its midpoint */
  const int axis = BB_widest_axis(&cb);
  const float mid = (cb.bmax[axis] + cb.bmin[axis]) * 0.5f;

  if (isnan(mid)) {
    printf("NAN ERROR! %s\n", __func__);
  }

  /* Add two new child nodes */
  const int children = pbvh->totnode;
  n->children_offset = children;
  pbvh_grow_nodes(pbvh, pbvh->totnode + 2);

  /* Array reallocated, update current node pointer */
  n = &pbvh->nodes[node_index];

  /* Initialize children */
  PBVHNode *c1 = &pbvh->nodes[children], *c2 = &pbvh->nodes[children + 1];

  c1->draw_batches = c2->draw_batches = nullptr;
  c1->depth = c2->depth = n->depth + 1;

  c1->flag |= PBVH_Leaf;
  c2->flag |= PBVH_Leaf;

  c1->bm_faces = BLI_table_gset_new_ex("bm_faces", BLI_table_gset_len(n->bm_faces) / 2);
  c2->bm_faces = BLI_table_gset_new_ex("bm_faces", BLI_table_gset_len(n->bm_faces) / 2);

  c1->bm_unique_verts = BLI_table_gset_new("bm_unique_verts");
  c2->bm_unique_verts = BLI_table_gset_new("bm_unique_verts");

  c1->bm_other_verts = c2->bm_other_verts = nullptr;

  /* Partition the parent node's faces between the two children */
  TGSET_ITER (f, n->bm_faces) {
    const BBC *bbc = &bbc_array[BM_elem_index_get(f)];

    if (bbc->bcentroid[axis] < mid) {
      BLI_table_gset_insert(c1->bm_faces, f);
    }
    else {
      BLI_table_gset_insert(c2->bm_faces, f);
    }
  }
  TGSET_ITER_END
#if 0
    /* Enforce at least one primitive in each node */
    TableGSet *empty = nullptr,*other;
  if (BLI_table_gset_len(c1->bm_faces) == 0) {
    empty = c1->bm_faces;
    other = c2->bm_faces;
  } else if (BLI_table_gset_len(c2->bm_faces) == 0) {
    empty = c2->bm_faces;
    other = c1->bm_faces;
  }

  if (empty) {
    void *key;
    TGSET_ITER (key,other) {
      BLI_table_gset_insert(empty,key);
      BLI_table_gset_remove(other,key,nullptr);
      break;
    }
    TGSET_ITER_END
  }
#endif
  /* Clear this node */

  BMVert *v;

  /* Assign verts to c1 and c2.  Note that the previous
     method of simply marking them as untaken and rebuilding
     unique verts later doesn't work, as it assumes that dyntopo
     never assigns verts to nodes that don't contain their
     faces.*/
  if (n->bm_unique_verts) {
    TGSET_ITER (v, n->bm_unique_verts) {
      if (v->co[axis] < mid) {
        BM_ELEM_CD_SET_INT(v, cd_vert_node_offset, (c1 - pbvh->nodes));
        BLI_table_gset_add(c1->bm_unique_verts, v);
      }
      else {
        BM_ELEM_CD_SET_INT(v, cd_vert_node_offset, (c2 - pbvh->nodes));
        BLI_table_gset_add(c2->bm_unique_verts, v);
      }
    }
    TGSET_ITER_END

    BLI_table_gset_free(n->bm_unique_verts, nullptr);
  }

  if (n->bm_faces) {
    /* Unclaim faces */
    TGSET_ITER (f, n->bm_faces) {
      BM_ELEM_CD_SET_INT(f, cd_face_node_offset, DYNTOPO_NODE_NONE);
    }
    TGSET_ITER_END

    BLI_table_gset_free(n->bm_faces, nullptr);
  }

  if (n->bm_other_verts) {
    BLI_table_gset_free(n->bm_other_verts, nullptr);
  }

  if (n->layer_disp) {
    MEM_freeN(n->layer_disp);
  }

  if (n->tribuf || n->tri_buffers) {
    BKE_pbvh_bmesh_free_tris(pbvh, n);
  }

  n->bm_faces = nullptr;
  n->bm_unique_verts = nullptr;
  n->bm_other_verts = nullptr;
  n->layer_disp = nullptr;

  if (n->draw_batches) {
    DRW_pbvh_node_free(n->draw_batches);
    n->draw_batches = nullptr;
  }
  n->flag &= ~PBVH_Leaf;

  /* Recurse */
  pbvh_bmesh_node_split(pbvh, bbc_array, children, add_orco, depth + 1);
  pbvh_bmesh_node_split(pbvh, bbc_array, children + 1, add_orco, depth + 1);

  /* Array maybe reallocated, update current node pointer */
  n = &pbvh->nodes[node_index];

  /* Update bounding box */
  BB_reset(&n->vb);
  BB_expand_with_bb(&n->vb, &pbvh->nodes[n->children_offset].vb);
  BB_expand_with_bb(&n->vb, &pbvh->nodes[n->children_offset + 1].vb);
  n->orig_vb = n->vb;
}

/* Recursively split the node if it exceeds the leaf_limit */
bool pbvh_bmesh_node_limit_ensure(PBVH *pbvh, int node_index)
{
  TableGSet *bm_faces = pbvh->nodes[node_index].bm_faces;
  const int bm_faces_size = BLI_table_gset_len(bm_faces);

  // pbvh_bmesh_check_nodes(pbvh);

  if (bm_faces_size <= pbvh->leaf_limit ||
      pbvh->nodes[node_index].depth >= PBVH_STACK_FIXED_DEPTH) {
    /* Node limit not exceeded */
    return false;
  }

  /* Trigger draw manager cache invalidation. */
  pbvh->draw_cache_invalid = true;

  /* For each BMFace, store the AABB and AABB centroid */
  BBC *bbc_array = MEM_cnew_array<BBC>(bm_faces_size, "BBC");

  BMFace *f;

  int i;

  /*
  TGSET_ITER_INDEX(f, bm_faces, i)
  {
  }
  TGSET_ITER_INDEX_END
  printf("size: %d %d\n", i + 1, bm_faces_size);
  */

  TGSET_ITER_INDEX(f, bm_faces, i)
  {
    BBC *bbc = &bbc_array[i];

    BB_reset((BB *)bbc);
    BMLoop *l_first = BM_FACE_FIRST_LOOP(f);
    BMLoop *l_iter = l_first;
    do {
      BB_expand((BB *)bbc, l_iter->v->co);
    } while ((l_iter = l_iter->next) != l_first);
    BBC_update_centroid(bbc);

    /* so we can do direct lookups on 'bbc_array' */
    BM_elem_index_set(f, i); /* set_dirty! */
  }
  TGSET_ITER_INDEX_END

  /* Likely this is already dirty. */
  pbvh->header.bm->elem_index_dirty |= BM_FACE;

  pbvh_bmesh_node_split(pbvh, bbc_array, node_index, false, 0);

  MEM_freeN(bbc_array);

  // pbvh_bmesh_check_nodes(pbvh);

  return true;
}

/**********************************************************************/

static bool point_in_node(const PBVHNode *node, const float co[3])
{
  return co[0] >= node->vb.bmin[0] && co[0] <= node->vb.bmax[0] && co[1] >= node->vb.bmin[1] &&
         co[1] <= node->vb.bmax[1] && co[2] >= node->vb.bmin[2] && co[2] <= node->vb.bmax[2];
}

void bke_pbvh_insert_face_finalize(PBVH *pbvh, BMFace *f, const int ni)
{
  PBVHNode *node = pbvh->nodes + ni;
  BM_ELEM_CD_SET_INT(f, pbvh->cd_face_node_offset, ni);

  if (!(node->flag & PBVH_Leaf)) {
    printf("%s: major pbvh corruption error\n", __func__);
    return;
  }

  BLI_table_gset_add(node->bm_faces, f);

  PBVHNodeFlags updateflag = PBVH_UpdateTris | PBVH_UpdateBB | PBVH_UpdateDrawBuffers |
                             PBVH_UpdateCurvatureDir | PBVH_UpdateOtherVerts;
  updateflag |= PBVH_UpdateColor | PBVH_UpdateMask | PBVH_UpdateNormals | PBVH_UpdateOriginalBB;
  updateflag |= PBVH_UpdateVisibility | PBVH_UpdateRedraw | PBVH_RebuildDrawBuffers;

  node->flag |= updateflag;

  // ensure verts are in pbvh
  BMLoop *l = f->l_first;
  do {
    const int ni2 = BM_ELEM_CD_GET_INT(l->v, pbvh->cd_vert_node_offset);
    MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, l->v);

    BB_expand(&node->vb, l->v->co);
    BB_expand(&node->orig_vb, mv->origco);

    if (ni2 == DYNTOPO_NODE_NONE) {
      BM_ELEM_CD_SET_INT(l->v, pbvh->cd_vert_node_offset, ni);
      BLI_table_gset_add(node->bm_unique_verts, l->v);
    }
    else {
      PBVHNode *node2 = pbvh->nodes + ni2;

      if (ni != ni2) {
        BLI_table_gset_add(node->bm_other_verts, l->v);
      }

      node2->flag |= updateflag;

      BB_expand(&node2->vb, l->v->co);
      BB_expand(&node2->orig_vb, mv->origco);
    }
    l = l->next;
  } while (l != f->l_first);
}

void bke_pbvh_insert_face(PBVH *pbvh, struct BMFace *f)
{
  int i = 0;
  bool ok = false;
  int ni = -1;

  while (i < pbvh->totnode) {
    PBVHNode *node = pbvh->nodes + i;
    bool ok2 = false;

    if (node->flag & PBVH_Leaf) {
      ok = true;
      ni = i;
      break;
    }

    if (node->children_offset == 0) {
      continue;
    }

    for (int j = 0; j < 2; j++) {
      int ni2 = node->children_offset + j;
      if (ni2 == 0) {
        continue;
      }

      PBVHNode *node2 = pbvh->nodes + ni2;
      BMLoop *l = f->l_first;

      do {
        if (point_in_node(node2, l->v->co)) {
          i = ni2;
          ok2 = true;
          break;
        }

        l = l->next;
      } while (l != f->l_first);

      if (ok2) {
        break;
      }
    }

    if (!ok2) {
      break;
    }
  }

  if (!ok) {
    // find closest node
    float co[3];
    int tot = 0;
    BMLoop *l = f->l_first;

    zero_v3(co);

    do {
      add_v3_v3(co, l->v->co);
      l = l->next;
      tot++;
    } while (l != f->l_first);

    mul_v3_fl(co, 1.0f / (float)tot);
    float mindis = 1e17;

    for (int i = 0; i < pbvh->totnode; i++) {
      PBVHNode *node = pbvh->nodes + i;

      if (!(node->flag & PBVH_Leaf)) {
        continue;
      }

      float cent[3];
      add_v3_v3v3(cent, node->vb.bmin, node->vb.bmax);
      mul_v3_fl(cent, 0.5f);

      float dis = len_squared_v3v3(co, cent);
      if (dis < mindis) {
        mindis = dis;
        ni = i;
      }
    }
  }

  if (ni < 0 || !(pbvh->nodes[ni].flag & PBVH_Leaf)) {
    fprintf(stderr, "pbvh error! failed to find node to insert face into!\n");
    fflush(stderr);
    return;
  }

  bke_pbvh_insert_face_finalize(pbvh, f, ni);
}

static void pbvh_bmesh_regen_node_verts(PBVH *pbvh, PBVHNode *node)
{
  node->flag &= ~PBVH_RebuildNodeVerts;

  int usize = BLI_table_gset_len(node->bm_unique_verts);
  int osize = BLI_table_gset_len(node->bm_other_verts);

  TableGSet *old_unique_verts = node->bm_unique_verts;

  BLI_table_gset_free(node->bm_other_verts, nullptr);

  BMVert *v;
  TGSET_ITER (v, old_unique_verts) {
    BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, -1);
  }
  TGSET_ITER_END;

  node->bm_unique_verts = BLI_table_gset_new("bm_unique_verts");
  node->bm_other_verts = BLI_table_gset_new("bm_other_verts");

  const int cd_vert_node = pbvh->cd_vert_node_offset;
  const int ni = (int)(node - pbvh->nodes);

  bool update = false;

  BMFace *f;
  TGSET_ITER (f, node->bm_faces) {
    BMLoop *l = f->l_first;
    do {
      int ni2 = BM_ELEM_CD_GET_INT(l->v, cd_vert_node);

      if (ni2 == DYNTOPO_NODE_NONE) {
        BM_ELEM_CD_SET_INT(l->v, cd_vert_node, ni);
        ni2 = ni;
        update = true;
      }

      if (ni2 == ni) {
        BLI_table_gset_add(node->bm_unique_verts, l->v);
      }
      else {
        BLI_table_gset_add(node->bm_other_verts, l->v);
      }
    } while ((l = l->next) != f->l_first);
  }
  TGSET_ITER_END;

  TGSET_ITER (v, old_unique_verts) {
    if (BM_ELEM_CD_GET_INT(v, pbvh->cd_vert_node_offset) == -1) {
      // try to find node to insert into
      BMIter iter2;
      BMFace *f2;
      bool ok = false;

      BM_ITER_ELEM (f2, &iter2, v, BM_FACES_OF_VERT) {
        int ni2 = BM_ELEM_CD_GET_INT(f2, pbvh->cd_face_node_offset);

        if (ni2 >= 0) {
          BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, ni2);
          PBVHNode *node = pbvh->nodes + ni2;

          BLI_table_gset_add(node->bm_unique_verts, v);
          BLI_table_gset_remove(node->bm_other_verts, v, nullptr);

          ok = true;
          break;
        }
      }

      if (!ok) {
        printf("pbvh error: orphaned vert node reference\n");
      }
    }
  }
  TGSET_ITER_END;

  if (usize != BLI_table_gset_len(node->bm_unique_verts)) {
    update = true;
#if 0
    printf("possible pbvh error: bm_unique_verts might have had bad data. old: %d, new: %d\n",
      usize,
      BLI_table_gset_len(node->bm_unique_verts));
#endif
  }

  if (osize != BLI_table_gset_len(node->bm_other_verts)) {
    update = true;
#if 0
    printf("possible pbvh error: bm_other_verts might have had bad data. old: %d, new: %d\n",
      osize,
      BLI_table_gset_len(node->bm_other_verts));
#endif
  }

  if (update) {
    node->flag |= PBVH_UpdateNormals | PBVH_UpdateDrawBuffers | PBVH_RebuildDrawBuffers |
                  PBVH_UpdateBB;
    node->flag |= PBVH_UpdateOriginalBB | PBVH_UpdateRedraw | PBVH_UpdateColor | PBVH_UpdateTris |
                  PBVH_UpdateVisibility;
  }

  BLI_table_gset_free(old_unique_verts, nullptr);
}

void BKE_pbvh_bmesh_mark_node_regen(PBVH *pbvh, PBVHNode *node)
{
  node->flag |= PBVH_RebuildNodeVerts;
}

PBVHNode *BKE_pbvh_get_node_leaf_safe(PBVH *pbvh, int i)
{
  if (i >= 0 && i < pbvh->totnode) {
    PBVHNode *node = pbvh->nodes + i;
    if ((node->flag & PBVH_Leaf) && !(node->flag & PBVH_Delete)) {
      return node;
    }
  }

  return nullptr;
}

void BKE_pbvh_bmesh_regen_node_verts(PBVH *pbvh)
{
  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (!(node->flag & PBVH_Leaf) || !(node->flag & PBVH_RebuildNodeVerts)) {
      continue;
    }

    pbvh_bmesh_regen_node_verts(pbvh, node);
  }
}

/************************* Called from pbvh.c *************************/

static bool pbvh_poly_hidden(PBVH *pbvh, BMFace *f)
{
  return BM_elem_flag_test(f, BM_ELEM_HIDDEN);
}

bool BKE_pbvh_bmesh_check_origdata(PBVH *pbvh, BMVert *v, int stroke_id)
{
  PBVHVertRef vertex = {(intptr_t)v};

  return BKE_pbvh_get_origvert(pbvh, vertex, nullptr, nullptr, nullptr);
}

bool pbvh_bmesh_node_raycast(PBVH *pbvh,
                             PBVHNode *node,
                             const float ray_start[3],
                             const float ray_normal[3],
                             struct IsectRayPrecalc *isect_precalc,
                             int *hit_count,
                             float *depth,
                             float *back_depth,
                             bool use_original,
                             PBVHVertRef *r_active_vertex,
                             PBVHFaceRef *r_active_face,
                             float *r_face_normal,
                             int stroke_id)
{
  bool hit = false;
  float nearest_vertex_co[3] = {0.0f};

  GSetIterator gs_iter;

  BKE_pbvh_bmesh_check_tris(pbvh, node);

  for (int i = 0; i < node->tribuf->tottri; i++) {
    PBVHTri *tri = node->tribuf->tris + i;
    BMVert *verts[3] = {
        (BMVert *)node->tribuf->verts[tri->v[0]].i,
        (BMVert *)node->tribuf->verts[tri->v[1]].i,
        (BMVert *)node->tribuf->verts[tri->v[2]].i,
    };

    float *cos[3];
    float *nos[3];

    if (use_original) {
      BKE_pbvh_bmesh_check_origdata(pbvh, verts[0], stroke_id);
      BKE_pbvh_bmesh_check_origdata(pbvh, verts[1], stroke_id);
      BKE_pbvh_bmesh_check_origdata(pbvh, verts[2], stroke_id);

      MSculptVert *mv1 = BM_ELEM_CD_PTR<MSculptVert *>(verts[0], pbvh->cd_sculpt_vert);
      MSculptVert *mv2 = BM_ELEM_CD_PTR<MSculptVert *>(verts[1], pbvh->cd_sculpt_vert);
      MSculptVert *mv3 = BM_ELEM_CD_PTR<MSculptVert *>(verts[2], pbvh->cd_sculpt_vert);

      cos[0] = mv1->origco;
      cos[1] = mv2->origco;
      cos[2] = mv3->origco;

      nos[0] = mv1->origno;
      nos[1] = mv2->origno;
      nos[2] = mv3->origno;
    }
    else {
      for (int j = 0; j < 3; j++) {
        cos[j] = verts[j]->co;
        nos[j] = verts[j]->no;
      }
    }

    if (ray_face_intersection_depth_tri(
            ray_start, isect_precalc, cos[0], cos[1], cos[2], depth, back_depth, hit_count)) {
      hit = true;

      if (r_face_normal) {
        normal_tri_v3(r_face_normal, cos[0], cos[1], cos[2]);
      }

      if (r_active_vertex) {
        float location[3] = {0.0f};
        madd_v3_v3v3fl(location, ray_start, ray_normal, *depth);
        for (int j = 0; j < 3; j++) {
          if (j == 0 ||
              len_squared_v3v3(location, cos[j]) < len_squared_v3v3(location, nearest_vertex_co)) {
            copy_v3_v3(nearest_vertex_co, cos[j]);
            r_active_vertex->i = (intptr_t)verts[j];
          }
        }
      }

      if (r_active_face) {
        *r_active_face = tri->f;
      }
    }
  }

  return hit;
}

bool BKE_pbvh_bmesh_node_raycast_detail(PBVH *pbvh,
                                        PBVHNode *node,
                                        const float ray_start[3],
                                        struct IsectRayPrecalc *isect_precalc,
                                        float *depth,
                                        float *r_edge_length)
{
  if (node->flag & PBVH_FullyHidden) {
    return false;
  }

  BKE_pbvh_bmesh_check_tris(pbvh, node);
  for (int i = 0; i < node->tribuf->tottri; i++) {
    PBVHTri *tri = node->tribuf->tris + i;
    BMVert *v1 = (BMVert *)node->tribuf->verts[tri->v[0]].i;
    BMVert *v2 = (BMVert *)node->tribuf->verts[tri->v[1]].i;
    BMVert *v3 = (BMVert *)node->tribuf->verts[tri->v[2]].i;
    BMFace *f = (BMFace *)tri->f.i;

    if (pbvh_poly_hidden(pbvh, f)) {
      continue;
    }

    bool hit_local = ray_face_intersection_tri(
        ray_start, isect_precalc, v1->co, v2->co, v3->co, depth);

    if (hit_local) {
      float len1 = len_squared_v3v3(v1->co, v2->co);
      float len2 = len_squared_v3v3(v2->co, v3->co);
      float len3 = len_squared_v3v3(v3->co, v1->co);

      /* detail returned will be set to the maximum allowed size, so take max here */
      *r_edge_length = sqrtf(max_fff(len1, len2, len3));

      return true;
    }
  }

  return false;
}

bool pbvh_bmesh_node_nearest_to_ray(PBVH *pbvh,
                                    PBVHNode *node,
                                    const float ray_start[3],
                                    const float ray_normal[3],
                                    float *depth,
                                    float *dist_sq,
                                    bool use_original,
                                    int stroke_id)
{
  bool hit = false;

  BKE_pbvh_bmesh_check_tris(pbvh, node);
  PBVHTriBuf *tribuf = node->tribuf;
  const int cd_sculpt_vert = pbvh->cd_sculpt_vert;

  for (int i = 0; i < tribuf->tottri; i++) {
    PBVHTri *tri = tribuf->tris + i;
    BMFace *f = (BMFace *)tri->f.i;

    if (pbvh_poly_hidden(pbvh, f)) {
      continue;
    }

    BMVert *v1 = (BMVert *)tribuf->verts[tri->v[0]].i;
    BMVert *v2 = (BMVert *)tribuf->verts[tri->v[1]].i;
    BMVert *v3 = (BMVert *)tribuf->verts[tri->v[2]].i;

    float *co1, *co2, *co3;

    if (use_original) {
      BKE_pbvh_bmesh_check_origdata(pbvh, v1, stroke_id);
      BKE_pbvh_bmesh_check_origdata(pbvh, v2, stroke_id);
      BKE_pbvh_bmesh_check_origdata(pbvh, v3, stroke_id);

      co1 = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, v1)->origco;
      co2 = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, v2)->origco;
      co3 = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, v3)->origco;
    }
    else {
      co1 = v1->co;
      co2 = v2->co;
      co3 = v3->co;
    }

    hit |= ray_face_nearest_tri(ray_start, ray_normal, co1, co2, co3, depth, dist_sq);
  }

  return hit;
}

struct UpdateNormalsTaskData {
  PBVHNode *node;
  Vector<BMVert *> border_verts;
  int cd_sculpt_vert;
  int cd_vert_node_offset;
  int cd_face_node_offset;
  int node_nr;
};

static void pbvh_update_normals_task_cb(void *__restrict userdata,
                                        const int n,
                                        const TaskParallelTLS *__restrict /* tls */)
{
  BMVert *v;
  BMFace *f;
  UpdateNormalsTaskData *data = ((UpdateNormalsTaskData *)userdata) + n;
  PBVHNode *node = data->node;
  const int node_nr = data->node_nr;

  const int cd_face_node_offset = data->cd_face_node_offset;
  const int cd_vert_node_offset = data->cd_vert_node_offset;

  node->flag |= PBVH_UpdateCurvatureDir;

#ifdef NORMAL_VERT_BAD
#  undef NORMAL_VERT_BAD
#endif
#define NORMAL_VERT_BAD(v) \
  (!v->e || BM_ELEM_CD_GET_INT((v), cd_vert_node_offset) != node_nr || \
   (BKE_PBVH_SCULPTVERT(data->cd_sculpt_vert, (v))->flag & SCULPTVERT_PBVH_BOUNDARY))

  const char tag = BM_ELEM_TAG_ALT;

  TGSET_ITER (v, node->bm_unique_verts) {
    PBVH_CHECK_NAN(v->no);

    if (NORMAL_VERT_BAD(v)) {
      v->head.hflag |= tag;
      data->border_verts.append(v);
      continue;
    }

    v->head.hflag &= ~tag;

    BMEdge *e = v->e;
    do {
      BMLoop *l = e->l;

      if (!l) {
        continue;
      }

      do {
        if (BM_ELEM_CD_GET_INT(l->f, cd_face_node_offset) != node_nr) {
          v->head.hflag |= tag;
          goto loop_exit;
        }
      } while ((l = l->radial_next) != e->l);
    } while ((e = BM_DISK_EDGE_NEXT(e, v)) != v->e);
  loop_exit:

    if (v->head.hflag & tag) {
      data->border_verts.append(v);
      continue;
    }

    zero_v3(v->no);
  }
  TGSET_ITER_END

  TGSET_ITER (f, node->bm_faces) {
    BM_face_normal_update(f);

    PBVH_CHECK_NAN(f->no);

    BMLoop *l = f->l_first;
    do {
      PBVH_CHECK_NAN(l->v->no);

      if (BM_ELEM_CD_GET_INT(l->v, cd_vert_node_offset) == node_nr && !(l->v->head.hflag & tag)) {
        add_v3_v3(l->v->no, f->no);
      }
    } while ((l = l->next) != f->l_first);
  }
  TGSET_ITER_END

  TGSET_ITER (v, node->bm_unique_verts) {
    PBVH_CHECK_NAN(v->no);

    if (dot_v3v3(v->no, v->no) == 0.0f) {
      data->border_verts.append(v);

      continue;
    }

    if (!(v->head.hflag & tag)) {
      normalize_v3(v->no);
    }
  }
  TGSET_ITER_END

  node->flag &= ~PBVH_UpdateNormals;
}

void pbvh_bmesh_normals_update(PBVH *pbvh, PBVHNode **nodes, int totnode)
{
  TaskParallelSettings settings;
  Vector<UpdateNormalsTaskData> datas;
  datas.resize(totnode);

  for (int i : IndexRange(totnode)) {
    datas[i].node = nodes[i];
    datas[i].cd_sculpt_vert = pbvh->cd_sculpt_vert;
    datas[i].cd_vert_node_offset = pbvh->cd_vert_node_offset;
    datas[i].cd_face_node_offset = pbvh->cd_face_node_offset;
    datas[i].node_nr = nodes[i] - pbvh->nodes;

    BKE_pbvh_bmesh_check_tris(pbvh, nodes[i]);
  }

  BKE_pbvh_parallel_range_settings(&settings, true, totnode);
  BLI_task_parallel_range(0, totnode, datas.data(), pbvh_update_normals_task_cb, &settings);

  /* not sure it's worth calling BM_mesh_elem_index_ensure here */
#if 0
  BLI_bitmap *visit = BLI_BITMAP_NEW(bm->totvert, "visit");
  BM_mesh_elem_index_ensure(bm, BM_VERT);
#endif

  for (int i = 0; i < totnode; i++) {
    UpdateNormalsTaskData *data = &datas[i];

#if 0
    printf("%.2f%% : %d %d\n",
      100.0f * (float)data->tot_border_verts / (float)data->node->bm_unique_verts->length,
      data->tot_border_verts,
      data->node->bm_unique_verts->length);
#endif

    for (int j = 0; j < data->border_verts.size(); j++) {
      BMVert *v = data->border_verts[j];

      if (BM_elem_is_free((BMElem *)v, BM_VERT)) {
        printf("%s: error, v was freed!\n", __func__);
        continue;
      }

#if 0
      if (v->head.index < 0 || v->head.index >= bm->totvert) {
        printf("%s: error, v->head.index was out of bounds!\n", __func__);
        continue;
      }

      if (BLI_BITMAP_TEST(visit, v->head.index)) {
        continue;
      }

      BLI_BITMAP_ENABLE(visit, v->head.index);
#endif

      // manual iteration
      BMEdge *e = v->e;

      if (!e) {
        continue;
      }

      zero_v3(v->no);

      do {
        if (e->l) {
          add_v3_v3(v->no, e->l->f->no);
        }
        e = BM_DISK_EDGE_NEXT(e, v);
      } while (e != v->e);

      normalize_v3(v->no);
    }
  }

#if 0
  MEM_SAFE_FREE(visit);
#endif
}

struct FastNodeBuildInfo {
  int totface; /* number of faces */
  int start;   /* start of faces in array */
  int depth;
  int node_index;
  struct FastNodeBuildInfo *child1;
  struct FastNodeBuildInfo *child2;
  float cent[3], no[3];
  int tag;
};

/**
 * Recursively split the node if it exceeds the leaf_limit.
 * This function is multi-thread-able since each invocation applies
 * to a sub part of the arrays.
 */
static void pbvh_bmesh_node_limit_ensure_fast(PBVH *pbvh,
                                              BMFace **nodeinfo,
                                              BBC *bbc_array,
                                              FastNodeBuildInfo *node,
                                              Vector<FastNodeBuildInfo *> &leaves,
                                              MemArena *arena)
{
  FastNodeBuildInfo *child1, *child2;

  if (node->totface <= pbvh->leaf_limit || node->depth >= PBVH_STACK_FIXED_DEPTH) {
    return;
  }

  /* Calculate bounding box around primitive centroids */
  BB cb;
  BB_reset(&cb);
  for (int i = 0; i < node->totface; i++) {
    BMFace *f = nodeinfo[i + node->start];
    BBC *bbc = &bbc_array[BM_elem_index_get(f)];

    BB_expand(&cb, bbc->bcentroid);
  }

  /* initialize the children */

  /* Find widest axis and its midpoint */
  const int axis = BB_widest_axis(&cb);
  const float mid = (cb.bmax[axis] + cb.bmin[axis]) * 0.5f;

  int num_child1 = 0, num_child2 = 0;

  /* split vertices along the middle line */
  const int end = node->start + node->totface;
  for (int i = node->start; i < end - num_child2; i++) {
    BMFace *f = nodeinfo[i];
    BBC *bbc = &bbc_array[BM_elem_index_get(f)];

    if (bbc->bcentroid[axis] > mid) {
      int i_iter = end - num_child2 - 1;
      int candidate = -1;
      /* found a face that should be part of another node, look for a face to substitute with */

      for (; i_iter > i; i_iter--) {
        BMFace *f_iter = nodeinfo[i_iter];
        const BBC *bbc_iter = &bbc_array[BM_elem_index_get(f_iter)];
        if (bbc_iter->bcentroid[axis] <= mid) {
          candidate = i_iter;
          break;
        }

        num_child2++;
      }

      if (candidate != -1) {
        BMFace *tmp = nodeinfo[i];
        nodeinfo[i] = nodeinfo[candidate];
        nodeinfo[candidate] = tmp;
        /* increase both counts */
        num_child1++;
        num_child2++;
      }
      else {
        /* not finding candidate means second half of array part is full of
         * second node parts, just increase the number of child nodes for it */
        num_child2++;
      }
    }
    else {
      num_child1++;
    }
  }

  /* ensure at least one child in each node */
  if (num_child2 == 0) {
    num_child2++;
    num_child1--;
  }
  else if (num_child1 == 0) {
    num_child1++;
    num_child2--;
  }

  /* at this point, faces should have been split along the array range sequentially,
   * each sequential part belonging to one node only */
  BLI_assert((num_child1 + num_child2) == node->totface);

  node->child1 = child1 = (FastNodeBuildInfo *)BLI_memarena_alloc(arena,
                                                                  sizeof(FastNodeBuildInfo));
  node->child2 = child2 = (FastNodeBuildInfo *)BLI_memarena_alloc(arena,
                                                                  sizeof(FastNodeBuildInfo));

  child1->totface = num_child1;
  child1->start = node->start;
  child1->depth = node->depth + 1;

  child2->totface = num_child2;
  child2->start = node->start + num_child1;
  child2->depth = node->depth + 2;

  child1->child1 = child1->child2 = child2->child1 = child2->child2 = nullptr;

  pbvh_bmesh_node_limit_ensure_fast(pbvh, nodeinfo, bbc_array, child1, leaves, arena);
  pbvh_bmesh_node_limit_ensure_fast(pbvh, nodeinfo, bbc_array, child2, leaves, arena);

  if (!child1->child1 && !child1->child2) {
    leaves.append(child1);
  }

  if (!child2->child1 && !child2->child2) {
    leaves.append(child2);
  }
}

struct LeafBuilderThreadData {
  PBVH *pbvh;
  BMFace **nodeinfo;
  BBC *bbc_array;
  Vector<FastNodeBuildInfo *> leaves;
};

static void pbvh_bmesh_create_leaf_fast_task_cb(void *__restrict userdata,
                                                const int i,
                                                const TaskParallelTLS *__restrict /* tls */)
{
  LeafBuilderThreadData *data = (LeafBuilderThreadData *)userdata;
  PBVH *pbvh = data->pbvh;
  BMFace **nodeinfo = data->nodeinfo;
  BBC *bbc_array = data->bbc_array;
  struct FastNodeBuildInfo *node = data->leaves[i];

  /* node does not have children so it's a leaf node, populate with faces and tag accordingly
   * this is an expensive part but it's not so easily thread-able due to vertex node indices */
  // const int cd_vert_node_offset = pbvh->cd_vert_node_offset;
  const int cd_face_node_offset = pbvh->cd_face_node_offset;

  PBVHNode *n = pbvh->nodes + node->node_index;
  const int node_index = node->node_index;

  bool has_visible = false;

  /* Build GPU buffers for new node */

  n->flag = PBVH_Leaf | PBVH_UpdateTris | PBVH_UpdateBB | PBVH_UpdateOriginalBB |
            PBVH_UpdateTriAreas | PBVH_UpdateColor | PBVH_UpdateVisibility |
            PBVH_UpdateDrawBuffers | PBVH_RebuildDrawBuffers | PBVH_UpdateCurvatureDir |
            PBVH_UpdateTriAreas | PBVH_UpdateMask | PBVH_UpdateRedraw;

  n->bm_faces = BLI_table_gset_new_ex("bm_faces", node->totface);

  /* Create vert hash sets */
  n->bm_unique_verts = BLI_table_gset_new_ex("bm_unique_verts", node->totface * 3);
  n->bm_other_verts = BLI_table_gset_new_ex("bm_other_verts", node->totface * 3);

  BB_reset(&n->vb);

  const int end = node->start + node->totface;

  for (int i = node->start; i < end; i++) {
    BMFace *f = nodeinfo[i];
    BBC *bbc = &bbc_array[BM_elem_index_get(f)];

    /* Update ownership of faces */
    BLI_table_gset_insert(n->bm_faces, f);
    BM_ELEM_CD_SET_INT(f, cd_face_node_offset, node_index);

    /* Update vertices */
    BMLoop *l_first = BM_FACE_FIRST_LOOP(f);
    BMLoop *l_iter = l_first;
    do {
      BMVert *v = l_iter->v;

      int old = BM_ELEM_CD_GET_INT(v, pbvh->cd_vert_node_offset);

      char *ptr = (char *)v->head.data;
      ptr += pbvh->cd_vert_node_offset;

      if (old == DYNTOPO_NODE_NONE &&
          atomic_cas_int32((int32_t *)ptr, DYNTOPO_NODE_NONE, node_index) == DYNTOPO_NODE_NONE) {
        BLI_table_gset_insert(n->bm_unique_verts, v);
      }
      else {
        BLI_table_gset_add(n->bm_other_verts, v);
      }
    } while ((l_iter = l_iter->next) != l_first);

    /* Update node bounding box */
    if (!pbvh_poly_hidden(pbvh, f)) {
      has_visible = true;
    }

    BB_expand_with_bb(&n->vb, (BB *)bbc);
  }

  BLI_assert(n->vb.bmin[0] <= n->vb.bmax[0] && n->vb.bmin[1] <= n->vb.bmax[1] &&
             n->vb.bmin[2] <= n->vb.bmax[2]);

  n->orig_vb = n->vb;

  BKE_pbvh_node_fully_hidden_set(n, !has_visible);
  n->flag |= PBVH_UpdateNormals | PBVH_UpdateCurvatureDir;
}

static void pbvh_bmesh_create_nodes_fast_recursive_create(PBVH *pbvh,
                                                          BMFace **nodeinfo,
                                                          BBC *bbc_array,
                                                          struct FastNodeBuildInfo *node)
{
  /* two cases, node does not have children or does have children */
  if (node->child1) {
    int children_offset = pbvh->totnode;
    pbvh_grow_nodes(pbvh, pbvh->totnode + 2);

    PBVHNode *n = pbvh->nodes + node->node_index;
    n->children_offset = children_offset;

    n->depth = node->depth;
    (n + 1)->depth = node->child1->depth;
    (n + 2)->depth = node->child2->depth;

    node->child1->node_index = children_offset;
    node->child2->node_index = children_offset + 1;

    pbvh_bmesh_create_nodes_fast_recursive_create(pbvh, nodeinfo, bbc_array, node->child1);
    pbvh_bmesh_create_nodes_fast_recursive_create(pbvh, nodeinfo, bbc_array, node->child2);
  }
}

static void pbvh_bmesh_create_nodes_fast_recursive_final(PBVH *pbvh,
                                                         BMFace **nodeinfo,
                                                         BBC *bbc_array,
                                                         struct FastNodeBuildInfo *node)
{
  /* two cases, node does not have children or does have children */
  if (node->child1) {
    pbvh_bmesh_create_nodes_fast_recursive_final(pbvh, nodeinfo, bbc_array, node->child1);
    pbvh_bmesh_create_nodes_fast_recursive_final(pbvh, nodeinfo, bbc_array, node->child2);

    PBVHNode *n = pbvh->nodes + node->node_index;

    /* Update bounding box */
    BB_reset(&n->vb);
    BB_expand_with_bb(&n->vb, &pbvh->nodes[node->child1->node_index].vb);
    BB_expand_with_bb(&n->vb, &pbvh->nodes[node->child2->node_index].vb);
    n->orig_vb = n->vb;
  }
}

/***************************** Public API *****************************/

struct FSetTemp {
  BMVert *v;
  int fset;
  bool boundary;
};

int BKE_pbvh_do_fset_symmetry(int fset, const int symflag, const float *co)
{
  fset = abs(fset);

  // don't affect base face set
  if (fset == 1) {
    return 1;
  }

  /*
    flag symmetry by shifting by 24 bits;
    surely we don't need more then 8 million face sets?
  */
  if (co[0] < 0.0f) {
    fset |= (symflag & 1) << 24;
  }

  if (co[1] < 0.0f) {
    fset |= (symflag & 2) << 24;
  }

  if (co[2] < 0.0f) {
    fset |= (symflag & 4) << 24;
  }

  return fset;
}

//#define MV_COLOR_BOUNDARY

#ifdef MV_COLOR_BOUNDARY
static int color_boundary_key(float col[4])
{
  const float steps = 2.0f;
  float hsv[3];

  rgb_to_hsv(col[0], col[1], col[2], hsv, hsv + 1, hsv + 2);

  int x = (int)((hsv[0] * 0.5f + 0.5f) * steps + 0.5f);
  int y = (int)(hsv[1] * steps + 0.5f);
  int z = (int)(hsv[2] * steps + 0.5f);

  return z * steps * steps + y * steps + x;
}
#endif

/* calls atomic_cas_uint32 on two adjacent (and int aligned) shorts */
BLI_INLINE void atomic_cas_short2(ushort *base, ushort olda, ushort oldb, ushort newa, ushort newb)
{
  uint oldi, newi;

  ((ushort *)&oldi)[0] = olda;
  ((ushort *)&oldi)[1] = oldb;

  ((ushort *)&newi)[0] = newa;
  ((ushort *)&newi)[1] = newb;

  atomic_cas_uint32((uint32_t *)base, oldi, newi);
}

void bke_pbvh_update_vert_boundary(int cd_sculpt_vert,
                                   int cd_faceset_offset,
                                   int cd_vert_node_offset,
                                   int cd_face_node_offset,
                                   int cd_vcol,
                                   int cd_boundary_flag,
                                   BMVert *v,
                                   int bound_symmetry,
                                   const CustomData *ldata,
                                   const int totuv)
{
  MSculptVert *mv = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, v);

  float curv = 0.0f, totcurv = 0.0f;

  int newflag = mv->flag;
  int oldflag = newflag;
  int oldval = mv->valence;
  int boundflag = 0;

  BMEdge *e = v->e;
  newflag &= ~(SCULPTVERT_VERT_FSET_HIDDEN | SCULPTVERT_PBVH_BOUNDARY);

  ushort stroke_id = (ushort)mv->stroke_id;

  if (!e) {
    boundflag |= SCULPT_BOUNDARY_MESH;

    int oldboundflag = BM_ELEM_CD_GET_INT(v, cd_boundary_flag);

    atomic_cas_int32(&mv->flag, oldflag, newflag);
    atomic_cas_int32(BM_ELEM_CD_PTR<int *>(v, cd_boundary_flag), oldboundflag, boundflag);

    atomic_cas_short2(&mv->valence, (ushort)oldval, stroke_id, 0, stroke_id);

    return;
  }

  int val = 0;

  int ni = BM_ELEM_CD_GET_INT(v, cd_vert_node_offset);

  int sharpcount = 0;
  int seamcount = 0;
  int quadcount = 0;

#ifdef MV_COLOR_BOUNDARY
  int last_key = -1;
#endif

#if 0
  struct FaceSetRef {
    int fset;
    BMVert *v2;
    BMEdge *e;
  } *fsets = nullptr;
#endif
  Vector<int, 16> fsets;

  float(*lastuv)[2] = (float(*)[2])BLI_array_alloca(lastuv, totuv);
  float(*lastuv2)[2] = (float(*)[2])BLI_array_alloca(lastuv2, totuv);

  int *disjount_uv_count = (int *)BLI_array_alloca(disjount_uv_count, totuv);
  int *cd_uvs = (int *)BLI_array_alloca(cd_uvs, totuv);
  int base_uv_idx = ldata->typemap[CD_PROP_FLOAT2];
  bool uv_first = true;

  for (int i = 0; i < totuv; i++) {
    CustomDataLayer *layer = ldata->layers + base_uv_idx + i;
    cd_uvs[i] = layer->offset;
    disjount_uv_count[i] = 0;
  }

  do {
    BMVert *v2 = v == e->v1 ? e->v2 : e->v1;

#if 0
    float tmp[3];
    sub_v3_v3v3(tmp, v2->co, v->co);
    madd_v3_v3fl(avg, v->no, -dot_v3v3(v->no, tmp));
    // madd_v3_v3fl(tmp, v->no, -dot_v3v3(v->no, tmp));
    add_v3_v3(avg, tmp);

    avg_len += len_squared_v3(tmp);
    totcurv += 1.0f;
#endif

    if (BM_ELEM_CD_GET_INT(v2, cd_vert_node_offset) != ni) {
      newflag |= SCULPTVERT_PBVH_BOUNDARY;
    }

    if (e->head.hflag & BM_ELEM_SEAM) {
      boundflag |= SCULPT_BOUNDARY_SEAM;
      seamcount++;

      if (seamcount > 2) {
        boundflag |= SCULPT_CORNER_SEAM;
      }
    }

#ifdef MV_COLOR_BOUNDARY
#  error "fix me! need to handle loops too"
#endif

#ifdef MV_COLOR_BOUNDARY  // CURRENTLY BROKEN
    if (cd_vcol >= 0) {
      float *color1 = BM_ELEM_CD_PTR<float *>(v, cd_vcol);
      float *color2 = BM_ELEM_CD_PTR<float *>(v2, cd_vcol);

      int colorkey1 = color_boundary_key(color1);
      int colorkey2 = color_boundary_key(color2);

      if (colorkey1 != colorkey2) {
        boundflag |= SCUKPT_BOUNDARY_FACE_SETS;
      }
    }
#endif

    if (!(e->head.hflag & BM_ELEM_SMOOTH)) {
      boundflag |= SCULPT_BOUNDARY_SHARP;
      sharpcount++;

      if (sharpcount > 2) {
        boundflag |= SCULPT_CORNER_SHARP;
      }
    }

    if (e->l) {
      /* detect uv island boundaries */
      if (totuv) {
        BMLoop *l_iter = e->l;
        do {
          BMLoop *l = l_iter->v != v ? l_iter->next : l_iter;

          for (int i = 0; i < totuv; i++) {
            float *luv = BM_ELEM_CD_PTR<float *>(l, cd_uvs[i]);

            if (uv_first) {
              copy_v2_v2(lastuv[i], luv);
              copy_v2_v2(lastuv2[i], luv);

              continue;
            }

            const float uv_snap_limit = 0.01f * 0.01f;

            float dist = len_squared_v2v2(luv, lastuv[i]);
            bool same = dist <= uv_snap_limit;

            bool corner = len_squared_v2v2(lastuv[i], lastuv2[i]) > uv_snap_limit &&
                          len_squared_v2v2(lastuv[i], luv) > uv_snap_limit &&
                          len_squared_v2v2(lastuv2[i], luv) > uv_snap_limit;

            if (!same) {
              boundflag |= SCULPT_BOUNDARY_UV;
            }

            if (corner) {
              boundflag |= SCULPT_CORNER_UV;
            }

            if (!same) {
              copy_v2_v2(lastuv2[i], lastuv[i]);
              copy_v2_v2(lastuv[i], luv);
            }
          }

          uv_first = false;
        } while ((l_iter = l_iter->radial_next) != e->l);
      }

      if (BM_ELEM_CD_GET_INT(e->l->f, cd_face_node_offset) != ni) {
        newflag |= SCULPTVERT_PBVH_BOUNDARY;
      }

      if (e->l != e->l->radial_next) {
        if (e->l->f->len > 3) {
          quadcount++;
        }

        if (e->l->radial_next->f->len > 3) {
          quadcount++;
        }

        if (BM_ELEM_CD_GET_INT(e->l->radial_next->f, cd_face_node_offset) != ni) {
          newflag |= SCULPTVERT_PBVH_BOUNDARY;
        }
      }

      if (e->l->f->len > 3) {
        newflag |= SCULPTVERT_NEED_TRIANGULATE;
      }

      if (cd_faceset_offset != -1) {
        int fset = BKE_pbvh_do_fset_symmetry(
            BM_ELEM_CD_GET_INT(e->l->f, cd_faceset_offset), bound_symmetry, v2->co);

        bool ok = true;
        for (int i = 0; i < fsets.size(); i++) {
          if (fsets[i] == fset) {
            ok = false;
          }
        }

        if (ok) {
          fsets.append(fset);
        }
      }

      // also check e->l->radial_next, in case we are not manifold
      // which can mess up the loop order
      if (e->l->radial_next != e->l) {
        float th = saacos(dot_v3v3(e->l->f->no, e->l->radial_next->f->no));

        th *= M_1_PI * 0.25f;
        // th = th * 0.5 + 0.5;
        curv += th;
        totcurv += 1.0f;

        if (cd_faceset_offset != -1) {
          // fset = abs(BM_ELEM_CD_GET_INT(e->l->radial_next->f, cd_faceset_offset));
          int fset2 = BKE_pbvh_do_fset_symmetry(
              BM_ELEM_CD_GET_INT(e->l->radial_next->f, cd_faceset_offset), bound_symmetry, v2->co);

          bool ok2 = true;
          for (int i = 0; i < fsets.size(); i++) {
            if (fsets[i] == fset2) {
              ok2 = false;
            }
          }

          if (ok2) {
            fsets.append(fset2);
          }
        }

        if (e->l->radial_next->f->len > 3) {
          newflag |= SCULPTVERT_NEED_TRIANGULATE;
        }
      }
    }

    if (!e->l || e->l->radial_next == e->l) {
      boundflag |= SCULPT_BOUNDARY_MESH;
    }

    val++;
  } while ((e = BM_DISK_EDGE_NEXT(e, v)) != v->e);

  if (fsets.size() > 1) {
    boundflag |= SCULPT_BOUNDARY_FACE_SET;
  }

  if (fsets.size() > 2) {
    boundflag |= SCULPT_CORNER_FACE_SET;
  }

  if (sharpcount == 1) {
    boundflag |= SCULPT_CORNER_SHARP;
  }

  if (seamcount == 1) {
    boundflag |= SCULPT_CORNER_SEAM;
  }

  if ((boundflag & SCULPT_BOUNDARY_MESH) && quadcount >= 3) {
    boundflag |= SCULPT_CORNER_MESH;
  }

#if 0
  if (totcurv > 0.0f) {
    mul_v3_fl(avg, 1.0f / totcurv);
    avg_len /= totcurv;
  }

  if (avg_len > 0.0f) {
    curv = len_squared_v3(avg) / avg_len;
  }
  else {
    curv = 0.0f;
  }
#else
  if (totcurv > 0.0f) {
    curv /= totcurv;
  }
#endif

  int oldboundflag = BM_ELEM_CD_GET_INT(v, cd_boundary_flag);

  atomic_cas_int32(&mv->flag, oldflag, newflag);
  atomic_cas_int32(BM_ELEM_CD_PTR<int *>(v, cd_boundary_flag), oldboundflag, boundflag);

  atomic_cas_short2(&mv->valence, (ushort)oldval, stroke_id, (ushort)val, stroke_id);

  /* no atomic_cas_int16, so do origmask and curv at once */

  ushort newcurv = (unsigned short)(min_ff(fabsf(curv), 1.0f) * 65535.0f);
  ushort oldcurv = mv->curv;
  ushort origmask = mv->origmask;

  atomic_cas_short2(&mv->origmask, origmask, oldcurv, origmask, newcurv);
}

bool BKE_pbvh_check_vert_boundary(PBVH *pbvh, BMVert *v)
{
  return pbvh_check_vert_boundary(pbvh, v);
}

void BKE_pbvh_update_vert_boundary(int cd_sculpt_vert,
                                   int cd_faceset_offset,
                                   int cd_vert_node_offset,
                                   int cd_face_node_offset,
                                   int cd_vcol,
                                   int cd_boundary_flag,
                                   BMVert *v,
                                   int bound_symmetry,
                                   const CustomData *ldata,
                                   const int totuv,
                                   const bool do_uvs)
{
  bke_pbvh_update_vert_boundary(cd_sculpt_vert,
                                cd_faceset_offset,
                                cd_vert_node_offset,
                                cd_face_node_offset,
                                cd_vcol,
                                cd_boundary_flag,
                                v,
                                bound_symmetry,
                                ldata,
                                do_uvs ? totuv : 0);
}

/*Used by symmetrize to update boundary flags*/
void BKE_pbvh_recalc_bmesh_boundary(PBVH *pbvh)
{
  BMVert *v;
  BMIter iter;

  BM_ITER_MESH (v, &iter, pbvh->header.bm, BM_VERTS_OF_MESH) {
    bke_pbvh_update_vert_boundary(pbvh->cd_sculpt_vert,
                                  pbvh->cd_faceset_offset,
                                  pbvh->cd_vert_node_offset,
                                  pbvh->cd_face_node_offset,
                                  pbvh->cd_vcol_offset,
                                  pbvh->cd_boundary_flag,
                                  v,
                                  pbvh->boundary_symmetry,
                                  &pbvh->header.bm->ldata,
                                  pbvh->flags & PBVH_IGNORE_UVS ? 0 : pbvh->totuv);
  }
}

void BKE_pbvh_update_sculpt_verts(PBVH *pbvh)
{
  int totvert = pbvh->totvert;

  for (int i = 0; i < totvert; i++) {
    PBVHVertRef vertex = BKE_pbvh_index_to_vertex(pbvh, i);

    BKE_pbvh_get_origvert(pbvh, vertex, nullptr, nullptr, nullptr);
  }
}

void BKE_pbvh_set_idmap(PBVH *pbvh, BMIdMap *idmap)
{
  pbvh->bm_idmap = idmap;
}

/* Build a PBVH from a BMesh */
void BKE_pbvh_build_bmesh(PBVH *pbvh,
                          Mesh *me,
                          BMesh *bm,
                          bool smooth_shading,
                          BMLog *log,
                          BMIdMap *idmap,
                          const int cd_vert_node_offset,
                          const int cd_face_node_offset,
                          const int cd_sculpt_vert,
                          const int cd_face_areas,
                          const int cd_boundary_flag,
                          bool fast_draw,
                          bool update_sculptverts)
{
  // coalese_pbvh(pbvh, bm);

  pbvh->bm_idmap = idmap;

  pbvh->cd_face_area = cd_face_areas;
  pbvh->cd_vert_node_offset = cd_vert_node_offset;
  pbvh->cd_face_node_offset = cd_face_node_offset;
  pbvh->cd_vert_mask_offset = CustomData_get_offset(&bm->vdata, CD_PAINT_MASK);
  pbvh->cd_sculpt_vert = cd_sculpt_vert;
  pbvh->cd_boundary_flag = cd_boundary_flag;

  pbvh->mesh = me;

  smooth_shading |= fast_draw;

  pbvh->header.bm = bm;

  BKE_pbvh_bmesh_detail_size_set(pbvh, 0.75f, 0.4f);

  pbvh->header.type = PBVH_BMESH;
  pbvh->bm_log = log;
  pbvh->cd_faceset_offset = CustomData_get_offset_named(
      &pbvh->header.bm->pdata, CD_PROP_INT32, ".sculpt_face_set");

  int tottri = poly_to_tri_count(bm->totface, bm->totloop);

  /* TODO: choose leaf limit better */
  if (tottri > 4 * 1000 * 1000) {
    pbvh->leaf_limit = 10000;
  }
  else {
    pbvh->leaf_limit = 1000;
  }

  BMIter iter;
  BMVert *v;

  if (smooth_shading) {
    pbvh->flags |= PBVH_DYNTOPO_SMOOTH_SHADING;
  }

  if (fast_draw) {
    pbvh->flags |= PBVH_FAST_DRAW;
  }

  /* bounding box array of all faces, no need to recalculate every time */
  BBC *bbc_array = MEM_cnew_array<BBC>(bm->totface, "BBC");
  BMFace **nodeinfo = MEM_cnew_array<BMFace *>(bm->totface, "nodeinfo");
  MemArena *arena = BLI_memarena_new(BLI_MEMARENA_STD_BUFSIZE, "fast PBVH node storage");

  BMFace *f;
  int i;
  BM_ITER_MESH_INDEX (f, &iter, bm, BM_FACES_OF_MESH, i) {
    BBC *bbc = &bbc_array[i];
    BMLoop *l_first = BM_FACE_FIRST_LOOP(f);
    BMLoop *l_iter = l_first;

    // check for currupted faceset
    if (pbvh->cd_faceset_offset != -1 && BM_ELEM_CD_GET_INT(f, pbvh->cd_faceset_offset) == 0) {
      BM_ELEM_CD_SET_INT(f, pbvh->cd_faceset_offset, 1);
    }

    BB_reset((BB *)bbc);
    do {
      BB_expand((BB *)bbc, l_iter->v->co);
    } while ((l_iter = l_iter->next) != l_first);
    BBC_update_centroid(bbc);

    /* so we can do direct lookups on 'bbc_array' */
    BM_elem_index_set(f, i); /* set_dirty! */
    nodeinfo[i] = f;
    BM_ELEM_CD_SET_INT(f, cd_face_node_offset, DYNTOPO_NODE_NONE);
  }
  /* Likely this is already dirty. */
  bm->elem_index_dirty |= BM_FACE;

  BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
    MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, v);

    mv->flag |= SCULPTVERT_NEED_VALENCE | SCULPTVERT_NEED_TRIANGULATE | SCULPTVERT_NEED_DISK_SORT;
    copy_v3_v3(mv->origco, v->co);
    copy_v3_v3(mv->origno, v->no);

    int *flags = BM_ELEM_CD_PTR<int *>(v, pbvh->cd_boundary_flag);
    *flags |= SCULPT_BOUNDARY_NEEDS_UPDATE;

    BM_ELEM_CD_SET_INT(v, cd_vert_node_offset, DYNTOPO_NODE_NONE);
  }

  /* setup root node */
  struct FastNodeBuildInfo rootnode = {0};

  Vector<FastNodeBuildInfo *> leaves;

  rootnode.totface = bm->totface;
  int totleaf = 0;

  /* start recursion, assign faces to nodes accordingly */
  pbvh_bmesh_node_limit_ensure_fast(pbvh, nodeinfo, bbc_array, &rootnode, leaves, arena);

  pbvh_grow_nodes(pbvh, 1);
  rootnode.node_index = 0;

  pbvh_bmesh_create_nodes_fast_recursive_create(pbvh, nodeinfo, bbc_array, &rootnode);
  totleaf = leaves.size();

  if (!totleaf) {
    leaves.append(&rootnode);
    totleaf = 1;
  }

  /* build leaf nodes */
  LeafBuilderThreadData tdata;
  tdata.pbvh = pbvh;
  tdata.nodeinfo = nodeinfo;
  tdata.bbc_array = bbc_array;
  tdata.leaves = leaves;

  TaskParallelSettings settings;
  BKE_pbvh_parallel_range_settings(&settings, true, totleaf);
  BLI_task_parallel_range(0, totleaf, &tdata, pbvh_bmesh_create_leaf_fast_task_cb, &settings);

  /* take root node and visit and populate children recursively */
  pbvh_bmesh_create_nodes_fast_recursive_final(pbvh, nodeinfo, bbc_array, &rootnode);

  BLI_memarena_free(arena);
  MEM_freeN(bbc_array);
  MEM_freeN(nodeinfo);

  if (me) { /* Ensure pbvh->vcol_type, vcol_domain and cd_vcol_offset are up to date. */
    CustomDataLayer *cl;
    eAttrDomain domain;

    BKE_pbvh_get_color_layer(pbvh, me, &cl, &domain);
  }

  /*final check that nodes are sufficiently subdivided*/
  int totnode = pbvh->totnode;

  for (int i = 0; i < totnode; i++) {
    PBVHNode *n = pbvh->nodes + i;

    if (totnode != pbvh->totnode) {
#ifdef PROXY_ADVANCED
      BKE_pbvh_free_proxyarray(pbvh, n);
#endif
    }

    if (n->flag & PBVH_Leaf) {
      /* Recursively split nodes that have gotten too many
       * elements */
      pbvh_bmesh_node_limit_ensure(pbvh, i);
    }
  }

  if (update_sculptverts) {
    BKE_pbvh_update_sculpt_verts(pbvh);
  }

  pbvh_print_mem_size(pbvh);

  /* update face areas */
  const int cd_face_area = pbvh->cd_face_area;
  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (!(node->flag & PBVH_Leaf)) {
      continue;
    }

    BKE_pbvh_bmesh_check_tris(pbvh, node);

    node->flag |= PBVH_UpdateTriAreas;
    BKE_pbvh_check_tri_areas(pbvh, node);

    int area_src_i = pbvh->face_area_i ^ 1;
    int area_dst_i = pbvh->face_area_i;

    /* make sure read side of double buffer is set too */
    TGSET_ITER (f, node->bm_faces) {
      float *areabuf = BM_ELEM_CD_PTR<float *>(f, cd_face_area);
      areabuf[area_dst_i] = areabuf[area_src_i];
    }
    TGSET_ITER_END;
  }
}

void BKE_pbvh_set_bm_log(PBVH *pbvh, BMLog *log)
{
  pbvh->bm_log = log;
  BM_log_set_idmap(log, pbvh->bm_idmap);
}

bool BKE_pbvh_bmesh_update_topology_nodes(PBVH *pbvh,
                                          bool (*searchcb)(PBVHNode *node, void *data),
                                          void (*undopush)(PBVHNode *node, void *data),
                                          void *searchdata,
                                          PBVHTopologyUpdateMode mode,
                                          const float center[3],
                                          const float view_normal[3],
                                          float radius,
                                          const bool use_frontface,
                                          const bool use_projected,
                                          int sym_axis,
                                          bool updatePBVH,
                                          DyntopoMaskCB mask_cb,
                                          void *mask_cb_data,
                                          bool disable_surface_relax,
                                          bool is_snake_hook)
{
  bool modified = false;
  Vector<PBVHNode *> nodes;

  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (!(node->flag & PBVH_Leaf) || !searchcb(node, searchdata)) {
      continue;
    }

    if (node->flag & PBVH_Leaf) {
      undopush(node, searchdata);

      nodes.append(node);
    }
  }

  for (PBVHNode *node : nodes) {
    node->flag |= PBVH_UpdateCurvatureDir;
    BKE_pbvh_node_mark_topology_update(node);
  }

  modified = BKE_pbvh_bmesh_update_topology(pbvh,
                                            mode,
                                            center,
                                            view_normal,
                                            radius,
                                            use_frontface,
                                            use_projected,
                                            sym_axis,
                                            updatePBVH,
                                            mask_cb,
                                            mask_cb_data,
                                            0,  // is_snake_hook ? 40960 : 0,
                                            disable_surface_relax,
                                            is_snake_hook);

  return modified;
}

static void pbvh_free_tribuf(PBVHTriBuf *tribuf)
{
  MEM_SAFE_FREE(tribuf->verts);
  MEM_SAFE_FREE(tribuf->tris);
  MEM_SAFE_FREE(tribuf->loops);
  MEM_SAFE_FREE(tribuf->edges);

  BLI_smallhash_release(&tribuf->vertmap);

  tribuf->verts = nullptr;
  tribuf->tris = nullptr;
  tribuf->loops = nullptr;
  tribuf->edges = nullptr;

  tribuf->totloop = tribuf->tottri = tribuf->totedge = tribuf->totvert = 0;

  tribuf->verts_size = 0;
  tribuf->tris_size = 0;
  tribuf->edges_size = 0;
}

PBVHTriBuf *BKE_pbvh_bmesh_get_tris(PBVH *pbvh, PBVHNode *node)
{
  BKE_pbvh_bmesh_check_tris(pbvh, node);

  return node->tribuf;
}

void BKE_pbvh_bmesh_free_tris(PBVH *pbvh, PBVHNode *node)
{
  if (node->tribuf) {
    pbvh_free_tribuf(node->tribuf);
    MEM_freeN(node->tribuf);
    node->tribuf = nullptr;
  }

  if (node->tri_buffers) {
    for (int i = 0; i < node->tot_tri_buffers; i++) {
      pbvh_free_tribuf(node->tri_buffers + i);
    }

    MEM_SAFE_FREE(node->tri_buffers);

    node->tri_buffers = nullptr;
    node->tot_tri_buffers = 0;
  }
}

BLI_INLINE PBVHTri *pbvh_tribuf_add_tri(PBVHTriBuf *tribuf)
{
  tribuf->tottri++;

  if (tribuf->tottri >= tribuf->tris_size) {
    size_t newsize = (size_t)32 + (size_t)tribuf->tris_size + (size_t)(tribuf->tris_size >> 1);

    if (!tribuf->tris) {
      tribuf->tris = MEM_cnew_array<PBVHTri>(newsize, "tribuf tris");
    }
    else {
      tribuf->tris = (PBVHTri *)MEM_reallocN_id(
          (void *)tribuf->tris, sizeof(*tribuf->tris) * newsize, "tribuf tris");
    }

    tribuf->tris_size = newsize;
  }

  return tribuf->tris + tribuf->tottri - 1;
}

BLI_INLINE void pbvh_tribuf_add_vert(PBVHTriBuf *tribuf, PBVHVertRef vertex, BMLoop *l)
{
  tribuf->totvert++;
  tribuf->totloop++;

  if (tribuf->totvert >= tribuf->verts_size) {
    size_t newsize = (size_t)32 + (size_t)(tribuf->verts_size << 1);

    if (!tribuf->verts) {
      tribuf->verts = MEM_cnew_array<PBVHVertRef>(newsize, "tribuf verts");
      tribuf->loops = MEM_cnew_array<intptr_t>(newsize, "tribuf loops");
    }
    else {
      tribuf->verts = (PBVHVertRef *)MEM_reallocN_id(
          (void *)tribuf->verts, sizeof(*tribuf->verts) * newsize, "tribuf verts");
      tribuf->loops = (intptr_t *)MEM_reallocN_id(
          (void *)tribuf->loops, sizeof(*tribuf->loops) * newsize, "tribuf loops");
    }

    tribuf->verts_size = newsize;
  }

  tribuf->verts[tribuf->totvert - 1] = vertex;
  tribuf->loops[tribuf->totloop - 1] = (uintptr_t)l;
}

BLI_INLINE void pbvh_tribuf_add_edge(PBVHTriBuf *tribuf, int v1, int v2)
{
  tribuf->totedge++;

  if (tribuf->totedge >= tribuf->edges_size) {
    size_t newsize = (size_t)32 + (size_t)(tribuf->edges_size << 1);

    if (!tribuf->edges) {
      tribuf->edges = MEM_cnew_array<int>(2ULL * newsize, "tribuf edges");
    }
    else {
      tribuf->edges = (int *)MEM_reallocN_id(
          (void *)tribuf->edges, sizeof(*tribuf->edges) * 2ULL * newsize, "tribuf edges");
    }

    tribuf->edges_size = newsize;
  }

  int i = (tribuf->totedge - 1) * 2;

  tribuf->edges[i] = v1;
  tribuf->edges[i + 1] = v2;
}

void pbvh_bmesh_check_other_verts(PBVHNode *node)
{
  if (!(node->flag & PBVH_UpdateOtherVerts)) {
    return;
  }

  node->flag &= ~PBVH_UpdateOtherVerts;

  if (node->bm_other_verts) {
    BLI_table_gset_free(node->bm_other_verts, nullptr);
  }

  node->bm_other_verts = BLI_table_gset_new("bm_other_verts");
  BMFace *f;

  TGSET_ITER (f, node->bm_faces) {
    BMLoop *l = f->l_first;

    do {
      if (!BLI_table_gset_haskey(node->bm_unique_verts, l->v)) {
        BLI_table_gset_add(node->bm_other_verts, l->v);
      }
    } while ((l = l->next) != f->l_first);
  }
  TGSET_ITER_END;
}

static void pbvh_init_tribuf(PBVHNode *node, PBVHTriBuf *tribuf)
{
  tribuf->tottri = 0;
  tribuf->tris_size = 0;
  tribuf->verts_size = 0;
  tribuf->mat_nr = 0;
  tribuf->tottri = 0;
  tribuf->totvert = 0;
  tribuf->totloop = 0;
  tribuf->totedge = 0;

  tribuf->edges = nullptr;
  tribuf->verts = nullptr;
  tribuf->tris = nullptr;
  tribuf->loops = nullptr;

  BLI_smallhash_init_ex(&tribuf->vertmap, node->bm_unique_verts->length);
}

static uintptr_t tri_loopkey(BMLoop *l, int mat_nr, int cd_fset, int cd_uvs[], int totuv)
{
  uintptr_t key = (uintptr_t)mat_nr;

  key ^= (uintptr_t)l->v;

  if (cd_fset >= 0) {
    // key ^= (uintptr_t)BLI_hash_int(BM_ELEM_CD_GET_INT(l->f, cd_fset));
    key ^= (uintptr_t)BM_ELEM_CD_GET_INT(l->f, cd_fset);
  }

  for (int i = 0; i < totuv; i++) {
    float *luv = BM_ELEM_CD_PTR<float *>(l, cd_uvs[i]);
    float snap = 4196.0f;

    uintptr_t x = (uintptr_t)(luv[0] * snap);
    uintptr_t y = (uintptr_t)(luv[1] * snap);

    uintptr_t key2 = y * snap + x;
    key ^= key2;
  }

  return key;
}
/* In order to perform operations on the original node coordinates
 * (currently just raycast), store the node's triangles and vertices.
 *
 * Skips triangles that are hidden. */
ATTR_NO_OPT bool BKE_pbvh_bmesh_check_tris(PBVH *pbvh, PBVHNode *node)
{
  BMesh *bm = pbvh->header.bm;

  if (!(node->flag & PBVH_UpdateTris) && node->tribuf) {
    return false;
  }

  int totuv = CustomData_number_of_layers(&bm->ldata, CD_PROP_FLOAT2);
  int *cd_uvs = (int *)BLI_array_alloca(cd_uvs, totuv);

  for (int i = 0; i < totuv; i++) {
    int idx = CustomData_get_layer_index_n(&bm->ldata, CD_PROP_FLOAT2, i);
    cd_uvs[i] = bm->ldata.layers[idx].offset;
  }

  node->flag |= PBVH_UpdateOtherVerts;

  int mat_map[MAXMAT];

  for (int i = 0; i < MAXMAT; i++) {
    mat_map[i] = -1;
  }

  if (node->tribuf || node->tri_buffers) {
    BKE_pbvh_bmesh_free_tris(pbvh, node);
  }

  node->tribuf = MEM_cnew<PBVHTriBuf>("node->tribuf");
  pbvh_init_tribuf(node, node->tribuf);

  Vector<BMLoop *, 128> loops;
  Vector<blender::uint3, 128> loops_idx;
  Vector<PBVHTriBuf> tribufs;

  node->flag &= ~PBVH_UpdateTris;

  const int edgeflag = BM_ELEM_TAG_ALT;

  BMFace *f;

  float min[3], max[3];

  INIT_MINMAX(min, max);

  TGSET_ITER (f, node->bm_faces) {
    if (pbvh_poly_hidden(pbvh, f)) {
      continue;
    }

    /* Clear edgeflag for building edge indices later. */
    BMLoop *l = f->l_first;
    do {
      l->e->head.hflag &= ~edgeflag;
    } while ((l = l->next) != f->l_first);

    const int mat_nr = f->mat_nr;

    if (mat_map[mat_nr] == -1) {
      PBVHTriBuf _tribuf = {0};

      mat_map[mat_nr] = tribufs.size();

      pbvh_init_tribuf(node, &_tribuf);
      _tribuf.mat_nr = mat_nr;
      tribufs.append(_tribuf);
    }

#ifdef DYNTOPO_DYNAMIC_TESS
    int tottri = (f->len - 2);

    loops.resize(f->len);
    loops_idx.resize(tottri);

    BM_face_calc_tessellation(f, true, loops.data(), (uint(*)[3])loops_idx.data());

    for (int i = 0; i < tottri; i++) {
      PBVHTri *tri = pbvh_tribuf_add_tri(node->tribuf);
      PBVHTriBuf *mat_tribuf = &tribufs[mat_map[mat_nr]];
      PBVHTri *mat_tri = pbvh_tribuf_add_tri(mat_tribuf);

      tri->eflag = mat_tri->eflag = 0;

      for (int j = 0; j < 3; j++) {
        // BMLoop *l0 = loops[loops_idx[i][(j + 2) % 3]];
        BMLoop *l = loops[loops_idx[i][j]];
        BMLoop *l2 = loops[loops_idx[i][(j + 1) % 3]];

        void **val = nullptr;
        BMEdge *e = BM_edge_exists(l->v, l2->v);

        if (e) {
          tri->eflag |= 1 << j;
          mat_tri->eflag |= 1 << j;
        }

        uintptr_t loopkey = tri_loopkey(l, mat_nr, pbvh->cd_faceset_offset, cd_uvs, totuv);

        if (!BLI_smallhash_ensure_p(&node->tribuf->vertmap, loopkey, &val)) {
          PBVHVertRef sv = {(intptr_t)l->v};

          minmax_v3v3_v3(min, max, l->v->co);

          *val = POINTER_FROM_INT(node->tribuf->totvert);
          pbvh_tribuf_add_vert(node->tribuf, sv, l);
        }

        tri->v[j] = (intptr_t)val[0];
        tri->l[j] = (intptr_t)l;

        val = nullptr;
        if (!BLI_smallhash_ensure_p(&mat_tribuf->vertmap, loopkey, &val)) {
          PBVHVertRef sv = {(intptr_t)l->v};

          minmax_v3v3_v3(min, max, l->v->co);

          *val = POINTER_FROM_INT(mat_tribuf->totvert);
          pbvh_tribuf_add_vert(mat_tribuf, sv, l);
        }

        mat_tri->v[j] = (intptr_t)val[0];
        mat_tri->l[j] = (intptr_t)l;
      }

      copy_v3_v3(tri->no, f->no);
      copy_v3_v3(mat_tri->no, f->no);
      tri->f.i = (intptr_t)f;
      mat_tri->f.i = (intptr_t)f;
    }
#else
    PBVHTri *tri = pbvh_tribuf_add_tri(node->tribuf);
    PBVHTriBuf *mat_tribuf = tribufs + mat_map[mat_nr];
    PBVHTri *mat_tri = pbvh_tribuf_add_tri(mat_tribuf);

    BMLoop *l = f->l_first;
    int j = 0;

    do {
      void **val = nullptr;

      if (!BLI_ghash_ensure_p(vmap, l->v, &val)) {
        PBVHVertRef sv = {(intptr_t)l->v};

        minmax_v3v3_v3(min, max, l->v->co);

        *val = (void *)node->tribuf->totvert;
        pbvh_tribuf_add_vert(node->tribuf, sv);
      }

      tri->v[j] = (intptr_t)val[0];
      tri->l[j] = (intptr_t)l;

      val = nullptr;
      if (!BLI_ghash_ensure_p(mat_vmaps[mat_nr], l->v, &val)) {
        PBVHVertRef sv = {(intptr_t)l->v};

        minmax_v3v3_v3(min, max, l->v->co);

        *val = (void *)mat_tribuf->totvert;
        pbvh_tribuf_add_vert(mat_tribuf, sv);
      }

      mat_tri->v[j] = (intptr_t)val[0];
      mat_tri->l[j] = (intptr_t)l;

      j++;

      if (j >= 3) {
        break;
      }

      l = l->next;
    } while (l != f->l_first);

    copy_v3_v3(tri->no, f->no);
    tri->f.i = (intptr_t)f;
#endif
  }
  TGSET_ITER_END

  TGSET_ITER (f, node->bm_faces) {
    if (pbvh_poly_hidden(pbvh, f)) {
      continue;
    }

    int mat_nr = f->mat_nr;
    PBVHTriBuf *mat_tribuf = &tribufs[mat_map[mat_nr]];

    BMLoop *l = f->l_first;
    do {
      if (l->e->head.hflag & edgeflag) {
        continue;
      }

      l->e->head.hflag |= edgeflag;

      int v1 = POINTER_AS_INT(BLI_smallhash_lookup(&node->tribuf->vertmap, (uintptr_t)l->e->v1));
      int v2 = POINTER_AS_INT(BLI_smallhash_lookup(&node->tribuf->vertmap, (uintptr_t)l->e->v2));

      pbvh_tribuf_add_edge(node->tribuf, v1, v2);

      v1 = POINTER_AS_INT(BLI_smallhash_lookup(&mat_tribuf->vertmap, (uintptr_t)l->e->v1));
      v2 = POINTER_AS_INT(BLI_smallhash_lookup(&mat_tribuf->vertmap, (uintptr_t)l->e->v2));

      pbvh_tribuf_add_edge(mat_tribuf, v1, v2);
    } while ((l = l->next) != f->l_first);
  }
  TGSET_ITER_END

  bm->elem_index_dirty |= BM_VERT;

  node->tri_buffers = c_array_from_vector<PBVHTriBuf>(tribufs);
  node->tot_tri_buffers = tribufs.size();

  if (node->tribuf->totvert) {
    copy_v3_v3(node->tribuf->min, min);
    copy_v3_v3(node->tribuf->max, max);
  }
  else {
    zero_v3(node->tribuf->min);
    zero_v3(node->tribuf->max);
  }

  return true;
}

static int pbvh_count_subtree_verts(PBVH *pbvh, PBVHNode *n)
{
  if (n->flag & PBVH_Leaf) {
    n->subtree_tottri = BLI_table_gset_len(
        n->bm_faces);  // n->tm_unique_verts->length + n->tm_other_verts->length;
    return n->subtree_tottri;
  }

  int ni = n->children_offset;

  int ret = pbvh_count_subtree_verts(pbvh, pbvh->nodes + ni);
  ret += pbvh_count_subtree_verts(pbvh, pbvh->nodes + ni + 1);

  n->subtree_tottri = ret;

  return ret;
}

void BKE_pbvh_bmesh_flag_all_disk_sort(PBVH *pbvh)
{
  BMVert *v;
  BMIter iter;

  BM_ITER_MESH (v, &iter, pbvh->header.bm, BM_VERTS_OF_MESH) {
    MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, v);
    mv->flag |= SCULPTVERT_NEED_DISK_SORT;
  }
}

void BKE_pbvh_bmesh_update_all_valence(PBVH *pbvh)
{
  BMIter iter;
  BMVert *v;

  BM_ITER_MESH (v, &iter, pbvh->header.bm, BM_VERTS_OF_MESH) {
    BKE_pbvh_bmesh_update_valence(pbvh->cd_sculpt_vert, BKE_pbvh_make_vref((intptr_t)v));
  }
}

void BKE_pbvh_bmesh_on_mesh_change(PBVH *pbvh)
{
  BMIter iter;
  BMVert *v;

  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (node->flag & PBVH_Leaf) {
      node->flag |= PBVH_UpdateTriAreas;
    }
  }

  const int cd_sculpt_vert = pbvh->cd_sculpt_vert;

  BM_ITER_MESH (v, &iter, pbvh->header.bm, BM_VERTS_OF_MESH) {
    MSculptVert *mv = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, v);
    int *flags = BM_ELEM_CD_PTR<int *>(v, pbvh->cd_boundary_flag);
    *flags |= SCULPT_BOUNDARY_NEEDS_UPDATE;

    MV_ADD_FLAG(mv, SCULPTVERT_NEED_DISK_SORT | SCULPTVERT_NEED_TRIANGULATE);
    BKE_pbvh_bmesh_update_valence(pbvh->cd_sculpt_vert, BKE_pbvh_make_vref((intptr_t)v));
  }
}

bool BKE_pbvh_bmesh_mark_update_valence(PBVH *pbvh, PBVHVertRef vertex)
{
  BMVert *v = (BMVert *)vertex.i;
  MSculptVert *mv = BM_ELEM_CD_PTR<MSculptVert *>(v, pbvh->cd_sculpt_vert);

  bool ret = mv->flag & SCULPTVERT_NEED_VALENCE;

  mv->flag |= SCULPTVERT_NEED_VALENCE;

  return ret;
}

bool BKE_pbvh_bmesh_check_valence(PBVH *pbvh, PBVHVertRef vertex)
{
  BMVert *v = (BMVert *)vertex.i;
  MSculptVert *mv = BM_ELEM_CD_PTR<MSculptVert *>(v, pbvh->cd_sculpt_vert);

  if (mv->flag & SCULPTVERT_NEED_VALENCE) {
    BKE_pbvh_bmesh_update_valence(pbvh->cd_sculpt_vert, vertex);
    return true;
  }

  return false;
}

void BKE_pbvh_bmesh_update_valence(int cd_sculpt_vert, PBVHVertRef vertex)
{
  BMVert *v = (BMVert *)vertex.i;
  BMEdge *e;

  MSculptVert *mv = BM_ELEM_CD_PTR<MSculptVert *>(v, cd_sculpt_vert);

  mv->flag &= ~SCULPTVERT_NEED_VALENCE;

  if (!v->e) {
    mv->valence = 0;
    return;
  }

  mv->valence = 0;

  e = v->e;

  if (!e) {
    return;
  }

  do {
    mv->valence++;

    e = v == e->v1 ? e->v1_disk_link.next : e->v2_disk_link.next;

    if (!e) {
      printf("bmesh error!\n");
      break;
    }
  } while (e != v->e);
}

static void pbvh_bmesh_join_subnodes(PBVH *pbvh, PBVHNode *node, PBVHNode *parent)
{
  if (!(node->flag & PBVH_Leaf)) {
    int ni = node->children_offset;

    if (ni > 0 && ni < pbvh->totnode - 1) {
      pbvh_bmesh_join_subnodes(pbvh, pbvh->nodes + ni, parent);
      pbvh_bmesh_join_subnodes(pbvh, pbvh->nodes + ni + 1, parent);
    }
    else {
      printf("node corruption: %d\n", ni);
      return;
    }
    if (node != parent) {
      node->flag |= PBVH_Delete; /* Mark for deletion. */
    }

    return;
  }

  if (node != parent) {
    node->flag |= PBVH_Delete; /* Mark for deletion. */
  }

  BMVert *v;

  TGSET_ITER (v, node->bm_unique_verts) {
    BLI_table_gset_add(parent->bm_unique_verts, v);

    int *flags = BM_ELEM_CD_PTR<int *>(v, pbvh->cd_boundary_flag);
    *flags |= SCULPT_BOUNDARY_NEEDS_UPDATE;

    BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, DYNTOPO_NODE_NONE);
  }
  TGSET_ITER_END

  // printf("  subtotface: %d\n", BLI_table_gset_len(node->bm_faces));

  BMFace *f;
  TGSET_ITER (f, node->bm_faces) {
    BLI_table_gset_add(parent->bm_faces, f);
    BM_ELEM_CD_SET_INT(f, pbvh->cd_face_node_offset, DYNTOPO_NODE_NONE);
  }
  TGSET_ITER_END
}

static void BKE_pbvh_bmesh_correct_tree(PBVH *pbvh, PBVHNode *node, PBVHNode *parent)
{
  const int size_lower = pbvh->leaf_limit - (pbvh->leaf_limit >> 1);

  if (node->flag & PBVH_Leaf) {
    return;
  }

  if (node->subtree_tottri < size_lower && node != pbvh->nodes) {
    node->bm_unique_verts = BLI_table_gset_new("bm_unique_verts");
    node->bm_other_verts = BLI_table_gset_new("bm_other_verts");
    node->bm_faces = BLI_table_gset_new("bm_faces");

    pbvh_bmesh_join_subnodes(pbvh, pbvh->nodes + node->children_offset, node);
    pbvh_bmesh_join_subnodes(pbvh, pbvh->nodes + node->children_offset + 1, node);

    node->children_offset = 0;
    node->flag |= PBVH_Leaf | PBVH_UpdateRedraw | PBVH_UpdateBB | PBVH_UpdateDrawBuffers |
                  PBVH_RebuildDrawBuffers | PBVH_UpdateOriginalBB | PBVH_UpdateMask |
                  PBVH_UpdateVisibility | PBVH_UpdateColor | PBVH_UpdateNormals | PBVH_UpdateTris;

    TableGSet *other = BLI_table_gset_new(__func__);
    BMVert *v;

    node->children_offset = 0;
    node->draw_batches = nullptr;

    // rebuild bm_other_verts
    BMFace *f;
    TGSET_ITER (f, node->bm_faces) {
      BMLoop *l = f->l_first;

      BM_ELEM_CD_SET_INT(f, pbvh->cd_face_node_offset, DYNTOPO_NODE_NONE);

      do {
        if (!BLI_table_gset_haskey(node->bm_unique_verts, l->v)) {
          BLI_table_gset_add(other, l->v);
        }
        l = l->next;
      } while (l != f->l_first);
    }
    TGSET_ITER_END

    BLI_table_gset_free(node->bm_other_verts, nullptr);
    node->bm_other_verts = other;

    BB_reset(&node->vb);

#if 1
    TGSET_ITER (v, node->bm_unique_verts) {
      BB_expand(&node->vb, v->co);
    }
    TGSET_ITER_END

    TGSET_ITER (v, node->bm_other_verts) {
      BB_expand(&node->vb, v->co);
    }
    TGSET_ITER_END
#endif

    // printf("totface: %d\n", BLI_table_gset_len(node->bm_faces));
    node->orig_vb = node->vb;

    return;
  }

  int ni = node->children_offset;

  for (int i = 0; i < 2; i++, ni++) {
    PBVHNode *child = pbvh->nodes + ni;
    BKE_pbvh_bmesh_correct_tree(pbvh, child, node);
  }
}

/* Deletes PBVH_Delete marked nodes. */
static void pbvh_bmesh_compact_tree(PBVH *bvh)
{
  // compact nodes
  int totnode = 0;
  for (int i = 0; i < bvh->totnode; i++) {
    PBVHNode *n = bvh->nodes + i;

    if (!(n->flag & PBVH_Delete)) {
      if (!(n->flag & PBVH_Leaf)) {
        PBVHNode *n1 = bvh->nodes + n->children_offset;
        PBVHNode *n2 = bvh->nodes + n->children_offset + 1;

        if ((n1->flag & PBVH_Delete) != (n2->flag & PBVH_Delete)) {
          printf("un-deleting an empty node\n");
          PBVHNode *n3 = n1->flag & PBVH_Delete ? n1 : n2;

          n3->flag = PBVH_Leaf | PBVH_UpdateTris;
          n3->bm_unique_verts = BLI_table_gset_new("bm_unique_verts");
          n3->bm_other_verts = BLI_table_gset_new("bm_other_verts");
          n3->bm_faces = BLI_table_gset_new("bm_faces");
          n3->tribuf = nullptr;
          n3->draw_batches = nullptr;
        }
        else if ((n1->flag & PBVH_Delete) && (n2->flag & PBVH_Delete)) {
          n->children_offset = 0;
          n->flag |= PBVH_Leaf | PBVH_UpdateTris;

          if (!n->bm_unique_verts) {
            // should not happen
            n->bm_unique_verts = BLI_table_gset_new("bm_unique_verts");
            n->bm_other_verts = BLI_table_gset_new("bm_other_verts");
            n->bm_faces = BLI_table_gset_new("bm_faces");
            n->tribuf = nullptr;
            n->draw_batches = nullptr;
          }
        }
      }

      totnode++;
    }
  }

  int *map = MEM_cnew_array<int>(bvh->totnode, "bmesh map temp");

  // build idx map for child offsets
  int j = 0;
  for (int i = 0; i < bvh->totnode; i++) {
    PBVHNode *n = bvh->nodes + i;

    if (!(n->flag & PBVH_Delete)) {
      map[i] = j++;
    }
    else if (1) {
      if (n->layer_disp) {
        MEM_freeN(n->layer_disp);
        n->layer_disp = nullptr;
      }

      pbvh_free_draw_buffers(bvh, n);

      if (n->vert_indices) {
        MEM_freeN((void *)n->vert_indices);
        n->vert_indices = nullptr;
      }
      if (n->face_vert_indices) {
        MEM_freeN((void *)n->face_vert_indices);
        n->face_vert_indices = nullptr;
      }

      if (n->tribuf || n->tri_buffers) {
        BKE_pbvh_bmesh_free_tris(bvh, n);
      }

      if (n->bm_unique_verts) {
        BLI_table_gset_free(n->bm_unique_verts, nullptr);
        n->bm_unique_verts = nullptr;
      }

      if (n->bm_other_verts) {
        BLI_table_gset_free(n->bm_other_verts, nullptr);
        n->bm_other_verts = nullptr;
      }

      if (n->bm_faces) {
        BLI_table_gset_free(n->bm_faces, nullptr);
        n->bm_faces = nullptr;
      }

#ifdef PROXY_ADVANCED
      BKE_pbvh_free_proxyarray(bvh, n);
#endif
    }
  }

  // compact node array
  j = 0;
  for (int i = 0; i < bvh->totnode; i++) {
    if (!(bvh->nodes[i].flag & PBVH_Delete)) {
      if (bvh->nodes[i].children_offset >= bvh->totnode - 1) {
        printf("error %i %i\n", i, bvh->nodes[i].children_offset);
        continue;
      }

      int i1 = map[bvh->nodes[i].children_offset];
      int i2 = map[bvh->nodes[i].children_offset + 1];

      if (bvh->nodes[i].children_offset >= bvh->totnode) {
        printf("bad child node reference %d->%d, totnode: %d\n",
               i,
               bvh->nodes[i].children_offset,
               bvh->totnode);
        continue;
      }

      if (bvh->nodes[i].children_offset && i2 != i1 + 1) {
        printf("      pbvh corruption during node join %d %d\n", i1, i2);
      }

      bvh->nodes[j] = bvh->nodes[i];
      bvh->nodes[j].children_offset = i1;

      j++;
    }
  }

  if (j != totnode) {
    printf("pbvh error: %s", __func__);
  }

  if (bvh->totnode != j) {
    memset(bvh->nodes + j, 0, sizeof(*bvh->nodes) * (bvh->totnode - j));
    bvh->node_mem_count = j;
  }

  bvh->totnode = j;

  // set vert/face node indices again
  for (int i = 0; i < bvh->totnode; i++) {
    PBVHNode *n = bvh->nodes + i;

    if (!(n->flag & PBVH_Leaf)) {
      continue;
    }

    if (!n->bm_unique_verts) {
      printf("ERROR!\n");
      n->bm_unique_verts = BLI_table_gset_new("bleh");
      n->bm_other_verts = BLI_table_gset_new("bleh");
      n->bm_faces = BLI_table_gset_new("bleh");
    }

    BMVert *v;

    TGSET_ITER (v, n->bm_unique_verts) {
      BM_ELEM_CD_SET_INT(v, bvh->cd_vert_node_offset, i);
    }
    TGSET_ITER_END

    BMFace *f;

    TGSET_ITER (f, n->bm_faces) {
      BM_ELEM_CD_SET_INT(f, bvh->cd_face_node_offset, i);
    }
    TGSET_ITER_END
  }

  Vector<BMVert *> scratch;

  for (int i = 0; i < bvh->totnode; i++) {
    PBVHNode *n = bvh->nodes + i;

    if (!(n->flag & PBVH_Leaf)) {
      continue;
    }

    scratch.clear();
    BMVert *v;

    TGSET_ITER (v, n->bm_other_verts) {
      int ni = BM_ELEM_CD_GET_INT(v, bvh->cd_vert_node_offset);
      if (ni == DYNTOPO_NODE_NONE) {
        scratch.append(v);
      }
      // BM_ELEM_CD_SET_INT(v, bvh->cd_vert_node_offset, i);
    }
    TGSET_ITER_END

    int slen = scratch.size();
    for (int j = 0; j < slen; j++) {
      BMVert *v = scratch[j];

      BLI_table_gset_remove(n->bm_other_verts, v, nullptr);
      BLI_table_gset_add(n->bm_unique_verts, v);
      BM_ELEM_CD_SET_INT(v, bvh->cd_vert_node_offset, i);
    }
  }

  MEM_freeN(map);
}

static void recursive_delete_nodes(PBVH *pbvh, int ni)
{
  PBVHNode *node = pbvh->nodes + ni;

  node->flag |= PBVH_Delete;

  if (!(node->flag & PBVH_Leaf) && node->children_offset) {
    if (node->children_offset < pbvh->totnode) {
      recursive_delete_nodes(pbvh, node->children_offset);
    }

    if (node->children_offset + 1 < pbvh->totnode) {
      recursive_delete_nodes(pbvh, node->children_offset + 1);
    }
  }
}

/* Prunes leaf nodes that are too small or degenerate. */
static void pbvh_bmesh_balance_tree(PBVH *pbvh)
{
  float *overlaps = MEM_cnew_array<float>(pbvh->totnode, "overlaps");
  PBVHNode **parentmap = MEM_cnew_array<PBVHNode *>(pbvh->totnode, "parentmap");
  int *depthmap = MEM_cnew_array<int>(pbvh->totnode, "depthmap");

  Vector<PBVHNode *> stack;
  Vector<BMFace *> faces;
  Vector<PBVHNode *> substack;

  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if ((node->flag & PBVH_Leaf) || node->children_offset == 0) {
      continue;
    }

    if (node->children_offset < pbvh->totnode) {
      parentmap[node->children_offset] = node;
    }

    if (node->children_offset + 1 < pbvh->totnode) {
      parentmap[node->children_offset + 1] = node;
    }
  }

  const int cd_vert_node = pbvh->cd_vert_node_offset;
  const int cd_face_node = pbvh->cd_face_node_offset;

  bool modified = false;

  stack.append(pbvh->nodes);

  while (stack.size() > 0) {
    PBVHNode *node = stack.pop_last();
    BB clip;

    if (!(node->flag & PBVH_Leaf) && node->children_offset > 0) {
      PBVHNode *child1 = pbvh->nodes + node->children_offset;
      PBVHNode *child2 = pbvh->nodes + node->children_offset + 1;

      float volume = BB_volume(&child1->vb) + BB_volume(&child2->vb);

      /* dissolve nodes whose children overlap by more then a percentage
        of the total volume.  we use a simple huerstic to calculate the
        cutoff threshold.*/

      BB_intersect(&clip, &child1->vb, &child2->vb);
      float overlap = BB_volume(&clip);
      float factor;

      /* use higher threshold for the root node and its immediate children */
      switch (stack.size()) {
        case 0:
          factor = 0.5;
          break;
        case 1:
        case 2:
          factor = 0.2;
          break;
        default:
          factor = 0.2;
          break;
      }

#if 0
      for (int k = 0; k < stack.size(); k++) {
        printf(" ");
      }

      printf("factor: %.3f\n", factor);
#endif

      bool bad = overlap > volume * factor;

      bad |= child1->bm_faces && !BLI_table_gset_len(child1->bm_faces);
      bad |= child2->bm_faces && !BLI_table_gset_len(child2->bm_faces);

      if (bad) {
        modified = true;

        substack.clear();
        substack.append(child1);
        substack.append(child2);

        while (substack.size() > 0) {
          PBVHNode *node2 = substack.pop_last();

          node2->flag |= PBVH_Delete;

          if (node2->flag & PBVH_Leaf) {
            BMFace *f;
            BMVert *v;

            TGSET_ITER (f, node2->bm_faces) {
              if (BM_ELEM_CD_GET_INT(f, cd_face_node) == -1) {
                // eek!
                continue;
              }

              BM_ELEM_CD_SET_INT(f, cd_face_node, DYNTOPO_NODE_NONE);
              faces.append(f);
            }
            TGSET_ITER_END;

            TGSET_ITER (v, node2->bm_unique_verts) {
              int *flags = BM_ELEM_CD_PTR<int *>(v, pbvh->cd_boundary_flag);
              *flags |= SCULPT_BOUNDARY_NEEDS_UPDATE;

              BM_ELEM_CD_SET_INT(v, cd_vert_node, DYNTOPO_NODE_NONE);
            }
            TGSET_ITER_END;
          }
          else if (node2->children_offset > 0 && node2->children_offset < pbvh->totnode) {
            substack.append(pbvh->nodes + node2->children_offset);

            if (node2->children_offset + 1 < pbvh->totnode) {
              substack.append(pbvh->nodes + node2->children_offset + 1);
            }
          }
        }
      }

      if (node->children_offset < pbvh->totnode) {
        stack.append(child1);
      }

      if (node->children_offset + 1 < pbvh->totnode) {
        stack.append(child2);
      }
    }
  }

  if (modified) {
    pbvh_bmesh_compact_tree(pbvh);

    printf("joined nodes; %d faces\n", (int)faces.size());

    for (int i : IndexRange(faces.size())) {
      if (BM_elem_is_free((BMElem *)faces[i], BM_FACE)) {
        printf("corrupted face in pbvh tree; faces[i]: %p\n", faces[i]);
        continue;
      }

      if (BM_ELEM_CD_GET_INT(faces[i], cd_face_node) != DYNTOPO_NODE_NONE) {
        // printf("duplicate faces in pbvh_bmesh_balance_tree!\n");
        continue;
      }

      bke_pbvh_insert_face(pbvh, faces[i]);
    }
  }

  MEM_SAFE_FREE(parentmap);
  MEM_SAFE_FREE(overlaps);
  MEM_SAFE_FREE(depthmap);
}

/* Fix any orphaned empty leaves that survived other stages of culling.*/
ATTR_NO_OPT static void pbvh_fix_orphan_leaves(PBVH *pbvh)
{
  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (!(node->flag & PBVH_Leaf) || (node->flag & PBVH_Delete) ||
        BLI_table_gset_len(node->bm_faces) != 0) {
      continue;
    }

    /* Find parent node. */
    PBVHNode *parent = nullptr;

    for (int j = 0; j < pbvh->totnode; j++) {
      if (pbvh->nodes[j].children_offset == i || pbvh->nodes[j].children_offset + 1 == i) {
        parent = pbvh->nodes + j;
        break;
      }
    }

    if (!parent) {
      printf("%s: node corruption, could not find parent node!\n", __func__);
      continue;
    }

    PBVHNode *other = pbvh->nodes + (parent->children_offset == i ? i + 1 : i - 1);

    if (other->flag & PBVH_Delete) {
      printf("%s: error!\n", __func__);
      continue;
    }

    while (!(other->flag & PBVH_Leaf)) {
      PBVHNode *a = pbvh->nodes + other->children_offset;
      PBVHNode *b = pbvh->nodes + other->children_offset + 1;

      if (!(a->flag & PBVH_Delete) && (a->flag & PBVH_Leaf) &&
          BLI_table_gset_len(a->bm_faces) > 1) {
        other = a;
      }
      else if (!(b->flag & PBVH_Delete) && (b->flag & PBVH_Leaf) &&
               BLI_table_gset_len(b->bm_faces) > 1) {
        other = b;
      }
      else {
        other = nullptr;
        break;
      }
    }

    if (other == nullptr || BLI_table_gset_len(other->bm_faces) < 1) {
      printf("%s: other was nullptr\n", __func__);
      continue;
    }

    /* Steal a single face from other */
    BMFace *f;
    PBVHNodeFlags updateflag = PBVH_UpdateOtherVerts | PBVH_UpdateBB | PBVH_UpdateOriginalBB |
                               PBVH_UpdateTris | PBVH_UpdateTriAreas | PBVH_RebuildDrawBuffers |
                               PBVH_RebuildNodeVerts | PBVH_RebuildPixels | PBVH_UpdateNormals |
                               PBVH_UpdateCurvatureDir | PBVH_UpdateRedraw | PBVH_UpdateVisibility;

    TGSET_ITER (f, other->bm_faces) {
      BLI_table_gset_remove(other->bm_faces, static_cast<void *>(f), nullptr);
      BLI_table_gset_add(node->bm_faces, static_cast<void *>(f));
      BM_ELEM_CD_SET_INT(f, pbvh->cd_face_node_offset, i);

      BMVert *v = f->l_first->v;
      int node_i = BM_ELEM_CD_GET_INT(v, pbvh->cd_vert_node_offset);
      if (node_i != DYNTOPO_NODE_NONE) {
        PBVHNode *node = pbvh->nodes + node_i;
        if (BLI_table_gset_haskey(node->bm_unique_verts, static_cast<void *>(v))) {
          BLI_table_gset_remove(node->bm_unique_verts, static_cast<void *>(v), nullptr);
        }
      }

      BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, i);
      BLI_table_gset_add(node->bm_unique_verts, v);

      node->flag |= updateflag;
      other->flag |= updateflag;

      printf("%s: Patched empty leaf node.\n", __func__);
      break;
    }
    TGSET_ITER_END;
  }
}

static void pbvh_bmesh_join_nodes(PBVH *pbvh)
{
  if (pbvh->totnode < 2) {
    return;
  }

  pbvh_count_subtree_verts(pbvh, pbvh->nodes);
  BKE_pbvh_bmesh_correct_tree(pbvh, pbvh->nodes, nullptr);
  pbvh_fix_orphan_leaves(pbvh);

  /* Compact nodes. */
  int totnode = 0;
  int *map = MEM_cnew_array<int>(pbvh->totnode, "bmesh map temp");

  for (int i = 0; i < pbvh->totnode; i++) {
    for (int j = 0; j < pbvh->totnode; j++) {
      if (i == j || !pbvh->nodes[i].draw_batches) {
        continue;
      }

      if (pbvh->nodes[i].draw_batches == pbvh->nodes[j].draw_batches) {
        printf("%s: error %d %d\n", __func__, i, j);

        pbvh->nodes[j].draw_batches = nullptr;
      }
    }
  }

  /* Build index map for child offsets. */
  int j = 0;
  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *n = pbvh->nodes + i;

    if (!(n->flag & PBVH_Delete)) {
      map[i] = j++;
      totnode++;
    }
    else {
      if (n->layer_disp) {
        MEM_freeN(n->layer_disp);
        n->layer_disp = nullptr;
      }

      pbvh_free_draw_buffers(pbvh, n);

      if (n->vert_indices) {
        MEM_freeN((void *)n->vert_indices);
        n->vert_indices = nullptr;
      }
      if (n->face_vert_indices) {
        MEM_freeN((void *)n->face_vert_indices);
        n->face_vert_indices = nullptr;
      }

      if (n->tribuf || n->tri_buffers) {
        BKE_pbvh_bmesh_free_tris(pbvh, n);
      }

      if (n->bm_unique_verts) {
        BLI_table_gset_free(n->bm_unique_verts, nullptr);
        n->bm_unique_verts = nullptr;
      }

      if (n->bm_other_verts) {
        BLI_table_gset_free(n->bm_other_verts, nullptr);
        n->bm_other_verts = nullptr;
      }

      if (n->bm_faces) {
        BLI_table_gset_free(n->bm_faces, nullptr);
        n->bm_faces = nullptr;
      }

#ifdef PROXY_ADVANCED
      BKE_pbvh_free_proxyarray(pbvh, n);
#endif
    }
  }

  // compact node array
  j = 0;
  for (int i = 0; i < pbvh->totnode; i++) {
    if (!(pbvh->nodes[i].flag & PBVH_Delete)) {
      if (pbvh->nodes[i].children_offset >= pbvh->totnode - 1) {
        printf("%s: error %i %i\n", __func__, i, pbvh->nodes[i].children_offset);
        continue;
      }

      int i1 = map[pbvh->nodes[i].children_offset];
      int i2 = map[pbvh->nodes[i].children_offset + 1];

      if (pbvh->nodes[i].children_offset >= pbvh->totnode) {
        printf("%s: Bad child node reference %d->%d, totnode: %d\n",
               __func__,
               i,
               pbvh->nodes[i].children_offset,
               pbvh->totnode);
        continue;
      }

      if (pbvh->nodes[i].children_offset && i2 != i1 + 1) {
        printf("      pbvh corruption during node join %d %d\n", i1, i2);
      }

      pbvh->nodes[j] = pbvh->nodes[i];
      pbvh->nodes[j].children_offset = i1;

      j++;
    }
  }

  if (j != totnode) {
    printf("%s: pbvh error.\n", __func__);
  }

  if (pbvh->totnode != j) {
    memset(pbvh->nodes + j, 0, sizeof(*pbvh->nodes) * (pbvh->totnode - j));
    pbvh->node_mem_count = j;
  }

  pbvh->totnode = j;

  // set vert/face node indices again
  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *n = pbvh->nodes + i;

    if (!(n->flag & PBVH_Leaf)) {
      continue;
    }

    if (!n->bm_unique_verts) {
      printf("%s: ERROR!\n", __func__);
      n->bm_unique_verts = BLI_table_gset_new("bleh");
      n->bm_other_verts = BLI_table_gset_new("bleh");
      n->bm_faces = BLI_table_gset_new("bleh");
    }

    BMVert *v;

    TGSET_ITER (v, n->bm_unique_verts) {
      BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, i);
    }
    TGSET_ITER_END

    BMFace *f;

    TGSET_ITER (f, n->bm_faces) {
      BM_ELEM_CD_SET_INT(f, pbvh->cd_face_node_offset, i);
    }
    TGSET_ITER_END
  }

  Vector<BMVert *> scratch;

  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *n = pbvh->nodes + i;

    if (!(n->flag & PBVH_Leaf)) {
      continue;
    }

    scratch.clear();
    BMVert *v;

    TGSET_ITER (v, n->bm_other_verts) {
      int ni = BM_ELEM_CD_GET_INT(v, pbvh->cd_vert_node_offset);
      if (ni == DYNTOPO_NODE_NONE) {
        scratch.append(v);
      }
      // BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, i);
    }
    TGSET_ITER_END

    for (int j : IndexRange(scratch.size())) {
      BMVert *v = scratch[j];

      BLI_table_gset_remove(n->bm_other_verts, v, nullptr);
      BLI_table_gset_add(n->bm_unique_verts, v);
      BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, i);
    }
  }

  MEM_freeN(map);
}

void BKE_pbvh_bmesh_after_stroke(PBVH *pbvh, bool force_balance)
{
  int totnode = pbvh->totnode;

  BKE_pbvh_update_bounds(pbvh, (PBVH_UpdateBB | PBVH_UpdateOriginalBB | PBVH_UpdateRedraw));

  pbvh_bmesh_check_nodes(pbvh);
  pbvh_bmesh_join_nodes(pbvh);
  pbvh_bmesh_check_nodes(pbvh);

  BKE_pbvh_update_bounds(pbvh, (PBVH_UpdateBB | PBVH_UpdateOriginalBB | PBVH_UpdateRedraw));

  if (force_balance || pbvh->balance_counter++ == 10) {
    pbvh_bmesh_balance_tree(pbvh);
    pbvh_bmesh_check_nodes(pbvh);
    pbvh->balance_counter = 0;

    totnode = pbvh->totnode;

    for (int i = 0; i < totnode; i++) {
      PBVHNode *n = pbvh->nodes + i;

      if (totnode != pbvh->totnode) {
#ifdef PROXY_ADVANCED
        BKE_pbvh_free_proxyarray(pbvh, n);
#endif
      }

      if (n->flag & PBVH_Leaf) {
        /* Recursively split nodes that have gotten too many
         * elements */
        pbvh_bmesh_node_limit_ensure(pbvh, i);
      }
    }
  }

  pbvh_print_mem_size(pbvh);
}

void BKE_pbvh_bmesh_detail_size_set(PBVH *pbvh, float detail_size, float detail_range)
{
  pbvh->bm_max_edge_len = detail_size;
  pbvh->bm_min_edge_len = pbvh->bm_max_edge_len * detail_range;
  pbvh->bm_detail_range = detail_range;
}

void BKE_pbvh_node_mark_topology_update(PBVHNode *node)
{
  node->flag |= PBVH_UpdateTopology;
}

TableGSet *BKE_pbvh_bmesh_node_unique_verts(PBVHNode *node)
{
  return node->bm_unique_verts;
}

TableGSet *BKE_pbvh_bmesh_node_other_verts(PBVHNode *node)
{
  pbvh_bmesh_check_other_verts(node);
  return node->bm_other_verts;
}

struct TableGSet *BKE_pbvh_bmesh_node_faces(PBVHNode *node)
{
  return node->bm_faces;
}

/****************************** Debugging *****************************/

void BKE_pbvh_update_offsets(PBVH *pbvh,
                             const int cd_vert_node_offset,
                             const int cd_face_node_offset,
                             const int cd_sculpt_vert,
                             const int cd_face_areas,
                             const int cd_boundary_flag)
{
  pbvh->cd_face_node_offset = cd_face_node_offset;
  pbvh->cd_vert_node_offset = cd_vert_node_offset;
  pbvh->cd_face_area = cd_face_areas;
  pbvh->cd_vert_mask_offset = CustomData_get_offset(&pbvh->header.bm->vdata, CD_PAINT_MASK);
  pbvh->cd_sculpt_vert = cd_sculpt_vert;
  pbvh->cd_faceset_offset = CustomData_get_offset_named(
      &pbvh->header.bm->pdata, CD_PROP_INT32, ".sculpt_face_set");

  pbvh->totuv = CustomData_number_of_layers(&pbvh->header.bm->ldata, CD_PROP_FLOAT2);
  pbvh->cd_boundary_flag = cd_boundary_flag;

  if (pbvh->bm_idmap) {
    BM_idmap_check_attributes(pbvh->bm_idmap);
  }
}

static void scan_edge_split(BMesh *bm, BMEdge **edges, int totedge)
{
  Vector<BMFace *> faces;
  Vector<BMEdge *> newedges;
  Vector<BMVert *> newverts;
  Vector<BMVert *> fmap;
  Vector<int> emap;

  // remove e from radial list of e->v2
  for (int i = 0; i < totedge; i++) {
    BMEdge *e = edges[i];

    BMDiskLink *prev;
    BMDiskLink *next;

    if (e->v2_disk_link.prev->v1 == e->v2) {
      prev = &e->v2_disk_link.prev->v1_disk_link;
    }
    else {
      prev = &e->v2_disk_link.prev->v2_disk_link;
    }

    if (e->v2_disk_link.next->v1 == e->v2) {
      next = &e->v2_disk_link.next->v1_disk_link;
    }
    else {
      next = &e->v2_disk_link.next->v2_disk_link;
    }

    prev->next = e->v2_disk_link.next;
    next->prev = e->v2_disk_link.prev;
  }

  for (int i = 0; i < totedge; i++) {
    BMEdge *e = edges[i];

    BMVert *v2 = static_cast<BMVert *>(BLI_mempool_alloc(bm->vpool));
    memset(v2, 0, sizeof(*v2));
    v2->head.data = BLI_mempool_alloc(bm->vdata.pool);

    newverts.append(v2);

    BMEdge *e2 = static_cast<BMEdge *>(BLI_mempool_alloc(bm->epool));
    newedges.append(e2);

    memset(e2, 0, sizeof(*e2));
    if (bm->edata.pool) {
      e2->head.data = BLI_mempool_alloc(bm->edata.pool);
    }

    BMLoop *l = e->l;

    if (!l) {
      continue;
    }

    do {
      faces.append(l->f);
      BMFace *f2 = static_cast<BMFace *>(BLI_mempool_alloc(bm->fpool));

      faces.append(l->f);
      fmap.append(v2);
      emap.append(i);

      faces.append(f2);
      fmap.append(v2);
      emap.append(i);

      memset(f2, 0, sizeof(*f2));
      f2->head.data = BLI_mempool_alloc(bm->ldata.pool);

      BMLoop *prev = nullptr;
      BMLoop *l2 = nullptr;

      for (int j = 0; j < 3; j++) {
        l2 = (BMLoop *)BLI_mempool_alloc(bm->lpool);
        memset(l2, 0, sizeof(*l2));
        l2->head.data = BLI_mempool_alloc(bm->ldata.pool);

        l2->prev = prev;

        if (prev) {
          prev->next = l2;
        }
        else {
          f2->l_first = l2;
        }
      }

      f2->l_first->prev = l2;
      l2->next = f2->l_first;

      faces.append(f2);
      l = l->radial_next;
    } while (l != e->l);
  }

  for (int i : IndexRange(newedges.size())) {
    BMEdge *e1 = edges[i];
    BMEdge *e2 = newedges[i];
    BMVert *v = newverts[i];

    add_v3_v3v3(v->co, e1->v1->co, e1->v2->co);
    mul_v3_fl(v->co, 0.5f);

    e2->v1 = v;
    e2->v2 = e1->v2;
    e1->v2 = v;

    v->e = e1;

    e1->v2_disk_link.next = e1->v2_disk_link.prev = e2;
    e2->v1_disk_link.next = e2->v1_disk_link.prev = e1;
  }

  for (int i : IndexRange(faces.size())) {
    BMFace *f1 = faces[i], *f2 = faces[i + 1];
    BMEdge *e1 = edges[emap[i]];
    BMEdge *e2 = newedges[emap[i]];
    BMVert *nv = fmap[i];

    // make sure first loop points to e1->v1
    BMLoop *l = f1->l_first;
    do {
      if (l->v == e1->v1) {
        break;
      }
      l = l->next;
    } while (l != f1->l_first);

    f1->l_first = l;

    BMLoop *l2 = f2->l_first;

    l2->f = l2->next->f = l2->prev->f = f2;
    l2->v = nv;
    l2->next->v = l->next->v;
    l2->prev->v = l->prev->v;
    l2->e = e2;
    l2->next->e = l->next->e;
    l2->prev->e = l->prev->e;

    l->next->v = nv;
    l->next->e = e2;
  }
}

#define MAX_RE_CHILD 3
struct ReVertNode {
  int totvert, totchild;
  struct ReVertNode *parent;
  struct ReVertNode *children[MAX_RE_CHILD];
  BMVert *verts[];
};

BMesh *BKE_pbvh_reorder_bmesh(PBVH *pbvh)
{
  /*try to compute size of verts per node*/
  int vsize = sizeof(BMVert);
  vsize += pbvh->header.bm->vdata.totsize;

  // perhaps aim for l2 cache?
  const int limit = 1024;
  int leaf_limit = MAX2(limit / vsize, 4);

  BLI_mempool *pool = BLI_mempool_create(sizeof(ReVertNode) + sizeof(void *) * vsize, 0, 8192, 0);
  ReVertNode **vnodemap = (ReVertNode **)MEM_calloc_arrayN(
      pbvh->header.bm->totvert, sizeof(void *), "vnodemap");

  printf("leaf_limit: %d\n", leaf_limit);

  BMIter iter;
  BMVert *v;
  const char flag = BM_ELEM_TAG_ALT;
  int i = 0;

  BM_ITER_MESH (v, &iter, pbvh->header.bm, BM_VERTS_OF_MESH) {
    v->head.hflag &= ~flag;
    v->head.index = i++;
  }

  Vector<BMVert *> stack;

  BM_ITER_MESH (v, &iter, pbvh->header.bm, BM_VERTS_OF_MESH) {
    if (v->head.hflag & flag) {
      continue;
    }

    ReVertNode *node = (ReVertNode *)BLI_mempool_calloc(pool);

    stack.clear();
    stack.append(v);

    v->head.hflag |= flag;

    vnodemap[v->head.index] = node;
    node->verts[node->totvert++] = v;

    while (stack.size() > 0) {
      BMVert *v2 = stack.pop_last();
      BMEdge *e;

      if (node->totvert >= leaf_limit) {
        break;
      }

      if (!v2->e) {
        continue;
      }

      int len = node->totvert;

      e = v2->e;
      do {
        BMVert *v3 = BM_edge_other_vert(e, v2);

        if (!BM_elem_flag_test(v3, flag) && len < leaf_limit) {
          v3->head.hflag |= flag;

          vnodemap[v3->head.index] = node;
          node->verts[node->totvert++] = v3;

          len++;

          stack.append(v3);
        }

        e = e->v1 == v2 ? e->v1_disk_link.next : e->v2_disk_link.next;
      } while (e != v2->e);
    }
  }

  const int steps = 4;
  Vector<ReVertNode *> roots;

  for (int step = 0; step < steps; step++) {
    const bool last_step = step == steps - 1;

    BM_ITER_MESH_INDEX (v, &iter, pbvh->header.bm, BM_VERTS_OF_MESH, i) {
      BMEdge *e = v->e;

      if (!e) {
        continue;
      }

      ReVertNode *node = vnodemap[v->head.index];
      if (node->parent) {
        continue;
      }

      ReVertNode *parent = (ReVertNode *)BLI_mempool_calloc(pool);
      parent->children[0] = node;
      parent->totchild = 1;

      do {
        BMVert *v2 = BM_edge_other_vert(e, v);

        ReVertNode *node2 = vnodemap[v2->head.index];

        bool ok = node != node2 && !node2->parent;
        ok = ok && parent->totchild < MAX_RE_CHILD;

        for (int j = 0; j < parent->totchild; j++) {
          if (parent->children[j] == node2) {
            ok = false;
            break;
          }
        }

        if (ok) {
          parent->children[parent->totchild++] = node2;
          node2->parent = parent;
          break;
        }

        e = e->v1 == v ? e->v1_disk_link.next : e->v2_disk_link.next;
      } while (e != v->e);

      if (last_step) {
        roots.append(parent);
      }

      for (int j = 0; j < parent->totchild; j++) {
        parent->children[j]->parent = parent;
      }
    }

    BM_ITER_MESH_INDEX (v, &iter, pbvh->header.bm, BM_VERTS_OF_MESH, i) {
      while (vnodemap[i]->parent) {
        vnodemap[i] = vnodemap[i]->parent;
      }
    }
  }

  BLI_mempool_iter loopiter;
  BLI_mempool_iternew(pbvh->header.bm->lpool, &loopiter);
  BMLoop *l = (BMLoop *)BLI_mempool_iterstep(&loopiter);
  BMEdge *e;
  BMFace *f;

  for (i = 0; l; l = (BMLoop *)BLI_mempool_iterstep(&loopiter), i++) {
    l->head.hflag &= ~flag;
  }
  BM_ITER_MESH (e, &iter, pbvh->header.bm, BM_EDGES_OF_MESH) {
    e->head.hflag &= ~flag;
  }

  BM_ITER_MESH (f, &iter, pbvh->header.bm, BM_FACES_OF_MESH) {
    f->head.hflag &= ~flag;
  }

  int totroot = roots.size();
  Vector<ReVertNode *> nstack;
  int vorder = 0, eorder = 0, lorder = 0, forder = 0;

  for (i = 0; i < totroot; i++) {
    nstack.clear();

    ReVertNode *node = roots[i];
    nstack.append(node);

    while (nstack.size() > 0) {
      ReVertNode *node2 = nstack.pop_last();

      if (node2->totchild == 0) {
        for (int j = 0; j < node2->totvert; j++) {
          v = node2->verts[j];

#if 0
          const int cd_vcol = CustomData_get_offset(&pbvh->header.bm->vdata,CD_PROP_COLOR);

          if (cd_vcol >= 0) {
            MPropCol *col = BM_ELEM_CD_PTR<MPropCol*>(node2->verts[j],cd_vcol);

            float r = 0.0f,g = 0.0f,b = 0.0f;

            ReVertNode *parent = node2->parent;
            for (int j = 0; parent->parent && j < 2; j++) {
              parent = parent->parent;
            }

            unsigned int p = (unsigned int)node2->parent;
            p = p % 65535;

            unsigned int p2 = (unsigned int)parent;
            p2 = p2 % 65535;

            r = ((float)vorder) * 0.01;
            g = ((float)p2) / 65535.0f;
            b = ((float)p2) / 65535.0f;

            r = cosf(r * 17.2343) * 0.5 + 0.5;
            g = cosf(g * 11.2343) * 0.5 + 0.5;
            b = cosf(b * 19.2343) * 0.5 + 0.5;

            col->color[0] = r;
            col->color[1] = g;
            col->color[2] = b;
            col->color[3] = 1.0f;
          }
#endif
          v->head.index = vorder++;

          BMEdge *e = v->e;
          if (!e) {
            continue;
          }

          do {
            if (!(e->head.hflag & flag)) {
              e->head.hflag |= flag;
              e->head.index = eorder++;
            }

            if (e->l) {
              BMLoop *l = e->l;

              do {
                if (!(l->head.hflag & flag)) {
                  l->head.hflag |= flag;
                  l->head.index = lorder++;
                }

                if (!(l->f->head.hflag & flag)) {
                  l->f->head.hflag |= flag;
                  l->f->head.index = forder++;
                }

                l = l->radial_next;
              } while (l != e->l);
            }
            e = e->v1 == v ? e->v1_disk_link.next : e->v2_disk_link.next;
          } while (e != v->e);
        }
      }
      else {
        for (int j = 0; j < node2->totchild; j++) {
          nstack.append(node2->children[j]);
        }
      }
    }
  }

  uint *vidx, *eidx, *lidx, *fidx;

  vidx = MEM_cnew_array<uint>(pbvh->header.bm->totvert, "vorder");
  eidx = MEM_cnew_array<uint>(pbvh->header.bm->totedge, "eorder");
  lidx = MEM_cnew_array<uint>(pbvh->header.bm->totloop, "lorder");
  fidx = MEM_cnew_array<uint>(pbvh->header.bm->totface, "forder");

  printf("v %d %d\n", vorder, pbvh->header.bm->totvert);
  printf("e %d %d\n", eorder, pbvh->header.bm->totedge);
  printf("l %d %d\n", lorder, pbvh->header.bm->totloop);
  printf("f %d %d\n", forder, pbvh->header.bm->totface);

  BM_ITER_MESH_INDEX (v, &iter, pbvh->header.bm, BM_VERTS_OF_MESH, i) {
    vidx[i] = (uint)v->head.index;
  }

  BM_ITER_MESH_INDEX (e, &iter, pbvh->header.bm, BM_EDGES_OF_MESH, i) {
    eidx[i] = (uint)e->head.index;
  }
  BM_ITER_MESH_INDEX (f, &iter, pbvh->header.bm, BM_FACES_OF_MESH, i) {
    fidx[i] = (uint)f->head.index;
  }

  BLI_mempool_iternew(pbvh->header.bm->lpool, &loopiter);
  l = (BMLoop *)BLI_mempool_iterstep(&loopiter);

  for (i = 0; l; l = (BMLoop *)BLI_mempool_iterstep(&loopiter), i++) {
    // handle orphaned loops
    if (!(l->head.hflag & flag)) {
      printf("warning in %s: orphaned loop!\n", __func__);
      l->head.index = lorder++;
    }

    lidx[i] = (uint)l->head.index;
  }

  printf("roots: %d\n", (int)roots.size());

  BM_mesh_remap(pbvh->header.bm, vidx, eidx, lidx, fidx);

  MEM_SAFE_FREE(vidx);
  MEM_SAFE_FREE(eidx);
  MEM_SAFE_FREE(lidx);
  MEM_SAFE_FREE(fidx);

  BLI_mempool_destroy(pool);
  MEM_SAFE_FREE(vnodemap);

  return pbvh->header.bm;
}

BMesh *BKE_pbvh_reorder_bmesh2(PBVH *pbvh)
{
  if (BKE_pbvh_type(pbvh) != PBVH_BMESH || pbvh->totnode == 0) {
    return pbvh->header.bm;
  }

  // try to group memory allocations by node
  struct TempNodeData {
    Vector<BMVert *> verts;
    Vector<BMEdge *> edges;
    Vector<BMFace *> faces;
  };

  Vector<TempNodeData> nodedata;
  nodedata.resize(pbvh->totnode);

  BMIter iter;
  int types[3] = {BM_VERTS_OF_MESH, BM_EDGES_OF_MESH, BM_FACES_OF_MESH};

#define VISIT_TAG BM_ELEM_TAG

  BM_mesh_elem_index_ensure(pbvh->header.bm, BM_VERT | BM_EDGE | BM_FACE);
  BM_mesh_elem_table_ensure(pbvh->header.bm, BM_VERT | BM_EDGE | BM_FACE);

  for (int i = 0; i < 3; i++) {
    BMHeader *elem;

    BM_ITER_MESH (elem, &iter, pbvh->header.bm, types[i]) {
      elem->hflag &= ~VISIT_TAG;
    }
  }

  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (!(node->flag & PBVH_Leaf)) {
      continue;
    }

    Vector<BMVert *> &verts = nodedata[i].verts;
    Vector<BMEdge *> &edges = nodedata[i].edges;
    Vector<BMFace *> &faces = nodedata[i].faces;

    BMVert *v;
    BMFace *f;

    TGSET_ITER (v, node->bm_unique_verts) {
      if (v->head.hflag & VISIT_TAG) {
        continue;
      }

      v->head.hflag |= VISIT_TAG;
      verts.append(v);

      BMEdge *e = v->e;
      do {
        if (!(e->head.hflag & VISIT_TAG)) {
          e->head.hflag |= VISIT_TAG;
          edges.append(e);
        }
        e = v == e->v1 ? e->v1_disk_link.next : e->v2_disk_link.next;
      } while (e != v->e);
    }
    TGSET_ITER_END;

    TGSET_ITER (f, node->bm_faces) {
      if (f->head.hflag & VISIT_TAG) {
        continue;
      }

      faces.append(f);
      f->head.hflag |= VISIT_TAG;
    }
    TGSET_ITER_END;
  }

  BMAllocTemplate templ = {pbvh->header.bm->totvert,
                           pbvh->header.bm->totedge,
                           pbvh->header.bm->totloop,
                           pbvh->header.bm->totface};
  struct BMeshCreateParams params = {0};

  BMesh *bm2 = BM_mesh_create(&templ, &params);

  CustomData_copy_all_layout(&pbvh->header.bm->vdata, &bm2->vdata);
  CustomData_copy_all_layout(&pbvh->header.bm->edata, &bm2->edata);
  CustomData_copy_all_layout(&pbvh->header.bm->ldata, &bm2->ldata);
  CustomData_copy_all_layout(&pbvh->header.bm->pdata, &bm2->pdata);

  CustomData_bmesh_init_pool(&bm2->vdata, pbvh->header.bm->totvert, BM_VERT);
  CustomData_bmesh_init_pool(&bm2->edata, pbvh->header.bm->totedge, BM_EDGE);
  CustomData_bmesh_init_pool(&bm2->ldata, pbvh->header.bm->totloop, BM_LOOP);
  CustomData_bmesh_init_pool(&bm2->pdata, pbvh->header.bm->totface, BM_FACE);

  Vector<BMVert *> verts;
  Vector<BMEdge *> edges;
  Vector<BMFace *> faces;

  for (int i = 0; i < pbvh->totnode; i++) {
    for (int j = 0; j < nodedata[i].verts.size(); j++) {
      BMVert *v1 = nodedata[i].verts[j];
      BMVert *v2 = BM_vert_create(bm2, v1->co, nullptr, BM_CREATE_NOP);
      BM_elem_attrs_copy_ex(pbvh->header.bm, bm2, v1, v2, 0, 0L);

      v2->head.index = v1->head.index = verts.size();
      verts.append(v2);
    }
  }

  for (int i = 0; i < pbvh->totnode; i++) {
    for (int j = 0; j < nodedata[i].edges.size(); j++) {
      BMEdge *e1 = nodedata[i].edges[j];
      BMEdge *e2 = BM_edge_create(
          bm2, verts[e1->v1->head.index], verts[e1->v2->head.index], nullptr, BM_CREATE_NOP);
      BM_elem_attrs_copy_ex(pbvh->header.bm, bm2, e1, e2, 0, 0L);

      e2->head.index = e1->head.index = edges.size();
      edges.append(e2);
    }
  }

  Vector<BMVert *> fvs;
  Vector<BMEdge *> fes;

  for (int i = 0; i < pbvh->totnode; i++) {
    for (int j = 0; j < nodedata[i].faces.size(); j++) {
      BMFace *f1 = nodedata[i].faces[j];

      fvs.clear();
      fes.clear();

      int totloop = 0;
      BMLoop *l1 = f1->l_first;
      do {
        fvs.append(verts[l1->v->head.index]);
        fes.append(edges[l1->e->head.index]);
        l1 = l1->next;
        totloop++;
      } while (l1 != f1->l_first);

      BMFace *f2 = BM_face_create(bm2, fvs.data(), fes.data(), totloop, nullptr, BM_CREATE_NOP);
      f1->head.index = f2->head.index = faces.size();
      faces.append(f2);

      // CustomData_bmesh_copy_data(&pbvh->header.bm->pdata, &bm2->pdata, f1->head.data,
      // &f2->head.data);
      BM_elem_attrs_copy_ex(pbvh->header.bm, bm2, f1, f2, 0, 0L);

      BMLoop *l2 = f2->l_first;
      do {
        BM_elem_attrs_copy_ex(pbvh->header.bm, bm2, l1, l2, 0, 0L);

        l1 = l1->next;
        l2 = l2->next;
      } while (l2 != f2->l_first);
    }
  }

  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (!(node->flag & PBVH_Leaf)) {
      continue;
    }

    int totunique = node->bm_unique_verts->length;
    int totother = node->bm_other_verts->length;
    int totface = node->bm_faces->length;

    TableGSet *bm_faces = BLI_table_gset_new_ex("bm_faces", totface);
    TableGSet *bm_other_verts = BLI_table_gset_new_ex("bm_other_verts", totunique);
    TableGSet *bm_unique_verts = BLI_table_gset_new_ex("bm_unique_verts", totother);

    BMVert *v;
    BMFace *f;

    TGSET_ITER (v, node->bm_unique_verts) {
      BLI_table_gset_insert(bm_unique_verts, verts[v->head.index]);
    }
    TGSET_ITER_END;
    TGSET_ITER (v, node->bm_other_verts) {
      BLI_table_gset_insert(bm_other_verts, verts[v->head.index]);
    }
    TGSET_ITER_END;
    TGSET_ITER (f, node->bm_faces) {
      BLI_table_gset_insert(bm_faces, faces[f->head.index]);
    }
    TGSET_ITER_END;

    BLI_table_gset_free(node->bm_faces, nullptr);
    BLI_table_gset_free(node->bm_other_verts, nullptr);
    BLI_table_gset_free(node->bm_unique_verts, nullptr);

    node->bm_faces = bm_faces;
    node->bm_other_verts = bm_other_verts;
    node->bm_unique_verts = bm_unique_verts;

    node->flag |= PBVH_UpdateTris | PBVH_UpdateRedraw;
  }

  BM_mesh_free(pbvh->header.bm);
  pbvh->header.bm = bm2;

  return bm2;
}

struct SortElem {
  BMElem *elem;
  int index;
  int cd_node_off;
};

static int sort_verts_faces(const void *va, const void *vb)
{
  SortElem *a = (SortElem *)va;
  SortElem *b = (SortElem *)vb;
  int ni1 = BM_ELEM_CD_GET_INT(a->elem, a->cd_node_off);
  int ni2 = BM_ELEM_CD_GET_INT(b->elem, b->cd_node_off);

  return ni1 - ni2;
}

static int sort_edges(const void *va, const void *vb)
{
  SortElem *a = (SortElem *)va;
  SortElem *b = (SortElem *)vb;

  BMEdge *e1 = (BMEdge *)a->elem;
  BMEdge *e2 = (BMEdge *)b->elem;

  int ni1 = BM_ELEM_CD_GET_INT(e1->v1, a->cd_node_off);
  int ni2 = BM_ELEM_CD_GET_INT(e1->v2, a->cd_node_off);
  int ni3 = BM_ELEM_CD_GET_INT(e2->v1, b->cd_node_off);
  int ni4 = BM_ELEM_CD_GET_INT(e2->v2, b->cd_node_off);

  return (ni1 + ni2) - (ni3 + ni4);
}

BMesh *BKE_pbvh_reorder_bmesh1(PBVH *pbvh)
{
  BMesh *bm = pbvh->header.bm;

  Vector<Vector<int>> save_other_vs;
  Vector<Vector<int>> save_unique_vs;
  Vector<Vector<int>> save_fs;

  save_other_vs.resize(pbvh->totnode);
  save_unique_vs.resize(pbvh->totnode);
  save_fs.resize(pbvh->totnode);

  SortElem *verts = MEM_cnew_array<SortElem>(bm->totvert, __func__);
  SortElem *edges = MEM_cnew_array<SortElem>(bm->totedge, __func__);
  SortElem *faces = MEM_cnew_array<SortElem>(bm->totface, __func__);

  BMIter iter;
  BMVert *v;
  BMEdge *e;
  BMFace *f;

  int i = 0;

  BM_ITER_MESH_INDEX (v, &iter, bm, BM_VERTS_OF_MESH, i) {
    verts[i].elem = (BMElem *)v;
    verts[i].cd_node_off = pbvh->cd_vert_node_offset;
    verts[i].index = i;
    v->head.index = i;
  }
  BM_ITER_MESH_INDEX (e, &iter, bm, BM_EDGES_OF_MESH, i) {
    edges[i].elem = (BMElem *)e;
    edges[i].cd_node_off = pbvh->cd_vert_node_offset;
    edges[i].index = i;
    e->head.index = i;
  }
  BM_ITER_MESH_INDEX (f, &iter, bm, BM_FACES_OF_MESH, i) {
    faces[i].elem = (BMElem *)f;
    faces[i].cd_node_off = pbvh->cd_face_node_offset;
    faces[i].index = i;
    f->head.index = i;
  }

  for (i = 0; i < pbvh->totnode; i++) {
    Vector<int> other_vs;
    Vector<int> unique_vs;
    Vector<int> fs;

    PBVHNode *node = pbvh->nodes + i;
    if (!(node->flag & PBVH_Leaf)) {
      continue;
    }

    BMVert *v;
    BMFace *f;

    TGSET_ITER (v, node->bm_unique_verts) {
      unique_vs.append(v->head.index);
    }
    TGSET_ITER_END;
    TGSET_ITER (v, node->bm_other_verts) {
      other_vs.append(v->head.index);
    }
    TGSET_ITER_END;
    TGSET_ITER (f, node->bm_faces) {
      fs.append(f->head.index);
    }
    TGSET_ITER_END;

    save_unique_vs[i] = unique_vs;
    save_other_vs[i] = other_vs;
    save_fs[i] = fs;
  }

  qsort(verts, bm->totvert, sizeof(SortElem), sort_verts_faces);
  qsort(edges, bm->totedge, sizeof(SortElem), sort_edges);
  qsort(faces, bm->totface, sizeof(SortElem), sort_verts_faces);

  uint *vs = MEM_cnew_array<uint>(bm->totvert, __func__);
  uint *es = MEM_cnew_array<uint>(bm->totedge, __func__);
  uint *fs = MEM_cnew_array<uint>(bm->totface, __func__);

  for (i = 0; i < bm->totvert; i++) {
    vs[i] = (uint)verts[i].index;
    verts[i].elem->head.index = verts[i].index;
  }
  for (i = 0; i < bm->totedge; i++) {
    es[i] = (uint)edges[i].index;
    edges[i].elem->head.index = edges[i].index;
  }
  for (i = 0; i < bm->totface; i++) {
    fs[i] = (uint)faces[i].index;
    faces[i].elem->head.index = faces[i].index;
  }

  BM_mesh_remap(bm, vs, es, nullptr, fs);

  // create new mappings
  BMVert **mapvs = MEM_cnew_array<BMVert *>(bm->totvert, __func__);
  BMEdge **mapes = MEM_cnew_array<BMEdge *>(bm->totedge, __func__);
  BMFace **mapfs = MEM_cnew_array<BMFace *>(bm->totface, __func__);

  BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
    mapvs[v->head.index] = v;
  }
  BM_ITER_MESH (e, &iter, bm, BM_EDGES_OF_MESH) {
    mapes[e->head.index] = e;
  }
  BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH) {
    mapfs[f->head.index] = f;
  }

  // rebuild bm_unique_verts bm_other_verts and bm_faces in pbvh nodes
  for (i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (!(node->flag & PBVH_Leaf)) {
      continue;
    }

    int tot_unique_vs = BLI_table_gset_len(node->bm_unique_verts);
    int tot_other_vs = BLI_table_gset_len(node->bm_other_verts);
    int tot_fs = BLI_table_gset_len(node->bm_faces);

    BLI_table_gset_free(node->bm_unique_verts, nullptr);
    BLI_table_gset_free(node->bm_other_verts, nullptr);
    BLI_table_gset_free(node->bm_faces, nullptr);

    node->bm_unique_verts = BLI_table_gset_new("bm_unique_verts");
    node->bm_other_verts = BLI_table_gset_new("bm_other_verts");
    node->bm_faces = BLI_table_gset_new("bm_faces");

    Vector<int> &unique_vs = save_unique_vs[i];
    Vector<int> &other_vs = save_other_vs[i];
    Vector<int> &fs = save_fs[i];

    for (int j = 0; j < tot_unique_vs; j++) {
      BLI_table_gset_add(node->bm_unique_verts, mapvs[unique_vs[j]]);
    }
    for (int j = 0; j < tot_other_vs; j++) {
      BLI_table_gset_add(node->bm_other_verts, mapvs[other_vs[j]]);
    }

    for (int j = 0; j < tot_fs; j++) {
      BLI_table_gset_add(node->bm_faces, mapfs[fs[j]]);
    }

    node->flag |= PBVH_UpdateTris;
  }

  MEM_SAFE_FREE(mapvs);
  MEM_SAFE_FREE(mapes);
  MEM_SAFE_FREE(mapfs);

  bm->elem_index_dirty |= BM_VERT | BM_EDGE | BM_FACE;
  bm->elem_table_dirty |= BM_VERT | BM_EDGE | BM_FACE;

  MEM_SAFE_FREE(vs);
  MEM_SAFE_FREE(es);
  MEM_SAFE_FREE(fs);

  MEM_SAFE_FREE(verts);
  MEM_SAFE_FREE(edges);
  MEM_SAFE_FREE(faces);

  return pbvh->header.bm;
}

// only floats! and 8 byte aligned!
struct CacheParams {
  float vchunk, echunk, lchunk, pchunk;
  int cluster_steps, cluster_size;
};

struct CacheParamDef {
  char name[32];
  float defvalue, min, max;
};

static CacheParamDef pbvh_bmesh_cache_param_def[] = {
    {"vchunk", 512.0f, 256.0f, 1024.0f * 12.0f},
    {"echunk", 512.0f, 256.0f, 1024.0f * 12.0f},
    {"lchunk", 512.0f, 256.0f, 1024.0f * 12.0f},
    {"pchunk", 512.0f, 256.0f, 1024.0f * 12.0f},
    {"cluster_steps", 512.0f, 1.0f, 256.0f},
    {"cluster_size", 512.0f, 1.0f, 8192.0f * 32.0f}};

int pbvh_bmesh_cache_test_totparams()
{
  return sizeof(pbvh_bmesh_cache_param_def) / sizeof(*pbvh_bmesh_cache_param_def);
}

void pbvh_bmesh_cache_test_default_params(CacheParams *params)
{
  float *fparams = (float *)params;
  int totparam = pbvh_bmesh_cache_test_totparams();

  for (int i = 0; i < totparam; i++) {
    fparams[i] = pbvh_bmesh_cache_param_def[i].defvalue;
  }
}

static void *hashco(float fx, float fy, float fz, float fdimen)
{
  double x = (double)fx;
  double y = (double)fy;
  double z = (double)fz;
  double dimen = (double)fdimen;

  return (void *)((intptr_t)(z * dimen * dimen * dimen + y * dimen * dimen + x * dimen));
}

struct MeshTest {
  float (*v_co)[3];
  float (*v_no)[3];
  int *v_e;
  int *v_index;
  int *v_flag;

  int *e_v1;
  int *e_v2;
  int *e_v1_next;
  int *e_v2_next;
  int *e_l;
  int *e_flag;
  int *e_index;

  int *l_v;
  int *l_e;
  int *l_f;
  int *l_next;
  int *l_prev;
  int *l_radial_next;
  int *l_radial_prev;

  int *f_l;
  int *f_index;
  int *f_flag;

  int totvert, totedge, totloop, totface;
  MemArena *arena;
};

struct ElemHeader {
  short type, hflag;
  int index;
  void *data;
};

struct MeshVert2 {
  ElemHeader head;
  float co[3];
  float no[3];
  int e;
};

struct MeshEdge2 {
  ElemHeader head;
  int v1, v2;
  int v1_next, v2_next;
  int l;
};

struct MeshLoop2 {
  ElemHeader head;
  int v, e, f, next, prev;
  int radial_next, radial_prev;
};

struct MeshFace2 {
  ElemHeader head;
  int l, len;
  float no[3];
};

struct MeshTest2 {
  MeshVert2 *verts;
  MeshEdge2 *edges;
  MeshLoop2 *loops;
  MeshFace2 *faces;

  int totvert, totedge, totloop, totface;
  MemArena *arena;
};

static MeshTest2 *meshtest2_from_bm(BMesh *bm)
{
  MeshTest2 *m2 = (MeshTest2 *)MEM_callocN(sizeof(MeshTest2), "MeshTest2");
  m2->arena = BLI_memarena_new(1024 * 32, "MeshTest2 arena");

  m2->totvert = bm->totvert;
  m2->totedge = bm->totedge;
  m2->totloop = bm->totloop;
  m2->totface = bm->totface;

  BMVert *v;
  BMEdge *e;
  BMFace *f;
  BMIter iter;

  int lindex = 0;

  BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH) {
    BMLoop *l = f->l_first;
    do {
      l->head.index = lindex++;
    } while ((l = l->next) != f->l_first);
  }

  m2->totloop = lindex;

  m2->verts = (MeshVert2 *)MEM_calloc_arrayN(bm->totvert, sizeof(MeshVert2), "MeshVert2s");
  m2->edges = (MeshEdge2 *)MEM_calloc_arrayN(bm->totedge, sizeof(MeshEdge2), "MeshEdge2s");
  m2->loops = (MeshLoop2 *)MEM_calloc_arrayN(m2->totloop, sizeof(MeshLoop2), "MeshLoop2s");
  m2->faces = (MeshFace2 *)MEM_calloc_arrayN(bm->totface, sizeof(MeshFace2), "MeshFace2s");

  bm->elem_index_dirty |= BM_VERT | BM_EDGE | BM_FACE;
  BM_mesh_elem_index_ensure(bm, BM_VERT | BM_EDGE | BM_FACE);

  BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
    const int vi = v->head.index;

    copy_v3_v3(m2->verts[vi].co, v->co);
    copy_v3_v3(m2->verts[vi].no, v->no);

    m2->verts[vi].e = v->e ? v->e->head.index : -1;
  }

  BM_ITER_MESH (e, &iter, bm, BM_EDGES_OF_MESH) {
    const int ei = e->head.index;

    m2->edges[ei].v1 = e->v1->head.index;
    m2->edges[ei].v2 = e->v2->head.index;
    m2->edges[ei].l = e->l ? e->l->head.index : -1;

    m2->edges[ei].v1_next = e->v1_disk_link.next->head.index;
    m2->edges[ei].v2_next = e->v2_disk_link.next->head.index;
  }

  BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH) {
    const int fi = f->head.index;

    m2->faces[fi].len = f->len;
    copy_v3_v3(m2->faces[fi].no, f->no);

    BMLoop *l = f->l_first;
    do {
      int li = l->head.index;

      m2->loops[li].v = l->v->head.index;
      m2->loops[li].e = l->e->head.index;
      m2->loops[li].f = l->f->head.index;

      m2->loops[li].radial_next = l->radial_next->head.index;
      m2->loops[li].radial_prev = l->radial_prev->head.index;

      m2->loops[li].next = l->next->head.index;
      m2->loops[li].prev = l->prev->head.index;
    } while ((l = l->next) != f->l_first);
  }

  return m2;
}

static void free_meshtest2(MeshTest2 *m2)
{
  BLI_memarena_free(m2->arena);
  MEM_freeN(m2);
}

static MeshTest *meshtest_from_bm(BMesh *bm)
{
  MeshTest *m = MEM_cnew<MeshTest>("MeshTest");
  m->arena = BLI_memarena_new(1024 * 32, "m->arena");

  m->v_co = (float(*)[3])BLI_memarena_alloc(m->arena, bm->totvert * sizeof(*m->v_co));
  m->v_no = (float(*)[3])BLI_memarena_alloc(m->arena, bm->totvert * sizeof(*m->v_no));
  m->v_e = (int *)BLI_memarena_alloc(m->arena, bm->totvert * sizeof(*m->v_e));
  m->v_flag = (int *)BLI_memarena_alloc(m->arena, bm->totvert * sizeof(*m->v_flag));
  m->v_index = (int *)BLI_memarena_alloc(m->arena, bm->totvert * sizeof(*m->v_index));

  m->e_v1 = (int *)BLI_memarena_alloc(m->arena, bm->totedge * sizeof(*m->e_v1));
  m->e_v1_next = (int *)BLI_memarena_alloc(m->arena, bm->totedge * sizeof(*m->e_v1));
  m->e_v2 = (int *)BLI_memarena_alloc(m->arena, bm->totedge * sizeof(*m->e_v1));
  m->e_v2_next = (int *)BLI_memarena_alloc(m->arena, bm->totedge * sizeof(*m->e_v1));
  m->e_l = (int *)BLI_memarena_alloc(m->arena, bm->totedge * sizeof(*m->e_v1));
  m->e_index = (int *)BLI_memarena_alloc(m->arena, bm->totedge * sizeof(*m->e_v1));
  m->e_flag = (int *)BLI_memarena_alloc(m->arena, bm->totedge * sizeof(*m->e_v1));

  m->l_v = (int *)BLI_memarena_alloc(m->arena, bm->totloop * sizeof(*m->l_e));
  m->l_e = (int *)BLI_memarena_alloc(m->arena, bm->totloop * sizeof(*m->l_e));
  m->l_f = (int *)BLI_memarena_alloc(m->arena, bm->totloop * sizeof(*m->l_e));
  m->l_next = (int *)BLI_memarena_alloc(m->arena, bm->totloop * sizeof(*m->l_e));
  m->l_prev = (int *)BLI_memarena_alloc(m->arena, bm->totloop * sizeof(*m->l_e));
  m->l_radial_next = (int *)BLI_memarena_alloc(m->arena, bm->totloop * sizeof(*m->l_e));
  m->l_radial_prev = (int *)BLI_memarena_alloc(m->arena, bm->totloop * sizeof(*m->l_e));

  m->f_l = (int *)BLI_memarena_alloc(m->arena, bm->totface * sizeof(*m->f_l));

  m->totvert = bm->totvert;
  m->totedge = bm->totedge;
  m->totface = bm->totface;
  m->totloop = bm->totloop;

  BMVert *v;
  BMEdge *e;
  BMFace *f;
  BMIter iter;

  int lindex = 0;

  BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH) {
    BMLoop *l = f->l_first;
    do {
      l->head.index = lindex++;
    } while ((l = l->next) != f->l_first);
  }

  bm->elem_index_dirty |= BM_VERT | BM_EDGE | BM_FACE;
  BM_mesh_elem_index_ensure(bm, BM_VERT | BM_EDGE | BM_FACE);

  BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
    copy_v3_v3(m->v_co[v->head.index], v->co);
    copy_v3_v3(m->v_no[v->head.index], v->no);

    m->v_e[v->head.index] = v->e ? v->e->head.index : -1;
  }

  BM_ITER_MESH (e, &iter, bm, BM_EDGES_OF_MESH) {
    m->e_v1[e->head.index] = e->v1->head.index;
    m->e_v2[e->head.index] = e->v2->head.index;

    m->e_v1_next[e->head.index] = e->v1_disk_link.next->head.index;
    m->e_v2_next[e->head.index] = e->v2_disk_link.next->head.index;

    m->e_l[e->head.index] = e->l ? e->l->head.index : -1;
  }

  BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH) {
    m->f_l[f->head.index] = f->l_first->head.index;

    BMLoop *l = f->l_first;
    do {
      const int li = l->head.index;

      m->l_e[li] = l->e->head.index;
      m->l_v[li] = l->v->head.index;
      m->l_f[li] = f->head.index;
      m->l_next[li] = l->next->head.index;
      m->l_prev[li] = l->prev->head.index;
      m->l_radial_next[li] = l->radial_next->head.index;
      m->l_radial_prev[li] = l->radial_prev->head.index;
    } while ((l = l->next) != f->l_first);
  }

  return m;
}

static void free_meshtest(MeshTest *m)
{
  BLI_memarena_free(m->arena);
  MEM_freeN(m);
}

#define SMOOTH_TEST_STEPS 20

double pbvh_bmesh_smooth_test(BMesh *bm, PBVH *pbvh)
{
  double average = 0.0f;
  double average_tot = 0.0f;

  for (int iter = 0; iter < SMOOTH_TEST_STEPS; iter++) {
    RNG *rng = BLI_rng_new(0);

    double time1 = PIL_check_seconds_timer();

    for (int step = 0; step < 5; step++) {
      for (int i = 0; i < bm->totvert; i++) {
        int vi = BLI_rng_get_int(rng) % bm->totvert;
        BMVert *v = bm->vtable[vi];
        BMEdge *e = v->e;
        float co[3];

        zero_v3(co);
        int tot = 0.0;

        if (!e) {
          continue;
        }

        do {
          BMVert *v2 = BM_edge_other_vert(e, v);
          float co2[3];

          sub_v3_v3v3(co2, v2->co, v->co);
          madd_v3_v3fl(co2, v->no, -dot_v3v3(v->no, co2) * 0.9f);
          add_v3_v3(co, co2);

          tot++;

          e = e->v1 == v ? e->v1_disk_link.next : e->v2_disk_link.next;
        } while (e != v->e);

        if (tot == 0.0) {
          continue;
        }

        mul_v3_fl(co, 1.0f / (float)tot);
        madd_v3_v3fl(v->co, co, 0.5f);
      }
    }

    double time2 = PIL_check_seconds_timer();

    double time = time2 - time1;

    printf("  time: %.5f, %d of %d\n", time, iter, SMOOTH_TEST_STEPS);

    // skip first five
    if (iter >= 5) {
      average += time;
      average_tot += 1.0f;
    }

    BLI_rng_free(rng);
  }

  printf("time: %.5f\n", average / average_tot);
  return average / average_tot;
}

double pbvh_meshtest2_smooth_test(MeshTest2 *m2, PBVH *pbvh)
{
  double average = 0.0f;
  double average_tot = 0.0f;

  for (int iter = 0; iter < SMOOTH_TEST_STEPS; iter++) {
    RNG *rng = BLI_rng_new(0);

    double time1 = PIL_check_seconds_timer();

    for (int step = 0; step < 5; step++) {
      for (int i = 0; i < m2->totvert; i++) {
        int vi = BLI_rng_get_int(rng) % m2->totvert;
        MeshVert2 *v = m2->verts + vi;
        MeshEdge2 *e = v->e != -1 ? m2->edges + v->e : nullptr;
        float co[3];

        zero_v3(co);
        int tot = 0.0;

        if (!e) {
          continue;
        }

        int enext = -1;

        do {
          MeshVert2 *v2 = vi == e->v1 ? m2->verts + e->v2 : m2->verts + e->v1;
          float co2[3];

          sub_v3_v3v3(co2, v2->co, v->co);
          madd_v3_v3fl(co2, v->no, -dot_v3v3(v->no, co2) * 0.9f);
          add_v3_v3(co, co2);

          tot++;

          enext = e->v1 == vi ? e->v1_next : e->v2_next;
          e = m2->edges + enext;
        } while (enext != v->e);

        if (tot == 0.0) {
          continue;
        }

        mul_v3_fl(co, 1.0f / (float)tot);
        madd_v3_v3fl(v->co, co, 0.5f);
      }
    }

    double time2 = PIL_check_seconds_timer();

    double time = time2 - time1;

    printf("  time: %.5f, %d of %d\n", time, iter, SMOOTH_TEST_STEPS);

    // skip first five
    if (iter >= 5) {
      average += time;
      average_tot += 1.0f;
    }

    BLI_rng_free(rng);
  }

  printf("time: %.5f\n", average / average_tot);
  return average / average_tot;
}

double pbvh_meshtest_smooth_test(MeshTest *m, PBVH *pbvh)
{
  double average = 0.0f;
  double average_tot = 0.0f;

  for (int iter = 0; iter < SMOOTH_TEST_STEPS; iter++) {
    RNG *rng = BLI_rng_new(0);

    double time1 = PIL_check_seconds_timer();

    for (int step = 0; step < 5; step++) {
      for (int i = 0; i < m->totvert; i++) {
        int vi = BLI_rng_get_int(rng) % m->totvert;
        // BMVert *v = bm->vtable[vi];
        const int startei = m->v_e[vi];
        int ei = startei;

        float co[3];

        zero_v3(co);
        int tot = 0.0;

        if (ei == -1) {
          continue;
        }

        const float *no = m->v_no[vi];
        const float *vco = m->v_co[vi];

        do {
          int ev1 = m->e_v1[ei];
          int ev2 = m->e_v2[ei];

          int v2i = ev1 == vi ? ev2 : ev1;

          float co2[3];

          sub_v3_v3v3(co2, m->v_co[v2i], vco);
          madd_v3_v3fl(co2, no, -dot_v3v3(no, co2) * 0.9f);
          add_v3_v3(co, co2);

          tot++;

          ei = ev1 == vi ? m->e_v1_next[ei] : m->e_v2_next[ei];
        } while (ei != startei);

        if (tot == 0.0) {
          continue;
        }

        mul_v3_fl(co, 1.0f / (float)tot);
        madd_v3_v3fl(m->v_co[vi], co, 0.5f);
      }
    }

    double time2 = PIL_check_seconds_timer();

    double time = time2 - time1;

    printf("  time: %.5f, %d of %d\n", time, iter + 1, SMOOTH_TEST_STEPS);

    // skip first five
    if (iter >= 5) {
      average += time;
      average_tot += 1.0f;
    }

    BLI_rng_free(rng);
  }

  printf("time: %.5f\n", average / average_tot);
  return average / average_tot;
}

/*
test results from blenderartists thread:

random, cluster,  percent,  data, data_perc, indices, ind_perc, mem (gb)
  [1.22,    1.04,   14.42,  0.73,     67,     0.94,    29,      0],
  [1.49,    1.46,   2.35,   1.10,     36,     1.17,    27,      0],
  [1.29,    1.13,     14,   0.75,  71.54,     0.89,    45.08 ,  0],
  [1.58,    1.40,   12.3,   1.09,   44.7,     1.11,    42.42,   16],
  [1.53,    1.36,   12.77,  1.08,  41.6,      1.07,    42.91,   0],
  [1.56,    1.39,   12.47,  1.09,  42.65,     1.10,    42.15,   16],
  [1.22,    1.06,   15.05,  0.75,   63.85,    0.82,    49.67,   32]

[random]           average: 1.41 variange: 0.15 median: 1.49
[cluster]          average: 1.26 variange: 0.17 median: 1.36
[cluster-percent]  average: 11.91 variange: 4.02 median: 12.77
[data]             average: 0.94 variange: 0.17 median: 1.08
[data-percent]     average: 52.48 variange: 13.37 median: 44.70
[indices]          average: 1.01 variange: 0.12 median: 1.07
[indices-percent]  average: 39.75 variange: 7.82 median: 42.42

So looks like the biggest gain is from replacing pointers with indices
(which lessens total memory bandwidth).  The pure data-oriented version
is a tad bit faster then the index-replacement one, but not by that much.
*/

void pbvh_bmesh_cache_test(CacheParams *params, BMesh **r_bm, PBVH **r_pbvh_out)
{
  // build mesh
  const int steps = 325;

  printf("== Starting Test ==\n");

  printf("building test mesh. . .\n");

  BMAllocTemplate templ = {0, 0, 0, 0};
  BMeshCreateParams bmparams = {};

  bmparams.id_elem_mask = BM_VERT | BM_EDGE | BM_FACE;
  bmparams.id_map = true;
  bmparams.create_unique_ids = true;
  bmparams.temporary_ids = false;
  bmparams.no_reuse_ids = false;

  BMesh *bm = BM_mesh_create(&templ, &bmparams);

  // reinit pools
  BLI_mempool_destroy(bm->vpool);
  BLI_mempool_destroy(bm->epool);
  BLI_mempool_destroy(bm->lpool);
  BLI_mempool_destroy(bm->fpool);

  bm->vpool = BLI_mempool_create(sizeof(BMVert), 0, (int)params->vchunk, BLI_MEMPOOL_ALLOW_ITER);
  bm->epool = BLI_mempool_create(sizeof(BMEdge), 0, (int)params->echunk, BLI_MEMPOOL_ALLOW_ITER);
  bm->lpool = BLI_mempool_create(sizeof(BMLoop), 0, (int)params->lchunk, BLI_MEMPOOL_ALLOW_ITER);
  bm->fpool = BLI_mempool_create(sizeof(BMFace), 0, (int)params->pchunk, BLI_MEMPOOL_ALLOW_ITER);

  GHash *vhash = BLI_ghash_ptr_new("vhash");

  float df = 1.0f / (float)steps;

  int hashdimen = steps * 8;

  BMVert **grid = MEM_cnew_array<BMVert *>(steps * steps, "bmvert grid");

  BM_data_layer_add_named(bm, &bm->vdata, CD_PROP_INT32, "__dyntopo_vert_node");
  BM_data_layer_add_named(bm, &bm->pdata, CD_PROP_INT32, "__dyntopo_face_node");
  BM_data_layer_add_named(bm, &bm->pdata, CD_PROP_INT32, ".sculpt_face_set");
  BM_data_layer_add_named(bm, &bm->pdata, CD_PROP_INT32, ".sculpt_boundary_flags");

  BM_data_layer_add(bm, &bm->vdata, CD_PAINT_MASK);
  BM_data_layer_add(bm, &bm->vdata, CD_DYNTOPO_VERT);
  BM_data_layer_add(bm, &bm->vdata, CD_PROP_COLOR);

  BMIdMap *idmap = BM_idmap_new(bm, BM_VERT | BM_EDGE | BM_FACE);

  for (int side = 0; side < 6; side++) {
    int axis = side >= 3 ? side - 3 : side;
    float sign = side >= 3 ? -1.0f : 1.0f;

    printf("AXIS: %d\n", axis);

    float u = 0.0f;

    for (int i = 0; i < steps; i++, u += df) {
      float v = 0.0f;

      for (int j = 0; j < steps; j++, v += df) {
        float co[3];

        co[axis] = u;
        co[(axis + 1) % 3] = v;
        co[(axis + 2) % 3] = sign;

        // turn into sphere
        normalize_v3(co);

        void *key = hashco(co[0], co[1], co[2], hashdimen);

#if 0
        printf("%.3f %.3f %.3f, key: %p i: %d j: %d df: %f, u: %f v: %f\n",
          co[0],
          co[1],
          co[2],
          key,
          i,
          j,
          df,
          u,
          v);
#endif

        void **val = nullptr;

        if (!BLI_ghash_ensure_p(vhash, key, &val)) {
          BMVert *v2 = BM_vert_create(bm, co, nullptr, BM_CREATE_NOP);

          *val = (void *)v2;
        }

        BMVert *v2 = (BMVert *)*val;
        int idx = j * steps + i;

        grid[idx] = v2;
      }
    }

    for (int i = 0; i < steps - 1; i++) {
      for (int j = 0; j < steps - 1; j++) {
        int idx1 = j * steps + i;
        int idx2 = (j + 1) * steps + i;
        int idx3 = (j + 1) * steps + i + 1;
        int idx4 = j * steps + i + 1;

        BMVert *v1 = grid[idx1];
        BMVert *v2 = grid[idx2];
        BMVert *v3 = grid[idx3];
        BMVert *v4 = grid[idx4];

        if (v1 == v2 || v1 == v3 || v1 == v4 || v2 == v3 || v2 == v4 || v3 == v4) {
          printf("ERROR!\n");
          continue;
        }

        if (sign < 0) {
          BMVert *vs[4] = {v4, v3, v2, v1};
          BM_face_create_verts(bm, vs, 4, nullptr, BM_CREATE_NOP, true);
        }
        else {
          BMVert *vs[4] = {v1, v2, v3, v4};
          BM_face_create_verts(bm, vs, 4, nullptr, BM_CREATE_NOP, true);
        }
      }
    }
  }

  // randomize
  uint *rands[4];
  uint tots[4] = {(uint)bm->totvert, (uint)bm->totedge, (uint)bm->totloop, (uint)bm->totface};

  RNG *rng = BLI_rng_new(0);

  for (uint i = 0; i < 4; i++) {
    rands[i] = MEM_cnew_array<uint>(tots[i], "rands[i]");

    for (uint j = 0; j < tots[i]; j++) {
      rands[i][j] = j;
    }

    for (uint j = 0; j < tots[i] >> 1; j++) {
      int j2 = BLI_rng_get_int(rng) % tots[i];
      SWAP(uint, rands[i][j], rands[i][j2]);
    }
  }

  BM_mesh_remap(bm, rands[0], rands[1], rands[2], rands[3]);

  for (int i = 0; i < 4; i++) {
    MEM_SAFE_FREE(rands[i]);
  }

  BLI_rng_free(rng);
  BLI_ghash_free(vhash, nullptr, nullptr);
  MEM_SAFE_FREE(grid);

  printf("totvert: %d, totface: %d, tottri: %d\n", bm->totvert, bm->totface, bm->totface * 2);

  int cd_vert_node = CustomData_get_named_layer_index(
      &bm->vdata, CD_PROP_INT32, "__dyntopo_vert_node");
  int cd_face_node = CustomData_get_named_layer_index(
      &bm->pdata, CD_PROP_INT32, "__dyntopo_face_node");
  int cd_face_area = CustomData_get_named_layer_index(
      &bm->pdata, CD_PROP_FLOAT2, "__dyntopo_face_areas");
  int cd_boundary_flag = CustomData_get_named_layer_index(
      &bm->vdata, CD_PROP_INT32, ".sculpt_boundary_flags");

  cd_vert_node = bm->vdata.layers[cd_vert_node].offset;
  cd_boundary_flag = bm->vdata.layers[cd_boundary_flag].offset;
  cd_face_node = bm->pdata.layers[cd_face_node].offset;
  cd_face_area = bm->pdata.layers[cd_face_area].offset;

  const int cd_sculpt_vert = CustomData_get_offset(&bm->vdata, CD_DYNTOPO_VERT);
  BMLog *bmlog = BM_log_create(bm, idmap);

  PBVH *pbvh = BKE_pbvh_new(PBVH_BMESH);

  bm->elem_table_dirty |= BM_VERT | BM_EDGE | BM_FACE;
  bm->elem_index_dirty |= BM_VERT | BM_EDGE | BM_FACE;
  BM_mesh_elem_table_ensure(bm, BM_VERT | BM_FACE);
  BM_mesh_elem_index_ensure(bm, BM_VERT | BM_EDGE | BM_FACE);

  BKE_pbvh_build_bmesh(pbvh,
                       nullptr,
                       bm,
                       false,
                       bmlog,
                       idmap,
                       cd_vert_node,
                       cd_face_node,
                       cd_sculpt_vert,
                       cd_face_area,
                       cd_boundary_flag,
                       false,
                       true);

  int loop_size = sizeof(BMLoop) - sizeof(void *) * 4;

  size_t s1 = 0, s2 = 0, s3 = 0;
  s1 = sizeof(BMVert) * (size_t)bm->totvert + sizeof(BMEdge) * (size_t)bm->totedge +
       sizeof(BMLoop) * (size_t)bm->totloop + sizeof(BMFace) * (size_t)bm->totface;
  s2 = sizeof(MeshVert2) * (size_t)bm->totvert + sizeof(MeshEdge2) * (size_t)bm->totedge +
       sizeof(MeshLoop2) * (size_t)bm->totloop + sizeof(MeshFace2) * (size_t)bm->totface;
  s3 = (size_t)loop_size * (size_t)bm->totvert + sizeof(BMEdge) * (size_t)bm->totedge +
       sizeof(BMLoop) * (size_t)bm->totloop + sizeof(BMFace) * (size_t)bm->totface;

  double times[4];
  const char *names[4];

  int cd_overhead = 0;
  CustomData *cdatas[4] = {&bm->vdata, &bm->edata, &bm->ldata, &bm->pdata};
  int ctots[4] = {bm->totvert, bm->totedge, bm->totloop, bm->totface};
  for (int i = 0; i < 4; i++) {
    cd_overhead += cdatas[i]->totsize * ctots[i];
  }

  s1 += cd_overhead;
  s2 += cd_overhead;

  printf("    bmesh mem size: %.2fmb %.2fmb\n",
         (float)s1 / 1024.0f / 1024.0f,
         (float)s3 / 1024.0f / 1024.0f);
  printf("meshtest2 mem size: %.2fmb\n", (float)s2 / 1024.0f / 1024.0f);

  printf("= BMesh random order\n");
  times[0] = pbvh_bmesh_smooth_test(bm, pbvh);
  names[0] = "random order";

  BMesh *bm2 = BKE_pbvh_reorder_bmesh(pbvh);

  printf("= BMesh vertex cluster order\n");

  bm2->elem_table_dirty |= BM_VERT | BM_EDGE | BM_FACE;
  bm2->elem_index_dirty |= BM_VERT | BM_EDGE | BM_FACE;
  BM_mesh_elem_table_ensure(bm2, BM_VERT | BM_FACE);
  BM_mesh_elem_index_ensure(bm2, BM_VERT | BM_EDGE | BM_FACE);

  times[1] = pbvh_bmesh_smooth_test(bm2, pbvh);
  names[1] = "vertex cluser";

  printf("= Pure data-oriented (struct of arrays)\n");
  MeshTest *m = meshtest_from_bm(bm2);

  times[2] = pbvh_meshtest_smooth_test(m, pbvh);
  names[2] = "data-oriented";

  free_meshtest(m);

  printf("= Object-oriented but with integer indices instead of pointers\n");
  MeshTest2 *m2 = meshtest2_from_bm(bm2);

  times[3] = pbvh_meshtest2_smooth_test(m2, pbvh);
  names[3] = "integer indices";

  free_meshtest2(m2);

  if (bm2 && bm2 != bm) {
    BM_mesh_free(bm2);
  }

  if (r_bm) {
    *r_bm = bm;
  }
  else {
    BM_mesh_free(bm);
  }

  if (r_pbvh_out) {
    *r_pbvh_out = pbvh;
  }
  else {
    BKE_pbvh_free(pbvh);
  }

  printf("\n== Times ==\n");

  for (int i = 0; i < ARRAY_SIZE(times); i++) {
    if (i > 0) {
      double perc = (times[0] / times[i] - 1.0) * 100.0;
      printf("  %s : %.2f (%.2f%% improvement)\n", names[i], times[i], perc);
    }
    else {
      printf("  %s : %.2f\n", names[i], times[i]);
    }
  }

  printf("== Test Finished ==\n");
}

#include "BLI_smallhash.h"

static void hash_test()
{
  const int count = 1024 * 1024 * 4;

  int *data = MEM_cnew_array<int>(count, "test data");

  TableGSet *gs = BLI_table_gset_new("test");
  GHash *gh = BLI_ghash_ptr_new("test");
  SmallHash sh;

  BLI_smallhash_init(&sh);
  RNG *rng = BLI_rng_new(0);

  for (int i = 0; i < count; i++) {
    data[i] = i;
  }

  printf("== creation: table_gset ==\n");
  double t = PIL_check_seconds_timer();

  for (int i = 0; i < count; i++) {
    int ri = BLI_rng_get_int(rng) % count;
    int *ptr = (int *)POINTER_FROM_INT(data[ri]);

    BLI_table_gset_add(gs, ptr);
  }
  printf("  %.3f\n", PIL_check_seconds_timer() - t);

  printf("== creation: ghash ==\n");
  t = PIL_check_seconds_timer();

  for (int i = 0; i < count; i++) {
    int ri = BLI_rng_get_int(rng) % count;
    int *ptr = (int *)POINTER_FROM_INT(data[ri]);

    BLI_ghash_insert(gh, ptr, POINTER_FROM_INT(i));
  }
  printf("  %.3f\n", PIL_check_seconds_timer() - t);

  printf("== creation: small hash ==\n");
  t = PIL_check_seconds_timer();

  for (int i = 0; i < count; i++) {
    int ri = BLI_rng_get_int(rng) % count;
    int *ptr = (int *)POINTER_FROM_INT(data[ri]);

    BLI_smallhash_insert(&sh, (uintptr_t)ptr, POINTER_FROM_INT(i));
  }
  printf("  %.3f\n", PIL_check_seconds_timer() - t);

  printf("== lookup: g hash ==\n");
  t = PIL_check_seconds_timer();

  for (int i = 0; i < count; i++) {
    int ri = BLI_rng_get_int(rng) % count;
    int *ptr = (int *)POINTER_FROM_INT(data[ri]);

    BLI_ghash_lookup(gh, ptr);
  }
  printf("  %.3f\n", PIL_check_seconds_timer() - t);

  printf("== lookup: small hash ==\n");
  t = PIL_check_seconds_timer();

  for (int i = 0; i < count; i++) {
    int ri = BLI_rng_get_int(rng) % count;
    int *ptr = (int *)POINTER_FROM_INT(data[ri]);

    BLI_smallhash_lookup(&sh, (uintptr_t)ptr);
  }
  printf("  %.3f\n", PIL_check_seconds_timer() - t);

  BLI_rng_free(rng);
  BLI_ghash_free(gh, nullptr, nullptr);
  BLI_smallhash_release(&sh);
  BLI_table_gset_free(gs, nullptr);

  MEM_freeN(data);
}

void pbvh_bmesh_do_cache_test()
{
  for (int i = 0; i < 15; i++) {
    printf("\n\n====== %d of %d =====\n", i + 1, 15);
    hash_test();
  }
  // pbvh_bmesh_cache_test_default_params(&params);
  // pbvh_bmesh_cache_test(&params, &bm, &pbvh);
}

/* saves all bmesh references to internal indices, to be restored later */
void BKE_pbvh_bmesh_save_indices(PBVH *pbvh)
{
  BM_mesh_elem_index_ensure(pbvh->header.bm, BM_VERT | BM_EDGE | BM_FACE);

  BMFace *f;
  BMVert *v;
  BMIter iter;

  int j = 0;

  BM_ITER_MESH (f, &iter, pbvh->header.bm, BM_FACES_OF_MESH) {
    BMLoop *l = f->l_first;

    do {
      l->head.index = j++;
    } while ((l = l->next) != f->l_first);
  }

  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (!(node->flag & PBVH_Leaf)) {
      continue;
    }

    node->prim_indices = (int *)MEM_calloc_arrayN(1 + node->bm_faces->length +
                                                      node->bm_unique_verts->length,
                                                  sizeof(int),
                                                  "saved bmesh indices");

    int j = 0;

    TGSET_ITER (f, node->bm_faces) {
      node->prim_indices[j++] = f->head.index;
    }
    TGSET_ITER_END;

    // flag start of vertex array
    node->prim_indices[j++] = -1;

    TGSET_ITER (v, node->bm_unique_verts) {
      node->prim_indices[j++] = v->head.index;
    }
    TGSET_ITER_END;

    node->totprim = j;

    // don't try to save invalid triangulation
    if (node->flag & PBVH_UpdateTris) {
      continue;
    }

    // now do tribufs
    for (j = 0; j < node->tot_tri_buffers + 1; j++) {
      PBVHTriBuf *tribuf = j == node->tot_tri_buffers ? node->tribuf : node->tri_buffers + j;

      if (!tribuf) {
        break;
      }

      for (int k = 0; k < tribuf->totvert; k++) {
        if (i == 35 && k == 12) {
          printf("eek!");
        }

        tribuf->verts[k].i = ((BMVert *)tribuf->verts[k].i)->head.index;
      }

      for (int k = 0; k < tribuf->totloop; k++) {
        tribuf->loops[k] = ((BMLoop *)tribuf->loops[k])->head.index;
      }

      for (int k = 0; k < tribuf->tottri; k++) {
        PBVHTri *tri = tribuf->tris + k;
        BMFace *f = (BMFace *)tri->f.i;

        tri->f.i = f->head.index;

        for (int l = 0; l < 3; l++) {
          tri->l[l] = ((BMLoop *)tri->l[l])->head.index;
        }
      }
    }
  }
}

/* restore bmesh references from previously indices saved by BKE_pbvh_bmesh_save_indices */
void BKE_pbvh_bmesh_from_saved_indices(PBVH *pbvh)
{
  BM_mesh_elem_table_ensure(pbvh->header.bm, BM_VERT | BM_EDGE | BM_FACE);
  BM_mesh_elem_index_ensure(pbvh->header.bm, BM_VERT | BM_EDGE | BM_FACE);

  Vector<BMLoop *> ltable;

  BMFace *f;
  BMIter iter;
  int i = 0;

  BM_ITER_MESH (f, &iter, pbvh->header.bm, BM_FACES_OF_MESH) {
    BMLoop *l = f->l_first;

    do {
      l->head.index = i++;
      ltable.append(l);
    } while ((l = l->next) != f->l_first);
  }

  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (!(node->flag & PBVH_Leaf)) {
      continue;
    }

    BLI_table_gset_free(node->bm_unique_verts, nullptr);
    BLI_table_gset_free(node->bm_faces, nullptr);

    if (node->bm_other_verts) {
      BLI_table_gset_free(node->bm_other_verts, nullptr);
    }

    node->bm_other_verts = BLI_table_gset_new("bm_other_verts");
    node->flag |= PBVH_UpdateOtherVerts;

    node->bm_faces = BLI_table_gset_new("bm_faces");
    node->bm_unique_verts = BLI_table_gset_new("bm_verts");

    int j = 0;
    int *data = node->prim_indices;

    while (data[j] != -1 && j < node->totprim) {
      BMFace *f = pbvh->header.bm->ftable[data[j]];
      BM_ELEM_CD_SET_INT(f, pbvh->cd_face_node_offset, i);

      BLI_table_gset_insert(node->bm_faces, f);
      j++;
    }

    j++;

    while (j < node->totprim) {
      if (data[j] < 0 || data[j] >= pbvh->header.bm->totvert) {
        printf("%s: bad vertex at index %d!\n", __func__, data[j]);
        continue;
      }
      BMVert *v = pbvh->header.bm->vtable[data[j]];
      BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, i);

      BLI_table_gset_insert(node->bm_unique_verts, v);
      j++;
    }

    MEM_SAFE_FREE(node->prim_indices);

    // don't try to load invalid triangulation
    if (node->flag & PBVH_UpdateTris) {
      continue;
    }

    for (j = 0; j < node->tot_tri_buffers + 1; j++) {
      PBVHTriBuf *tribuf = j == node->tot_tri_buffers ? node->tribuf : node->tri_buffers + j;

      if (!tribuf) {
        break;
      }

      for (int k = 0; k < tribuf->totvert; k++) {
        tribuf->verts[k].i = (intptr_t)pbvh->header.bm->vtable[tribuf->verts[k].i];
      }

      for (int k = 0; k < tribuf->totloop; k++) {
        tribuf->loops[k] = (uintptr_t)ltable[tribuf->loops[k]];
      }

      for (int k = 0; k < tribuf->tottri; k++) {
        PBVHTri *tri = tribuf->tris + k;

        for (int l = 0; l < 3; l++) {
          tri->l[l] = (uintptr_t)ltable[tri->l[l]];
        }

        tri->f.i = (intptr_t)pbvh->header.bm->ftable[tri->f.i];
      }
    }

    node->prim_indices = nullptr;
    node->totprim = 0;
  }
}

static void pbvh_bmesh_fetch_cdrefs(PBVH *pbvh)
{
  BMesh *bm = pbvh->header.bm;

  int idx = CustomData_get_named_layer_index(
      &bm->vdata, CD_PROP_INT32, SCULPT_ATTRIBUTE_NAME(dyntopo_node_id_vertex));
  pbvh->cd_vert_node_offset = bm->vdata.layers[idx].offset;

  idx = CustomData_get_named_layer_index(
      &bm->pdata, CD_PROP_INT32, SCULPT_ATTRIBUTE_NAME(dyntopo_node_id_face));
  pbvh->cd_face_node_offset = bm->pdata.layers[idx].offset;

  idx = CustomData_get_named_layer_index(
      &bm->pdata, CD_PROP_FLOAT2, SCULPT_ATTRIBUTE_NAME(face_areas));
  pbvh->cd_face_area = bm->pdata.layers[idx].offset;

  pbvh->cd_vert_mask_offset = CustomData_get_offset(&bm->vdata, CD_PAINT_MASK);
  pbvh->cd_faceset_offset = CustomData_get_offset_named(
      &pbvh->header.bm->pdata, CD_PROP_INT32, ".sculpt_face_set");
  pbvh->cd_sculpt_vert = CustomData_get_offset(&bm->vdata, CD_DYNTOPO_VERT);
}

void BKE_pbvh_bmesh_set_toolflags(PBVH *pbvh, bool use_toolflags)
{
  if (use_toolflags == pbvh->header.bm->use_toolflags) {
    return;
  }

  // BKE_pbvh_bmesh_save_indices(pbvh);
  BM_mesh_toolflags_set(pbvh->header.bm, use_toolflags);

  // customdata layout might've changed
  pbvh_bmesh_fetch_cdrefs(pbvh);

  // BKE_pbvh_bmesh_from_saved_indices(pbvh);
}
