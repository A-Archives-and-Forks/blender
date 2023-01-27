/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2020 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup edsculpt
 */

#include "MEM_guardedalloc.h"

#include "BLI_alloca.h"
#include "BLI_array.h"
#include "BLI_blenlib.h"
#include "BLI_compiler_attrs.h"
#include "BLI_compiler_compat.h"
#include "BLI_hash.h"

#include "BLI_math.h"
#include "BLI_rand.h"
#include "BLI_task.h"
#include "BLI_threads.h"
#include "BLI_utildefines.h"

#include "DNA_brush_types.h"
#include "DNA_meshdata_types.h"

#include "BKE_attribute.h"
#include "BKE_brush.h"
#include "BKE_context.h"
#include "BKE_global.h"
#include "BKE_mesh.h"
#include "BKE_mesh_mapping.h"
#include "BKE_object.h"
#include "BKE_paint.h"
#include "BKE_pbvh.h"

#include "sculpt_intern.h"

#include "RNA_access.h"
#include "RNA_define.h"

#include "atomic_ops.h"
#include "bmesh.h"

#include <math.h>
#include <stdlib.h>

#include "BLI_compiler_compat.h"
#include "BLI_utildefines.h"

void SCULPT_reproject_cdata(SculptSession *ss,
                            PBVHVertRef vertex,
                            float origco[3],
                            float origno[3])
{
  BMVert *v = (BMVert *)vertex.i;

  if (!ss->bm || !v->e) {
    return;
  }

  MSculptVert *mv = BKE_PBVH_SCULPTVERT(ss->cd_sculpt_vert, v);

  // int totuv = CustomData_number_of_layers(&ss->bm->ldata, CD_PROP_FLOAT2);
  CustomData *ldata = &ss->bm->ldata;

  int totuv = 0;
  CustomDataLayer *uvlayer = NULL;

  if (ldata->typemap[CD_PROP_FLOAT2] != -1) {
    for (int i = ldata->typemap[CD_PROP_FLOAT2];
         i < ldata->totlayer && ldata->layers[i].type == CD_PROP_FLOAT2;
         i++) {
      totuv++;
    }

    uvlayer = ldata->layers + ldata->typemap[CD_PROP_FLOAT2];
  }

  BMEdge *e;
  int tag = BM_ELEM_TAG_ALT;

  float origin[3];
  float ray[3];

  copy_v3_v3(origin, v->co);
  copy_v3_v3(ray, v->no);
  negate_v3(ray);

  struct IsectRayPrecalc precalc;
  isect_ray_tri_watertight_v3_precalc(&precalc, ray);

  float *lastuvs = BLI_array_alloca(lastuvs, totuv * 2);
  bool *snapuvs = BLI_array_alloca(snapuvs, totuv);

  e = v->e;

  /* first clear some flags */
  do {
    e->head.api_flag &= ~tag;

    if (!e->l) {
      continue;
    }

    BMLoop *l = e->l;
    do {
      l->head.hflag &= ~tag;
      l->next->head.hflag &= ~tag;
      l->prev->head.hflag &= ~tag;
    } while ((l = l->radial_next) != e->l);
  } while ((e = BM_DISK_EDGE_NEXT(e, v)) != v->e);

  BMLoop **ls = NULL;
  BLI_array_staticdeclare(ls, 32);

  bool first = true;
  bool bad = false;

  for (int i = 0; i < totuv; i++) {
    snapuvs[i] = true;  //!(mv->flag & SCULPTVERT_UV_BOUNDARY);
  }

  do {
    BMLoop *l = e->l;

    if (!l) {
      continue;
    }
#if 0
    bool bound = l == l->radial_next;

    // check for faceset boundaries
    bound = bound || (BM_ELEM_CD_GET_INT(l->f,ss->cd_faceset_offset) !=
      BM_ELEM_CD_GET_INT(l->radial_next->f,ss->cd_faceset_offset));

    // check for seam and sharp edges
    bound = bound || (e->head.hflag & BM_ELEM_SEAM) || !(e->head.hflag & BM_ELEM_SMOOTH);

    if (bound) {
      continue;
    }
#endif
    do {
      BMLoop *l2 = l->v != v ? l->next : l;

      if (l2->head.hflag & tag) {
        continue;
      }

      l2->head.hflag |= tag;
      BLI_array_append(ls, l2);

      for (int i = 0; i < totuv; i++) {
        const int cd_uv = uvlayer[i].offset;
        float *luv = BM_ELEM_CD_GET_VOID_P(l2, cd_uv);

        // check that we are not part of a uv seam
        if (!first) {
          const float dx = lastuvs[i * 2] - luv[0];
          const float dy = lastuvs[i * 2 + 1] - luv[1];
          const float eps = 0.00001f;

          if (dx * dx + dy * dy > eps) {
            bad = true;
            snapuvs[i] = false;
          }
        }

        lastuvs[i * 2] = luv[0];
        lastuvs[i * 2 + 1] = luv[1];
      }

      first = false;

      if (bad) {
        break;
      }
    } while ((l = l->radial_next) != e->l);

    if (bad) {
      break;
    }
  } while ((e = BM_DISK_EDGE_NEXT(e, v)) != v->e);

  if (bad || !BLI_array_len(ls)) {
    return;
  }

  int totloop = BLI_array_len(ls);

  const float *v_proj_axis = v->no;
  float v_proj[3][3];

  project_plane_normalized_v3_v3v3(v_proj[1], mv->origco, v_proj_axis);

  /* original (l->prev, l, l->next) projections for each loop ('l' remains unchanged) */

  char *_blocks = alloca(ldata->totsize * totloop);
  void **blocks = BLI_array_alloca(blocks, totloop);

  for (int i = 0; i < totloop; i++, _blocks += ldata->totsize) {
    blocks[i] = (void *)_blocks;
  }

  float vco[3], vno[3];

  copy_v3_v3(vco, v->co);
  copy_v3_v3(vno, v->no);

  BMFace _fakef, *fakef = &_fakef;

#if 0
  BMFace *projf = NULL;
  // find face vertex projects into
  for (int i = 0; i < totloop; i++) {
    BMLoop *l = ls[i];

    copy_v3_v3(ray,l->f->no);
    negate_v3(ray);

    float t,uv[2];

    //*
    bool hit = isect_ray_tri_v3(origin,ray,l->prev->v->co,origco,l->next->v->co,&t,uv);
    if (hit) {
      projf = l->f;
      break;
    }  //*/
  }

  if (!projf) {
    return;
  }
#endif

  // build fake f with original coordinates
  for (int i = 0; i < totloop; i++) {
    // create fake face
    BMLoop *l = ls[i];
    float no[3] = {0.0f, 0.0f, 0.0f};

    BMLoop *fakels = BLI_array_alloca(fakels, l->f->len);
    BMVert *fakevs = BLI_array_alloca(fakevs, l->f->len);
    BMLoop *l2 = l->f->l_first;
    BMLoop *fakel = fakels;
    BMVert *fakev = fakevs;
    int j = 0;

    do {
      *fakel = *l2;
      fakel->next = fakels + ((j + 1) % l->f->len);
      fakel->prev = fakels + ((j + l->f->len - 1) % l->f->len);

      *fakev = *l2->v;
      fakel->v = fakev;

      SCULPT_vertex_check_origdata(ss, (PBVHVertRef){.i = (intptr_t)l2->v});

      if (l2->v == v) {
        copy_v3_v3(fakev->co, origco);
        copy_v3_v3(fakev->no, origno);
        add_v3_v3(no, origno);
      }
      else {
        add_v3_v3(no, l2->v->no);
      }

      fakel++;
      fakev++;
      j++;
    } while ((l2 = l2->next) != l->f->l_first);

    *fakef = *l->f;
    fakef->l_first = fakels;

    // set original face normal
    // normalize_v3(no);
    // copy_v3_v3(fakef->no, no);

    // interpolate
    BMLoop _interpl, *interpl = &_interpl;

    MSculptVert saved = *mv;

    *interpl = *l;
    interpl->head.data = blocks[i];
    // memcpy(interpl->head.data, l2->head.data, ldata->totsize);

    BM_loop_interp_from_face(ss->bm, interpl, fakef, false, false);

    *mv = saved;

    CustomData_bmesh_copy_data(&ss->bm->ldata, &ss->bm->ldata, interpl->head.data, &l->head.data);
  }

  int *tots = BLI_array_alloca(tots, totuv);

  for (int i = 0; i < totuv; i++) {
    lastuvs[i * 2] = lastuvs[i * 2 + 1] = 0.0f;
    tots[i] = 0;
  }

  // re-snap uvs
  v = (BMVert *)vertex.i;

  e = v->e;
  do {
    if (!e->l) {
      continue;
    }

    BMLoop *l_iter = e->l;
    do {
      BMLoop *l = l_iter->v != v ? l_iter->next : l_iter;

      for (int i = 0; i < totuv; i++) {
        const int cd_uv = uvlayer[i].offset;
        float *luv = BM_ELEM_CD_GET_VOID_P(l, cd_uv);

        add_v2_v2(lastuvs + i * 2, luv);
        tots[i]++;
      }
    } while ((l_iter = l_iter->radial_next) != e->l);
  } while ((e = BM_DISK_EDGE_NEXT(e, v)) != v->e);

  for (int i = 0; i < totuv; i++) {
    if (tots[i]) {
      mul_v2_fl(lastuvs + i * 2, 1.0f / (float)tots[i]);
    }
  }

  e = v->e;
  do {
    if (!e->l) {
      continue;
    }

    BMLoop *l_iter = e->l;
    do {
      BMLoop *l = l_iter->v != v ? l_iter->next : l_iter;

      for (int i = 0; i < totuv; i++) {
        const int cd_uv = uvlayer[i].offset;
        float *luv = BM_ELEM_CD_GET_VOID_P(l, cd_uv);

        if (snapuvs[i]) {
          copy_v2_v2(luv, lastuvs + i * 2);
        }
      }
    } while ((l_iter = l_iter->radial_next) != e->l);
  } while ((e = BM_DISK_EDGE_NEXT(e, v)) != v->e);

  BLI_array_free(ls);
}

MINLINE float safe_shell_angle_to_dist(const float angle)
{
  float th = cosf(angle);

  if (th == 0.0f) {
    return 10.0f;
  }

  return (UNLIKELY(angle < 1.e-8f)) ? 1.0f : fabsf(1.0f / th);
}

ATTR_NO_OPT static void SCULPT_neighbor_coords_average_interior_boundary(SculptSession *ss,
                                                                         float result[3],
                                                                         PBVHVertRef vertex,
                                                                         SculptSmoothArgs *args)
{
  float avg[3] = {0.0f, 0.0f, 0.0f};

  const float bevel_smooth_factor = 1.0f - args->bevel_smooth_factor;
  float projection = args->projection;
  float slide_fset = args->slide_fset;
  float bound_smooth = args->bound_smooth;
  bool do_origco = args->do_origco;
  SculptAttribute *bound_scl = args->bound_scl;

  MSculptVert *mv = SCULPT_vertex_get_sculptvert(ss, vertex);

  float bound1[3], bound2[3];
  int totbound = 0;

  if (do_origco) {
    SCULPT_vertex_check_origdata(ss, vertex);
  }

  float total = 0.0f;
  int neighbor_count = 0;
  bool check_fsets = args->preserve_fset_boundaries;

  int bflag = SCULPT_BOUNDARY_MESH | SCULPT_BOUNDARY_SHARP;

  slide_fset = MAX2(slide_fset, bound_smooth);

  if (check_fsets) {
    bflag |= SCULPT_BOUNDARY_FACE_SET | SCULPT_BOUNDARY_SEAM | SCULPT_BOUNDARY_UV;
  }

  const eSculptBoundary is_boundary = SCULPT_vertex_is_boundary(ss, vertex, bflag);

  const float *co = do_origco ? mv->origco : SCULPT_vertex_co_get(ss, vertex);
  float no[3];

  PBVH_CHECK_NAN(co);

  if (true || projection > 0.0f) {
    if (do_origco) {
      copy_v3_v3(no, mv->origno);
    }
    else {
      SCULPT_vertex_normal_get(ss, vertex, no);
    }
  }

  float startco[3], startno[3];
  copy_v3_v3(startco, co);
  copy_v3_v3(startno, no);

  const bool weighted = args->do_weighted_smooth && !is_boundary;
  float *areas = NULL;

  eSculptCorner ctype = SCULPT_CORNER_MESH | SCULPT_CORNER_SHARP;
  if (check_fsets) {
    ctype |= SCULPT_CORNER_FACE_SET | SCULPT_CORNER_SEAM | SCULPT_CORNER_UV;
  }

  // bool have_bmesh = ss->bm;

  int val = SCULPT_vertex_valence_get(ss, vertex);
  areas = BLI_array_alloca(areas, val);

  BKE_pbvh_get_vert_face_areas(ss->pbvh, vertex, areas, val);

  /* normalize areas, then apply a 0.25/val floor */

  float totarea = 0.0f;

  for (int i = 0; i < val; i++) {
    totarea += areas[i];
  }

  totarea = totarea != 0.0f ? 1.0f / totarea : 0.0f;

  float df = 0.05f / (float)val;

  for (int i = 0; i < val; i++) {
    areas[i] = (areas[i] * totarea) + df;
  }

  float *b1 = NULL, btot = 0.0f, b1_orig;

  b1 = SCULPT_vertex_attr_get(vertex, bound_scl);
  b1_orig = *b1;
  *b1 = 0.0f;

  float vel[3] = {0.0f, 0.0f, 0.0f};
  int totvel = 0;

  SculptVertexNeighborIter ni;
  SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vertex, ni) {
    MSculptVert *mv2 = SCULPT_vertex_get_sculptvert(ss, ni.vertex);
    const float *co2;

    if (!do_origco || mv2->stroke_id != ss->stroke_id) {
      co2 = SCULPT_vertex_co_get(ss, ni.vertex);
    }
    else {
      co2 = mv2->origco;
    }

    neighbor_count++;

    float tmp[3], w;
    bool ok = false;

    if (weighted) {
      w = areas[ni.i];
    }
    else {
      w = 1.0f;
    }

    bool do_diffuse = false;

    /*use the new edge api if edges are available, if not estimate boundary
      from verts*/

    eSculptBoundary final_boundary = 0;

    if (ni.has_edge) {
      final_boundary = SCULPT_edge_is_boundary(ss, ni.edge, bflag);

#ifdef SCULPT_DIAGONAL_EDGE_MARKS
      if (ss->bm) {
        BMEdge *e = (BMEdge *)ni.edge.i;
        if (!(e->head.hflag & BM_ELEM_DRAW)) {
          neighbor_count--;
          continue;
        }
      }
#endif
    }
    else {
      final_boundary = is_boundary & SCULPT_vertex_is_boundary(ss, ni.vertex, bflag);
    }

    do_diffuse = bound_scl != NULL;

    if (final_boundary) {
      copy_v3_v3(!totbound ? bound1 : bound2, co2);

      totbound++;
    }

    if (is_boundary) {
      /* Boundary vertices use only other boundary vertices.

      This if statement needs to be refactored a bit, it's confusing.

      */

      bool slide = (slide_fset > 0.0f &&
                    (is_boundary &
                     (SCULPT_BOUNDARY_FACE_SET | SCULPT_BOUNDARY_SEAM | SCULPT_BOUNDARY_UV))) ||
                   bound_smooth > 0.0f;
      slide = slide && !final_boundary;

      if (slide) {
        // project non-boundary offset onto boundary normal
        float t[3];

        w *= slide_fset;

        sub_v3_v3v3(t, co2, co);
        madd_v3_v3v3fl(tmp, co, no, dot_v3v3(t, no));
        ok = true;
      }
      else if (final_boundary & is_boundary) {
        copy_v3_v3(tmp, co2);
        ok = true;
        do_diffuse = false;
      }
      else {
        ok = false;
      }
    }
    else {
      copy_v3_v3(tmp, co2);
      ok = true;
    }

    if (do_diffuse && bound_scl && !is_boundary) {
      /*
      simple boundary inflator using an ad-hoc diffusion-based pseudo-geodesic field

      makes more rounded edges.
      */
      copy_v3_v3(tmp, co2);
      ok = true;

      float len = len_v3v3(co, tmp);
      float w2 = 1.0f;

      float *b2 = SCULPT_vertex_attr_get(ni.vertex, bound_scl);
      float b2_val = *b2 + len;

      if (SCULPT_vertex_is_boundary(ss, ni.vertex, bflag)) {
        w2 = 1000.0f;
        b2_val = len;
      }

      *b1 += b2_val * w2;
      btot += w2;

      float no2[3];

      if (!do_origco || mv2->stroke_id != ss->stroke_id) {
        SCULPT_vertex_normal_get(ss, ni.vertex, no2);
      }
      else {
        copy_v3_v3(no2, mv2->origno);
      }

      float radius;
      if (!args->bound_smooth_radius && ss->cache) {
        radius = ss->cache->radius;
      }
      else {
        radius = args->bound_smooth_radius;
      }

      radius = radius == 0.0f ? 0.0001f : radius;

#if 0
      float color[4];
      SCULPT_vertex_color_get(ss,ni.vertex, color);

      color[0] = color[1] = color[2] = th;
      color[3] = 1.0f;

      SCULPT_vertex_color_set(ss, ni.vertex, color);
#endif

      float th = min_ff(b1_orig / radius, bevel_smooth_factor);

      /*smooth bevel edges slightly to avoid artifacts.
        not entire sure why this works.*/
      float shellth = saacos(dot_v3v3(no, no2));
      shellth = safe_shell_angle_to_dist(shellth * 2.0f);
      th /= 0.00001f + shellth;

      sub_v3_v3v3(tmp, co2, co);
      madd_v3_v3fl(tmp, no, -dot_v3v3(no, tmp) * th);
      add_v3_v3(tmp, co);
    }

    if (!ok) {
      continue;
    }

    if (projection > 0.0f) {
      sub_v3_v3(tmp, co);
      float fac = dot_v3v3(tmp, no);
      madd_v3_v3fl(tmp, no, -fac * projection);
      madd_v3_v3fl(avg, tmp, w);
    }
    else {
      madd_v3_v3fl(avg, tmp, w);
    }

    total += w;
  }
  SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

  PBVH_CHECK_NAN(co);
  PBVH_CHECK_NAN(avg);

  if (btot != 0.0f) {
    *b1 /= btot;
  }
  else if (b1) {
    *b1 = b1_orig;
  }

  /*try to prevent shrinkage of smooth closed boundaries like circles*/
  if (totbound == 2) {
    sub_v3_v3(bound1, co);
    sub_v3_v3(bound2, co);
    negate_v3(bound2);

    float bound[3];
    add_v3_v3v3(bound, bound1, bound2);
    float tan[3];
    cross_v3_v3v3(tan, bound, no);
    normalize_v3(tan);

    // project to plane, remember we negated bound2 earlier
    madd_v3_v3fl(avg, tan, -dot_v3v3(bound1, tan) * 0.75);
    madd_v3_v3fl(avg, tan, dot_v3v3(bound2, tan) * 0.75);
  }

  if (args->vel_scl && totvel > 1) {
    float *final_vel = SCULPT_vertex_attr_get(vertex, args->vel_scl);
    mul_v3_fl(vel, 1.0f / (float)totvel);

    interp_v3_v3v3(final_vel, final_vel, vel, args->vel_smooth_fac);
  }

  /* Do not modify corner vertices. */
  if (neighbor_count <= 2 && is_boundary) {
    copy_v3_v3(result, co);
    return;
  }

  /* Avoid division by 0 when there are no neighbors. */
  if (total == 0.0f) {
    copy_v3_v3(result, co);
    return;
  }

  mul_v3_v3fl(result, avg, 1.0f / total);

  PBVH_CHECK_NAN(co);

  if (projection > 0.0f) {
    add_v3_v3(result, co);
  }

  PBVH_CHECK_NAN(co);

  eSculptCorner c = SCULPT_vertex_is_corner(ss, vertex, ctype);
  float corner_smooth;

  if (c == 0) {
    return;
  }

  if (c & (SCULPT_CORNER_FACE_SET | SCULPT_CORNER_SEAM | SCULPT_CORNER_UV)) {
    corner_smooth = MAX2(slide_fset, 2.0f * bound_smooth);
  }
  else {
    corner_smooth = 2.0f * bound_smooth;
  }

  interp_v3_v3v3(result, result, co, 1.0f - corner_smooth);
  PBVH_CHECK_NAN(co);
}

ATTR_NO_OPT void SCULPT_neighbor_coords_average_interior(SculptSession *ss,
                                                         float result[3],
                                                         PBVHVertRef vertex,
                                                         SculptSmoothArgs *args)
{
  if (args->bound_smooth > 0.0f && args->bound_scl) {
    SCULPT_neighbor_coords_average_interior_boundary(ss, result, vertex, args);
    return;
  }

  float avg[3] = {0.0f, 0.0f, 0.0f};

  float projection = args->projection;
  float slide_fset = args->slide_fset;
  float bound_smooth = args->bound_smooth;
  bool do_origco = args->do_origco;

  MSculptVert *mv = SCULPT_vertex_get_sculptvert(ss, vertex);

  float bound1[3], bound2[3];
  int totbound = 0;

  if (do_origco) {
    SCULPT_vertex_check_origdata(ss, vertex);
  }

  float total = 0.0f;
  int neighbor_count = 0;
  bool check_fsets = args->preserve_fset_boundaries;

  int bflag = SCULPT_BOUNDARY_MESH | SCULPT_BOUNDARY_SHARP;

  slide_fset = MAX2(slide_fset, bound_smooth);

  if (check_fsets) {
    bflag |= SCULPT_BOUNDARY_FACE_SET | SCULPT_BOUNDARY_SEAM | SCULPT_BOUNDARY_UV;
  }

  const eSculptBoundary is_boundary = SCULPT_vertex_is_boundary(ss, vertex, bflag);

  const float *co = do_origco ? mv->origco : SCULPT_vertex_co_get(ss, vertex);
  float no[3];

  PBVH_CHECK_NAN(co);

  if (do_origco) {
    copy_v3_v3(no, mv->origno);
  }
  else {
    SCULPT_vertex_normal_get(ss, vertex, no);
  }

  float startco[3], startno[3];
  copy_v3_v3(startco, co);
  copy_v3_v3(startno, no);

  const bool weighted = args->do_weighted_smooth && !is_boundary;
  float *areas = NULL;

  eSculptCorner ctype = SCULPT_CORNER_MESH | SCULPT_CORNER_SHARP;
  if (check_fsets) {
    ctype |= SCULPT_CORNER_FACE_SET | SCULPT_CORNER_SEAM | SCULPT_CORNER_UV;
  }

  if (weighted) {
    int val = SCULPT_vertex_valence_get(ss, vertex);
    areas = BLI_array_alloca(areas, val);

    BKE_pbvh_get_vert_face_areas(ss->pbvh, vertex, areas, val);

    /* normalize areas, then apply a 0.25/val floor */

    float totarea = 0.0f;
    for (int i = 0; i < val; i++) {
      totarea += areas[i];
    }

    totarea = totarea > 0.0f ? 1.0f / totarea : 0.0f;

    float df = 0.25f / (float)val;

    for (int i = 0; i < val; i++) {
      areas[i] = (areas[i] * totarea) + df;
    }
  }

  float vel[3] = {0.0f, 0.0f, 0.0f};
  float totvel = 0.0f;

  SculptVertexNeighborIter ni;
  SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vertex, ni) {
    MSculptVert *mv2 = SCULPT_vertex_get_sculptvert(ss, ni.vertex);
    const float *co2;
    float w;

    if (weighted) {
      w = areas[ni.i];
    }
    else {
      w = 1.0f;
    }

    if (args->vel_scl) {
      /* propagate velocities */
      float *vel2 = SCULPT_vertex_attr_get(ni.vertex, args->vel_scl);
      madd_v3_v3fl(vel, vel2, w);
      totvel += w;
    }

    if (!do_origco || mv2->stroke_id != ss->stroke_id) {
      co2 = SCULPT_vertex_co_get(ss, ni.vertex);
    }
    else {
      co2 = mv2->origco;
    }

    neighbor_count++;

    float tmp[3];
    bool ok = false;

    /* use the new edge api if edges are available, if not estimate boundary
       from verts
     */
    eSculptBoundary final_boundary = 0;

    if (ni.has_edge) {
      final_boundary = SCULPT_edge_is_boundary(ss, ni.edge, bflag);

#ifdef SCULPT_DIAGONAL_EDGE_MARKS
      if (ss->bm) {
        BMEdge *e = (BMEdge *)ni.edge.i;
        if (!(e->head.hflag & BM_ELEM_DRAW)) {
          neighbor_count--;
          continue;
        }
      }
#endif
    }
    else {
      final_boundary = is_boundary & SCULPT_vertex_is_boundary(ss, ni.vertex, bflag);
    }

    if (final_boundary) {
      copy_v3_v3(!totbound ? bound1 : bound2, co2);
      totbound++;
    }

    if (is_boundary) {
      /*
      Boundary rules:

      Hard edges: Boundary vertices use only other boundary vertices.
      Slide: Boundary vertices use normal component of non-boundary vertices
      */

      bool slide = slide_fset > 0.0f &&
                   (is_boundary &
                    (SCULPT_BOUNDARY_FACE_SET | SCULPT_BOUNDARY_SEAM | SCULPT_BOUNDARY_UV));
      slide = slide && !final_boundary;

      if (slide) {
        /* project non-boundary offset onto boundary normal*/
        float t[3];

        w *= slide_fset;

        sub_v3_v3v3(t, co2, co);
        madd_v3_v3v3fl(tmp, co, no, dot_v3v3(t, no));
        ok = true;
      }
      else if (final_boundary & is_boundary) {
        copy_v3_v3(tmp, co2);
        ok = true;
      }
      else {
        ok = false;
      }
    }
    else {
      copy_v3_v3(tmp, co2);
      ok = true;
    }

    if (!ok) {
      continue;
    }

    if (projection > 0.0f) {
      sub_v3_v3(tmp, co);
      float fac = dot_v3v3(tmp, no);
      madd_v3_v3fl(tmp, no, -fac * projection);
      madd_v3_v3fl(avg, tmp, w);
    }
    else {
      madd_v3_v3fl(avg, tmp, w);
    }

    total += w;
  }
  SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

  PBVH_CHECK_NAN(co);
  PBVH_CHECK_NAN(avg);

  /*try to prevent shrinkage of smooth closed boundaries like circles*/
  if (totbound == 2) {
    /* find tangent to boundary */
    sub_v3_v3(bound1, co);
    sub_v3_v3(bound2, co);

    negate_v3(bound2);

    float bound[3];
    add_v3_v3v3(bound, bound1, bound2);
    float tan[3];

    cross_v3_v3v3(tan, bound, no);
    normalize_v3(tan);

    float w = total / (float)neighbor_count;

    /* project to plane, remember we negated bound2 earlier */
    madd_v3_v3fl(avg, tan, -w * dot_v3v3(bound1, tan) * 0.75);
    madd_v3_v3fl(avg, tan, w * dot_v3v3(bound2, tan) * 0.75);
  }

  if (args->vel_scl && totvel != 0.0f) {
    float *final_vel = SCULPT_vertex_attr_get(vertex, args->vel_scl);
    mul_v3_fl(vel, 1.0f / totvel);

    interp_v3_v3v3(final_vel, final_vel, vel, args->vel_smooth_fac);
  }

  /* Do not modify corner vertices. */
  if (total == 0.0f || (neighbor_count <= 2 && is_boundary)) {
    copy_v3_v3(result, co);
    return;
  }

  mul_v3_v3fl(result, avg, 1.0f / total);

  PBVH_CHECK_NAN(co);

  if (projection > 0.0f) {
    add_v3_v3(result, co);
  }

  PBVH_CHECK_NAN(co);

  eSculptCorner c = SCULPT_vertex_is_corner(ss, vertex, ctype);

  if (!c) {
    return;
  }

  if (c & (SCULPT_CORNER_FACE_SET | SCULPT_CORNER_SEAM | SCULPT_CORNER_UV)) {
    interp_v3_v3v3(result, result, co, 1.0f - slide_fset);
  }

  PBVH_CHECK_NAN(result);
}

int closest_vec_to_perp(float dir[3], float r_dir2[3], float no[3], float *buckets, float w)
{
  int bits = 0;

  if (dot_v3v3(r_dir2, dir) < 0.0f) {
    negate_v3(r_dir2);
    bits |= 1;
  }

  float dir4[3];
  cross_v3_v3v3(dir4, r_dir2, no);
  normalize_v3(dir4);

  if (dot_v3v3(dir4, dir) < 0.0f) {
    negate_v3(dir4);
    bits |= 2;
  }

  if (dot_v3v3(dir4, dir) > dot_v3v3(r_dir2, dir)) {
    copy_v3_v3(r_dir2, dir4);
    bits |= 4;
  }

  buckets[bits] += w;

  return bits;
}

void vec_transform(float r_dir2[3], float no[3], int bits)
{
  if (bits & 4) {
    float dir4[3];

    copy_v3_v3(dir4, r_dir2);

    if (bits & 2) {
      negate_v3(dir4);
    }

    float dir5[3];

    cross_v3_v3v3(dir5, no, dir4);
    normalize_v3(dir5);

    copy_v3_v3(r_dir2, dir5);
  }

  if (bits & 1) {
    negate_v3(r_dir2);
  }
}

/* For bmesh: Average surrounding verts based on an orthogonality measure.
 * Naturally converges to a quad-like structure. Builds a simple cross
 * field combining either mesh curvature principle direction or
   stroke direction, and edge boundaries. */
void SCULPT_bmesh_four_neighbor_average(SculptSession *ss,
                                        float avg[3],
                                        float direction[3],
                                        BMVert *v,
                                        float projection,
                                        bool check_fsets,
                                        int cd_temp,
                                        int cd_sculpt_vert,
                                        bool do_origco)
{
  float avg_co[3] = {0.0f, 0.0f, 0.0f};
  float tot_co = 0.0f;

  float buckets[8] = {0};

  // zero_v3(direction);

  MSculptVert *mv = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, v);

  float *col = BM_ELEM_CD_GET_VOID_P(v, cd_temp);
  float dir[3];
  float dir3[3] = {0.0f, 0.0f, 0.0f};

  const bool weighted = (ss->cache->brush->flag2 & BRUSH_SMOOTH_USE_AREA_WEIGHT);
  float *areas;

  SCULPT_vertex_check_origdata(ss, (PBVHVertRef){.i = (intptr_t)v});

  if (do_origco) {
    // SCULPT_vertex_check_origdata(ss, (PBVHVertRef){.i = (intptr_t)v});
    madd_v3_v3fl(direction, mv->origno, -dot_v3v3(mv->origno, direction));
    normalize_v3(direction);
  }

  float *co1 = do_origco ? mv->origco : v->co;
  float *no1 = do_origco ? mv->origno : v->no;

  if (weighted) {
    PBVHVertRef vertex = {.i = (intptr_t)v};

    int val = SCULPT_vertex_valence_get(ss, vertex);
    areas = BLI_array_alloca(areas, val * 2);

    BKE_pbvh_get_vert_face_areas(ss->pbvh, vertex, areas, val);
    float totarea = 0.0f;

    for (int i = 0; i < val; i++) {
      totarea += areas[i];
    }

    totarea = totarea != 0.0f ? 1.0f / totarea : 0.0f;

    for (int i = 0; i < val; i++) {
      areas[i] *= totarea;
    }
  }

  copy_v3_v3(dir, col);

  if (dot_v3v3(dir, dir) == 0.0f) {
    copy_v3_v3(dir, direction);
  }
  else {
    closest_vec_to_perp(dir, direction, no1, buckets, 1.0f);  // col[3]);
  }

  float totdir3 = 0.0f;

  const float selfw = (float)mv->valence * 0.0025f;
  madd_v3_v3fl(dir3, direction, selfw);

  totdir3 += selfw;

  BMIter eiter;
  BMEdge *e;
  bool had_bound = false;
  int area_i = 0;

  BM_ITER_ELEM_INDEX (e, &eiter, v, BM_EDGES_OF_VERT, area_i) {
    BMVert *v_other = (e->v1 == v) ? e->v2 : e->v1;

    float dir2[3];
    float *col2 = BM_ELEM_CD_GET_VOID_P(v_other, cd_temp);

    float bucketw = 1.0f;

    MSculptVert *mv2 = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, v_other);
    float *co2;

    if (!do_origco || mv2->stroke_id != ss->stroke_id) {
      co2 = v_other->co;
    }
    else {
      co2 = mv2->origco;
    }

    eSculptBoundary bflag = SCULPT_BOUNDARY_FACE_SET | SCULPT_BOUNDARY_MESH |
                            SCULPT_BOUNDARY_SHARP | SCULPT_BOUNDARY_SEAM | SCULPT_BOUNDARY_UV;

    int bound = SCULPT_edge_is_boundary(ss, (PBVHEdgeRef){.i = (intptr_t)e}, bflag);
    float dirw = 1.0f;

    if (bound) {
      had_bound = true;

      sub_v3_v3v3(dir2, co2, co1);
      madd_v3_v3fl(dir2, no1, -dot_v3v3(no1, dir2));
      normalize_v3(dir2);
      dirw = 100000.0f;
    }
    else {
      dirw = col2[3];

      copy_v3_v3(dir2, col2);
      if (dot_v3v3(dir2, dir2) == 0.0f) {
        copy_v3_v3(dir2, dir);
      }
    }

    closest_vec_to_perp(dir, dir2, no1, buckets, bucketw);  // col2[3]);

    madd_v3_v3fl(dir3, dir2, dirw);
    totdir3 += dirw;

    if (had_bound) {
      tot_co = 0.0f;
      continue;
    }
    float vec[3];
    sub_v3_v3v3(vec, co2, co1);

    madd_v3_v3fl(vec, no1, -dot_v3v3(vec, no1) * projection);
    normalize_v3(vec);

    /* fac is a measure of how orthogonal or parallel the edge is
     * relative to the direction. */
    float fac = dot_v3v3(vec, dir);
#ifdef SCULPT_DIAGONAL_EDGE_MARKS
    float th = fabsf(saacos(fac)) / M_PI + 0.5f;
    th -= floorf(th);

    const float limit = 0.045;

    if (fabsf(th - 0.25) < limit || fabsf(th - 0.75) < limit) {
      BMEdge enew = *e, eold = *e;

      enew.head.hflag &= ~BM_ELEM_DRAW;
      // enew.head.hflag |= BM_ELEM_SEAM;  // XXX debug

      atomic_cas_int64((intptr_t *)(&e->head.index),
                       *(intptr_t *)(&eold.head.index),
                       *(intptr_t *)(&enew.head.index));
    }
#endif

    fac = fac * fac - 0.5f;
    fac *= fac;

    if (weighted) {
      fac *= areas[area_i];
    }

    madd_v3_v3fl(avg_co, co2, fac);
    tot_co += fac;
  }

  /* In case vert has no Edge s. */
  if (tot_co > 0.0f) {
    mul_v3_v3fl(avg, avg_co, 1.0f / tot_co);

    /* Preserve volume. */
    float vec[3];
    sub_v3_v3(avg, co1);
    mul_v3_v3fl(vec, no1, dot_v3v3(avg, no1) * projection);
    sub_v3_v3(avg, vec);
    add_v3_v3(avg, co1);
  }
  else {
    // zero_v3(avg);
    copy_v3_v3(avg, co1);
  }

  PBVH_CHECK_NAN(avg);

  // do not update in do_origco
  if (do_origco) {
    return;
  }

  if (totdir3 > 0.0f) {
    float outdir = totdir3 / (float)mv->valence;

    // mul_v3_fl(dir3, 1.0 / totdir3);
    normalize_v3(dir3);
    if (had_bound) {
      copy_v3_v3(col, dir3);
      col[3] = 1000.0f;
    }
    else {

      mul_v3_fl(col, col[3]);
      madd_v3_v3fl(col, dir3, outdir);

      col[3] = (col[3] + outdir) * 0.4;
      normalize_v3(col);
    }

    float maxb = 0.0f;
    int bi = 0;
    for (int i = 0; i < 8; i++) {
      if (buckets[i] > maxb) {
        maxb = buckets[i];
        bi = i;
      }
    }

    // negate_v3(col);
    vec_transform(col, no1, bi);
    // negate_v3(col);
  }
}

static void sculpt_neighbor_coords_average_fset(
    SculptSession *ss, float result[3], PBVHVertRef vertex, float projection, bool weighted)
{
  float avg[3] = {0.0f, 0.0f, 0.0f};
  float *co, no[3];
  float total = 0.0f;

  bool boundary = !SCULPT_vertex_has_unique_face_set(ss, vertex);

  if (projection > 0.0f) {
    co = (float *)SCULPT_vertex_co_get(ss, vertex);
    SCULPT_vertex_normal_get(ss, vertex, no);
  }

  float *areas;

  if (weighted) {
    int val = SCULPT_vertex_valence_get(ss, vertex);
    areas = BLI_array_alloca(areas, val);

    BKE_pbvh_get_vert_face_areas(ss->pbvh, vertex, areas, val);
  }

  SculptVertexNeighborIter ni;
  SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vertex, ni) {
    const float *co2 = SCULPT_vertex_co_get(ss, ni.vertex);
    float w;

    if (weighted) {
      w = areas[ni.i];
    }
    else {
      w = 1.0f;
    }

    if (boundary && SCULPT_vertex_has_unique_face_set(ss, ni.vertex)) {
      continue;
    }

    if (projection > 0.0f) {
      float tmp[3];

      sub_v3_v3v3(tmp, co2, co);
      float fac = dot_v3v3(tmp, no);
      madd_v3_v3fl(tmp, no, -fac * projection);

      madd_v3_v3fl(avg, tmp, w);
    }
    else {
      madd_v3_v3fl(avg, co2, w);
    }
    total += w;
  }
  SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

  if (total > (boundary ? 1.0f : 0.0f)) {
    mul_v3_v3fl(result, avg, 1.0f / total);

    if (projection > 0.0) {
      add_v3_v3(result, co);
    }
  }
  else {
    copy_v3_v3(result, SCULPT_vertex_co_get(ss, vertex));
  }
}
/* Generic functions for laplacian smoothing. These functions do not take boundary vertices into
 * account. */

void SCULPT_neighbor_coords_average(SculptSession *ss,
                                    float result[3],
                                    PBVHVertRef vertex,
                                    float projection,
                                    bool check_fsets,
                                    bool weighted)
{
  if (check_fsets) {
    // sculpt_neighbor_coords_average_fset(ss, result, vertex, projection, weighted);
    // return;
  }

  float avg[3] = {0.0f, 0.0f, 0.0f};
  float *co, no[3];
  float total = 0.0f;
  int bound_mask = SCULPT_BOUNDARY_ALL;

  if (!check_fsets) {
    bound_mask &= ~SCULPT_BOUNDARY_FACE_SET;
  }

  bool boundary = SCULPT_vertex_is_boundary(ss, vertex, bound_mask);

  co = (float *)SCULPT_vertex_co_get(ss, vertex);
  SCULPT_vertex_normal_get(ss, vertex, no);

  float *areas;

  if (weighted) {
    int val = SCULPT_vertex_valence_get(ss, vertex);
    areas = BLI_array_alloca(areas, val);

    BKE_pbvh_get_vert_face_areas(ss->pbvh, vertex, areas, val);
  }

  SculptVertexNeighborIter ni;
  SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vertex, ni) {
    const float *co2 = SCULPT_vertex_co_get(ss, ni.vertex);
    float w;

    if (weighted) {
      w = areas[ni.i];
    }
    else {
      w = 1.0f;
    }

    if (boundary) {
      bool boundary2 = SCULPT_vertex_is_boundary(ss, ni.vertex, bound_mask);

      if (!boundary2 && (boundary & SCULPT_BOUNDARY_FACE_SET)) {
        float tmp[3];

        sub_v3_v3v3(tmp, co2, co);
        madd_v3_v3fl(avg, no, dot_v3v3(tmp, no) * w);
        total += w;
        continue;
      }
      else if (!boundary2) {
        continue;
      }
    }

    if (projection > 0.0f) {
      float tmp[3];

      sub_v3_v3v3(tmp, co2, co);
      float fac = dot_v3v3(tmp, no);
      madd_v3_v3fl(tmp, no, -fac * projection);

      madd_v3_v3fl(avg, tmp, w);
    }
    else {
      madd_v3_v3fl(avg, co2, w);
    }
    total += w;
  }
  SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

  if (total > 0.0f) {
    mul_v3_v3fl(result, avg, 1.0f / total);

    if (projection > 0.0) {
      add_v3_v3(result, co);
    }
  }
  else {
    copy_v3_v3(result, SCULPT_vertex_co_get(ss, vertex));
  }
}

float SCULPT_neighbor_mask_average(SculptSession *ss, PBVHVertRef vertex)
{
  float avg = 0.0f;
  int total = 0;

  SculptVertexNeighborIter ni;
  SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vertex, ni) {
    avg += SCULPT_vertex_mask_get(ss, ni.vertex);
    total++;
  }
  SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

  if (total > 0) {
    return avg / total;
  }
  return SCULPT_vertex_mask_get(ss, vertex);
}

void SCULPT_neighbor_color_average(SculptSession *ss, float result[4], PBVHVertRef vertex)
{
  float avg[4] = {0.0f, 0.0f, 0.0f, 0.0f};
  float total = 0.0f;

  AutomaskingCache *automasking = SCULPT_automasking_active_cache_get(ss);

  AutomaskingNodeData automask_data = {0};
  automask_data.have_orig_data = true;
  automask_data.orig_data.co = SCULPT_vertex_origco_get(ss, vertex);
  automask_data.orig_data.no = SCULPT_vertex_origno_get(ss, vertex);

  SculptVertexNeighborIter ni;
  SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vertex, ni) {
    float tmp[4] = {0};

    float w = automasking ?
                  SCULPT_automasking_factor_get(automasking, ss, ni.vertex, &automask_data) :
                  1.0f;

    SCULPT_vertex_color_get(ss, ni.vertex, tmp);

    madd_v4_v4fl(avg, tmp, w);
    total += w;
  }
  SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

  if (total > 0.0f) {
    mul_v4_v4fl(result, avg, 1.0f / total);
  }
  else {
    SCULPT_vertex_color_get(ss, vertex, result);
  }
}

static void do_enhance_details_brush_task_cb_ex(void *__restrict userdata,
                                                const int n,
                                                const TaskParallelTLS *__restrict tls)
{
  SculptThreadedTaskData *data = userdata;
  SculptSession *ss = data->ob->sculpt;
  Sculpt *sd = data->sd;
  const Brush *brush = data->brush;

  PBVHVertexIter vd;

  float bstrength = ss->cache->bstrength;
  CLAMP(bstrength, -1.0f, 1.0f);

  SculptBrushTest test;
  SculptBrushTestFn sculpt_brush_test_sq_fn = SCULPT_brush_test_init(
      ss, &test, data->brush->falloff_shape);

  const int thread_id = BLI_task_parallel_thread_id(tls);
  AutomaskingNodeData automask_data;
  SCULPT_automasking_node_begin(
      data->ob, ss, ss->cache->automasking, &automask_data, data->nodes[n]);

  BKE_pbvh_vertex_iter_begin (ss->pbvh, data->nodes[n], vd, PBVH_ITER_UNIQUE) {
    if (!sculpt_brush_test_sq_fn(&test, vd.co)) {
      continue;
    }

    SCULPT_automasking_node_update(ss, &automask_data, &vd);

    const float fade = bstrength * SCULPT_brush_strength_factor(ss,
                                                                brush,
                                                                vd.co,
                                                                sqrtf(test.dist),
                                                                vd.no,
                                                                vd.fno,
                                                                vd.mask ? *vd.mask : 0.0f,
                                                                vd.vertex,
                                                                thread_id,
                                                                &automask_data);

    float disp[3];
    float *dir = SCULPT_vertex_attr_get(vd.vertex, data->scl);

    madd_v3_v3v3fl(disp, vd.co, dir, fade);
    SCULPT_clip(sd, ss, vd.co, disp);

    if (vd.is_mesh) {
      BKE_pbvh_vert_tag_update_normal(ss->pbvh, vd.vertex);
    }
  }
  BKE_pbvh_vertex_iter_end;
}

static void do_enhance_details_brush_dir_task_cb_ex(void *__restrict userdata,
                                                    const int n,
                                                    const TaskParallelTLS *__restrict tls)
{
  SculptThreadedTaskData *data = userdata;
  SculptSession *ss = data->ob->sculpt;

  PBVHVertexIter vd;

  SculptBrushTest test;
  SculptBrushTestFn sculpt_brush_test_sq_fn = SCULPT_brush_test_init(
      ss, &test, data->brush->falloff_shape);

  SculptAttribute *scl = data->scl;

  bool modified = false;

  BKE_pbvh_vertex_iter_begin (ss->pbvh, data->nodes[n], vd, PBVH_ITER_UNIQUE) {
    if (!sculpt_brush_test_sq_fn(&test, vd.co)) {
      continue;
    }

    if (SCULPT_stroke_id_test(ss, vd.vertex, STROKEID_USER_SMOOTH)) {
      modified = true;

      float avg[3];
      float *dir = SCULPT_vertex_attr_get(vd.vertex, scl);
      float no[3];

      SCULPT_vertex_normal_get(ss, vd.vertex, no);

      SCULPT_neighbor_coords_average(ss, avg, vd.vertex, 0.0f, false, data->use_area_cos);
      sub_v3_v3v3(dir, avg, SCULPT_vertex_co_get(ss, vd.vertex));

      /* get rid of tangential displacement */
      float fac = dot_v3v3(dir, no);
      copy_v3_v3(dir, no);
      mul_v3_fl(dir, -fac);
    }
  }
  BKE_pbvh_vertex_iter_end;

  /* tell main thread we modified something */
  if (modified) {
    atomic_add_and_fetch_uint32((uint32_t *)(&data->cd_temp), 1UL);
  }
}

/* Diffuses detail directions */
static void do_enhance_details_brush_dir2_task_cb_ex(void *__restrict userdata,
                                                     const int n,
                                                     const TaskParallelTLS *__restrict tls)
{
  SculptThreadedTaskData *data = userdata;
  SculptSession *ss = data->ob->sculpt;

  PBVHVertexIter vd;

  SculptBrushTest test;
  SculptBrushTestFn sculpt_brush_test_sq_fn = SCULPT_brush_test_init(
      ss, &test, data->brush->falloff_shape);

  // SculptAttribute *strokeid_scl = data->scl2;
  SculptAttribute *scl = data->scl;
  bool use_area_weights = data->use_area_cos;

  int lastvalence = 8;
  float *areas = MEM_malloc_arrayN(lastvalence, sizeof(float), __func__);

  BKE_pbvh_vertex_iter_begin (ss->pbvh, data->nodes[n], vd, PBVH_ITER_UNIQUE) {
    if (!sculpt_brush_test_sq_fn(&test, vd.co)) {
      continue;
    }

    int valence;

    if (use_area_weights) {
      valence = SCULPT_vertex_valence_get(ss, vd.vertex);

      if (valence > lastvalence) {
        areas = MEM_reallocN(areas, sizeof(float) * valence);
      }

      BKE_pbvh_get_vert_face_areas(ss->pbvh, vd.vertex, areas, valence);
    }

    /* this check here is overly restrictive,
       we already get filtered by whether stage
       1 did anything */
    if (1) {  //*strokeid == current_stroke_id << 1UL) {
      // if (data->cd_temp2) {
      //(*strokeid)++;
      //}

      float avg[3] = {0.0f, 0.0f, 0.0f};
      float tot = 0.0f;

      SculptVertexNeighborIter ni;
      SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vd.vertex, ni) {
        float *dir2 = SCULPT_vertex_attr_get(ni.vertex, scl);

        float w = 1.0f;

        if (use_area_weights) {
          w = areas[ni.i];
        }

        madd_v3_v3fl(avg, dir2, w);
        tot += w;
      }
      SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

      if (tot == 0.0f) {
        continue;
      }

      mul_v3_fl(avg, 1.0f / (float)tot);

      float *dir = SCULPT_vertex_attr_get(vd.vertex, scl);
      interp_v3_v3v3(dir, dir, avg, 0.5f);  // valence == 1 ? 0.5f : 0.75f);
    }
  }
  BKE_pbvh_vertex_iter_end;

  MEM_freeN(areas);
}

void SCULPT_enhance_details_brush(
    Sculpt *sd, Object *ob, PBVHNode **nodes, const int totnode, int presteps)
{
  SculptSession *ss = ob->sculpt;
  Brush *brush = BKE_paint_brush(&sd->paint);
  SculptAttribute *strokeid_scl;

  bool use_area_weights = (ss->cache->brush->flag2 & BRUSH_SMOOTH_USE_AREA_WEIGHT);

  if (SCULPT_stroke_is_first_brush_step(ss->cache)) {
    SCULPT_vertex_random_access_ensure(ss);
  }

  if (use_area_weights) {
    if (SCULPT_stroke_is_first_brush_step(ss->cache)) {
      BKE_pbvh_update_all_tri_areas(ss->pbvh);
      PBVHNode **nodes;
      int totnode;

      BKE_pbvh_get_nodes(ss->pbvh, PBVH_Leaf, &nodes, &totnode);
      for (int i = 0; i < totnode; i++) {
        BKE_pbvh_check_tri_areas(ss->pbvh, nodes[i]);
      }
    }
    else {
      BKE_pbvh_face_areas_begin(ss->pbvh);
    }
  }

  SCULPT_boundary_info_ensure(ob);

  SculptAttribute *scl;
  SculptAttributeParams params = {.permanent = false, .simple_array = false, .stroke_only = true};
  bool weighted = SCULPT_get_int(ss, use_weighted_smooth, sd, brush);

  scl = BKE_sculpt_attribute_ensure(
      ob, ATTR_DOMAIN_POINT, CD_PROP_FLOAT3, "__dyntopo_detail_dir", &params);
  strokeid_scl = SCULPT_stroke_id_attribute_ensure(ob);

  if (SCULPT_stroke_is_first_brush_step(ss->cache)) {
    SCULPT_vertex_random_access_ensure(ss);

#if 0
    const int totvert = SCULPT_vertex_count_get(ss);
    ss->cache->detail_directions = MEM_malloc_arrayN(
        totvert, sizeof(float[3]), "details directions");

    for (int i = 0; i < totvert; i++) {
      PBVHVertRef vertex = BKE_pbvh_index_to_vertex(ss->pbvh, i);

      float avg[3];
      PBVHVertRef vertex = BKE_pbvh_index_to_vertex(ss->pbvh, i);
      float *dir = SCULPT_vertex_attr_get(vertex, &scl);
      float no[3];

      SCULPT_vertex_normal_get(ss, vertex, no);

      SCULPT_neighbor_coords_average(ss, avg, vertex, 0.0f, false, weighted);
      sub_v3_v3v3(dir, avg, SCULPT_vertex_co_get(ss, vertex));

      /* get rid of tangential displacement */
      float fac = dot_v3v3(dir, no);
      copy_v3_v3(dir, no);
      mul_v3_fl(dir, -fac);
    }

    int lastvalence = 5;
    float *areas = MEM_malloc_arrayN(lastvalence, sizeof(float), __func__);

    /* smooth offsets */
    for (int step = 0; step < presteps; step++) {
      for (int i = 0; i < totvert; i++) {
        float avg[3] = {0.0f, 0.0f, 0.0f};

        PBVHVertRef vertex = BKE_pbvh_index_to_vertex(ss->pbvh, i);
        float *dir = SCULPT_vertex_attr_get(vertex, &scl);
        float tot = 0.0f;
        float *areas = NULL;
        int valence;

        if (use_area_weights) {
          valence = SCULPT_vertex_valence_get(ss, vertex);

          if (valence > lastvalence) {
            areas = MEM_reallocN(areas, sizeof(float) * valence);
          }

          areas = MEM_callocN(sizeof(float) * valence * 4, "sdaf");
          // BLI_array_alloca(areas, valence * 4);

          BKE_pbvh_get_vert_face_areas(ss->pbvh, vertex, areas, valence);
        }

        SculptVertexNeighborIter ni;
        SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vertex, ni) {
          float *dir2 = SCULPT_vertex_attr_get(ni.vertex, &scl);

          float w = 1.0f;

          if (use_area_weights) {
            w = areas[ni.i];
          }

          madd_v3_v3fl(avg, dir2, w);
          tot += w;
        }
        SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

        if (tot > 0.0f) {
          mul_v3_fl(avg, 1.0f / tot);
          interp_v3_v3v3(dir, dir, avg, 0.75f);
        }
      }
    }
    MEM_freeN(areas);
#endif
  }

  SculptThreadedTaskData data = {.sd = sd,
                                 .ob = ob,
                                 .brush = brush,
                                 .use_area_cos = weighted,
                                 .cd_temp = 0,
                                 .cd_temp2 = 0,
                                 .nodes = nodes,
                                 .scl = scl};

  TaskParallelSettings settings;
  BKE_pbvh_parallel_range_settings(&settings, true, totnode);

  BLI_task_parallel_range(0, totnode, &data, do_enhance_details_brush_dir_task_cb_ex, &settings);
  if (data.cd_temp) { /* did something change? if so run diffuse task*/
    for (int i = 0; i < presteps; i++) {
      if (i == presteps - 1) {
        data.cd_temp2 = 1;  // flag tasks to update to final stroke id
      }

      BLI_task_parallel_range(
          0, totnode, &data, do_enhance_details_brush_dir2_task_cb_ex, &settings);
    }
  }

  BLI_task_parallel_range(0, totnode, &data, do_enhance_details_brush_task_cb_ex, &settings);
}

//#  define SMOOTH_ITER_IN_THREADS
static void do_smooth_brush_task_cb_ex(void *__restrict userdata,
                                       const int n,
                                       const TaskParallelTLS *__restrict tls)
{
  SculptThreadedTaskData *data = userdata;
  SculptSession *ss = data->ob->sculpt;
  Sculpt *sd = data->sd;
  const float bound_smooth = data->bound_smooth;
  const float fset_slide = data->fset_slide;
  const Brush *brush = data->brush;
  const bool smooth_mask = data->smooth_mask;
  float bstrength = data->strength;
  float projection = data->smooth_projection;

  if (!data->nodes[n]) {
    return;
  }

  PBVHVertexIter vd;

  CLAMP(bstrength, 0.0f, 1.0f);

  SculptBrushTest test;
  SculptBrushTestFn sculpt_brush_test_sq_fn = SCULPT_brush_test_init(
      ss, &test, data->brush->falloff_shape);

  const int thread_id = BLI_task_parallel_thread_id(tls);
  const bool weighted = ss->cache->brush->flag2 & BRUSH_SMOOTH_USE_AREA_WEIGHT;
  const bool check_fsets = ss->cache->brush->flag2 & BRUSH_SMOOTH_PRESERVE_FACE_SETS;

  AutomaskingNodeData automask_data;
  SCULPT_automasking_node_begin(
      data->ob, ss, ss->cache->automasking, &automask_data, data->nodes[n]);

  if (weighted || ss->cache->brush->boundary_smooth_factor > 0.0f) {
    BKE_pbvh_check_tri_areas(ss->pbvh, data->nodes[n]);
  }

  bool modified = false;
  // const float bound_smooth = powf(ss->cache->brush->boundary_smooth_factor,
  // BOUNDARY_SMOOTH_EXP);
  // const float slide_fset = BKE_brush_fset_slide_get(ss->scene, ss->cache->brush);

  SculptAttribute *bound_scl = data->scl2;
  SculptAttribute *vel_scl = data->scl;
#ifdef SMOOTH_ITER_IN_THREADS
  for (int iteration = 0; iteration < data->iterations; iteration++) {
    if (!data->nodes[n]) {
      break;
    }

    if (iteration == data->iterations - 1) {
      bstrength = data->strength;
    }
    else {
      bstrength = 1.0f;
    }
#endif

    BKE_pbvh_vertex_iter_begin (ss->pbvh, data->nodes[n], vd, PBVH_ITER_UNIQUE) {
      if (!sculpt_brush_test_sq_fn(&test, vd.co)) {
        continue;
      }

      modified = true;

      SCULPT_automasking_node_update(ss, &automask_data, &vd);

      // check origdata to be sure we don't mess it up
      SCULPT_vertex_check_origdata(ss, vd.vertex);

      const float fade = bstrength * SCULPT_brush_strength_factor(
                                         ss,
                                         brush,
                                         vd.co,
                                         sqrtf(test.dist),
                                         vd.no,
                                         vd.fno,
                                         smooth_mask ? 0.0f : (vd.mask ? *vd.mask : 0.0f),
                                         vd.vertex,
                                         thread_id,
                                         &automask_data);
      if (smooth_mask) {
        float val = SCULPT_neighbor_mask_average(ss, vd.vertex) - *vd.mask;
        val *= fade * bstrength;
        *vd.mask += val;
        CLAMP(*vd.mask, 0.0f, 1.0f);
      }
      else {
        float avg[3], off[3];

        // if (SCULPT_vertex_is_corner(ss, vd.vertex, ctype) & ~SCULPT_CORNER_FACE_SET) {
        // continue;
        //}

        float startco[3], startno[3];

        copy_v3_v3(startco, vd.co);
        SCULPT_vertex_normal_get(ss, vd.vertex, startno);

        /* Ensure we have a valid normal. */
        if (fabsf(dot_v3v3(startno, startno) - 1.0f) > 0.001) {
          normalize_v3(startno);
        }

        int steps = data->do_origco ? 2 : 1;

        for (int step = 0; step < steps; step++) {
          float *co = step ? (float *)SCULPT_vertex_origco_get(ss, vd.vertex) : vd.co;

          float startvel[3];

          if (vel_scl) {
            float *vel = SCULPT_vertex_attr_get(vd.vertex, vel_scl);
#if 1
            if (isnan(dot_v3v3(vel, vel)) || !isfinite(dot_v3v3(vel, vel))) {
              printf("NaN!");
              zero_v3(vel);
            }
#endif

            copy_v3_v3(startvel, vel);
          }

          SCULPT_neighbor_coords_average_interior(
              ss,
              avg,
              vd.vertex,
              &((SculptSmoothArgs){.projection = projection,
                                   .slide_fset = fset_slide,
                                   .bound_smooth = bound_smooth,
                                   .bound_scl = bound_scl,
                                   .do_origco = step,
                                   .vel_scl = vel_scl,
                                   .vel_smooth_fac = data->vel_smooth_fac,
                                   .do_weighted_smooth = ss->cache->brush->flag2 &
                                                         BRUSH_SMOOTH_USE_AREA_WEIGHT,
                                   .preserve_fset_boundaries = ss->cache->brush->flag2 &
                                                               BRUSH_SMOOTH_PRESERVE_FACE_SETS}));

          sub_v3_v3v3(off, avg, co);

          /* Apply velocity smooth.  The point of this is to
             improve convergence for very high levels of smoothing*/
          if (vel_scl) {
            float *vel = SCULPT_vertex_attr_get(vd.vertex, vel_scl);

            float veltmp[3];
            copy_v3_v3(veltmp, vel);

            /* Remove tangental component. */
            float fac = dot_v3v3(vel, startno);
            mul_v3_v3fl(vel, startno, fac);

            madd_v3_v3fl(vel, off, 0.2f);
            mul_v3_fl(vel, 0.8f);

            /* Apply velocity. */
            add_v3_v3(off, veltmp);
          }

          madd_v3_v3v3fl(off, co, off, fade);
          SCULPT_clip(sd, ss, co, off);

          if (step == 0) {
            SCULPT_reproject_cdata(ss, vd.vertex, startco, startno);
          }
        }
      }
      if (vd.is_mesh) {
        BKE_pbvh_vert_tag_update_normal(ss->pbvh, vd.vertex);
      }
    }
    BKE_pbvh_vertex_iter_end;

    if (modified) {
      if (weighted) {
        BKE_pbvh_vert_tag_update_normal_tri_area(data->nodes[n]);
      }

      BKE_pbvh_node_mark_update(data->nodes[n]);
    }
    else {
      // not modified? remove from future iterations
      data->nodes[n] = NULL;
    }
#ifdef SMOOTH_ITER_IN_THREADS
  }
#endif
}

void SCULPT_bound_smooth_ensure(SculptSession *ss, Object *ob)
{
  SculptAttributeParams params = {0};

  if (!ss->attrs.smooth_bdist) {
    ss->attrs.smooth_bdist = BKE_sculpt_attribute_ensure(
        ob, ATTR_DOMAIN_POINT, CD_PROP_COLOR, SCULPT_ATTRIBUTE_NAME(smooth_bdist), &params);
  }
}

ATTR_NO_OPT void SCULPT_smooth(Sculpt *sd,
                               Object *ob,
                               PBVHNode **nodes,
                               const int totnode,
                               float bstrength,
                               const bool smooth_mask,
                               float projection,
                               bool do_origco)
{
  SculptSession *ss = ob->sculpt;
  Brush *brush = ss->cache && ss->cache->brush ? ss->cache->brush : BKE_paint_brush(&sd->paint);

  const float vel_smooth_cutoff = 0.5;
  const bool do_vel_smooth = bstrength > vel_smooth_cutoff && G.debug_value != 895;

  const int max_iterations = MAX2((int)(2.5f * ceilf(bstrength)), 1);
  const float fract = 1.0f / max_iterations;
  PBVHType type = BKE_pbvh_type(ss->pbvh);
  int iteration, count;
  float last;

  if (bstrength == 0.0f) {
    return;
  }

  if ((ss->cache->brush->flag2 & BRUSH_SMOOTH_USE_AREA_WEIGHT) ||
      ss->cache->brush->boundary_smooth_factor > 0.0f) {
    if (SCULPT_stroke_is_first_brush_step(ss->cache)) {
      BKE_pbvh_update_all_tri_areas(ss->pbvh);
    }
    else {
      BKE_pbvh_face_areas_begin(ss->pbvh);
    }
  }

  SculptAttributeParams params = {.permanent = false, .simple_array = false};

  if (do_vel_smooth && !ss->attrs.smooth_vel) {
    ss->attrs.smooth_vel = BKE_sculpt_attribute_ensure(
        ob, ATTR_DOMAIN_POINT, CD_PROP_FLOAT3, SCULPT_ATTRIBUTE_NAME(smooth_vel), &params);
  }

  float bstrength2 = bstrength;
  CLAMP(bstrength2, 0.0f, 1.0f);

  count = (int)(bstrength2 * max_iterations);
  last = max_iterations * (bstrength2 - count * fract);

  if (last == 0.0f) {
    count--;
    last = 1.0f;

    if (do_vel_smooth) {
      count = MAX2(count, 1);
    }
    else {
      count = MAX2(count, 0);
    }
  }

  // increase strength of last to compensate for velocity smooth in previous iterations
  if (do_vel_smooth && count > 1) {
    // last = sqrtf(last);
    last = 1.0f;
  }

  // printf("smooth iterations: %d, last: %.4f\n", count + 1, last);

  if (type == PBVH_FACES && !ss->pmap) {
    BLI_assert_msg(0, "sculpt smooth: pmap missing");
    return;
  }

  SCULPT_boundary_info_ensure(ob);

  float bound_smooth = SCULPT_get_float(ss, boundary_smooth, sd, brush);
  float fset_slide = SCULPT_get_float(ss, fset_slide, sd, brush);

  /* create temp layer for psuedo-geodesic field */
  if (bound_smooth > 0.0f) {
    bound_smooth = powf(ss->cache->brush->boundary_smooth_factor, BOUNDARY_SMOOTH_EXP);

    SCULPT_bound_smooth_ensure(ss, ob);
  }

#ifndef SMOOTH_ITER_IN_THREADS
  for (iteration = 0; iteration <= count; iteration++) {
    const float strength = (iteration != count) ? 1.0f : last;
#else
  const float strength = last;
#endif

    /* turn off velocity smooth for final iteration or two to smooth out ripples */
    bool do_vel = do_vel_smooth && iteration != count;

    if (count > 1) {
      do_vel = do_vel && iteration != count - 1;
    }

    float vel_fac = 1.0f;
    if (do_vel) {
      vel_fac = 0.25f + (bstrength - vel_smooth_cutoff) / 2.0f;
      vel_fac = min_ff(vel_fac, 1.0f);
    }

    SculptThreadedTaskData data = {
        .sd = sd,
        .ob = ob,
        .brush = brush,
        .nodes = nodes,
        .smooth_mask = smooth_mask,
        .strength = strength,
        .smooth_projection = projection,
        .fset_slide = fset_slide,
        .bound_smooth = bound_smooth,
        .scl = do_vel ? ss->attrs.smooth_vel : NULL,
        .scl2 = bound_smooth > 0.0f ? ss->attrs.smooth_bdist : NULL,
        .vel_smooth_fac = vel_fac,
        .do_origco = do_origco,
        .iterations = count + 1,
    };

    TaskParallelSettings settings;
    BKE_pbvh_parallel_range_settings(&settings, false /* XXX */, totnode);
    BLI_task_parallel_range(0, totnode, &data, do_smooth_brush_task_cb_ex, &settings);

#ifndef SMOOTH_ITER_IN_THREADS
  }
#endif
}

void SCULPT_do_smooth_brush(
    Sculpt *sd, Object *ob, PBVHNode **nodes, int totnode, float projection, bool do_origco)
{
  SculptSession *ss = ob->sculpt;

  if (SCULPT_stroke_is_first_brush_step(ss->cache) &&
      ((ss->cache->brush->flag2 & BRUSH_SMOOTH_USE_AREA_WEIGHT) ||
       ss->cache->brush->boundary_smooth_factor > 0.0f)) {
    BKE_pbvh_update_all_tri_areas(ss->pbvh);
  }

  /* NOTE: The enhance brush needs to initialize its state on the first brush step. The stroke
   * strength can become 0 during the stroke, but it can not change sign (the sign is determined
   * in the beginning of the stroke. So here it is important to not switch to enhance brush in
   * the middle of the stroke. */
  if (ss->cache->bstrength < 0.0f) {
    /* Invert mode, intensify details. */
    ss->cache->bstrength = -ss->cache->bstrength;
    SCULPT_enhance_details_brush(
        sd,
        ob,
        nodes,
        totnode,
        SCULPT_get_int(ss, enhance_detail_presteps, sd, BKE_paint_brush(&sd->paint)));
    ss->cache->bstrength = -ss->cache->bstrength;
  }
  else {
    /* Regular mode, smooth. */
    SCULPT_smooth(sd, ob, nodes, totnode, ss->cache->bstrength, false, projection, do_origco);
  }
}

/* HC Smooth Algorithm. */
/* From: Improved Laplacian Smoothing of Noisy Surface Meshes */

void SCULPT_surface_smooth_laplacian_step(SculptSession *ss,
                                          float *disp,
                                          const float co[3],
                                          SculptAttribute *scl,
                                          const PBVHVertRef v_index,
                                          const float origco[3],
                                          const float alpha,
                                          const float projection,
                                          bool check_fsets,
                                          bool weighted)
{
  float laplacian_smooth_co[3];
  float weigthed_o[3], weigthed_q[3], d[3];
  SCULPT_neighbor_coords_average(
      ss, laplacian_smooth_co, v_index, projection, check_fsets, weighted);

  // int index = BKE_pbvh_vertex_to_index(ss->pbvh, v_index);

  mul_v3_v3fl(weigthed_o, origco, alpha);
  mul_v3_v3fl(weigthed_q, co, 1.0f - alpha);
  add_v3_v3v3(d, weigthed_o, weigthed_q);
  sub_v3_v3v3((float *)SCULPT_vertex_attr_get(v_index, scl), laplacian_smooth_co, d);

  sub_v3_v3v3(disp, laplacian_smooth_co, co);
}

void SCULPT_surface_smooth_displace_step(SculptSession *ss,
                                         float *co,
                                         SculptAttribute *scl,
                                         const PBVHVertRef vertex,
                                         const float beta,
                                         const float fade)
{
  float b_avg[3] = {0.0f, 0.0f, 0.0f};
  float b_current_vertex[3];
  int total = 0;
  // int index = BKE_pbvh_vertex_to_index(ss->pbvh, v_index);

  SculptVertexNeighborIter ni;
  SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vertex, ni) {
    add_v3_v3(b_avg, (float *)SCULPT_vertex_attr_get(ni.vertex, scl));
    total++;
  }

  SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

  if (total > 0) {
    mul_v3_v3fl(b_current_vertex, b_avg, (1.0f - beta) / total);
    madd_v3_v3fl(b_current_vertex, (float *)SCULPT_vertex_attr_get(vertex, scl), beta);
    mul_v3_fl(b_current_vertex, clamp_f(fade, 0.0f, 1.0f));
    sub_v3_v3(co, b_current_vertex);
  }
}

static void SCULPT_do_surface_smooth_brush_laplacian_task_cb_ex(
    void *__restrict userdata, const int n, const TaskParallelTLS *__restrict tls)
{
  SculptThreadedTaskData *data = userdata;
  SculptSession *ss = data->ob->sculpt;
  const Brush *brush = data->brush;
  const float bstrength = ss->cache->bstrength;
  float alpha = brush->surface_smooth_shape_preservation;

  PBVHVertexIter vd;
  SculptOrigVertData orig_data;

  SculptBrushTest test;
  SculptBrushTestFn sculpt_brush_test_sq_fn = SCULPT_brush_test_init(
      ss, &test, data->brush->falloff_shape);
  const int thread_id = BLI_task_parallel_thread_id(tls);

  const bool weighted = ss->cache->brush->flag2 & BRUSH_SMOOTH_USE_AREA_WEIGHT;

  if (weighted) {
    BKE_pbvh_check_tri_areas(ss->pbvh, data->nodes[n]);
  }

  bool modified = false;

  bool check_fsets = ss->cache->brush->flag2 & BRUSH_SMOOTH_PRESERVE_FACE_SETS;
  SCULPT_orig_vert_data_init(&orig_data, data->ob, data->nodes[n], SCULPT_UNDO_COORDS);
  AutomaskingNodeData automask_data;
  SCULPT_automasking_node_begin(
      data->ob, ss, ss->cache->automasking, &automask_data, data->nodes[n]);

  bool do_reproject = SCULPT_need_reproject(ss);

  BKE_pbvh_vertex_iter_begin (ss->pbvh, data->nodes[n], vd, PBVH_ITER_UNIQUE) {
    SCULPT_orig_vert_data_update(&orig_data, vd.vertex);

    if (!sculpt_brush_test_sq_fn(&test, vd.co)) {
      continue;
    }

    SCULPT_automasking_node_update(ss, &automask_data, &vd);

    const float fade = bstrength * SCULPT_brush_strength_factor(ss,
                                                                brush,
                                                                vd.co,
                                                                sqrtf(test.dist),
                                                                vd.no,
                                                                vd.fno,
                                                                vd.mask ? *vd.mask : 0.0f,
                                                                vd.vertex,
                                                                thread_id,
                                                                &automask_data);

    float disp[3];
    float oldco[3], oldno[3];

    copy_v3_v3(oldco, vd.co);
    SCULPT_vertex_normal_get(ss, vd.vertex, oldno);

    SCULPT_surface_smooth_laplacian_step(ss,
                                         disp,
                                         vd.co,
                                         data->scl,
                                         vd.vertex,
                                         orig_data.co,
                                         alpha,
                                         data->smooth_projection,
                                         check_fsets,
                                         weighted);
    madd_v3_v3fl(vd.co, disp, clamp_f(fade, 0.0f, 1.0f));
    if (vd.is_mesh) {
      BKE_pbvh_vert_tag_update_normal(ss->pbvh, vd.vertex);
    }

    if (do_reproject) {
      SCULPT_reproject_cdata(ss, vd.vertex, oldco, oldno);
    }

    modified = true;
  }
  BKE_pbvh_vertex_iter_end;

  if (modified && weighted) {
    BKE_pbvh_vert_tag_update_normal_tri_area(data->nodes[n]);
  }
}

static void SCULPT_do_surface_smooth_brush_displace_task_cb_ex(
    void *__restrict userdata, const int n, const TaskParallelTLS *__restrict tls)
{
  SculptThreadedTaskData *data = userdata;
  SculptSession *ss = data->ob->sculpt;
  const Brush *brush = data->brush;
  const float bstrength = ss->cache->bstrength;
  const float beta = brush->surface_smooth_current_vertex;

  PBVHVertexIter vd;

  SculptBrushTest test;
  SculptBrushTestFn sculpt_brush_test_sq_fn = SCULPT_brush_test_init(
      ss, &test, data->brush->falloff_shape);
  const int thread_id = BLI_task_parallel_thread_id(tls);
  AutomaskingNodeData automask_data;
  SCULPT_automasking_node_begin(
      data->ob, ss, ss->cache->automasking, &automask_data, data->nodes[n]);

  BKE_pbvh_vertex_iter_begin (ss->pbvh, data->nodes[n], vd, PBVH_ITER_UNIQUE) {
    if (!sculpt_brush_test_sq_fn(&test, vd.co)) {
      continue;
    }

    SCULPT_automasking_node_update(ss, &automask_data, &vd);

    const float fade = bstrength * SCULPT_brush_strength_factor(ss,
                                                                brush,
                                                                vd.co,
                                                                sqrtf(test.dist),
                                                                vd.no,
                                                                vd.fno,
                                                                vd.mask ? *vd.mask : 0.0f,
                                                                vd.vertex,
                                                                thread_id,
                                                                &automask_data);
    SCULPT_surface_smooth_displace_step(ss, vd.co, data->scl, vd.vertex, beta, fade);
  }
  BKE_pbvh_vertex_iter_end;
}

void SCULPT_do_surface_smooth_brush(Sculpt *sd, Object *ob, PBVHNode **nodes, int totnode)
{
  Brush *brush = BKE_paint_brush(&sd->paint);
  SculptSession *ss = ob->sculpt;

  SculptAttributeParams params = {.permanent = false, .simple_array = false, .stroke_only = true};
  SculptAttribute *scl = BKE_sculpt_attribute_ensure(
      ob, ATTR_DOMAIN_POINT, CD_PROP_FLOAT3, "__dyntopo_lapsmooth", &params);

  if (SCULPT_stroke_is_first_brush_step(ss->cache) &&
      (ss->cache->brush->flag2 & BRUSH_SMOOTH_USE_AREA_WEIGHT)) {
    BKE_pbvh_update_all_tri_areas(ss->pbvh);
  }

  /* Threaded loop over nodes. */
  SculptThreadedTaskData data = {.sd = sd,
                                 .ob = ob,
                                 .brush = brush,
                                 .nodes = nodes,
                                 .smooth_projection = brush->autosmooth_projection,
                                 .scl = scl};

  TaskParallelSettings settings;
  BKE_pbvh_parallel_range_settings(&settings, true, totnode);
  for (int i = 0; i < brush->surface_smooth_iterations; i++) {
    BLI_task_parallel_range(
        0, totnode, &data, SCULPT_do_surface_smooth_brush_laplacian_task_cb_ex, &settings);
    BLI_task_parallel_range(
        0, totnode, &data, SCULPT_do_surface_smooth_brush_displace_task_cb_ex, &settings);
  }
}

static void SCULPT_do_directional_smooth_task_cb_ex(void *__restrict userdata,
                                                    const int n,
                                                    const TaskParallelTLS *__restrict tls)
{
  SculptThreadedTaskData *data = userdata;
  SculptSession *ss = data->ob->sculpt;
  const Brush *brush = data->brush;
  const float bstrength = ss->cache->bstrength;

  PBVHVertexIter vd;

  SculptBrushTest test;
  SculptBrushTestFn sculpt_brush_test_sq_fn = SCULPT_brush_test_init(
      ss, &test, data->brush->falloff_shape);
  const int thread_id = BLI_task_parallel_thread_id(tls);

  AutomaskingNodeData automask_data;
  SCULPT_automasking_node_begin(
      data->ob, ss, ss->cache->automasking, &automask_data, data->nodes[n]);

  BKE_pbvh_vertex_iter_begin (ss->pbvh, data->nodes[n], vd, PBVH_ITER_UNIQUE) {
    if (!sculpt_brush_test_sq_fn(&test, vd.co)) {
      continue;
    }

    SCULPT_automasking_node_update(ss, &automask_data, &vd);
    const float fade = bstrength * SCULPT_brush_strength_factor(ss,
                                                                brush,
                                                                vd.co,
                                                                sqrtf(test.dist),
                                                                vd.no,
                                                                vd.fno,
                                                                vd.mask ? *vd.mask : 0.0f,
                                                                vd.vertex,
                                                                thread_id,
                                                                &automask_data);

    float stroke_disp[3];
    sub_v3_v3v3(stroke_disp, ss->cache->location, ss->cache->last_location);
    normalize_v3(stroke_disp);

    float avg[3] = {0.0f, 0.0f, 0.0f};
    int neighbor_count = 0;

    float oldco[3], oldno[3];

    copy_v3_v3(oldco, vd.co);
    SCULPT_vertex_normal_get(ss, vd.vertex, oldno);

    bool do_reproject = SCULPT_need_reproject(ss);

    SculptVertexNeighborIter ni;
    SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vd.vertex, ni) {
      float vertex_neighbor_disp[3];
      const float *neighbor_co = SCULPT_vertex_co_get(ss, ni.vertex);
      sub_v3_v3v3(vertex_neighbor_disp, neighbor_co, vd.co);
      normalize_v3(vertex_neighbor_disp);
      if (fabsf(dot_v3v3(stroke_disp, vertex_neighbor_disp)) > 0.6f) {
        neighbor_count++;
        add_v3_v3(avg, neighbor_co);
      }
    }
    SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

    /* Avoid division by 0 when there are no neighbors. */
    if (neighbor_count == 0) {
      continue;
    }

    float smooth_co[3];
    mul_v3_v3fl(smooth_co, avg, 1.0f / neighbor_count);

    float final_disp[3];
    sub_v3_v3v3(final_disp, smooth_co, vd.co);
    madd_v3_v3v3fl(final_disp, vd.co, final_disp, fade);
    SCULPT_clip(data->sd, ss, vd.co, final_disp);

    if (do_reproject) {
      SCULPT_reproject_cdata(ss, vd.vertex, oldco, oldno);
    }

    if (vd.is_mesh) {
      BKE_pbvh_vert_tag_update_normal(ss->pbvh, vd.vertex);
    }

    BKE_pbvh_vertex_iter_end;
  }
}

void SCULPT_do_directional_smooth_brush(Sculpt *sd, Object *ob, PBVHNode **nodes, int totnode)
{
  Brush *brush = BKE_paint_brush(&sd->paint);

  /* Threaded loop over nodes. */
  SculptThreadedTaskData data = {
      .sd = sd,
      .ob = ob,
      .brush = brush,
      .nodes = nodes,
  };

  TaskParallelSettings settings;
  BKE_pbvh_parallel_range_settings(&settings, true, totnode);
  for (int i = 0; i < brush->surface_smooth_iterations; i++) {
    BLI_task_parallel_range(0, totnode, &data, SCULPT_do_directional_smooth_task_cb_ex, &settings);
  }
}

static void SCULPT_do_uniform_weigths_smooth_task_cb_ex(void *__restrict userdata,
                                                        const int n,
                                                        const TaskParallelTLS *__restrict tls)
{
  SculptThreadedTaskData *data = userdata;
  SculptSession *ss = data->ob->sculpt;
  const Brush *brush = data->brush;
  const float bstrength = ss->cache->bstrength;

  PBVHVertexIter vd;

  SculptBrushTest test;
  SculptBrushTestFn sculpt_brush_test_sq_fn = SCULPT_brush_test_init(
      ss, &test, data->brush->falloff_shape);
  const int thread_id = BLI_task_parallel_thread_id(tls);

  bool do_reproject = SCULPT_need_reproject(ss);
  AutomaskingNodeData automask_data;
  SCULPT_automasking_node_begin(
      data->ob, ss, ss->cache->automasking, &automask_data, data->nodes[n]);

  BKE_pbvh_vertex_iter_begin (ss->pbvh, data->nodes[n], vd, PBVH_ITER_UNIQUE) {
    if (!sculpt_brush_test_sq_fn(&test, vd.co)) {
      continue;
    }
    SCULPT_automasking_node_update(ss, &automask_data, &vd);

    const float fade = bstrength * SCULPT_brush_strength_factor(ss,
                                                                brush,
                                                                vd.co,
                                                                sqrtf(test.dist),
                                                                vd.no,
                                                                vd.fno,
                                                                vd.mask ? *vd.mask : 0.0f,
                                                                vd.vertex,
                                                                thread_id,
                                                                &automask_data);

    float len_accum = 0;
    int tot_neighbors = 0;

    float oldco[3], oldno[3];

    copy_v3_v3(oldco, vd.co);
    SCULPT_vertex_normal_get(ss, vd.vertex, oldno);

    SculptVertexNeighborIter ni;
    SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vd.vertex, ni) {
      len_accum += len_v3v3(SCULPT_vertex_co_get(ss, vd.vertex),
                            SCULPT_vertex_co_get(ss, ni.vertex));
      tot_neighbors++;
    }
    SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

    /* Avoid division by 0 when there are no neighbors. */
    if (tot_neighbors == 0) {
      continue;
    }

    const float len_avg = bstrength * len_accum / tot_neighbors;

    float co_accum[3] = {0.0f};

    SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vd.vertex, ni) {
      float neighbor_co[3];
      float neighbor_disp[3];

      sub_v3_v3v3(
          neighbor_disp, SCULPT_vertex_co_get(ss, ni.vertex), SCULPT_vertex_co_get(ss, vd.vertex));
      normalize_v3(neighbor_disp);
      mul_v3_fl(neighbor_disp, len_avg);
      add_v3_v3v3(neighbor_co, SCULPT_vertex_co_get(ss, vd.vertex), neighbor_disp);
      add_v3_v3(co_accum, neighbor_co);
    }
    SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

    float smooth_co[3];
    mul_v3_v3fl(smooth_co, co_accum, 1.0f / tot_neighbors);

    float final_disp[3];
    sub_v3_v3v3(final_disp, smooth_co, vd.co);
    madd_v3_v3v3fl(final_disp, vd.co, final_disp, fade);
    SCULPT_clip(data->sd, ss, vd.co, final_disp);

    if (vd.is_mesh) {
      BKE_pbvh_vert_tag_update_normal(ss->pbvh, vd.vertex);
    }

    if (do_reproject) {
      SCULPT_reproject_cdata(ss, vd.vertex, oldco, oldno);
    }

    BKE_pbvh_vertex_iter_end;
  }
}

void SCULPT_do_uniform_weights_smooth_brush(Sculpt *sd, Object *ob, PBVHNode **nodes, int totnode)
{
  Brush *brush = BKE_paint_brush(&sd->paint);

  /* Threaded loop over nodes. */
  SculptThreadedTaskData data = {
      .sd = sd,
      .ob = ob,
      .brush = brush,
      .nodes = nodes,
  };

  TaskParallelSettings settings;
  BKE_pbvh_parallel_range_settings(&settings, true, totnode);
  for (int i = 0; i < brush->surface_smooth_iterations; i++) {
    BLI_task_parallel_range(
        0, totnode, &data, SCULPT_do_uniform_weigths_smooth_task_cb_ex, &settings);
  }
}

static void do_smooth_vcol_boundary_brush_task_cb_ex(void *__restrict userdata,
                                                     const int n,
                                                     const TaskParallelTLS *__restrict tls)
{
  SculptThreadedTaskData *data = userdata;
  SculptSession *ss = data->ob->sculpt;
  Sculpt *sd = data->sd;
  const Brush *brush = data->brush;
  const bool smooth_mask = data->smooth_mask;
  float bstrength = data->strength;

  PBVHVertexIter vd;

  CLAMP(bstrength, 0.0f, 1.0f);

  SculptBrushTest test;
  SculptBrushTestFn sculpt_brush_test_sq_fn = SCULPT_brush_test_init(
      ss, &test, data->brush->falloff_shape);

  const int thread_id = BLI_task_parallel_thread_id(tls);

  AutomaskingNodeData automask_data;
  SCULPT_automasking_node_begin(
      data->ob, ss, ss->cache->automasking, &automask_data, data->nodes[n]);

  float avg[4] = {0.0f, 0.0f, 0.0f, 0.0f};
  float tot = 0.0f;
  BKE_pbvh_vertex_iter_begin (ss->pbvh, data->nodes[n], vd, PBVH_ITER_UNIQUE) {
    float vcolor[4];

    SCULPT_vertex_color_get(ss, vd.vertex, vcolor);

    if (sculpt_brush_test_sq_fn(&test, vd.co)) {
      SCULPT_automasking_node_update(ss, &automask_data, &vd);

      const float fade = bstrength * SCULPT_brush_strength_factor(
                                         ss,
                                         brush,
                                         vd.co,
                                         sqrtf(test.dist),
                                         vd.no,
                                         vd.fno,
                                         smooth_mask ? 0.0f : (vd.mask ? *vd.mask : 0.0f),
                                         vd.vertex,
                                         thread_id,
                                         &automask_data);

      madd_v3_v3fl(avg, vcolor, fade);
      tot += fade;
    }
  }
  BKE_pbvh_vertex_iter_end;

  if (tot == 0.0f) {
    return;
  }
  tot = 1.0f / tot;

  mul_v3_fl(avg, tot);

  float exp = brush->vcol_boundary_exponent;
  // detect bad value

  if (exp == 0.0f) {
    exp = 1.0f;
  }

  SCULPT_automasking_node_begin(
      data->ob, ss, ss->cache->automasking, &automask_data, data->nodes[n]);

  BKE_pbvh_vertex_iter_begin (ss->pbvh, data->nodes[n], vd, PBVH_ITER_UNIQUE) {
    if (sculpt_brush_test_sq_fn(&test, vd.co)) {
      SCULPT_automasking_node_update(ss, &automask_data, &vd);

      const float fade = bstrength * SCULPT_brush_strength_factor(
                                         ss,
                                         brush,
                                         vd.co,
                                         sqrtf(test.dist),
                                         vd.no,
                                         vd.fno,
                                         smooth_mask ? 0.0f : (vd.mask ? *vd.mask : 0.0f),
                                         vd.vertex,
                                         thread_id,
                                         &automask_data);

      float vcolor[4];

      SCULPT_vertex_color_get(ss, vd.vertex, vcolor);

      float avg2[3], avg3[3], val[3];
      float tot2 = 0.0f, tot4 = 0.0f;

      copy_v4_v4(avg, vcolor);

      zero_v3(avg2);
      zero_v3(avg3);

      madd_v3_v3fl(avg2, vd.co, 0.5f);
      tot2 += 0.5f;

      SculptVertexNeighborIter ni;
      SCULPT_VERTEX_NEIGHBORS_ITER_BEGIN (ss, vd.vertex, ni) {
        float col[4];

        SCULPT_vertex_color_get(ss, ni.vertex, col);
        const float *co = SCULPT_vertex_co_get(ss, ni.vertex);

        // simple color metric.  TODO: plug in appropriate color space code?
        float dv[4];
        sub_v4_v4v4(dv, col, avg);
        float w = (fabs(dv[0]) + fabs(dv[1]) + fabs(dv[2]) + fabs(dv[3])) / 4.0;

        w = powf(w, exp);

        madd_v3_v3fl(avg3, co, 1.0f);
        tot4 += 1.0f;

        madd_v3_v3fl(avg2, co, w);
        tot2 += w;
      }
      SCULPT_VERTEX_NEIGHBORS_ITER_END(ni);

      if (tot2 == 0.0f) {
        continue;
      }

      if (tot4 > 0.0f) {
        mul_v3_fl(avg3, 1.0f / tot4);
      }

      /* try to avoid perfectly colinear triangles, and the normal discontinuities they create,
         by blending slightly with unweighted smoothed position */
      mul_v3_fl(avg2, 1.0f / tot2);
      interp_v3_v3v3(avg2, avg2, avg3, 0.025);

      sub_v3_v3v3(val, avg2, vd.co);
      madd_v3_v3v3fl(val, vd.co, val, fade);
      SCULPT_clip(sd, ss, vd.co, val);

      if (vd.is_mesh) {
        BKE_pbvh_vert_tag_update_normal(ss->pbvh, vd.vertex);
      }
    }
  }
  BKE_pbvh_vertex_iter_end;
}

void SCULPT_smooth_vcol_boundary(
    Sculpt *sd, Object *ob, PBVHNode **nodes, const int totnode, float bstrength)
{
  SculptSession *ss = ob->sculpt;

  Brush *brush = BKE_paint_brush(&sd->paint);

  const int max_iterations = 4;
  const float fract = 1.0f / max_iterations;
  PBVHType type = BKE_pbvh_type(ss->pbvh);
  int iteration, count;
  float last;

  // Did I do this? Why? -joeedh
  // CLAMP(bstrength, 0.0f, 1.0f);

  count = (int)(bstrength * max_iterations);
  last = max_iterations * (bstrength - count * fract);

  if (type == PBVH_FACES && !ss->pmap) {
    BLI_assert(!"sculpt smooth: pmap missing");
    return;
  }

  SCULPT_vertex_random_access_ensure(ss);
  SCULPT_boundary_info_ensure(ob);

  for (iteration = 0; iteration <= count; iteration++) {
    const float strength = (iteration != count) ? 1.0f : last;

    SculptThreadedTaskData data = {
        .sd = sd,
        .ob = ob,
        .brush = brush,
        .nodes = nodes,
        .smooth_mask = false,
        .strength = strength,
    };

    TaskParallelSettings settings;
    BKE_pbvh_parallel_range_settings(&settings, true, totnode);
    BLI_task_parallel_range(
        0, totnode, &data, do_smooth_vcol_boundary_brush_task_cb_ex, &settings);
  }
}
