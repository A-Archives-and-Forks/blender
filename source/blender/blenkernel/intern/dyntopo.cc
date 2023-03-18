#include "MEM_guardedalloc.h"

#include "DNA_customdata_types.h"
#include "DNA_mesh_types.h"
#include "DNA_meshdata_types.h"
#include "DNA_modifier_types.h"

#include "BLI_array.hh"
#include "BLI_index_range.hh"
#include "BLI_map.hh"
#include "BLI_set.hh"
#include "BLI_task.hh"
#include "BLI_vector.hh"

#include "BLI_alloca.h"
#include "BLI_array.h"
#include "BLI_asan.h"
#include "BLI_bitmap.h"
#include "BLI_buffer.h"
#include "BLI_compiler_attrs.h"
#include "BLI_compiler_compat.h"
#include "BLI_ghash.h"
#include "BLI_heap.h"
#include "BLI_heap_minmax.hh"
#include "BLI_heap_simple.h"
#include "BLI_linklist.h"
#include "BLI_math.h"
#include "BLI_memarena.h"
#include "BLI_rand.h"
#include "BLI_smallhash.h"
#include "BLI_task.h"
#include "BLI_utildefines.h"
#include "PIL_time.h"
#include "atomic_ops.h"

#include "BKE_customdata.h"
#include "BKE_dyntopo.h"
#include "BKE_paint.h"
#include "BKE_pbvh.h"

#include "bmesh.h"
#include "bmesh_log.h"

#include "dyntopo_intern.hh"
#include "pbvh_intern.hh"

#include <cstdio>

namespace blender::dyntopo {

using blender::float2;
using blender::float3;
using blender::float4;
using blender::IndexRange;
using blender::Map;
using blender::Set;
using blender::Vector;

static void pbvh_split_edges(struct EdgeQueueContext *eq_ctx,
                             PBVH *pbvh,
                             BMesh *bm,
                             BMEdge **edges,
                             int totedge,
                             bool ignore_isolated_edges);
extern "C" {
void bmesh_disk_edge_append(BMEdge *e, BMVert *v);
void bmesh_radial_loop_append(BMEdge *e, BMLoop *l);
void bm_kill_only_edge(BMesh *bm, BMEdge *e);
void bm_kill_only_loop(BMesh *bm, BMLoop *l);
void bm_kill_only_face(BMesh *bm, BMFace *f);
}

static bool edge_queue_test(struct EdgeQueueContext *eq_ctx, PBVH *pbvh, BMEdge *e);
static void edge_queue_create_local(struct EdgeQueueContext *eq_ctx,
                                    PBVH *pbvh,
                                    const float center[3],
                                    const float view_normal[3],
                                    float radius,
                                    const bool use_frontface,
                                    const bool use_projected,
                                    PBVHTopologyUpdateMode local_mode);

BLI_INLINE void surface_smooth_v_safe(PBVH *pbvh, BMVert *v, float fac)
{
  float co[3];
  float origco[3], origco1[3];
  float origno1[3];
  float tan[3];
  float tot = 0.0;

  MSculptVert *mv1 = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, v);

  if (mv1->stroke_id != pbvh->stroke_id) {
    copy_v3_v3(origco1, v->co);
    copy_v3_v3(origno1, v->no);
  }
  else {
    copy_v3_v3(origco1, mv1->origco);
    copy_v3_v3(origno1, dot_v3v3(mv1->origno, mv1->origno) == 0.0f ? v->no : mv1->origno);
  }

  zero_v3(co);
  zero_v3(origco);

  // this is a manual edge walk

  BMEdge *e = v->e;
  if (!e) {
    return;
  }

  if (pbvh_boundary_needs_update_bmesh(pbvh, v)) {
    pbvh_check_vert_boundary(pbvh, v);
  }

  // pbvh_check_vert_boundary(pbvh, v);

  const int cd_sculpt_vert = pbvh->cd_sculpt_vert;
  const int boundflag = BM_ELEM_CD_GET_INT(v, pbvh->cd_boundary_flag);

  const bool bound1 = boundflag & SCULPTVERT_SMOOTH_BOUNDARY;

  if (boundflag & SCULPTVERT_SMOOTH_CORNER) {
    return;
  }

  if (bound1) {
    fac *= 0.1;
  }

  do {
    BMVert *v2 = e->v1 == v ? e->v2 : e->v1;

    // can't check for boundary here, thread
    // pbvh_check_vert_boundary(pbvh, v2);

    MSculptVert *mv2 = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, v2);
    int boundflag2 = BM_ELEM_CD_GET_INT(v2, pbvh->cd_boundary_flag);

    const bool bound2 = boundflag2 & SCULPTVERT_SMOOTH_BOUNDARY;

    if (bound1 && !bound2) {
      continue;
    }

    sub_v3_v3v3(tan, v2->co, v->co);
    float d = dot_v3v3(tan, v->no);

    madd_v3_v3fl(tan, v->no, -d * 0.99f);
    add_v3_v3(co, tan);

    if (mv2->stroke_id == pbvh->stroke_id) {
      sub_v3_v3v3(tan, mv2->origco, origco1);
    }
    else {
      sub_v3_v3v3(tan, v2->co, origco1);
    }

    d = dot_v3v3(tan, origno1);
    madd_v3_v3fl(tan, origno1, -d * 0.99f);
    add_v3_v3(origco, tan);

    tot += 1.0f;

  } while ((e = BM_DISK_EDGE_NEXT(e, v)) != v->e);

  if (tot == 0.0f) {
    return;
  }

  mul_v3_fl(co, 1.0f / tot);
  mul_v3_fl(origco, 1.0f / tot);

  volatile float x = v->co[0], y = v->co[1], z = v->co[2];
  volatile float nx = x + co[0] * fac, ny = y + co[1] * fac, nz = z + co[2] * fac;

  /* Conflicts here should be pretty rare. */
  atomic_cas_float(&v->co[0], x, nx);
  atomic_cas_float(&v->co[1], y, ny);
  atomic_cas_float(&v->co[2], z, nz);

  /* Conflicts here should be pretty rare. */
  x = mv1->origco[0];
  y = mv1->origco[1];
  z = mv1->origco[2];

  nx = x + origco[0] * fac;
  ny = y + origco[1] * fac;
  nz = z + origco[2] * fac;

  atomic_cas_float(&mv1->origco[0], x, nx);
  atomic_cas_float(&mv1->origco[1], y, ny);
  atomic_cas_float(&mv1->origco[2], z, nz);

  PBVH_CHECK_NAN(mv1->origco);
  PBVH_CHECK_NAN(v->co);
  // atomic_cas_int32(&mv1->stroke_id, stroke_id, pbvh->stroke_id);
}

/****************************** EdgeQueue *****************************/

static float maskcb_get(EdgeQueueContext *eq_ctx, BMVert *v1, BMVert *v2)
{
  if (eq_ctx->mask_cb) {
    PBVHVertRef sv1 = {(intptr_t)v1};
    PBVHVertRef sv2 = {(intptr_t)v2};

    float w1 = eq_ctx->mask_cb(sv1, eq_ctx->mask_cb_data);
    float w2 = eq_ctx->mask_cb(sv2, eq_ctx->mask_cb_data);

    // float limit = 0.5;
    // if (w1 > limit || w2 > limit) {
    return min_ff(w1, w2);
    //}

    // return (w1 + w2) * 0.5f;
  }

  return 1.0f;
}

BLI_INLINE float calc_edge_length(EdgeQueueContext *eq_ctx, BMVert *v1, BMVert *v2)
{
  return len_squared_v3v3(v1->co, v2->co);
}

BLI_INLINE float calc_weighted_length(EdgeQueueContext *eq_ctx, BMVert *v1, BMVert *v2, float sign)
{
  float w = 1.0 - maskcb_get(eq_ctx, v1, v2);
  float len = len_squared_v3v3(v1->co, v2->co);

  w = 1.0 + w * sign;

  return len * w;
}

static void edge_queue_insert_unified(EdgeQueueContext *eq_ctx, BMEdge *e)
{
  if (/*BLI_mm_heap_len(eq_ctx->heap_mm) < eq_ctx->max_heap_mm &&*/ !(e->head.hflag &
                                                                      BM_ELEM_TAG)) {
    float lensqr = len_squared_v3v3(e->v1->co, e->v2->co);

    float len = sqrtf(lensqr);

    eq_ctx->avg_elen += len;
    eq_ctx->min_elen = min_ff(eq_ctx->min_elen, len);
    eq_ctx->max_elen = max_ff(eq_ctx->max_elen, len);
    eq_ctx->totedge++;

    // lensqr += (BLI_thread_frand(0) - 0.5f) * 0.1 * eq_ctx->limit_mid;

    eq_ctx->edge_heap.insert(lensqr, e);
    e->head.hflag |= BM_ELEM_TAG;
  }
}

static void edge_queue_insert_val34_vert(EdgeQueueContext *eq_ctx, BMVert *v)
{
  // BLI_table_gset_add(eq_ctx->used_verts, v);
  eq_ctx->tot_used_verts++;

  if (!eq_ctx->used_verts || eq_ctx->tot_used_verts >= eq_ctx->used_verts_size) {
    int newlen = (eq_ctx->tot_used_verts + 1);
    newlen += newlen >> 1;

    eq_ctx->used_verts_size = newlen;
    if (eq_ctx->used_verts) {
      eq_ctx->used_verts = (BMVert **)MEM_reallocN(eq_ctx->used_verts, newlen * sizeof(BMVert *));
    }
    else {
      eq_ctx->used_verts = (BMVert **)MEM_malloc_arrayN(
          newlen, sizeof(BMVert *), "eq_ctx->used_verts");
    }
  }

  eq_ctx->used_verts[eq_ctx->tot_used_verts - 1] = v;
}

BLI_INLINE float calc_curvature_weight(EdgeQueueContext *eq_ctx, BMVert *v1, BMVert *v2)
{
#ifdef WITH_ADAPTIVE_CURVATURE
  MSculptVert *mv1 = BKE_PBVH_SCULPTVERT(eq_ctx->cd_sculpt_vert, v1);
  MSculptVert *mv2 = BKE_PBVH_SCULPTVERT(eq_ctx->cd_sculpt_vert, v2);

  float c1 = (float)mv1->curv / 65535.0f;
  float c2 = (float)mv2->curv / 65535.0f;
  float fac = 1.0f + powf((c1 + c2) * 100.0, 4.0f);
  fac = min_ff(fac, 4.0f);

  return fac;
#else
  return 1.0;
#endif
}

static bool edge_queue_vert_in_sphere(const EdgeQueueContext *eq_ctx, BMVert *v)
{
  /* Check if triangle intersects the sphere */
  return len_squared_v3v3(eq_ctx->center, v->co) <= eq_ctx->radius_squared;
}

/*
  Profiling revealed the accurate distance to tri in blenlib was too slow,
  so we use a simpler version here
  */
/* reduce script

on factor;
off period;

load_package "avector";

comment: origin at p;

p := avec(0, 0, 0);
n :=- avec(nx, ny, nz);
v1 := avec(v1x, v1y, v1z);
v2 := avec(v2x, v2y, v2z);
v3 := avec(v3x, v3y, v3z);

comment: -((p - v1) dot n);simplified to this;
fac := v1 dot n;

co := fac*n;

a := co - v1;
b := co - v2;
c := co - v3;

a := v1;
b := v2;
c := v3;

t1 := a cross b;
t2 := b cross c;
t3 := c cross a;

on fort;
w1 := t1 dot n;
w2 := t2 dot n;
w3 := t3 dot n;
off fort;

inside := sign(a dot n) + sign(b dot n) + sign(c dot n);


*/
static bool point_in_tri_v3(float p[3], float v1[3], float v2[3], float v3[3], float n[3])
{
  float t1[3], t2[3], t3[3];
  sub_v3_v3v3(t1, v1, p);
  sub_v3_v3v3(t2, v2, p);
  sub_v3_v3v3(t3, v3, p);

  float c1[3], c2[3], c3[3];
  cross_v3_v3v3(c1, t1, t2);
  cross_v3_v3v3(c2, t2, t3);
  cross_v3_v3v3(c3, t3, t1);

  bool w1 = dot_v3v3(c1, n) >= 0.0f;
  bool w2 = dot_v3v3(c2, n) >= 0.0f;
  bool w3 = dot_v3v3(c3, n) >= 0.0f;

  return w1 == w2 && w2 == w3;

#if 0
  const float nx = n[0], ny = n[1], nz = n[2];

  float v1x = v1[0] - p[0], v1y = v1[1] - p[1], v1z = v1[2] - p[2];
  float v2x = v2[0] - p[0], v2y = v2[1] - p[1], v2z = v2[2] - p[2];
  float v3x = v3[0] - p[0], v3y = v3[1] - p[1], v3z = v3[2] - p[2];

  const float w1 = -(nx * v1y * v2z - nx * v1z * v2y - ny * v1x * v2z + ny * v1z * v2x +
                     nz * v1x * v2y - nz * v1y * v2x);
  const float w2 = -(nx * v2y * v3z - nx * v2z * v3y - ny * v2x * v3z + ny * v2z * v3x +
                     nz * v2x * v3y - nz * v2y * v3x);
  const float w3 = nx * v1y * v3z - nx * v1z * v3y - ny * v1x * v3z + ny * v1z * v3x +
                   nz * v1x * v3y - nz * v1y * v3x;
  return !((w1 >= 0.0f) && (w2 >= 0.0f) && (w3 >= 0.0f));
#endif
}

static float dist_to_tri_sphere_simple(
    float p[3], float v1[3], float v2[3], float v3[3], float n[3])
{
  float co[3];
  float t1[3], t2[3], t3[3];

  if (dot_v3v3(n, n) == 0.0f) {
    normal_tri_v3(n, v1, v2, v3);
  }

  if (point_in_tri_v3(p, v1, v2, v3, n)) {
    sub_v3_v3v3(co, p, v2);

    float dist = dot_v3v3(co, n);
    return dist * dist;
  }

  sub_v3_v3v3(co, p, v1);
  madd_v3_v3fl(co, n, -dot_v3v3(n, co));

  sub_v3_v3v3(t1, v1, co);
  sub_v3_v3v3(t2, v2, co);
  sub_v3_v3v3(t3, v3, co);

  float dis = len_squared_v3v3(p, v1);
  dis = fmin(dis, len_squared_v3v3(p, v2));
  dis = fmin(dis, len_squared_v3v3(p, v3));

  add_v3_v3v3(co, v1, v2);
  mul_v3_fl(co, 0.5f);
  dis = fmin(dis, len_squared_v3v3(p, co));

  add_v3_v3v3(co, v2, v3);
  mul_v3_fl(co, 0.5f);
  dis = fmin(dis, len_squared_v3v3(p, co));

  add_v3_v3v3(co, v3, v1);
  mul_v3_fl(co, 0.5f);
  dis = fmin(dis, len_squared_v3v3(p, co));

  add_v3_v3v3(co, v1, v2);
  add_v3_v3(co, v3);
  mul_v3_fl(co, 1.0f / 3.0f);
  dis = fmin(dis, len_squared_v3v3(p, co));

  return dis;
}

static bool edge_queue_tri_in_sphere(const EdgeQueueContext *q, BMVert *vs[3], float no[3])
{
#if 0
  float cent[3];

  zero_v3(cent);
  add_v3_v3(cent, l->v->co);
  add_v3_v3(cent, l->next->v->co);
  add_v3_v3(cent, l->prev->v->co);

  mul_v3_fl(cent, 1.0f / 3.0f);
  return len_squared_v3v3(cent, q->center) < q->radius_squared;
#endif

  /* Check if triangle intersects the sphere */
#if 1
  float dis = dist_to_tri_sphere_simple(
      (float *)q->center, (float *)vs[0]->co, (float *)vs[1]->co, (float *)vs[2]->co, (float *)no);
#else
  float dis = len_squared_v3v3(q->center, l->v->co);
#endif

  return dis <= q->radius_squared;
}

static bool edge_queue_tri_in_circle(const EdgeQueueContext *q, BMVert *v_tri[3], float no[3])
{
  float c[3];
  float tri_proj[3][3];

  project_plane_normalized_v3_v3v3(tri_proj[0], v_tri[0]->co, q->view_normal);
  project_plane_normalized_v3_v3v3(tri_proj[1], v_tri[1]->co, q->view_normal);
  project_plane_normalized_v3_v3v3(tri_proj[2], v_tri[2]->co, q->view_normal);

  closest_on_tri_to_point_v3(c, q->center_proj, tri_proj[0], tri_proj[1], tri_proj[2]);

  /* Check if triangle intersects the sphere */
  return len_squared_v3v3(q->center_proj, c) <= q->radius_squared;
}

typedef struct EdgeQueueThreadData {
  PBVH *pbvh;
  PBVHNode *node;
  BMEdge **edges;
  EdgeQueueContext *eq_ctx;
  int totedge;
  int size;
  bool is_collapse;
  int seed;
} EdgeQueueThreadData;

static void edge_thread_data_insert(EdgeQueueThreadData *tdata, BMEdge *e)
{
  if (tdata->size <= tdata->totedge) {
    tdata->size = (tdata->totedge + 1) << 1;
    if (!tdata->edges) {
      tdata->edges = (BMEdge **)MEM_mallocN(sizeof(void *) * tdata->size,
                                            "edge_thread_data_insert");
    }
    else {
      tdata->edges = (BMEdge **)MEM_reallocN(tdata->edges, sizeof(void *) * tdata->size);
    }
  }

  BMElem elem;
  memcpy(&elem, (BMElem *)e, sizeof(BMElem));

  elem.head.hflag = e->head.hflag | BM_ELEM_TAG;
  int64_t iold = *((int64_t *)&e->head.index);
  int64_t inew = *((int64_t *)&elem.head.index);

  atomic_cas_int64((int64_t *)&e->head.index, iold, inew);

  tdata->edges[tdata->totedge] = e;
  tdata->totedge++;
}

static bool edge_queue_vert_in_circle(const EdgeQueueContext *eq_ctx, BMVert *v)
{
  float c[3];

  project_plane_normalized_v3_v3v3(c, v->co, eq_ctx->view_normal);

  return len_squared_v3v3(eq_ctx->center_proj, c) <= eq_ctx->radius_squared;
}

static void long_edge_queue_edge_add_recursive(EdgeQueueContext *eq_ctx,
                                               BMLoop *l_edge,
                                               BMLoop *l_end,
                                               const float len_sq,
                                               float limit_len,
                                               int depth)
{
  BLI_assert(len_sq > square_f(limit_len));

  if (depth > DEPTH_START_LIMIT && eq_ctx->use_view_normal) {
    if (dot_v3v3(l_edge->f->no, eq_ctx->view_normal) < 0.0f) {
      return;
    }
  }

  edge_queue_insert_unified(eq_ctx, l_edge->e);

  if ((l_edge->radial_next != l_edge)) {
    const float len_sq_cmp = len_sq * EVEN_EDGELEN_THRESHOLD;

    limit_len *= EVEN_GENERATION_SCALE;
    const float limit_len_sq = square_f(limit_len);

    BMLoop *l_iter = l_edge;
    do {
      BMLoop *l_adjacent[2] = {l_iter->next, l_iter->prev};
      for (int i = 0; i < (int)ARRAY_SIZE(l_adjacent); i++) {
        float len_sq_other = calc_weighted_length(
            eq_ctx, l_adjacent[i]->e->v1, l_adjacent[i]->e->v2, -1.0f);

        bool insert_ok = len_sq_other > max_ff(len_sq_cmp, limit_len_sq);
        if (!insert_ok) {
          continue;
        }

        long_edge_queue_edge_add_recursive(
            eq_ctx, l_adjacent[i]->radial_next, l_adjacent[i], len_sq_other, limit_len, depth + 1);
      }
    } while ((l_iter = l_iter->radial_next) != l_end);
  }
}

static void long_edge_queue_edge_add_recursive_2(EdgeQueueThreadData *tdata,
                                                 BMLoop *l_edge,
                                                 BMLoop *l_end,
                                                 const float len_sq,
                                                 float limit_len,
                                                 int depth,
                                                 bool insert)
{
  BLI_assert(len_sq > square_f(limit_len));

  if (l_edge->e->head.hflag & BM_ELEM_TAG) {
    // return;
  }

#ifdef USE_EDGEQUEUE_FRONTFACE
  if (depth > DEPTH_START_LIMIT && tdata->eq_ctx->use_view_normal) {
    if (dot_v3v3(l_edge->f->no, tdata->eq_ctx->view_normal) < 0.0f) {
      return;
    }
  }
#endif

  if (insert) {
    edge_thread_data_insert(tdata, l_edge->e);
  }

  if ((l_edge->radial_next != l_edge)) {
    const float len_sq_cmp = len_sq * EVEN_EDGELEN_THRESHOLD;

    limit_len *= EVEN_GENERATION_SCALE;
    const float limit_len_sq = square_f(limit_len);

    BMLoop *l_iter = l_edge;
    do {
      BMLoop *l_adjacent[2] = {l_iter->next, l_iter->prev};
      for (int i = 0; i < (int)ARRAY_SIZE(l_adjacent); i++) {
        float len_sq_other = calc_edge_length(
            tdata->eq_ctx, l_adjacent[i]->e->v1, l_adjacent[i]->e->v2);

        bool insert_ok = len_sq_other > max_ff(len_sq_cmp, limit_len_sq);
#ifdef EVEN_NO_TEST_DEPTH_LIMIT
        if (!insert_ok && depth >= EVEN_NO_TEST_DEPTH_LIMIT) {
          continue;
        }
#else
        if (!insert_ok) {
          continue;
        }
#endif

        long_edge_queue_edge_add_recursive_2(tdata,
                                             l_adjacent[i]->radial_next,
                                             l_adjacent[i],
                                             len_sq_other,
                                             limit_len,
                                             depth + 1,
                                             insert_ok);
      }
    } while ((l_iter = l_iter->radial_next) != l_end);
  }
}

BLI_INLINE int dyntopo_thread_rand(int seed)
{
  // glibc
  const uint32_t multiplier = 1103515245;
  const uint32_t addend = 12345;
  const uint32_t mask = (1 << 30) - 1;

  return (seed * multiplier + addend) & mask;
}

static void unified_edge_queue_task_cb(void *__restrict userdata,
                                       const int n,
                                       const TaskParallelTLS *__restrict tls)
{
  EdgeQueueThreadData *tdata = ((EdgeQueueThreadData *)userdata) + n;
  PBVH *pbvh = tdata->pbvh;
  PBVHNode *node = tdata->node;
  EdgeQueueContext *eq_ctx = tdata->eq_ctx;
  RNG *rng = BLI_rng_new(POINTER_AS_UINT(tdata));

  int seed = tdata->seed + n;

  BMFace *f;
  bool do_smooth = eq_ctx->surface_smooth_fac > 0.0f;

  BKE_pbvh_bmesh_check_tris(pbvh, node);
  int ni = node - pbvh->nodes;

  const char facetag = BM_ELEM_TAG_ALT;

#if 1
#  if 0
  // try to be nice to branch predictor
  int off = (seed = dyntopo_thread_rand(seed));
  int stepi = off & 65535;
#  endif
  /*
  we care more about convergence to accurate results
  then accuracy in any individual runs.  profiling
  has shown this loop overwhelms the L3 cache,
  so randomly skip bits of it.
  */
  TGSET_ITER (f, node->bm_faces) {
    BMLoop *l = f->l_first;

    f->head.hflag &= ~facetag;

#  if 0
    if ((stepi++) & 3) {
      continue;
    }
#  else
    if ((seed = dyntopo_thread_rand(seed)) & 3) {
      continue;
    }
#  endif
    do {
      /* kind of tricky to atomicly update flags here. . . */
      BMEdge edge = *l->e;
      edge.head.hflag &= ~BM_ELEM_TAG;

      int64_t *t1 = (int64_t *)&edge.head.index;
      int64_t *t2 = (int64_t *)&l->e->head.index;

      atomic_cas_int64(t2, *t2, *t1);

      l->e->head.hflag &= ~BM_ELEM_TAG;
      l = l->next;
    } while (l != f->l_first);
  }
  TGSET_ITER_END
#endif

  PBVHTriBuf *tribuf = node->tribuf;
  for (int i = 0; i < node->tribuf->tottri; i++) {
    PBVHTri *tri = node->tribuf->tris + i;
    BMFace *f = (BMFace *)tri->f.i;

    if (f->head.hflag & facetag) {
      continue;
    }

#ifdef USE_EDGEQUEUE_FRONTFACE
    if (eq_ctx->use_view_normal) {
      if (dot_v3v3(f->no, eq_ctx->view_normal) < 0.0f) {
        continue;
      }
    }
#endif

    BMVert *vs[3] = {(BMVert *)tribuf->verts[tri->v[0]].i,
                     (BMVert *)tribuf->verts[tri->v[1]].i,
                     (BMVert *)tribuf->verts[tri->v[2]].i};
    if (eq_ctx->edge_queue_tri_in_range(eq_ctx, vs, f->no)) {
      f->head.hflag |= facetag;

      /* Check each edge of the face */
      BMLoop *l_first = BM_FACE_FIRST_LOOP(f);
      BMLoop *l_iter = l_first;
      do {
        /* are we owned by this node? if so, make sure origdata is up to date */
        if (BM_ELEM_CD_GET_INT(l_iter->v, pbvh->cd_vert_node_offset) == ni) {
          BKE_pbvh_bmesh_check_origdata(pbvh, l_iter->v, pbvh->stroke_id);
        }

        /* try to improve convergence by applying a small amount of smoothing to topology,
           but tangentially to surface.
         */
        int randval = (seed = dyntopo_thread_rand(seed)) & 255;

        if (do_smooth && randval > 127) {
          PBVHVertRef sv = {(intptr_t)l_iter->v};
          surface_smooth_v_safe(tdata->pbvh,
                                l_iter->v,
                                eq_ctx->surface_smooth_fac *
                                    eq_ctx->mask_cb(sv, eq_ctx->mask_cb_data));
        }

        float len_sq = calc_edge_length(eq_ctx, l_iter->e->v1, l_iter->e->v2);

        /* subdivide walks the mesh a bit for better transitions in the topology */

        if ((eq_ctx->mode & PBVH_Subdivide) && (len_sq > eq_ctx->limit_len_max_sqr)) {
          long_edge_queue_edge_add_recursive_2(
              tdata, l_iter->radial_next, l_iter, len_sq, eq_ctx->limit_len_max, 0, true);
        }

        if (edge_queue_test(eq_ctx, pbvh, l_iter->e)) {
          edge_thread_data_insert(tdata, l_iter->e);
        }
      } while ((l_iter = l_iter->next) != l_first);
    }
  }

  BLI_rng_free(rng);
}

ATTR_NO_OPT bool check_face_is_tri(PBVH *pbvh, BMFace *f)
{
#if DYNTOPO_DISABLE_FLAG & DYNTOPO_DISABLE_TRIANGULATOR
  return true;
#endif

  if (f->len == 3) {
    return true;
  }

  if (f->len < 3) {
    printf("pbvh had < 3 vert face!\n");
    BKE_pbvh_bmesh_remove_face(pbvh, f, false);
    return false;
  }

  bm_logstack_push();

  LinkNode *dbl = nullptr;

  Vector<BMFace *, 32> fs;
  Vector<BMEdge *, 32> es;

  BMLoop *l = f->l_first;
  do {
    validate_vert(pbvh, pbvh->header.bm, l->v, true, true);

    if (l->e->head.index == -1) {
      l->e->head.index = 0;
    }
  } while ((l = l->next) != f->l_first);

  // BKE_pbvh_bmesh_remove_face(pbvh, f, true);
  pbvh_bmesh_face_remove(pbvh, f, false, true, true);
  BM_log_face_removed(pbvh->header.bm, pbvh->bm_log, f);
  BM_idmap_release(pbvh->bm_idmap, (BMElem *)f, true);

  int len = (f->len - 2) * 3;

  fs.resize(len);
  es.resize(len);

  int totface = 0;
  int totedge = 0;
  MemArena *arena = nullptr;
  struct Heap *heap = nullptr;

  arena = BLI_memarena_new(512, "ngon arena");
  heap = BLI_heap_new();

  BM_face_triangulate(pbvh->header.bm,
                      f,
                      fs.data(),
                      &totface,
                      es.data(),
                      &totedge,
                      &dbl,
                      MOD_TRIANGULATE_QUAD_FIXED,
                      MOD_TRIANGULATE_NGON_BEAUTY,
                      false,
                      arena,
                      heap);

  while (totface && dbl) {
    BMFace *f2 = (BMFace *)dbl->link;
    LinkNode *next = dbl->next;

    for (int i = 0; i < totface; i++) {
      if (fs[i] == f2) {
        fs[i] = nullptr;
      }
    }

    if (f == f2) {
      BM_log_face_added(pbvh->header.bm, pbvh->bm_log, f);
      BM_log_face_removed_no_check(pbvh->header.bm, pbvh->bm_log, f);
      f = nullptr;
    }

    BM_idmap_release(pbvh->bm_idmap, (BMElem *)dbl->link, true);
    BM_face_kill(pbvh->header.bm, (BMFace *)dbl->link);

    MEM_freeN(dbl);
    dbl = next;
  }

  for (int i = 0; i < totface; i++) {
    BMFace *f2 = fs[i];

    if (!f2) {
      continue;
    }

    if (f == f2) {
      printf("%s: error\n", __func__);
      continue;
    }

    // detect new edges
    BMLoop *l = f2->l_first;
    do {
      if (l->e->head.index == -1) {
        BM_log_edge_added(pbvh->header.bm, pbvh->bm_log, l->e);
        l->e->head.index = 0;
      }
    } while ((l = l->next) != f2->l_first);

    validate_face(pbvh, pbvh->header.bm, f2, false, true);

    BKE_pbvh_bmesh_add_face(pbvh, f2, false, true);
    // BM_log_face_post(pbvh->bm_log, f2);
    BM_log_face_added(pbvh->header.bm, pbvh->bm_log, f2);
  }

  if (f) {
    BKE_pbvh_bmesh_add_face(pbvh, f, false, true);
    BM_log_face_added(pbvh->header.bm, pbvh->bm_log, f);
  }

  if (arena) {
    BLI_memarena_free(arena);
  }

  if (heap) {
    BLI_heap_free(heap, nullptr);
  }

  bm_logstack_pop();

  return false;
}

bool destroy_nonmanifold_fins(PBVH *pbvh, BMEdge *e_root)
{
#if !(DYNTOPO_DISABLE_FLAG & DYNTOPO_DISABLE_FIN_REMOVAL)
  bm_logstack_push();

  static int max_faces = 64;
  Vector<BMFace *, 32> stack;

  BMLoop *l = e_root->l;
  Vector<BMLoop *, 5> ls;
  Vector<BMFace *, 32> minfs;

  if (!l) {
    bm_logstack_pop();
    return false;
  }

  do {
    ls.append(l);
  } while ((l = l->radial_next) != e_root->l);

  for (int i = 0; i < ls.size(); i++) {
    Set<BMFace *, 64> visit;

    BMLoop *l = ls[i];
    BMFace *f = l->f;
    Vector<BMFace *, 32> fs2;
    stack.clear();

    stack.append(f);
    fs2.append(f);

    visit.add(f);

    bool bad = false;

    while (stack.size() > 0) {
      f = stack.pop_last();
      BMLoop *l = f->l_first;

      do {
        if (l->radial_next == l || l->radial_next->radial_next != l) {
          continue;
        }

        BMFace *f2 = l->radial_next->f;

        if (visit.add(f2)) {
          if (fs2.size() > max_faces) {
            bad = true;
            break;
          }

          stack.append(f2);
          fs2.append(f2);
        }
      } while ((l = l->next) != f->l_first);

      if (bad) {
        break;
      }
    }

    if (!bad && fs2.size() && (minfs.size() == 0 || fs2.size() < minfs.size())) {
      minfs = fs2;
    }
  }

  int nupdateflag = PBVH_UpdateOtherVerts | PBVH_UpdateDrawBuffers | PBVH_UpdateBB |
                    PBVH_UpdateTriAreas;
  nupdateflag = nupdateflag | PBVH_UpdateNormals | PBVH_UpdateTris | PBVH_RebuildDrawBuffers;

  if (!minfs.size()) {
    bm_logstack_pop();
    return false;
  }

  printf("manifold fin size: %d\n", (int)minfs.size());
  const int tag = BM_ELEM_TAG_ALT;

  for (int i = 0; i < minfs.size(); i++) {
    BMFace *f = minfs[i];

    BMLoop *l = f->l_first;
    do {
      BMLoop *l2 = l;
      do {
        BMLoop *l3 = l2;
        do {
          l3->v->head.hflag &= ~tag;
          l3->e->head.hflag &= ~tag;
        } while ((l3 = l3->next) != l2);
      } while ((l2 = l2->radial_next) != l);

      l->v->head.hflag &= ~tag;
      l->e->head.hflag &= ~tag;
    } while ((l = l->next) != f->l_first);
  }

  Vector<BMVert *, 32> vs;
  Vector<BMEdge *, 32> es;

  for (int i = 0; i < minfs.size(); i++) {
    BMFace *f = minfs[i];

    BMLoop *l = f->l_first;
    do {
      if (!(l->v->head.hflag & tag)) {
        l->v->head.hflag |= tag;
        vs.append(l->v);
      }

      if (!(l->e->head.hflag & tag)) {
        l->e->head.hflag |= tag;
        es.append(l->e);
      }
    } while ((l = l->next) != f->l_first);
  }

  for (int i = 0; i < minfs.size(); i++) {
    for (int j = 0; j < minfs.size(); j++) {
      if (i != j && minfs[i] == minfs[j]) {
        printf("%s: duplicate faces\n", __func__);
        continue;
      }
    }
  }

  for (int i = 0; i < minfs.size(); i++) {
    BMFace *f = minfs[i];

    if (f->head.htype != BM_FACE) {
      printf("%s: corruption!\n", __func__);
      continue;
    }

    int ni = BM_ELEM_CD_GET_INT(f, pbvh->cd_face_node_offset);
    if (ni >= 0 && ni < pbvh->totnode) {
      pbvh->nodes[ni].flag |= (PBVHNodeFlags)nupdateflag;
    }

    pbvh_bmesh_face_remove(pbvh, f, true, false, false);
    BM_idmap_release(pbvh->bm_idmap, (BMElem *)f, true);
    BM_face_kill(pbvh->header.bm, f);
  }

  const int mupdateflag = SCULPTVERT_NEED_DISK_SORT | SCULPTVERT_NEED_VALENCE;

  for (int i = 0; i < es.size(); i++) {
    BMEdge *e = es[i];

    if (!e->l) {
      BM_log_edge_removed(pbvh->header.bm, pbvh->bm_log, e);
      BM_idmap_release(pbvh->bm_idmap, (BMElem *)e, true);
      BM_edge_kill(pbvh->header.bm, e);
    }
  }

  for (int i = 0; i < vs.size(); i++) {
    BMVert *v = vs[i];

    if (!v->e) {
      pbvh_bmesh_vert_remove(pbvh, v);

      BM_log_vert_removed(pbvh->header.bm, pbvh->bm_log, v);
      BM_idmap_release(pbvh->bm_idmap, (BMElem *)v, true);
      BM_vert_kill(pbvh->header.bm, v);
    }
    else {
      pbvh_boundary_update_bmesh(pbvh, v);

      MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, v);
      mv->flag |= mupdateflag;
    }
  }

  bm_logstack_pop();
  return true;
#else
  return false;
#endif
}

bool check_for_fins(PBVH *pbvh, BMVert *v)
{
  BMEdge *e = v->e;
  if (!e) {
    return false;
  }

  bm_logstack_push();

  do {
    if (!e) {
      printf("%s: e was nullptr\n", __func__);
      break;
    }
    if (e->l) {
      BMLoop *l = e->l->f->l_first;

      do {
        if (l != l->radial_next && l != l->radial_next->radial_next) {
          if (destroy_nonmanifold_fins(pbvh, e)) {
            bm_logstack_pop();
            return true;
          }
        }
      } while ((l = l->next) != e->l->f->l_first);
    }
  } while ((e = BM_DISK_EDGE_NEXT(e, v)) != v->e);

  bm_logstack_pop();
  return false;
}

bool check_vert_fan_are_tris(PBVH *pbvh, BMVert *v)
{
  MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, v);

  if (!(mv->flag & SCULPTVERT_NEED_TRIANGULATE)) {
    return true;
  }

  bm_log_message("  == triangulate == ");

  Vector<BMFace *, 32> fs;

  validate_vert(pbvh, pbvh->header.bm, v, true, true);

  if (v->head.htype != BM_VERT) {
    printf("non-vert %p fed to %s\n", v, __func__);
    return false;
  }

  BMIter iter;
  BMFace *f;
  BM_ITER_ELEM (f, &iter, v, BM_FACES_OF_VERT) {
    BMLoop *l = f->l_first;

    do {
      MSculptVert *mv_l = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, l->v);
      pbvh_boundary_update_bmesh(pbvh, l->v);

      MV_ADD_FLAG(mv_l, SCULPTVERT_NEED_VALENCE | SCULPTVERT_NEED_DISK_SORT);
    } while ((l = l->next) != f->l_first);
    fs.append(f);
  }

  mv->flag &= ~SCULPTVERT_NEED_TRIANGULATE;

  for (int i = 0; i < fs.size(); i++) {
    check_face_is_tri(pbvh, fs[i]);
  }

  return false;
}

static void edge_queue_init(EdgeQueueContext *eq_ctx,
                            bool use_projected,
                            bool use_frontface,
                            const float center[3],
                            const float view_normal[3],
                            const float radius)
{
  if (use_projected) {
    eq_ctx->edge_queue_tri_in_range = edge_queue_tri_in_circle;
    eq_ctx->edge_queue_vert_in_range = edge_queue_vert_in_circle;
    project_plane_normalized_v3_v3v3(eq_ctx->center_proj, center, view_normal);
  }
  else {
    eq_ctx->edge_queue_tri_in_range = edge_queue_tri_in_sphere;
    eq_ctx->edge_queue_vert_in_range = edge_queue_vert_in_sphere;
  }

  eq_ctx->center = center;
  copy_v3_v3(eq_ctx->view_normal, view_normal);
  eq_ctx->radius_squared = radius * radius;
  eq_ctx->limit_len_min_sqr = eq_ctx->limit_len_min * eq_ctx->limit_len_min;
  eq_ctx->limit_len_max_sqr = eq_ctx->limit_len_max * eq_ctx->limit_len_max;

#ifdef USE_EDGEQUEUE_FRONTFACE
  eq_ctx->use_view_normal = use_frontface;
#else
  UNUSED_VARS(use_frontface);
#endif
}

static bool edge_queue_test(EdgeQueueContext *eq_ctx, PBVH *pbvh, BMEdge *e)
{
  float len = len_squared_v3v3(e->v1->co, e->v2->co);
  float min = eq_ctx->limit_len_min_sqr;
  float max = eq_ctx->limit_len_max_sqr;

  bool ret = false;

  if (eq_ctx->mode & PBVH_Subdivide) {
    ret |= len > max;
  }

  if (eq_ctx->mode & PBVH_Collapse) {
    ret |= len < min;
  }

  return ret;
}

/* Create a priority queue containing vertex pairs connected by a long
 * edge as defined by PBVH.bm_max_edge_len.
 *
 * Only nodes marked for topology update are checked, and in those
 * nodes only edges used by a face intersecting the (center, radius)
 * sphere are checked.
 *
 * The highest priority (lowest number) is given to the longest edge.
 */
static void unified_edge_queue_create(EdgeQueueContext *eq_ctx,
                                      PBVH *pbvh,
                                      const float center[3],
                                      const float view_normal[3],
                                      float radius,
                                      const bool use_frontface,
                                      const bool use_projected,
                                      PBVHTopologyUpdateMode local_mode)
{
  if (local_mode) {
    edge_queue_create_local(
        eq_ctx, pbvh, center, view_normal, radius, use_frontface, use_projected, local_mode);
    return;
  }

  eq_ctx->radius_squared = radius * radius;

  eq_ctx->limit_len_min = pbvh->bm_min_edge_len;
  eq_ctx->limit_len_max = pbvh->bm_max_edge_len;

  edge_queue_init(eq_ctx, use_projected, use_frontface, center, eq_ctx->view_normal, radius);

#ifdef USE_EDGEQUEUE_TAG_VERIFY
  pbvh_bmesh_edge_tag_verify(pbvh);
#endif

  Vector<EdgeQueueThreadData> tdata;

  bool push_subentry = false;

  for (int n = 0; n < pbvh->totnode; n++) {
    PBVHNode *node = &pbvh->nodes[n];

    /* Check leaf nodes marked for topology update */
    bool ok = ((node->flag & PBVH_Leaf) && (node->flag & PBVH_UpdateTopology) &&
               !(node->flag & PBVH_FullyHidden));

    if (!ok) {
      continue;
    }

    EdgeQueueThreadData td;

    memset(&td, 0, sizeof(td));

    td.seed = BLI_thread_rand(0);
    td.pbvh = pbvh;
    td.node = node;
    td.eq_ctx = eq_ctx;

    tdata.append(td);

    /* Check each face */
    /*
    BMFace *f;
    TGSET_ITER (f, node->bm_faces) {
      long_edge_queue_face_add(eq_ctx, f);
    }
    TGSET_ITER_END
    */
  }

  int count = tdata.size();

  TaskParallelSettings settings;

  BLI_parallel_range_settings_defaults(&settings);
#ifdef DYNTOPO_NO_THREADING
  settings.use_threading = false;
#endif

  BLI_task_parallel_range(0, count, (void *)tdata.data(), unified_edge_queue_task_cb, &settings);

  const int cd_sculpt_vert = pbvh->cd_sculpt_vert;

  for (int i = 0; i < count; i++) {
    EdgeQueueThreadData *td = &tdata[i];
    BMEdge **edges = td->edges;

    for (int j = 0; j < td->totedge; j++) {
      BMEdge *e = edges[j];

      e->head.hflag &= ~BM_ELEM_TAG;
    }
  }

  Vector<BMVert *> verts;

  for (int i = 0; i < count; i++) {
    EdgeQueueThreadData *td = &tdata[i];

    BMEdge **edges = td->edges;
    for (int j = 0; j < td->totedge; j++) {
      BMEdge *e = edges[j];

      if (bm_elem_is_free((BMElem *)e, BM_EDGE)) {
        continue;
      }

      e->head.hflag &= ~BM_ELEM_TAG;

      if (e->l && e->l != e->l->radial_next->radial_next) {
        // deal with non-manifold iffyness
        destroy_nonmanifold_fins(pbvh, e);
        push_subentry = true;

        if (bm_elem_is_free((BMElem *)e, BM_EDGE)) {
          continue;
        }
      }

      MSculptVert *mv1 = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, e->v1);
      MSculptVert *mv2 = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, e->v2);

      if (mv1->flag & SCULPTVERT_NEED_VALENCE) {
        BKE_pbvh_bmesh_update_valence(pbvh->cd_sculpt_vert, {(intptr_t)e->v1});
      }

      if (mv2->flag & SCULPTVERT_NEED_VALENCE) {
        BKE_pbvh_bmesh_update_valence(pbvh->cd_sculpt_vert, {(intptr_t)e->v2});
      }

      // check_vert_fan_are_tris(pbvh, e->v1);
      // check_vert_fan_are_tris(pbvh, e->v2);

      // float w = -calc_weighted_edge_split(eq_ctx, e->v1, e->v2);
      // float w2 = maskcb_get(eq_ctx, e);

      // w *= w2 * w2;

      if (eq_ctx->use_view_normal && (dot_v3v3(e->v1->no, eq_ctx->view_normal) < 0.0f &&
                                      dot_v3v3(e->v2->no, eq_ctx->view_normal) < 0.0f)) {
        return;
      }

      verts.append(e->v1);
      verts.append(e->v2);

      e->v1->head.hflag |= BM_ELEM_TAG;
      e->v2->head.hflag |= BM_ELEM_TAG;

      if (edge_queue_test(eq_ctx, pbvh, e)) {
        edge_queue_insert_unified(eq_ctx, e);
      }

      // edge_queue_insert(eq_ctx, e, w, eq_ctx->limit_len);
    }

    MEM_SAFE_FREE(td->edges);
  }

  for (int i = 0; i < verts.size(); i++) {
    BMVert *v = verts[i];

    if (v->head.hflag & BM_ELEM_TAG) {
      v->head.hflag &= ~BM_ELEM_TAG;

      edge_queue_insert_val34_vert(eq_ctx, v);
    }
  }

  if (push_subentry) {
    BM_log_entry_add_ex(pbvh->header.bm, pbvh->bm_log, true);
  }
}

static void short_edge_queue_task_cb_local(void *__restrict userdata,
                                           const int n,
                                           const TaskParallelTLS *__restrict tls)
{
  EdgeQueueThreadData *tdata = ((EdgeQueueThreadData *)userdata) + n;
  PBVHNode *node = tdata->node;
  EdgeQueueContext *eq_ctx = tdata->eq_ctx;

  BMFace *f;

  TGSET_ITER (f, node->bm_faces) {
#ifdef USE_EDGEQUEUE_FRONTFACE
    if (eq_ctx->use_view_normal) {
      if (dot_v3v3(f->no, eq_ctx->view_normal) < 0.0f) {
        continue;
      }
    }
#endif

    BMVert *vs[3] = {f->l_first->v, f->l_first->next->v, f->l_first->next->next->v};
    if (eq_ctx->edge_queue_tri_in_range(eq_ctx, vs, f->no)) {
      BMLoop *l = f->l_first;

      do {
        edge_thread_data_insert(tdata, l->e);

      } while ((l = l->next) != f->l_first);
    }
  }
  TGSET_ITER_END
}

static void edge_queue_create_local(EdgeQueueContext *eq_ctx,
                                    PBVH *pbvh,
                                    const float center[3],
                                    const float view_normal[3],
                                    float radius,
                                    const bool use_frontface,
                                    const bool use_projected,
                                    PBVHTopologyUpdateMode local_mode)
{
  eq_ctx->limit_len_min = pbvh->bm_min_edge_len;
  eq_ctx->limit_len_max = pbvh->bm_max_edge_len;
  eq_ctx->local_mode = true;

  edge_queue_init(eq_ctx, use_projected, use_frontface, center, eq_ctx->view_normal, radius);

  Vector<EdgeQueueThreadData> tdata;

  for (int n = 0; n < pbvh->totnode; n++) {
    PBVHNode *node = &pbvh->nodes[n];
    EdgeQueueThreadData td;

    if ((node->flag & PBVH_Leaf) && (node->flag & PBVH_UpdateTopology) &&
        !(node->flag & PBVH_FullyHidden)) {
      memset(&td, 0, sizeof(td));
      td.pbvh = pbvh;
      td.node = node;
      td.is_collapse = local_mode & PBVH_LocalCollapse;
      td.eq_ctx = eq_ctx;

      tdata.append(td);
    }
  }

  int count = tdata.size();

  TaskParallelSettings settings;

  BLI_parallel_range_settings_defaults(&settings);
#ifdef DYNTOPO_NO_THREADING
  settings.use_threading = false;
#endif

  BLI_task_parallel_range(
      0, count, (void *)tdata.data(), short_edge_queue_task_cb_local, &settings);

  const int cd_sculpt_vert = pbvh->cd_sculpt_vert;
  Vector<float> lens;
  Vector<BMEdge *> edges;

  for (int i = 0; i < count; i++) {
    EdgeQueueThreadData *td = &tdata[i];

    BMEdge **edges2 = td->edges;
    for (int j = 0; j < td->totedge; j++) {
      edges2[j]->head.hflag &= ~BM_ELEM_TAG;
    }
  }

  for (int i = 0; i < count; i++) {
    EdgeQueueThreadData *td = &tdata[i];

    BMEdge **edges2 = td->edges;
    for (int j = 0; j < td->totedge; j++) {
      BMEdge *e = edges2[j];

      e->v1->head.hflag &= ~BM_ELEM_TAG;
      e->v2->head.hflag &= ~BM_ELEM_TAG;

      if (!(e->head.hflag & BM_ELEM_TAG)) {
        edges.append(e);
        e->head.hflag |= BM_ELEM_TAG;
      }
    }
  }

  for (int i = 0; i < count; i++) {
    EdgeQueueThreadData *td = &tdata[i];
    MEM_SAFE_FREE(td->edges);
  }

  for (int i = 0; i < edges.size(); i++) {
    BMEdge *e = edges[i];
    float len = len_v3v3(e->v1->co, e->v2->co);

    for (int j = 0; j < 2; j++) {
      BMVert *v = j ? e->v2 : e->v1;

      if (!(local_mode & PBVH_LocalCollapse)) {
        if (!(v->head.hflag & BM_ELEM_TAG)) {
          v->head.hflag |= BM_ELEM_TAG;
          MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, v);

          if (mv->flag & SCULPTVERT_NEED_VALENCE) {
            BKE_pbvh_bmesh_update_valence(pbvh->cd_sculpt_vert, {(intptr_t)v});
          }

          edge_queue_insert_val34_vert(eq_ctx, v);
        }
      }
    }
    e->head.index = i;
    lens.append(len);
  }

  // make sure tags around border edges are unmarked
  for (int i = 0; i < edges.size(); i++) {
    BMEdge *e = edges[i];

    for (int j = 0; j < 2; j++) {
      BMVert *v1 = j ? e->v2 : e->v1;
      BMEdge *e1 = v1->e;

      do {
        e1->head.hflag &= ~BM_ELEM_TAG;

        e1 = BM_DISK_EDGE_NEXT(e1, v1);
      } while (e1 != v1->e);
    }
  }

  // re-tag edge list
  for (int i = 0; i < edges.size(); i++) {
    edges[i]->head.hflag |= BM_ELEM_TAG;
  }

  int totstep = 3;

  // blur lengths
  for (int step = 0; step < totstep; step++) {
    for (int i = 0; i < edges.size(); i++) {
      BMEdge *e = edges[i];

      float len = lens[i];
      float totlen = 0.0f;

      for (int j = 0; j < 2; j++) {
        BMVert *v1 = j ? e->v2 : e->v1;
        BMEdge *e1 = v1->e;

        do {
          if (e1->head.hflag & BM_ELEM_TAG) {
            len += lens[e1->head.index];
            totlen += 1.0f;
          }

          e1 = BM_DISK_EDGE_NEXT(e1, v1);
        } while (e1 != v1->e);
      }

      if (totlen != 0.0f) {
        len /= totlen;
        lens[i] += (len - lens[i]) * 0.5;
      }
    }
  }

  pbvh->header.bm->elem_index_dirty |= BM_EDGE;

  float limit = 0.0f;
  float tot = 0.0f;

  for (int i = 0; i < edges.size(); i++) {
    BMEdge *e = edges[i];
    MSculptVert *mv1, *mv2;

    e->head.hflag &= ~BM_ELEM_TAG;

    mv1 = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, e->v1);
    mv2 = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, e->v2);

    pbvh_check_vert_boundary(pbvh, e->v1);
    pbvh_check_vert_boundary(pbvh, e->v2);

    int boundflag1 = BM_ELEM_CD_GET_INT(e->v1, pbvh->cd_boundary_flag);
    int boundflag2 = BM_ELEM_CD_GET_INT(e->v2, pbvh->cd_boundary_flag);

    if ((boundflag1 & SCULPTVERT_ALL_CORNER) || (boundflag2 & SCULPTVERT_ALL_CORNER)) {
      continue;
    }

    if ((boundflag1 & SCULPTVERT_ALL_BOUNDARY) != (boundflag2 & SCULPTVERT_ALL_BOUNDARY)) {
      continue;
    }

    // check seam/sharp flags here
    // if (!(e->head.hflag & BM_ELEM_SMOOTH) || e->head.hflag & BM_ELEM_SEAM) {
    //  continue;
    // }

    limit += lens[i];
    tot += 1.0f;
  }

  if (tot > 0.0f) {
    limit /= tot;

    if (local_mode & PBVH_LocalCollapse) {
      eq_ctx->limit_len_min = limit * pbvh->bm_detail_range;
      eq_ctx->limit_len_min_sqr = eq_ctx->limit_len_min * eq_ctx->limit_len_min;
    }

    if (local_mode & PBVH_LocalSubdivide) {
      eq_ctx->limit_len_max = limit;
      eq_ctx->limit_len_max_sqr = eq_ctx->limit_len_max * eq_ctx->limit_len_max;
    }
  }

  for (int i = 0; i < edges.size(); i++) {
    BMEdge *e = edges[i];

    MSculptVert *mv1 = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, e->v1);
    MSculptVert *mv2 = BKE_PBVH_SCULPTVERT(cd_sculpt_vert, e->v2);

    int boundflag1 = BM_ELEM_CD_GET_INT(e->v1, pbvh->cd_boundary_flag);
    int boundflag2 = BM_ELEM_CD_GET_INT(e->v2, pbvh->cd_boundary_flag);

    if ((boundflag1 & SCULPTVERT_ALL_CORNER) || (boundflag2 & SCULPTVERT_ALL_CORNER)) {
      continue;
    }

    if ((boundflag1 & SCULPTVERT_ALL_BOUNDARY) != (boundflag2 & SCULPTVERT_ALL_BOUNDARY)) {
      continue;
    }

    bool ok = false;

    bool a = eq_ctx->mode & (PBVH_Subdivide | PBVH_LocalSubdivide);
    bool b = eq_ctx->mode & (PBVH_Collapse | PBVH_LocalCollapse);

    float len = len_squared_v3v3(e->v1->co, e->v2->co);

    if (a && b) {
      ok = len < eq_ctx->limit_len_min_sqr || len > eq_ctx->limit_len_max_sqr;
      ok = ok || (len < pbvh->bm_min_edge_len || len > pbvh->bm_max_edge_len);
    }
    else if (a) {
      ok = len > eq_ctx->limit_len_max || len > pbvh->bm_max_edge_len;
    }
    else if (b) {
      ok = len < eq_ctx->limit_len_min || len < pbvh->bm_min_edge_len;
    }

    if (!ok) {
      continue;
    }

    edge_queue_insert_unified(eq_ctx, e);
  }
}

// need to file a CLANG bug, getting weird behavior here
#ifdef __clang__
__attribute__((optnone))
#endif

static bool
cleanup_valence_3_4(EdgeQueueContext *ectx,
                    PBVH *pbvh,
                    const float center[3],
                    const float view_normal[3],
                    float radius,
                    const bool use_frontface,
                    const bool use_projected)
{
  return false;  // XXX not working with local collapse/subdivide mode
  bool modified = false;

  bm_logstack_push();

  bm_log_message("  == cleanup_valence_3_4 == ");

  // push log subentry
  BM_log_entry_add_ex(pbvh->header.bm, pbvh->bm_log, true);

  float radius2 = radius * 1.25;
  float rsqr = radius2 * radius2;

  const int cd_vert_node = pbvh->cd_vert_node_offset;

  int updateflag = SCULPTVERT_NEED_DISK_SORT | SCULPTVERT_NEED_VALENCE;

  for (int i = 0; i < ectx->tot_used_verts; i++) {
    BMVert *v = ectx->used_verts[i];

    // TGSET_ITER (v, ectx->used_verts) {
    if (bm_elem_is_free((BMElem *)v, BM_VERT)) {
      continue;
    }

    const int n = BM_ELEM_CD_GET_INT(v, cd_vert_node);

    if (n == DYNTOPO_NODE_NONE) {
      continue;
    }

    MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, v);

    BKE_pbvh_bmesh_check_valence(pbvh, {(intptr_t)v});

    int val = mv->valence;
    if (val != 4 && val != 3) {
      continue;
    }

    PBVHVertRef sv = {(intptr_t)v};

    if (len_squared_v3v3(v->co, center) >= rsqr || !v->e ||
        ectx->mask_cb(sv, ectx->mask_cb_data) < 0.5f) {
      continue;
    }

    validate_vert(pbvh, pbvh->header.bm, v, false, true);
    check_vert_fan_are_tris(pbvh, v);
    validate_vert(pbvh, pbvh->header.bm, v, true, true);

#if 0
    // check valence again
    val = BM_vert_edge_count(v);
    if (val != 4 && val != 3) {
      printf("pbvh valence error\n");
      continue;
    }
#endif

    pbvh_check_vert_boundary(pbvh, v);
    int boundflag = BM_ELEM_CD_GET_INT(v, pbvh->cd_boundary_flag);

    if (boundflag & SCULPTVERT_ALL_BOUNDARY) {
      continue;
    }

    BMIter iter;
    BMLoop *l;
    BMLoop *ls[4];
    BMVert *vs[4];

    l = v->e->l;

    if (!l) {
      continue;
    }

    if (l->v != v) {
      l = l->next;
    }

    bool bad = false;
    int ls_i = 0;

    for (int j = 0; j < val; j++) {
      ls[ls_i++] = l->v == v ? l->next : l;

      MSculptVert *mv_l;

      if (l->v == v) {
        mv_l = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, l->next->v);
        pbvh_boundary_update_bmesh(pbvh, l->next->v);
      }
      else {
        mv_l = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, l->v);
        pbvh_boundary_update_bmesh(pbvh, l->v);
      }

      MV_ADD_FLAG(mv_l, updateflag);

      l = l->prev->radial_next;

      if (l->v != v) {
        l = l->next;
      }

      /*ignore non-manifold edges along with ones flagged as sharp*/
      if (l->radial_next == l || l->radial_next->radial_next != l ||
          !(l->e->head.hflag & BM_ELEM_SMOOTH)) {
        bad = true;
        break;
      }

      if (l->radial_next != l && l->radial_next->v == l->v) {
        bad = true;  // bad normals
        break;
      }

      for (int k = 0; k < j; k++) {
        if (ls[k]->v == ls[j]->v) {
          if (ls[j]->next->v != v) {
            ls[j] = ls[j]->next;
          }
          else {
            bad = true;
            break;
          }
        }

        // check for non-manifold edges
        if (ls[k] != ls[k]->radial_next->radial_next) {
          bad = true;
          break;
        }

        if (ls[k]->f == ls[j]->f) {
          bad = true;
          break;
        }
      }
    }

    if (bad) {
      continue;
    }

    int ni = BM_ELEM_CD_GET_INT(v, pbvh->cd_vert_node_offset);

    if (ni < 0) {
      printf("cleanup_valence_3_4 error!\n");

      // attempt to recover

      BMFace *f;
      BM_ITER_ELEM (f, &iter, v, BM_FACES_OF_VERT) {
        int ni2 = BM_ELEM_CD_GET_INT(f, pbvh->cd_face_node_offset);

        if (ni2 != DYNTOPO_NODE_NONE) {
          PBVHNode *node2 = pbvh->nodes + ni2;

          BLI_table_gset_remove(node2->bm_unique_verts, v, nullptr);
        }
      }
    }

    pbvh_bmesh_vert_remove(pbvh, v);

    BMFace *f;
    BM_ITER_ELEM (f, &iter, v, BM_FACES_OF_VERT) {
      int ni2 = BM_ELEM_CD_GET_INT(f, pbvh->cd_face_node_offset);

      if (ni2 != DYNTOPO_NODE_NONE) {
        // PBVHNode *node2 = pbvh->nodes + ni2;
        // BLI_table_gset_remove(node2->bm_unique_verts, v, nullptr);

        pbvh_bmesh_face_remove(pbvh, f, true, true, true);
      }
      else {
        BM_log_face_removed(pbvh->header.bm, pbvh->bm_log, f);
      }
    }

    modified = true;

    if (!v->e) {
      printf("mesh error!\n");
      continue;
    }

    validate_vert(pbvh, pbvh->header.bm, v, false, true);

    l = v->e->l;

    if (val == 4) {
      // check which quad diagonal to use to split quad
      // try to preserve hard edges

      float n1[3], n2[3], th1, th2;
      normal_tri_v3(n1, ls[0]->v->co, ls[1]->v->co, ls[2]->v->co);
      normal_tri_v3(n2, ls[0]->v->co, ls[2]->v->co, ls[3]->v->co);

      th1 = dot_v3v3(n1, n2);

      normal_tri_v3(n1, ls[1]->v->co, ls[2]->v->co, ls[3]->v->co);
      normal_tri_v3(n2, ls[1]->v->co, ls[3]->v->co, ls[0]->v->co);

      th2 = dot_v3v3(n1, n2);

      if (th1 > th2) {
        BMLoop *ls2[4] = {ls[0], ls[1], ls[2], ls[3]};

        for (int j = 0; j < 4; j++) {
          ls[j] = ls2[(j + 1) % 4];
        }
      }
    }

    vs[0] = ls[0]->v;
    vs[1] = ls[1]->v;
    vs[2] = ls[2]->v;

    validate_vert(pbvh, pbvh->header.bm, v, false, false);

    MSculptVert *mv1 = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, vs[0]);
    MSculptVert *mv2 = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, vs[1]);
    MSculptVert *mv3 = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, vs[2]);

    pbvh_boundary_update_bmesh(pbvh, vs[0]);
    pbvh_boundary_update_bmesh(pbvh, vs[1]);
    pbvh_boundary_update_bmesh(pbvh, vs[2]);

    MV_ADD_FLAG(mv1, updateflag);
    MV_ADD_FLAG(mv2, updateflag);
    MV_ADD_FLAG(mv3, updateflag);

    BMFace *f1 = nullptr;
    bool ok1 = vs[0] != vs[1] && vs[1] != vs[2] && vs[0] != vs[2];
    ok1 = ok1 && !BM_face_exists(vs, 3);

    if (ok1) {
      f1 = pbvh_bmesh_face_create(pbvh, n, vs, nullptr, l->f, true, false);
      normal_tri_v3(
          f1->no, f1->l_first->v->co, f1->l_first->next->v->co, f1->l_first->prev->v->co);

      validate_face(pbvh, pbvh->header.bm, f1, false, false);
    }
    else {
      // printf("eek1!\n");
    }

    BMFace *f2 = nullptr;

    if (val == 4) {
      vs[0] = ls[0]->v;
      vs[1] = ls[2]->v;
      vs[2] = ls[3]->v;
    }

    bool ok2 = val == 4 && vs[0] != vs[2] && vs[2] != vs[3] && vs[0] != vs[3];
    ok2 = ok2 && !BM_face_exists(vs, 3);

    if (ok2) {
      MSculptVert *mv1 = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, vs[0]);
      MSculptVert *mv2 = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, vs[1]);
      MSculptVert *mv3 = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, vs[2]);

      pbvh_boundary_update_bmesh(pbvh, vs[0]);
      pbvh_boundary_update_bmesh(pbvh, vs[1]);
      pbvh_boundary_update_bmesh(pbvh, vs[2]);

      MV_ADD_FLAG(mv1, updateflag);
      MV_ADD_FLAG(mv2, updateflag);
      MV_ADD_FLAG(mv3, updateflag);

      BMFace *example = nullptr;
      if (v->e && v->e->l) {
        example = v->e->l->f;
      }

      f2 = pbvh_bmesh_face_create(pbvh, n, vs, nullptr, example, true, false);

      CustomData_bmesh_swap_data_simple(
          &pbvh->header.bm->ldata, &f2->l_first->prev->head.data, &ls[3]->head.data);
      CustomData_bmesh_copy_data(&pbvh->header.bm->ldata,
                                 &pbvh->header.bm->ldata,
                                 ls[0]->head.data,
                                 &f2->l_first->head.data);
      CustomData_bmesh_copy_data(&pbvh->header.bm->ldata,
                                 &pbvh->header.bm->ldata,
                                 ls[2]->head.data,
                                 &f2->l_first->next->head.data);

      normal_tri_v3(
          f2->no, f2->l_first->v->co, f2->l_first->next->v->co, f2->l_first->prev->v->co);
      BM_log_face_added(pbvh->header.bm, pbvh->bm_log, f2);

      validate_face(pbvh, pbvh->header.bm, f2, false, false);
    }

    if (f1) {
      CustomData_bmesh_swap_data_simple(
          &pbvh->header.bm->ldata, &f1->l_first->head.data, &ls[0]->head.data);
      CustomData_bmesh_swap_data_simple(
          &pbvh->header.bm->ldata, &f1->l_first->next->head.data, &ls[1]->head.data);
      CustomData_bmesh_swap_data_simple(
          &pbvh->header.bm->ldata, &f1->l_first->prev->head.data, &ls[2]->head.data);

      BM_log_face_added(pbvh->header.bm, pbvh->bm_log, f1);
    }

    validate_vert(pbvh, pbvh->header.bm, v, false, false);
    pbvh_kill_vert(pbvh, v, true, true);

    if (f1 && !bm_elem_is_free((BMElem *)f1, BM_FACE)) {
      if (!bm_elem_is_free((BMElem *)f1, BM_FACE)) {
        check_face_is_manifold(pbvh, pbvh->header.bm, f1);
      }
    }

    if (f2 && !bm_elem_is_free((BMElem *)f2, BM_FACE)) {
      if (!bm_elem_is_free((BMElem *)f2, BM_FACE)) {
        check_face_is_manifold(pbvh, pbvh->header.bm, f2);
      }
    }
  }
#ifdef DYNTOPO_USE_MINMAX_HEAP
  // TGSET_ITER_END;
#endif

  if (modified) {
    pbvh->header.bm->elem_index_dirty |= BM_VERT | BM_FACE | BM_EDGE;
    pbvh->header.bm->elem_table_dirty |= BM_VERT | BM_FACE | BM_EDGE;
  }

  bm_logstack_pop();

  return modified;
}

//#define DEFRAGMENT_MEMORY
bool BM_defragment_vertex(BMesh *bm,
                          BMVert *v,
                          RNG *rand,
                          void (*on_vert_swap)(BMVert *a, BMVert *b, void *userdata),
                          void *userdata);

typedef struct SwapData {
  PBVH *pbvh;
} SwapData;

void pbvh_tribuf_swap_verts(
    PBVH *pbvh, PBVHNode *node, PBVHNode *node2, PBVHTriBuf *tribuf, BMVert *v1, BMVert *v2)
{

  void *val;

  bool ok = BLI_smallhash_remove_p(&tribuf->vertmap, (uintptr_t)v1, &val);

  if (!ok) {
    // printf("eek, missing vertex!");
    return;
  }

  int idx = POINTER_AS_INT(val);

  tribuf->verts[idx].i = (intptr_t)v2;

  void **val2;
  if (BLI_smallhash_ensure_p(&tribuf->vertmap, (intptr_t)v2, &val2)) {
    // v2 was already in hash? add v1 back in then with v2's index

    int idx2 = POINTER_AS_INT(*val2);
    BLI_smallhash_insert(&tribuf->vertmap, (intptr_t)v1, POINTER_FROM_INT(idx2));
  }

  *val2 = POINTER_FROM_INT(idx);
}

void pbvh_node_tribuf_swap_verts(
    PBVH *pbvh, PBVHNode *node, PBVHNode *node2, BMVert *v1, BMVert *v2)
{
  if (node != node2) {
    BLI_table_gset_remove(node->bm_unique_verts, v1, nullptr);
    BLI_table_gset_add(node->bm_unique_verts, v2);
  }

  if (node->tribuf) {
    pbvh_tribuf_swap_verts(pbvh, node, node2, node->tribuf, v1, v2);
  }

  for (int i = 0; i < node->tot_tri_buffers; i++) {
    pbvh_tribuf_swap_verts(pbvh, node, node2, node->tri_buffers + i, v1, v2);
  }
}

static bool do_cleanup_3_4(EdgeQueueContext *eq_ctx,
                           PBVH *pbvh,
                           const float center[3],
                           const float view_normal[3],
                           float radius,
                           bool use_frontface,
                           bool use_projected)
{
  bool modified = false;

#if 1
  if (!(eq_ctx->mode &
        (PBVH_Subdivide | PBVH_Collapse | PBVH_LocalCollapse | PBVH_LocalSubdivide))) {
    for (int n = 0; n < pbvh->totnode; n++) {
      PBVHNode *node = pbvh->nodes + n;
      BMVert *v;

      if (!(node->flag & PBVH_Leaf) || !(node->flag & PBVH_UpdateTopology)) {
        continue;
      }

      TGSET_ITER (v, node->bm_unique_verts) {
        if (!eq_ctx->edge_queue_vert_in_range(eq_ctx, v)) {
          continue;
        }

        if (use_frontface && dot_v3v3(v->no, view_normal) < 0.0f) {
          continue;
        }

        MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, v);

        if (mv->flag & SCULPTVERT_NEED_VALENCE) {
          BKE_pbvh_bmesh_update_valence(pbvh->cd_sculpt_vert, {(intptr_t)v});
        }

        if (mv->valence < 5) {
          edge_queue_insert_val34_vert(eq_ctx, v);
        }
      }
      TGSET_ITER_END;
    }
  }
#endif

  BM_log_entry_add_ex(pbvh->header.bm, pbvh->bm_log, true);

  pbvh_bmesh_check_nodes(pbvh);

  modified |= cleanup_valence_3_4(
      eq_ctx, pbvh, center, view_normal, radius, use_frontface, use_projected);
  pbvh_bmesh_check_nodes(pbvh);

  return modified;
}

float mask_cb_nop(PBVHVertRef vertex, void *userdata)
{
  return 1.0f;
}

/* Collapse short edges, subdivide long edges */
extern "C" bool BKE_pbvh_bmesh_update_topology(PBVH *pbvh,
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
                                               int custom_max_steps,
                                               bool disable_surface_relax,
                                               bool is_snake_hook)
{
  /* Disable surface smooth if uv layers are present, to avoid expensive reprojection operation. */
  if (!is_snake_hook && CustomData_has_layer(&pbvh->header.bm->ldata, CD_PROP_FLOAT2)) {
    disable_surface_relax = true;
  }

  /* Push a subentry. */
  BM_log_entry_add_ex(pbvh->header.bm, pbvh->bm_log, true);

  /* 2 is enough for edge faces - manifold edge */
  BLI_buffer_declare_static(BMLoop *, edge_loops, BLI_BUFFER_NOP, 2);
  BLI_buffer_declare_static(BMFace *, deleted_faces, BLI_BUFFER_NOP, 32);

  const int cd_vert_mask_offset = CustomData_get_offset(&pbvh->header.bm->vdata, CD_PAINT_MASK);
  const int cd_vert_node_offset = pbvh->cd_vert_node_offset;
  const int cd_face_node_offset = pbvh->cd_face_node_offset;
  const int cd_sculpt_vert = pbvh->cd_sculpt_vert;
  // float ratio = 1.0f;

  bool modified = false;

  if (view_normal) {
    BLI_assert(len_squared_v3(view_normal) != 0.0f);
  }

#ifdef DYNTOPO_REPORT
  {
    BMesh *bm = pbvh->header.bm;

    int vmem = (int)((size_t)bm->totvert * (sizeof(BMVert) + bm->vdata.totsize));
    int emem = (int)((size_t)bm->totedge * (sizeof(BMEdge) + bm->edata.totsize));
    int lmem = (int)((size_t)bm->totloop * (sizeof(BMLoop) + bm->ldata.totsize));
    int fmem = (int)((size_t)bm->totface * (sizeof(BMFace) + bm->pdata.totsize));

    double fvmem = (double)vmem / 1024.0 / 1024.0;
    double femem = (double)emem / 1024.0 / 1024.0;
    double flmem = (double)lmem / 1024.0 / 1024.0;
    double ffmem = (double)fmem / 1024.0 / 1024.0;

    printf("totmem: %.2fmb\n", fvmem + femem + flmem + ffmem);
    printf("v: %.2f e: %.2f l: %.2f f: %.2f\n", fvmem, femem, flmem, ffmem);

    printf("custom attributes only:\n");
    vmem = (int)((size_t)bm->totvert * (bm->vdata.totsize));
    emem = (int)((size_t)bm->totedge * (bm->edata.totsize));
    lmem = (int)((size_t)bm->totloop * (bm->ldata.totsize));
    fmem = (int)((size_t)bm->totface * (bm->pdata.totsize));

    fvmem = (double)vmem / 1024.0 / 1024.0;
    femem = (double)emem / 1024.0 / 1024.0;
    flmem = (double)lmem / 1024.0 / 1024.0;
    ffmem = (double)fmem / 1024.0 / 1024.0;

    printf("v: %.2f e: %.2f l: %.2f f: %.2f\n", fvmem, femem, flmem, ffmem);
  }
#endif

  float safe_smooth;

  if (disable_surface_relax) {
    safe_smooth = 0.0f;
  }
  else if ((mode & PBVH_Subdivide) && (!(mode & PBVH_Collapse) || (mode & PBVH_LocalCollapse))) {
    safe_smooth = DYNTOPO_SAFE_SMOOTH_SUBD_ONLY_FAC;
  }
  else {
    safe_smooth = DYNTOPO_SAFE_SMOOTH_FAC;
  }

  EdgeQueueContext eq_ctx = {};

  eq_ctx.pool = nullptr;
  eq_ctx.bm = pbvh->header.bm;
  eq_ctx.mask_cb = mask_cb;
  eq_ctx.mask_cb_data = mask_cb_data;

  eq_ctx.cd_sculpt_vert = cd_sculpt_vert;
  eq_ctx.cd_vert_mask_offset = cd_vert_mask_offset;
  eq_ctx.cd_vert_node_offset = cd_vert_node_offset;
  eq_ctx.cd_face_node_offset = cd_face_node_offset;
  eq_ctx.avg_elen = 0.0f;
  eq_ctx.max_elen = -1e17;
  eq_ctx.min_elen = 1e17;
  eq_ctx.totedge = 0.0f;
  eq_ctx.local_mode = false;
  eq_ctx.surface_smooth_fac = safe_smooth;
  eq_ctx.mode = mode;

#ifdef DYNTOPO_USE_MINMAX_HEAP
  // eq_ctx.used_verts = BLI_table_gset_new(__func__);
  eq_ctx.max_heap_mm = DYNTOPO_MAX_ITER << 8;
  eq_ctx.limit_min = pbvh->bm_min_edge_len;
  eq_ctx.limit_max = pbvh->bm_max_edge_len;
  eq_ctx.limit_mid = eq_ctx.limit_max * 0.5f + eq_ctx.limit_min * 0.5f;

  eq_ctx.use_view_normal = use_frontface;
  if (view_normal) {
    copy_v3_v3(eq_ctx.view_normal, view_normal);
  }
  else {
    zero_v3(eq_ctx.view_normal);
    eq_ctx.view_normal[2] = 1.0f;
  }

#endif

  if (mode & PBVH_LocalSubdivide) {
    mode |= PBVH_Subdivide;
  }
  if (mode & PBVH_LocalSubdivide) {
    mode |= PBVH_Collapse;
  }

#if DYNTOPO_DISABLE_FLAG & DYNTOPO_DISABLE_COLLAPSE
  mode &= ~PBVH_Collapse;
#endif
#if DYNTOPO_DISABLE_FLAG & DYNTOPO_DISABLE_SPLIT_EDGES
  mode &= ~PBVH_Subdivide;
#endif

  if (mode & (PBVH_Subdivide | PBVH_Collapse)) {
    unified_edge_queue_create(&eq_ctx,
                              pbvh,
                              center,
                              view_normal,
                              radius,
                              use_frontface,
                              use_projected,
                              mode & (PBVH_LocalSubdivide | PBVH_LocalCollapse));
  }
  else {
    edge_queue_init(&eq_ctx, use_projected, use_frontface, center, eq_ctx.view_normal, radius);
  }

#ifdef SKINNY_EDGE_FIX
  // prevent remesher thrashing by throttling edge splitting in pathological case of skinny
  // edges
  float avg_elen = eq_ctx.avg_elen;
  if (eq_ctx.totedge > 0.0f) {
    avg_elen /= eq_ctx.totedge;

    float emax = eq_ctx.max_elen;
    if (emax == 0.0f) {
      emax = 0.0001f;
    }

    if (pbvh->bm_min_edge_len > 0.0f && avg_elen > 0.0f) {
      ratio = avg_elen / (pbvh->bm_min_edge_len * 0.5 + emax * 0.5);
      ratio = MAX2(ratio, 0.25f);
      ratio = MIN2(ratio, 5.0f);
    }
  }
#else
  // ratio = 1.0f;
#endif

  int steps[2] = {0, 0};

  if ((mode & PBVH_Subdivide) && (mode & PBVH_Collapse)) {
    steps[0] = 4096;
    steps[1] = 1024;
  }
  else if (mode & PBVH_Subdivide) {
    steps[0] = 4096;
  }
  else if (mode & PBVH_Collapse) {
    steps[0] = 4096;
  }

  int edges_size = steps[0];
  BMEdge **edges = (BMEdge **)MEM_malloc_arrayN(edges_size, sizeof(void *), __func__);
  int etot = 0;

  PBVHTopologyUpdateMode ops[2];
  int count = 0;
  int i = 0;
  int max_steps = max_ii(custom_max_steps, DYNTOPO_MAX_ITER) << 1;

  int totop = 0;

  if (mode & PBVH_Subdivide) {
    ops[totop++] = PBVH_Subdivide;
  }
  if (mode & PBVH_Collapse) {
    ops[totop++] = PBVH_Collapse;
  }

  int curop = 0;
  float limit_len_subd = eq_ctx.limit_len_max_sqr;
  float limit_len_cold = eq_ctx.limit_len_min_sqr;
  // limit_len_cold = limit_len_cold * limit_len_cold;

  // printf(" minmax queue size: %d\n", BLI_mm_heap_len(eq_ctx.heap_mm));

  SmallHash subd_edges;
  BLI_smallhash_init(&subd_edges);

  while (totop > 0 && !eq_ctx.edge_heap.empty() && i < max_steps) {
    BMEdge *e = nullptr;

    if (count >= steps[curop]) {
      if (ops[curop] == PBVH_Subdivide) {
        modified = true;
        BLI_smallhash_clear(&subd_edges, 0);
        pbvh_split_edges(&eq_ctx, pbvh, pbvh->header.bm, edges, etot, false);
        VALIDATE_LOG(pbvh->bm_log);
        etot = 0;
      }

      curop = (curop + 1) % totop;
      count = 0;
    }

    switch (ops[curop]) {
      case PBVH_Subdivide: {
        if (eq_ctx.edge_heap.max_weight() < limit_len_subd) {
          break;
        }

        e = eq_ctx.edge_heap.pop_max();
        while (!eq_ctx.edge_heap.empty() && e &&
               (bm_elem_is_free((BMElem *)e, BM_EDGE) ||
                calc_weighted_length(&eq_ctx, e->v1, e->v2, -1.0) < limit_len_subd)) {
          e = eq_ctx.edge_heap.pop_max();
        }

        if (!e || bm_elem_is_free((BMElem *)e, BM_EDGE)) {
          break;
        }

        /* add complete triangles */
        BMLoop *l = e->l;
        if (l) {
          do {
            BMLoop *l2 = l;
            do {
              if (etot >= edges_size) {
                break;
              }

              void **val = nullptr;

              if (!BLI_smallhash_ensure_p(&subd_edges, (uintptr_t)l->e, &val)) {
                *val = nullptr;
                edges[etot++] = e;
              }
            } while ((l2 = l2->next) != l);
          } while ((l = l->radial_next) != e->l);
        }

        // edges[etot++] = e;
        break;
      }
      case PBVH_Collapse: {
        if (eq_ctx.edge_heap.min_weight() > limit_len_cold) {
          break;
        }

        e = eq_ctx.edge_heap.pop_min();
        while (!eq_ctx.edge_heap.empty() && e &&
               (bm_elem_is_free((BMElem *)e, BM_EDGE) ||
                calc_weighted_length(&eq_ctx, e->v1, e->v2, 1.0) > limit_len_cold)) {
          e = eq_ctx.edge_heap.pop_min();
        }

        if (!e || bm_elem_is_free((BMElem *)e, BM_EDGE)) {
          break;
        }

        if (bm_elem_is_free((BMElem *)e->v1, BM_VERT) ||
            bm_elem_is_free((BMElem *)e->v2, BM_VERT)) {
          printf("%s: error! operated on freed bmesh elements! e: %p, e->v1: %p, e->v2: %p\n",
                 __func__,
                 e,
                 e->v1,
                 e->v2);
          break;
        }

        printf("\n");
        modified = true;
        pbvh_bmesh_collapse_edge(pbvh, e, e->v1, e->v2, nullptr, nullptr, &eq_ctx);
        VALIDATE_LOG(pbvh->bm_log);
        printf("\n");

        // XXX
        BM_log_entry_add_ex(pbvh->header.bm, pbvh->bm_log, true);
        break;
      }
      default:
        break;
    }

    count++;
    i++;
  }

  if (etot > 0) {
    modified = true;
    BLI_smallhash_clear(&subd_edges, 0);
    pbvh_split_edges(&eq_ctx, pbvh, pbvh->header.bm, edges, etot, false);
    VALIDATE_LOG(pbvh->bm_log);
    etot = 0;
  }

  MEM_SAFE_FREE(edges);
  BLI_smallhash_release(&subd_edges);

  if (mode & PBVH_Cleanup) {
    modified |= do_cleanup_3_4(
        &eq_ctx, pbvh, center, eq_ctx.view_normal, radius, use_frontface, use_projected);

    VALIDATE_LOG(pbvh->bm_log);
  }

  if (modified) {
    // avoid potential infinite loops
    const int totnode = pbvh->totnode;

    for (int i = 0; i < totnode; i++) {
      PBVHNode *node = pbvh->nodes + i;

      if ((node->flag & PBVH_Leaf) && (node->flag & PBVH_UpdateTopology) &&
          !(node->flag & PBVH_FullyHidden)) {

        /* do not clear PBVH_UpdateTopology here in case split messes with it */

        /* Recursively split nodes that have gotten too many
         * elements */
        if (updatePBVH) {
          pbvh_bmesh_node_limit_ensure(pbvh, i);
        }
      }
    }
  }

  /* clear PBVH_UpdateTopology flags */
  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (!(node->flag & PBVH_Leaf)) {
      continue;
    }

    node->flag &= ~PBVH_UpdateTopology;
  }

#ifdef DYNTOPO_USE_MINMAX_HEAP
  if (eq_ctx.used_verts) {
    MEM_SAFE_FREE(eq_ctx.used_verts);
    // BLI_table_gset_free(eq_ctx.used_verts, nullptr);
  }
#endif

  BLI_buffer_free(&edge_loops);
  BLI_buffer_free(&deleted_faces);

#ifdef USE_VERIFY
  pbvh_bmesh_verify(pbvh);
#endif

  // ensure triangulations are all up to date
  for (int i = 0; i < pbvh->totnode; i++) {
    PBVHNode *node = pbvh->nodes + i;

    if (node->flag & PBVH_Leaf) {
      //  BKE_pbvh_bmesh_check_tris(pbvh, node);
    }
  }

  /* Push a subentry. */
  BM_log_entry_add_ex(pbvh->header.bm, pbvh->bm_log, true);

  return modified;
}

#define SPLIT_TAG BM_ELEM_TAG

/*

#generate shifted and mirrored patterns
# [number of verts, vert_connections... ]
table = [
  [4,     -1,  3, -1, -1],
  [5,     -1,  3, -1,  0, -1],
  [6,     -1,  3, -1,  5, -1, 1]
]

table2 = {}

def getmask(row):
  mask = 0
  for i in range(len(row)):
    if row[i] >= 0:
      mask |= 1 << i
  return mask

ii = 0
for row in table:
  n = row[0]
  row = row[1:]

  mask = getmask(row)
  table2[mask] = [n] + row

  for step in range(2):
    for i in range(n):
      row2 = []
      for j in range(n):
        j2 = row[(j + i) % n]

        if j2 != -1:
          j2 = (j2 - i + n) % n

        row2.append(j2)

      if row2[0] != -1:
        continue

      mask2 = getmask(row2)
      if mask2 not in table2:
        table2[mask2] = [n] + row2

    #reverse row
    for i in range(n):
      if row[i] != -1:
        row[i] = n - row[i]

    row.reverse()

maxk = 0
for k in table2:
  maxk = max(maxk, k)

buf = 'static const int splitmap[%i][16] = {\n' % (maxk+1)
buf += '  //{numverts, vert_connections...}\n'

for k in range(maxk+1):
  line = ""

  if k not in table2:
    line += '  {-1},'
  else:
    line += '  {'
    row = table2[k]
    for j in range(len(row)):
      if j > 0:
        line += ", "
      line += str(row[j])
    line += '},'

  while len(line) < 35:
    line += " "
  line += "//" + str(k) + " "

  if k in table2:
    for i in range(table2[k][0]):
      ch = "1" if k & (1 << i) else "0"
      line += str(ch) + " "

  buf += line + "\n"
buf += '};\n'
print(buf)

*/
static const int splitmap[43][16] = {
    //{numverts, vert_connections...}
    {-1},                      // 0
    {-1},                      // 1
    {4, -1, 3, -1, -1},        // 2 0 1 0 0
    {-1},                      // 3
    {4, -1, -1, 0, -1},        // 4 0 0 1 0
    {-1},                      // 5
    {-1},                      // 6
    {-1},                      // 7
    {4, -1, -1, -1, 1},        // 8 0 0 0 1
    {-1},                      // 9
    {5, -1, 3, -1, 0, -1},     // 10 0 1 0 1 0
    {-1},                      // 11
    {-1},                      // 12
    {-1},                      // 13
    {-1},                      // 14
    {-1},                      // 15
    {-1},                      // 16
    {-1},                      // 17
    {5, -1, 3, -1, -1, 1},     // 18 0 1 0 0 1
    {-1},                      // 19
    {5, -1, -1, 4, -1, 1},     // 20 0 0 1 0 1
    {-1},                      // 21
    {-1},                      // 22
    {-1},                      // 23
    {-1},                      // 24
    {-1},                      // 25
    {-1},                      // 26
    {-1},                      // 27
    {-1},                      // 28
    {-1},                      // 29
    {-1},                      // 30
    {-1},                      // 31
    {-1},                      // 32
    {-1},                      // 33
    {-1},                      // 34
    {-1},                      // 35
    {-1},                      // 36
    {-1},                      // 37
    {-1},                      // 38
    {-1},                      // 39
    {-1},                      // 40
    {-1},                      // 41
    {6, -1, 3, -1, 5, -1, 1},  // 42 0 1 0 1 0 1
};

static void pbvh_split_edges(EdgeQueueContext *eq_ctx,
                             PBVH *pbvh,
                             BMesh *bm,
                             BMEdge **edges1,
                             int totedge,
                             bool ignore_isolated_edges)
{
  bm_logstack_push();
  bm_log_message("  == split edges == ");

  BMEdge **edges = edges1;
  Vector<BMFace *> faces;

  const int node_updateflag = PBVH_UpdateBB | PBVH_UpdateOriginalBB | PBVH_UpdateNormals |
                              PBVH_UpdateOtherVerts | PBVH_UpdateCurvatureDir |
                              PBVH_UpdateTriAreas | PBVH_UpdateDrawBuffers |
                              PBVH_RebuildDrawBuffers | PBVH_UpdateTris | PBVH_UpdateNormals;

  for (int i = 0; i < totedge; i++) {
    BMEdge *e = edges[i];

    check_vert_fan_are_tris(pbvh, e->v1);
    check_vert_fan_are_tris(pbvh, e->v2);
  }

  for (int i = 0; i < totedge; i++) {
    BMEdge *e = edges[i];
    BMLoop *l = e->l;

    /* Clear tags. */
    e->head.hflag &= ~SPLIT_TAG;
    e->v1->head.hflag &= ~SPLIT_TAG;
    e->v2->head.hflag &= ~SPLIT_TAG;

    if (!l) {
      continue;
    }

    /* Clear tags in wider neighborhood and flag valence/boundary for update. */
    do {
      BMLoop *l2 = l->f->l_first;

      do {
        l2->e->head.hflag &= ~SPLIT_TAG;
        l2->v->head.hflag &= ~SPLIT_TAG;

        pbvh_boundary_update_bmesh(pbvh, l2->v);
        MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, l2->v);
        MV_ADD_FLAG(
            mv, SCULPTVERT_NEED_VALENCE | SCULPTVERT_NEED_DISK_SORT | SCULPTVERT_NEED_TRIANGULATE);
      } while ((l2 = l2->next) != l->f->l_first);

      l->f->head.hflag &= ~SPLIT_TAG;
    } while ((l = l->radial_next) != e->l);
  }

  /* Tag edges and faces to split. */
  for (int i = 0; i < totedge; i++) {
    BMEdge *e = edges[i];
    BMLoop *l = e->l;

    e->head.index = 0;
    e->head.hflag |= SPLIT_TAG;

    if (!l) {
      continue;
    }

    do {
      if (!(l->f->head.hflag & SPLIT_TAG)) {
        BMLoop *l2 = l;
        do {
          l2->v->head.hflag &= ~BM_ELEM_TAG;
        } while ((l2 = l2->next) != l);
        l->f->head.hflag |= SPLIT_TAG;

        if (l->f->len == 3) {
          l->f->head.index = l->f->len;
          faces.append(l->f);
        }
      }

    } while ((l = l->radial_next) != e->l);
  }

  int totface = faces.size();
  for (int i = 0; i < totface; i++) {
    BMFace *f = faces[i];
    BMLoop *l = f->l_first;

    if (!ignore_isolated_edges) {
      f->head.hflag |= SPLIT_TAG;
      BM_log_face_removed(pbvh->header.bm, pbvh->bm_log, f);
      BM_idmap_release(pbvh->bm_idmap, (BMElem *)f, true);
    }
    else {
      f->head.hflag &= ~SPLIT_TAG;
    }

    int mask = 0;
    int j = 0;
    int count = 0;

    do {
      if (l->e->head.hflag & SPLIT_TAG) {
        mask |= 1 << j;
        count++;
      }

      j++;
    } while ((l = l->next) != f->l_first);

    if (ignore_isolated_edges) {
      do {
        l->e->head.index = MAX2(l->e->head.index, count);
      } while ((l = l->next) != f->l_first);
    }

    f->head.index = mask;
  }

  bm_log_message("  == split edges (edge split) == ");

  for (int i = 0; i < totedge; i++) {
    BMEdge *e = edges[i];

    if (ignore_isolated_edges && e->head.index < 2) {
      BMLoop *l = e->l;

      do {
        l->f->head.hflag &= ~SPLIT_TAG;
      } while ((l = l->radial_next) != e->l);

      edges[i] = nullptr;
      continue;
    }

    if (ignore_isolated_edges) {
      BMLoop *l = e->l;
      do {
        if (!(l->f->head.hflag & SPLIT_TAG)) {
          l->f->head.hflag |= SPLIT_TAG;
          BM_log_face_removed(pbvh->header.bm, pbvh->bm_log, l->f);
        }
      } while ((l = l->radial_next) != e->l);
    }
  }

  for (int i = 0; i < totedge; i++) {
    BMEdge *e = edges[i];
    BMVert *v1 = e->v1;
    BMVert *v2 = e->v2;
    BMEdge *newe = nullptr;

    if (!e || !(e->head.hflag & SPLIT_TAG)) {
      continue;
    }

    e->head.hflag &= ~SPLIT_TAG;

    MSculptVert *mv1 = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, e->v1);
    MSculptVert *mv2 = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, e->v2);

    if (mv1->stroke_id != pbvh->stroke_id) {
      BKE_pbvh_bmesh_check_origdata(pbvh, e->v1, pbvh->stroke_id);
    }
    if (mv2->stroke_id != pbvh->stroke_id) {
      BKE_pbvh_bmesh_check_origdata(pbvh, e->v2, pbvh->stroke_id);
    }

    if (mv1->stroke_id != mv2->stroke_id) {
      printf("stroke_id error\n");
    }

    validate_edge(pbvh, pbvh->header.bm, e, true, true);

    BM_idmap_check_assign(pbvh->bm_idmap, (BMElem *)e->v1);
    BM_idmap_check_assign(pbvh->bm_idmap, (BMElem *)e->v2);
    BM_log_edge_removed(pbvh->header.bm, pbvh->bm_log, e);

    BM_idmap_release(pbvh->bm_idmap, (BMElem *)e, true);

    // BMVert *newv = BM_log_edge_split_do(pbvh->bm_log, e, e->v1, &newe, 0.5f);
    BMVert *newv = BM_edge_split(pbvh->header.bm, e, e->v1, &newe, 0.5f);
    newe->head.hflag &= ~SPLIT_TAG;

    BM_log_vert_added(bm, pbvh->bm_log, newv);
    BM_log_edge_added(bm, pbvh->bm_log, e);
    BM_log_edge_added(bm, pbvh->bm_log, newe);

    edge_queue_insert_val34_vert(eq_ctx, newv);

#ifdef DYNTOPO_USE_MINMAX_HEAP
    const float elimit = eq_ctx->limit_len_max;

    if (0 && e->l) {
      e->head.hflag &= ~BM_ELEM_TAG;
      newe->head.hflag &= ~BM_ELEM_TAG;

      long_edge_queue_edge_add_recursive(
          eq_ctx, e->l->radial_next, e->l, len_squared_v3v3(e->v1->co, e->v2->co), elimit, 0);
      long_edge_queue_edge_add_recursive(eq_ctx,
                                         newe->l->radial_next,
                                         newe->l,
                                         len_squared_v3v3(newe->v1->co, newe->v2->co),
                                         elimit,
                                         0);
    }

    // edge_queue_insert_unified(eq_ctx, e);
    // edge_queue_insert_unified(eq_ctx, newe);

    // BLI_table_gset_add(eq_ctx->used_verts, newv);
#endif

    PBVH_CHECK_NAN(newv->co);

    validate_edge(pbvh, pbvh->header.bm, e, true, true);
    validate_edge(pbvh, pbvh->header.bm, newe, true, true);
    validate_vert(pbvh, pbvh->header.bm, newv, true, true);

    MSculptVert *mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, newv);

    newv->head.hflag |= SPLIT_TAG;

    pbvh_boundary_update_bmesh(pbvh, newv);
    MV_ADD_FLAG(mv,
                SCULPTVERT_NEED_VALENCE | SCULPTVERT_NEED_DISK_SORT | SCULPTVERT_NEED_TRIANGULATE);
    mv->stroke_id = pbvh->stroke_id;

    mv = BKE_PBVH_SCULPTVERT(pbvh->cd_sculpt_vert, e->v1 != newv ? e->v1 : e->v2);
    pbvh_boundary_update_bmesh(pbvh, e->v1 != newv ? e->v1 : e->v2);
    MV_ADD_FLAG(mv,
                SCULPTVERT_NEED_DISK_SORT | SCULPTVERT_NEED_VALENCE | SCULPTVERT_NEED_TRIANGULATE);

    BM_ELEM_CD_SET_INT(newv, pbvh->cd_vert_node_offset, DYNTOPO_NODE_NONE);

    int ni = BM_ELEM_CD_GET_INT(v1, pbvh->cd_vert_node_offset);

    if (ni == DYNTOPO_NODE_NONE) {
      ni = BM_ELEM_CD_GET_INT(v2, pbvh->cd_vert_node_offset);
    }

    if (ni >= pbvh->totnode || !(pbvh->nodes[ni].flag & PBVH_Leaf)) {
      printf("%s: error\n", __func__);
    }
    /* this should rarely happen */
    // if (ni < 0 || ni >= pbvh->totnode || !(pbvh->nodes[ni].flag & PBVH_Leaf)) {
    if (ni == DYNTOPO_NODE_NONE) {
      ni = DYNTOPO_NODE_NONE;

      for (int j = 0; j < 2; j++) {
        BMVert *v = nullptr;

        switch (j) {
          case 0:
            v = v1;
            break;
          case 1:
            v = v2;
            break;
        }

        if (!v->e) {
          continue;
        }

        BMEdge *e2 = v->e;
        do {
          if (!e2->l) {
            break;
          }

          BMLoop *l = e2->l;
          do {
            int ni2 = BM_ELEM_CD_GET_INT(l->f, pbvh->cd_face_node_offset);

            if (ni2 >= 0 && ni2 < pbvh->totnode && (pbvh->nodes[ni2].flag & PBVH_Leaf)) {
              ni = ni2;
              goto outerbreak;
            }
          } while ((l = l->radial_next) != e2->l);
        } while ((e2 = BM_DISK_EDGE_NEXT(e2, v)) != v->e);
      }
    outerbreak:;
    }

    if (ni != DYNTOPO_NODE_NONE) {
      PBVHNode *node = pbvh->nodes + ni;

      if (!(node->flag & PBVH_Leaf)) {
        printf("pbvh error in pbvh_split_edges!\n");

        BM_ELEM_CD_SET_INT(newv, pbvh->cd_vert_node_offset, DYNTOPO_NODE_NONE);

        continue;
      }

      node->flag |= (PBVHNodeFlags)node_updateflag;

      BLI_table_gset_add(node->bm_unique_verts, newv);

      BM_ELEM_CD_SET_INT(newv, pbvh->cd_vert_node_offset, ni);
      // BM_ELEM_CD_SET_INT(newv, pbvh->cd_vert_node_offset, -1);
    }
    else {
      BM_ELEM_CD_SET_INT(newv, pbvh->cd_vert_node_offset, DYNTOPO_NODE_NONE);
      printf("%s: error!\n", __func__);
    }
  }

  bm_log_message("  == split edges (triangulate) == ");

  /* Subdivide from template. */

  Vector<BMVert *, 32> vs;
  Vector<BMFace *> newfaces;

  for (int i = 0; i < totface; i++) {
    BMFace *f = faces[i];
    int mask = 0;

    if (!(f->head.hflag & SPLIT_TAG)) {
      continue;
    }

    int ni = BM_ELEM_CD_GET_INT(f, pbvh->cd_face_node_offset);

    if (ni < 0 || ni >= pbvh->totnode || !(pbvh->nodes[ni].flag & PBVH_Leaf)) {
      printf("%s: error!\n", __func__);
      ni = DYNTOPO_NODE_NONE;
    }

    BMLoop *l = f->l_first;
    int j = 0;
    int totmask = 0;

    do {
      if (l->v->head.hflag & SPLIT_TAG) {
        mask |= 1 << j;
        totmask++;
      }
      j++;
    } while ((l = l->next) != f->l_first);

    int flen = j;

    if (mask < 0 || mask >= (int)ARRAY_SIZE(splitmap)) {
      printf("splitmap error! flen: %d totmask: %d mask: %d\n", flen, totmask, mask);
      continue;
    }

    const int *pat = splitmap[mask];
    int n = pat[0];

    if (n < 0) {
      continue;
    }

    if (n != f->len || n != flen) {
      printf("%s: error! %d %d\n", __func__, n, flen);
      continue;
    }

    BMFace *f2 = f;

    vs.resize(n);

    l = f->l_first;
    j = 0;
    do {
      vs[j++] = l->v;
    } while ((l = l->next) != f->l_first);

    if (j != n) {
      printf("error! %s\n", __func__);
      continue;
    }

    newfaces.resize(newfaces.size() + n);

    int count = 0;

    for (j = 0; j < n; j++) {
      if (pat[j + 1] < 0) {
        continue;
      }

      BMVert *v1 = vs[j], *v2 = vs[pat[j + 1]];
      BMLoop *l1 = nullptr, *l2 = nullptr;
      BMLoop *rl = nullptr;

      BMLoop *l3 = f2->l_first;
      do {
        if (l3->v == v1) {
          l1 = l3;
        }
        else if (l3->v == v2) {
          l2 = l3;
        }
      } while ((l3 = l3->next) != f2->l_first);

      if (l1 == l2 || !l1 || !l2) {
        printf("errorl!\n");
        continue;
      }

      validate_face(pbvh, bm, f2, false, true);

      bool log_edge = true;
      BMFace *newf = nullptr;
      BMEdge *exist_e;

      if ((exist_e = BM_edge_exists(v1, v2))) {
        log_edge = false;

        BMLoop *l1 = exist_e->l;

        if (l1 && l1->f == f2) {
          l1 = l1->radial_next;
        }

        if (l1 && l1->f != f2) {
          // newf = l1->f;
        }
      }
      else {
        newf = BM_face_split(bm, f2, l1, l2, &rl, nullptr, false);
      }

      if (newf) {
        // rl->e->head.hflag &= ~BM_ELEM_TAG;
        // edge_queue_insert_unified(eq_ctx, rl->e);
        /*
        long_edge_queue_edge_add_recursive_3(eq_ctx,
                                             rl->radial_next,
                                             rl,
                                             len_squared_v3v3(rl->e->v1->co, rl->e->v2->co),
                                             pbvh->bm_max_edge_len,
                                             0);
                                             */

        check_face_is_manifold(pbvh, bm, newf);
        check_face_is_manifold(pbvh, bm, f2);
        check_face_is_manifold(pbvh, bm, f);

        validate_face(pbvh, bm, f2, false, true);
        validate_face(pbvh, bm, newf, false, true);

        if (log_edge) {
          BM_log_edge_added(bm, pbvh->bm_log, rl->e);
        }

        bool ok = ni != DYNTOPO_NODE_NONE;
        ok = ok && BM_ELEM_CD_GET_INT(v1, pbvh->cd_vert_node_offset) != DYNTOPO_NODE_NONE;
        ok = ok && BM_ELEM_CD_GET_INT(v2, pbvh->cd_vert_node_offset) != DYNTOPO_NODE_NONE;

        if (ok) {
          PBVHNode *node = pbvh->nodes + ni;

          node->flag |= (PBVHNodeFlags)node_updateflag;

          BLI_table_gset_add(node->bm_faces, newf);
          BM_ELEM_CD_SET_INT(newf, pbvh->cd_face_node_offset, ni);
        }
        else {
          BM_ELEM_CD_SET_INT(newf, pbvh->cd_face_node_offset, DYNTOPO_NODE_NONE);
        }

        if (count < n) {
          newfaces[count++] = newf;
        }
        else {
          printf("%s: error 4!\n", __func__);
        }
        f2 = newf;
      }
      else {
        printf("%s: split error 2!\n", __func__);
        continue;
      }
    }

    for (j = 0; j < count; j++) {
      if (BM_ELEM_CD_GET_INT(newfaces[j], pbvh->cd_face_node_offset) == DYNTOPO_NODE_NONE) {
        BKE_pbvh_bmesh_add_face(pbvh, newfaces[j], false, true);
      }

      if (newfaces[j] != f) {
        BM_log_face_added(bm, pbvh->bm_log, newfaces[j]);
      }
#if 1
      if (newfaces[j]->len != 3) {
        printf("%s: tesselation error!\n", __func__);
      }
#endif
    }

    if (f->len != 3) {
      printf("%s: tesselation error!\n", __func__);
    }

    if (BM_ELEM_CD_GET_INT(f, pbvh->cd_face_node_offset) == DYNTOPO_NODE_NONE) {
      BKE_pbvh_bmesh_add_face(pbvh, f, false, true);
    }

    BM_log_face_added(bm, pbvh->bm_log, f);
  }

  bm_logstack_pop();
}

typedef struct DynTopoState {
  PBVH *pbvh;
  bool is_fake_pbvh;
} DynTopoState;

/* existing_pbvh may be nullptr, if so a fake one will be created.
Note that all the sculpt customdata layers will be created
if they don't exist, so cd_vert/face_node_offset, cd_mask_offset,
cd_sculpt_vert, etc*/
DynTopoState *BKE_dyntopo_init(BMesh *bm, PBVH *existing_pbvh)
{
  PBVH *pbvh;

  if (!existing_pbvh) {
    pbvh = MEM_cnew<PBVH>("pbvh");

    pbvh->nodes = static_cast<PBVHNode *>(MEM_callocN(sizeof(PBVHNode), "PBVHNode"));
    pbvh->header.type = PBVH_BMESH;
    pbvh->totnode = 1;

    PBVHNode *node = pbvh->nodes;

    node->flag = PBVH_Leaf | PBVH_UpdateTris | PBVH_UpdateTriAreas;
    node->bm_faces = BLI_table_gset_new_ex("node->bm_faces", bm->totface);
    node->bm_unique_verts = BLI_table_gset_new_ex("node->bm_unique_verts", bm->totvert);
  }
  else {
    pbvh = existing_pbvh;
  }

  if (!pbvh->bm_idmap) {
    pbvh->bm_idmap = BM_idmap_new(bm, BM_VERT | BM_EDGE | BM_FACE);
  }

  const bool isfake = pbvh != existing_pbvh;

  BMCustomLayerReq vlayers[] = {{CD_PAINT_MASK, nullptr, 0},
                                {CD_DYNTOPO_VERT, nullptr, CD_FLAG_TEMPORARY | CD_FLAG_NOCOPY},
                                {CD_PROP_INT32,
                                 SCULPT_ATTRIBUTE_NAME(dyntopo_node_id_vertex),
                                 CD_FLAG_TEMPORARY | CD_FLAG_NOCOPY}};

  BMCustomLayerReq flayers[] = {
      {CD_PROP_FLOAT2, SCULPT_ATTRIBUTE_NAME(face_areas), CD_FLAG_TEMPORARY | CD_FLAG_NOCOPY},
      {CD_PROP_INT32, ".sculpt_face_set", 0},
      {CD_PROP_INT32,
       SCULPT_ATTRIBUTE_NAME(dyntopo_node_id_face),
       CD_FLAG_TEMPORARY | CD_FLAG_NOCOPY}};

  BM_data_layers_ensure(bm, &bm->vdata, vlayers, 3);
  BM_data_layers_ensure(bm, &bm->pdata, flayers, 3);

  pbvh->header.bm = bm;

  pbvh->cd_vert_node_offset = CustomData_get_named_offset(
      &bm->vdata, CD_PROP_INT32, SCULPT_ATTRIBUTE_NAME(dyntopo_node_id_vertex));

  pbvh->cd_face_node_offset = CustomData_get_named_offset(
      &bm->pdata, CD_PROP_INT32, SCULPT_ATTRIBUTE_NAME(dyntopo_node_id_face));

  pbvh->cd_face_area = CustomData_get_named_offset(
      &bm->pdata, CD_PROP_FLOAT2, SCULPT_ATTRIBUTE_NAME(face_areas));

  pbvh->cd_sculpt_vert = CustomData_get_offset(&bm->vdata, CD_DYNTOPO_VERT);
  pbvh->cd_vert_mask_offset = CustomData_get_offset(&bm->vdata, CD_PAINT_MASK);
  pbvh->cd_faceset_offset = CustomData_get_offset_named(
      &bm->pdata, CD_PROP_INT32, ".sculpt_face_set");
  pbvh->cd_vcol_offset = -1;

  if (isfake) {
    pbvh->bm_log = BM_log_create(bm, pbvh->bm_idmap);
  }

  BMVert *v;
  BMFace *f;
  BMIter iter;

  if (isfake) {
    BM_ITER_MESH (v, &iter, bm, BM_VERTS_OF_MESH) {
      BM_ELEM_CD_SET_INT(v, pbvh->cd_vert_node_offset, 0);
      BLI_table_gset_add(pbvh->nodes->bm_unique_verts, v);
    }

    BM_ITER_MESH (f, &iter, bm, BM_FACES_OF_MESH) {
      BM_ELEM_CD_SET_INT(f, pbvh->cd_face_node_offset, 0);
      BLI_table_gset_add(pbvh->nodes->bm_faces, f);
    }

    BKE_pbvh_bmesh_check_tris(pbvh, pbvh->nodes);
  }

  DynTopoState *ds = MEM_cnew<DynTopoState>("DynTopoState");

  ds->pbvh = pbvh;
  ds->is_fake_pbvh = isfake;

  return ds;
}

void BKE_dyntopo_default_params(DynRemeshParams *params, float edge_size)
{
  memset(params, 0, sizeof(*params));
  params->detail_range = 0.45f;
  params->edge_size = edge_size;
}

void BKE_dyntopo_free(DynTopoState *ds)
{
  if (ds->is_fake_pbvh) {
    BM_log_free(ds->pbvh->bm_log, false);
    BM_idmap_destroy(ds->pbvh->bm_idmap);

    PBVHNode *node = ds->pbvh->nodes;

    if (node->tribuf || node->tri_buffers) {
      BKE_pbvh_bmesh_free_tris(ds->pbvh, node);
    }

    BLI_table_gset_free(node->bm_faces, nullptr);
    BLI_table_gset_free(node->bm_unique_verts, nullptr);

    MEM_freeN(ds->pbvh->nodes);
    MEM_freeN(ds->pbvh);
  }

  MEM_freeN(ds);
}

/*
bool BKE_pbvh_bmesh_update_topology(PBVH *pbvh,
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
                                    int custom_max_steps,
                                    bool disable_surface_relax)
*/
void BKE_dyntopo_remesh(DynTopoState *ds,
                        DynRemeshParams *params,
                        int steps,
                        PBVHTopologyUpdateMode mode)
{
  float cent[3] = {0.0f, 0.0f, 0.0f};
  int totcent = 0;
  float view[3] = {0.0f, 0.0f, 1.0f};

  BMIter iter;
  BMVert *v;

  BM_ITER_MESH (v, &iter, ds->pbvh->header.bm, BM_VERTS_OF_MESH) {
    MSculptVert *mv = BKE_PBVH_SCULPTVERT(ds->pbvh->cd_sculpt_vert, v);

    pbvh_boundary_update_bmesh(ds->pbvh, v);
    mv->flag |= SCULPTVERT_NEED_TRIANGULATE;
    mv->valence = BM_vert_edge_count(v);

    pbvh_check_vert_boundary(ds->pbvh, v);

    add_v3_v3(cent, v->co);
    totcent++;
  }

  if (totcent) {
    mul_v3_fl(cent, 1.0f / (float)totcent);
  }

  ds->pbvh->bm_max_edge_len = params->edge_size;
  ds->pbvh->bm_min_edge_len = params->edge_size * params->detail_range;
  ds->pbvh->bm_detail_range = params->detail_range;

  /* subdivide once */
  if (mode & PBVH_Subdivide) {
    BKE_pbvh_bmesh_update_topology(ds->pbvh,
                                   PBVH_Subdivide,
                                   cent,
                                   view,
                                   1e17,
                                   false,
                                   false,
                                   0,
                                   false,
                                   mask_cb_nop,
                                   nullptr,
                                   ds->pbvh->header.bm->totedge,
                                   false,
                                   false);
  }

  for (int i = 0; i < steps; i++) {
    for (int j = 0; j < ds->pbvh->totnode; j++) {
      PBVHNode *node = ds->pbvh->nodes + j;

      if (node->flag & PBVH_Leaf) {
        node->flag |= PBVH_UpdateTopology;
      }
    }

    BKE_pbvh_bmesh_update_topology(ds->pbvh,
                                   mode,
                                   cent,
                                   view,
                                   1e17,
                                   false,
                                   false,
                                   0,
                                   false,
                                   mask_cb_nop,
                                   nullptr,
                                   ds->pbvh->header.bm->totedge * 5,
                                   true,
                                   false);

    BKE_pbvh_update_normals(ds->pbvh, nullptr);

    BM_ITER_MESH (v, &iter, ds->pbvh->header.bm, BM_VERTS_OF_MESH) {
      MSculptVert *mv = BKE_PBVH_SCULPTVERT(ds->pbvh->cd_sculpt_vert, v);

      pbvh_check_vert_boundary(ds->pbvh, v);
      int boundflag = BM_ELEM_CD_GET_INT(v, ds->pbvh->cd_boundary_flag);

      float avg[3] = {0.0f, 0.0f, 0.0f};
      float totw = 0.0f;

      bool bound1 = boundflag & SCULPTVERT_ALL_BOUNDARY;
      if (bound1) {
        continue;
      }

      if (boundflag & SCULPTVERT_ALL_CORNER) {
        continue;
      }

      if (!v->e) {
        continue;
      }

      BMEdge *e = v->e;
      do {
        BMVert *v2 = BM_edge_other_vert(e, v);
        MSculptVert *mv2 = BKE_PBVH_SCULPTVERT(ds->pbvh->cd_sculpt_vert, v2);

        pbvh_check_vert_boundary(ds->pbvh, v2);
        int boundflag2 = BM_ELEM_CD_GET_INT(v2, ds->pbvh->cd_boundary_flag);

        bool bound2 = boundflag2 & SCULPTVERT_ALL_BOUNDARY;

        if (bound1 && !bound2) {
          continue;
        }

        float tmp[3];
        float w = 1.0f;

        sub_v3_v3v3(tmp, v2->co, v->co);
        madd_v3_v3fl(tmp, v->no, -dot_v3v3(v->no, tmp) * 0.75);
        add_v3_v3(tmp, v->co);
        madd_v3_v3fl(avg, tmp, w);

        totw += w;
      } while ((e = BM_DISK_EDGE_NEXT(e, v)) != v->e);

      if (totw == 0.0f) {
        continue;
      }

      mul_v3_fl(avg, 1.0f / totw);
      interp_v3_v3v3(v->co, v->co, avg, 0.5f);
    }
  }
}
}  // namespace blender::dyntopo

void BKE_pbvh_bmesh_remove_face(PBVH *pbvh, BMFace *f, bool log_face)
{
  blender::dyntopo::pbvh_bmesh_face_remove(pbvh, f, log_face, true, true);
}

void BKE_pbvh_bmesh_remove_edge(PBVH *pbvh, BMEdge *e, bool log_edge)
{
  if (log_edge) {
    bm_logstack_push();
    BM_log_edge_removed(pbvh->header.bm, pbvh->bm_log, e);
    bm_logstack_pop();
  }
}

void BKE_pbvh_bmesh_remove_vertex(PBVH *pbvh, BMVert *v, bool log_vert)
{
  blender::dyntopo::pbvh_bmesh_vert_remove(pbvh, v);

  if (log_vert) {
    BM_log_vert_removed(pbvh->header.bm, pbvh->bm_log, v);
  }
}

void BKE_pbvh_bmesh_add_face(PBVH *pbvh, struct BMFace *f, bool log_face, bool force_tree_walk)
{
  bm_logstack_push();

  int ni = DYNTOPO_NODE_NONE;

  if (force_tree_walk) {
    bke_pbvh_insert_face(pbvh, f);

    if (log_face) {
      BM_log_face_added(pbvh->header.bm, pbvh->bm_log, f);
    }

    bm_logstack_pop();
    return;
  }

  // look for node in surrounding geometry
  BMLoop *l = f->l_first;
  do {
    ni = BM_ELEM_CD_GET_INT(l->radial_next->f, pbvh->cd_face_node_offset);

    if (ni >= 0 && (!(pbvh->nodes[ni].flag & PBVH_Leaf) || ni >= pbvh->totnode)) {
      printf("%s: error: ni: %d totnode: %d\n", __func__, ni, pbvh->totnode);
      l = l->next;
      continue;
    }

    if (ni >= 0 && (pbvh->nodes[ni].flag & PBVH_Leaf)) {
      break;
    }

    l = l->next;
  } while (l != f->l_first);

  if (ni < 0) {
    bke_pbvh_insert_face(pbvh, f);
  }
  else {
    BM_ELEM_CD_SET_INT(f, pbvh->cd_face_node_offset, ni);
    bke_pbvh_insert_face_finalize(pbvh, f, ni);
  }

  if (log_face) {
    BM_log_face_added(pbvh->header.bm, pbvh->bm_log, f);
  }

  bm_logstack_pop();
}
