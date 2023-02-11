/* SPDX-License-Identifier: GPL-2.0-or-later
 * Copyright 2008 Blender Foundation. All rights reserved. */

/** \file
 * \ingroup bli
 *
 * Simple, fast memory allocator for allocating many elements of the same size.
 *
 * Supports:
 *
 * - Freeing chunks.
 * - Iterating over allocated chunks
 *   (optionally when using the #BLI_MEMPOOL_ALLOW_ITER flag).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "atomic_ops.h"

#include "BLI_string.h"
#include "BLI_utildefines.h"
#
#include "BLI_asan.h"
#include "BLI_mempool.h"         /* own include */
#include "BLI_mempool_private.h" /* own include */

#ifndef WITH_DNA_GHASH
#  include "BLI_threads.h"
#endif

#include "MEM_guardedalloc.h"

#include "BLI_strict_flags.h" /* keep last */

#ifdef WITH_MEM_VALGRIND
#  include "valgrind/memcheck.h"
#endif

#if (defined(__SANITIZE_ADDRESS__) || __has_feature(address_sanitizer))
#  define POISON_REDZONE_SIZE 32
#  define HAVE_MEMPOOL_ASAN
#else
#  define POISON_REDZONE_SIZE 0
#endif

/* NOTE: copied from BLO_blend_defs.h, don't use here because we're in BLI. */
#ifdef __BIG_ENDIAN__
/* Big Endian */
#  define MAKE_ID(a, b, c, d) ((int)(a) << 24 | (int)(b) << 16 | (c) << 8 | (d))
#  define MAKE_ID_8(a, b, c, d, e, f, g, h) \
    ((int64_t)(a) << 56 | (int64_t)(b) << 48 | (int64_t)(c) << 40 | (int64_t)(d) << 32 | \
     (int64_t)(e) << 24 | (int64_t)(f) << 16 | (int64_t)(g) << 8 | (h))
#else
/* Little Endian */
#  define MAKE_ID(a, b, c, d) ((int)(d) << 24 | (int)(c) << 16 | (b) << 8 | (a))
#  define MAKE_ID_8(a, b, c, d, e, f, g, h) \
    ((int64_t)(h) << 56 | (int64_t)(g) << 48 | (int64_t)(f) << 40 | (int64_t)(e) << 32 | \
     (int64_t)(d) << 24 | (int64_t)(c) << 16 | (int64_t)(b) << 8 | (a))
#endif

/**
 * Important that this value is an is _not_  aligned with `sizeof(void *)`.
 * So having a pointer to 2/4/8... aligned memory is enough to ensure
 * the `freeword` will never be used.
 * To be safe, use a word that's the same in both directions.
 */
#define FREEWORD \
  ((sizeof(void *) > sizeof(int32_t)) ? MAKE_ID_8('e', 'e', 'r', 'f', 'f', 'r', 'e', 'e') : \
                                        MAKE_ID('e', 'f', 'f', 'e'))

/**
 * The 'used' word just needs to be set to something besides FREEWORD.
 */
#define USEDWORD MAKE_ID('u', 's', 'e', 'd')

/* Currently totalloc isn't used. */
// #define USE_TOTALLOC

/* optimize pool size */
#define USE_CHUNK_POW2

#ifndef NDEBUG
static bool mempool_debug_memset = false;
#endif

/**
 * A free element from #BLI_mempool_chunk. Data is cast to this type and stored in
 * #BLI_mempool.free as a single linked list, each item #BLI_mempool.esize large.
 *
 * Each element represents a block which BLI_mempool_alloc may return.
 */
typedef struct BLI_freenode {
  struct BLI_freenode *next;
  /** Used to identify this as a freed node. */
  intptr_t freeword;
} BLI_freenode;

/**
 * A chunk of memory in the mempool stored in
 * #BLI_mempool.chunks as a double linked list.
 */
typedef struct BLI_mempool_chunk {
  struct BLI_mempool_chunk *next;
} BLI_mempool_chunk;

/**
 * The mempool, stores and tracks memory \a chunks and elements within those chunks \a free.
 */
struct BLI_mempool {
#if !defined(WITH_DNA_GHASH) && (defined(__SANITIZE_ADDRESS__) || __has_feature(address_sanitizer))
  ThreadMutex mutex;
#endif
  /** Single linked list of allocated chunks. */
  BLI_mempool_chunk *chunks;
  /** Keep a pointer to the last, so we can append new chunks there
   * this is needed for iteration so we can loop over chunks in the order added. */
  BLI_mempool_chunk *chunk_tail;

  /* only used if BLI_MEMPOOL_RANDOM_ACCESS is true*/
  BLI_mempool_chunk **chunktable;

  /* only used if BLI_MEMPOOL_RANDOM_ACCESS is true*/
  int totchunk;

  /** Element size in bytes. */
  uint esize;
  /** Chunk size in bytes. */
  uint csize;
  /** Number of elements per chunk. */
  uint pchunk;
  uint flag;
  /* keeps aligned to 16 bits */

#ifdef HAVE_MEMPOOL_ASAN
  char poisoned[256];
#endif

  /** Free element list. Interleaved into chunk data. */
  BLI_freenode *free;
  /** Use to know how many chunks to keep for #BLI_mempool_clear. */
  uint maxchunks;
  /** Number of elements currently in use. */
  int totused;
#ifdef USE_TOTALLOC
  /** Number of elements allocated in total. */
  uint totalloc;
#endif

  char *memtag;
  char *memtag_chunk;
};

#define MEMPOOL_ELEM_SIZE_MIN (sizeof(void *) * 2)

#define CHUNK_DATA(chunk) (CHECK_TYPE_INLINE(chunk, BLI_mempool_chunk *), (void *)((chunk) + 1))

#define NODE_STEP_NEXT(node) ((void *)((char *)(node) + esize))
#define NODE_STEP_PREV(node) ((void *)((char *)(node)-esize))

/** Extra bytes implicitly used for every chunk alloc. */
#define CHUNK_OVERHEAD (uint)(MEM_SIZE_OVERHEAD + sizeof(BLI_mempool_chunk))

static void mempool_poison(BLI_mempool *pool)
{
#if (defined(__SANITIZE_ADDRESS__) || __has_feature(address_sanitizer))
#  ifndef WITH_DNA_GHASH
  BLI_mutex_unlock(&pool->mutex);
#  endif

#  ifndef WITH_DNA_GHASH
  char *ptr = (char *)pool;
  size_t offset = offsetof(BLI_mempool, mutex) + sizeof(ThreadMutex);
  BLI_asan_poison(ptr + offset, sizeof(*pool) - offset);
#  else
  BLI_asan_poison(pool, sizeof(*pool));
#  endif
#endif
}

static void mempool_poison_no_unlock(BLI_mempool *pool)
{
#if (defined(__SANITIZE_ADDRESS__) || __has_feature(address_sanitizer))

#  ifndef WITH_DNA_GHASH
  char *ptr = (char *)pool;
  size_t offset = offsetof(BLI_mempool, mutex) + sizeof(ThreadMutex);
  BLI_asan_poison(ptr + offset, sizeof(*pool) - offset);
#  else
  BLI_asan_poison(pool, sizeof(*pool));
#  endif
#endif
}

static void mempool_unpoison(BLI_mempool *pool)
{
#if (defined(__SANITIZE_ADDRESS__) || __has_feature(address_sanitizer))
#  ifndef WITH_DNA_GHASH
  BLI_mutex_lock(&pool->mutex);
#  endif

#  ifndef WITH_DNA_GHASH
  char *ptr = (char *)pool;
  size_t offset = offsetof(BLI_mempool, mutex) + sizeof(ThreadMutex);
  BLI_asan_unpoison(ptr + offset, sizeof(*pool) - offset);
#  else
  BLI_asan_unpoison(pool, sizeof(*pool));
#  endif
#endif
}

#ifdef USE_CHUNK_POW2
static uint power_of_2_max_u(uint x)
{
  x -= 1;
  x = x | (x >> 1);
  x = x | (x >> 2);
  x = x | (x >> 4);
  x = x | (x >> 8);
  x = x | (x >> 16);
  return x + 1;
}
#endif

static void mempool_update_chunktable(BLI_mempool *pool)
{
  if (!(pool->flag & BLI_MEMPOOL_RANDOM_ACCESS)) {
    return;
  }

  pool->totchunk = 0;
  BLI_mempool_chunk *chunk = pool->chunks;

  while (chunk) {
    pool->totchunk++;
    chunk = chunk->next;
  }

  MEM_SAFE_FREE(pool->chunktable);
  pool->chunktable = (BLI_mempool_chunk **)MEM_mallocN(
      sizeof(pool->chunktable) * (size_t)pool->totchunk, "mempool chunktable");

  int i = 0;
  chunk = pool->chunks;

  while (chunk) {
    pool->chunktable[i++] = chunk;
    chunk = chunk->next;
  }
}

BLI_INLINE BLI_mempool_chunk *mempool_chunk_find(BLI_mempool_chunk *head, uint index)
{
  while (index-- && head) {
    head = head->next;
  }
  return head;
}

/**
 * \return the number of chunks to allocate based on how many elements are needed.
 *
 * \note for small pools 1 is a good default, the elements need to be initialized,
 * adding overhead on creation which is redundant if they aren't used.
 */
BLI_INLINE uint mempool_maxchunks(const uint elem_num, const uint pchunk)
{
  return (elem_num <= pchunk) ? 1 : ((elem_num / pchunk) + 1);
}

static BLI_mempool_chunk *mempool_chunk_alloc(BLI_mempool *pool)
{
  return BLI_asan_safe_malloc(sizeof(BLI_mempool_chunk) + (size_t)pool->csize, pool->memtag);
}

/**
 * Initialize a chunk and add into \a pool->chunks
 *
 * \param pool: The pool to add the chunk into.
 * \param mpchunk: The new uninitialized chunk (can be malloc'd)
 * \param last_tail: The last element of the previous chunk
 * (used when building free chunks initially)
 * \return The last chunk,
 */
static BLI_freenode *mempool_chunk_add(BLI_mempool *pool,
                                       BLI_mempool_chunk *mpchunk,
                                       BLI_freenode *last_tail)
{
  const uint esize = pool->esize;
  BLI_freenode *curnode = CHUNK_DATA(mpchunk);
  uint j;

  if (pool->flag & BLI_MEMPOOL_RANDOM_ACCESS) {
    if (!pool->chunktable ||
        MEM_allocN_len(pool->chunktable) / sizeof(void *) <= (size_t)pool->totchunk) {
      void *old = pool->chunktable;

      int size = (int)pool->totchunk + 2;
      size += size >> 1;

      pool->chunktable = MEM_mallocN(sizeof(void *) * (size_t)size, "mempool chunktable");

      if (old) {
        memcpy(pool->chunktable, old, sizeof(void *) * (size_t)pool->totchunk);
      }

      MEM_SAFE_FREE(old);
    }

    pool->chunktable[pool->totchunk++] = mpchunk;
  }

  /* append */
  if (pool->chunk_tail) {
    pool->chunk_tail->next = mpchunk;
  }
  else {
    BLI_assert(pool->chunks == NULL);
    pool->chunks = mpchunk;
  }

  mpchunk->next = NULL;
  pool->chunk_tail = mpchunk;

  if (UNLIKELY(pool->free == NULL)) {
    pool->free = curnode;
  }

  /* loop through the allocated data, building the pointer structures */
  j = pool->pchunk;
  if (pool->flag & BLI_MEMPOOL_ALLOW_ITER) {
    while (j--) {
      BLI_freenode *next;

      BLI_asan_unpoison(curnode, pool->esize - POISON_REDZONE_SIZE);

      curnode->next = next = NODE_STEP_NEXT(curnode);
      curnode->freeword = FREEWORD;

      BLI_asan_poison(curnode, pool->esize);

      curnode = next;
    }
  }
  else {
    while (j--) {
      BLI_freenode *next;

      BLI_asan_unpoison(curnode, pool->esize - POISON_REDZONE_SIZE);
      curnode->next = next = NODE_STEP_NEXT(curnode);
      BLI_asan_poison(curnode, pool->esize);

      curnode = next;
    }
  }

  /* terminate the list (rewind one)
   * will be overwritten if 'curnode' gets passed in again as 'last_tail' */

  BLI_asan_unpoison(curnode, pool->esize - POISON_REDZONE_SIZE);
  BLI_asan_poison(curnode, pool->esize);

  curnode = NODE_STEP_PREV(curnode);

  BLI_asan_unpoison(curnode, pool->esize - POISON_REDZONE_SIZE);
  curnode->next = NULL;
  BLI_asan_poison(curnode, pool->esize);

#ifdef USE_TOTALLOC
  pool->totalloc += pool->pchunk;
#endif

  /* final pointer in the previously allocated chunk is wrong */
  if (last_tail) {
    BLI_asan_unpoison(last_tail, pool->esize - POISON_REDZONE_SIZE);
    last_tail->next = CHUNK_DATA(mpchunk);
    BLI_asan_poison(last_tail, pool->esize);
  }

  return curnode;
}

/**
 * Preallocates a mempool suitable for threading.  elem_num elements are preallocated
 * in chunks of size pchunk, and returned in r_chunks.  The idea is to pass these
 * to tasks.
 */
BLI_mempool *BLI_mempool_create_for_tasks(const unsigned int esize,
                                          int elem_num,
                                          const int pchunk,
                                          void ***r_chunks,
                                          int *r_totchunk,
                                          int *r_esize,
                                          int flag)
{
  BLI_mempool *pool = BLI_mempool_create(esize, 0, (uint)pchunk, (uint)flag);

  mempool_unpoison(pool);

  // override pchunk, may not be a power of 2
  pool->pchunk = (uint)pchunk;
  pool->csize = (uint)pchunk * pool->esize;

  if (elem_num % pchunk == 0) {
    pool->maxchunks = (uint)elem_num / (uint)pchunk;
  }
  else {
    pool->maxchunks = (uint)elem_num / (uint)pchunk + 1;
  }

  if (elem_num) {
    BLI_freenode *last_tail = NULL;

    /* Allocate the actual chunks. */
    for (uint i = 0; i < pool->maxchunks; i++) {
      BLI_mempool_chunk *mpchunk = mempool_chunk_alloc(pool);
      last_tail = mempool_chunk_add(pool, mpchunk, last_tail);
    }
  }

  void **chunks = MEM_callocN(sizeof(void *) * pool->maxchunks,
                              "BLI_mempool_create_for_tasks r_chunks");

  unsigned int totalloc = 0;
  *r_totchunk = 0;

  BLI_mempool_chunk *chunk = pool->chunks, *lastchunk = NULL;

  while (chunk) {
    lastchunk = chunk;
    totalloc += pool->pchunk;
    chunk = chunk->next;
  }

  pool->totused = (int)totalloc;
  pool->free = NULL;

  int i = (int)pool->pchunk - 1;

  while (lastchunk && totalloc > (uint)elem_num) {
    if (i < 0) {
      BLI_mempool_chunk *lastchunk2 = NULL;

      for (chunk = pool->chunks; chunk; chunk = chunk->next) {
        if (chunk == lastchunk) {
          lastchunk = lastchunk2;
        }
      }

      if (!lastchunk) {
        break;
      }

      i = (int)pool->pchunk - 1;
    }

    char *elem = CHUNK_DATA(lastchunk);
    elem += pool->esize * (unsigned int)i;

    BLI_mempool_free(pool, elem);

    totalloc--;
    i--;
  }

  int ci = 0;

  chunk = pool->chunks;
  while (chunk && chunk != lastchunk) {
    chunks[ci++] = CHUNK_DATA(chunk);
    chunk = chunk->next;
  }

  if (lastchunk && i >= 0) {
    chunks[ci++] = CHUNK_DATA(lastchunk);
  }

  *r_totchunk = ci;
  *r_chunks = (void **)chunks;
  *r_esize = (int)pool->esize;

  return pool;
}

static void mempool_chunk_free(BLI_mempool_chunk *mpchunk, BLI_mempool *pool)
{
  BLI_asan_unpoison(mpchunk, sizeof(BLI_mempool_chunk) + pool->esize * pool->csize);
  BLI_asan_safe_free(mpchunk);
}

ATTR_NO_OPT static void mempool_chunk_free_all(BLI_mempool_chunk *mpchunk, BLI_mempool *pool)
{
  BLI_mempool_chunk *mpchunk_next;

  for (; mpchunk; mpchunk = mpchunk_next) {
    mpchunk_next = mpchunk->next;
    mempool_chunk_free(mpchunk, pool);
  }
}

#ifdef BLI_mempool_create
#  undef BLI_mempool_create
#endif

BLI_mempool *BLI_mempool_create(uint esize, uint elem_num, uint pchunk, uint flag)
{
  return BLI_mempool_create_ex(esize, elem_num, pchunk, flag, "");
}

BLI_mempool *BLI_mempool_create_ex(
    uint esize, uint elem_num, uint pchunk, uint flag, const char *tag)
{
  BLI_mempool *pool;
  BLI_freenode *last_tail = NULL;
  uint i, maxchunks;

  char buf[512];
  if (tag) {
    sprintf(buf, "Mempool:%s", tag);
  }
  else {
    sprintf(buf, "memory pool");
  }

  char *memtag = strdup(buf);

  /* allocate the pool structure */
  pool = MEM_mallocN(sizeof(BLI_mempool), memtag);

  strcat(buf, " chunk");

#if !defined(WITH_DNA_GHASH) && (defined(__SANITIZE_ADDRESS__) || __has_feature(address_sanitizer))
  BLI_mutex_init(&pool->mutex);
#endif

  pool->memtag = memtag;
  pool->memtag_chunk = strdup(buf);

  pool->totchunk = 0;
  pool->chunktable = NULL;

  /* set the elem size */
  if (esize < (int)MEMPOOL_ELEM_SIZE_MIN) {
    esize = (int)MEMPOOL_ELEM_SIZE_MIN;
  }

  if (flag & BLI_MEMPOOL_ALLOW_ITER) {
    esize = MAX2(esize, (uint)sizeof(BLI_freenode));
  }

  esize += POISON_REDZONE_SIZE;

  maxchunks = mempool_maxchunks(elem_num, pchunk);

  pool->chunks = NULL;
  pool->chunk_tail = NULL;
  pool->esize = esize;

  /* Optimize chunk size to powers of 2, accounting for slop-space. */
#ifdef USE_CHUNK_POW2
  {
    BLI_assert(power_of_2_max_u(pchunk * esize) > CHUNK_OVERHEAD);
    pchunk = (power_of_2_max_u(pchunk * esize) - CHUNK_OVERHEAD) / esize;
  }
#endif

  pool->csize = esize * pchunk;

  /* Ensure this is a power of 2, minus the rounding by element size. */
#if defined(USE_CHUNK_POW2) && !defined(NDEBUG)
  {
    uint final_size = (uint)MEM_SIZE_OVERHEAD + (uint)sizeof(BLI_mempool_chunk) + pool->csize;
    BLI_assert(((uint)power_of_2_max_u(final_size) - final_size) < pool->esize);
  }
#endif

  pool->pchunk = pchunk;
  pool->flag = flag;
  pool->free = NULL; /* mempool_chunk_add assigns */
  pool->maxchunks = maxchunks;
#ifdef USE_TOTALLOC
  pool->totalloc = 0;
#endif
  pool->totused = 0;

  if (elem_num) {
    /* Allocate the actual chunks. */
    for (i = 0; i < maxchunks; i++) {
      BLI_mempool_chunk *mpchunk = mempool_chunk_alloc(pool);
      last_tail = mempool_chunk_add(pool, mpchunk, last_tail);
    }
  }

#ifdef WITH_MEM_VALGRIND
  VALGRIND_CREATE_MEMPOOL(pool, 0, false);
#endif

#ifdef HAVE_MEMPOOL_ASAN
  BLI_asan_poison(pool->poisoned, sizeof(pool->poisoned));
#endif

  mempool_poison_no_unlock(pool);

  return pool;
}

void *BLI_mempool_alloc(BLI_mempool *pool)
{
  mempool_unpoison(pool);

  BLI_freenode *free_pop;

  if (UNLIKELY(pool->free == NULL)) {
    /* Need to allocate a new chunk. */
    BLI_mempool_chunk *mpchunk = mempool_chunk_alloc(pool);
    mempool_chunk_add(pool, mpchunk, NULL);
  }

  free_pop = pool->free;

  BLI_asan_unpoison(free_pop, pool->esize - POISON_REDZONE_SIZE);

  BLI_assert(pool->chunk_tail->next == NULL);

  if (pool->flag & BLI_MEMPOOL_ALLOW_ITER) {
    free_pop->freeword = USEDWORD;
  }

  pool->free = free_pop->next;
  pool->totused++;

#ifdef WITH_MEM_VALGRIND
  VALGRIND_MEMPOOL_ALLOC(pool, free_pop, pool->esize);
#endif

  mempool_poison(pool);

  return (void *)free_pop;
}

void *BLI_mempool_calloc(BLI_mempool *pool)
{
  void *retval = BLI_mempool_alloc(pool);

  mempool_unpoison(pool);
  memset(retval, 0, (size_t)pool->esize - POISON_REDZONE_SIZE);
  mempool_poison(pool);

  return retval;
}

int BLI_mempool_find_real_index(BLI_mempool *pool, void *ptr)
{
  mempool_unpoison(pool);

  BLI_mempool_chunk *chunk = pool->chunks;
  uintptr_t uptr = (uintptr_t)ptr;
  uintptr_t cptr;
  int chunki = 0;

  while (chunk) {
    cptr = (uintptr_t)chunk;

    if (uptr >= cptr && uptr < cptr + pool->csize) {
      break;
    }

    chunk = chunk->next;
    chunki++;
  }

  if (!chunk) {
    mempool_poison(pool);
    return -1;  // failed
  }

  int ret = chunki * (int)pool->pchunk + ((int)(uptr - cptr)) / (int)pool->esize;

  mempool_poison(pool);

  return ret;
}

/*finds an element in pool that's roughly at idx, idx*/
int BLI_mempool_find_elems_fuzzy(
    BLI_mempool *pool, int idx, int range, void **r_elems, int r_elems_size)
{
  mempool_unpoison(pool);

  int istart = idx - range, iend = idx + range;
  istart = MAX2(istart, 0);

  int elem_num = 0;

  for (int i = istart; i < iend; i++) {
    int chunki = i / (int)pool->pchunk;
    if (chunki >= (int)pool->totchunk) {
      break;
    }

    int idx2 = i % (int)pool->pchunk;

    BLI_mempool_chunk *chunk = pool->chunktable[chunki];
    char *data = (char *)CHUNK_DATA(chunk);
    void *ptr = data + idx2 * (int)pool->esize;

    BLI_asan_unpoison(ptr, pool->esize - POISON_REDZONE_SIZE);

    BLI_freenode *fnode = (BLI_freenode *)ptr;
    if (fnode->freeword == FREEWORD) {
      BLI_asan_poison(ptr, pool->esize);
      continue;
    }

    r_elems[elem_num++] = ptr;

    if (elem_num == r_elems_size) {
      break;
    }
  }

  mempool_poison(pool);
  return elem_num;
}

int BLI_mempool_get_size(BLI_mempool *pool)
{
  mempool_unpoison(pool);

  BLI_mempool_chunk *chunk = pool->chunks;
  int ret = 0;

  while (chunk) {
    chunk = chunk->next;

    ret += (int)pool->pchunk;
  }

  mempool_poison(pool);
  return ret;
}

/**
 * Free an element from the mempool.
 *
 * \note doesn't protect against double frees, take care!
 */
void BLI_mempool_free(BLI_mempool *pool, void *addr)
{
  mempool_unpoison(pool);

  BLI_freenode *newhead = addr;

#ifndef NDEBUG
  {
    BLI_mempool_chunk *chunk;
    bool found = false;
    for (chunk = pool->chunks; chunk; chunk = chunk->next) {
      if (ARRAY_HAS_ITEM((char *)addr, (char *)CHUNK_DATA(chunk), pool->csize)) {
        found = true;
        break;
      }
    }
    if (!found) {
      BLI_assert_msg(0, "Attempt to free data which is not in pool.\n");
    }
  }

  /* Enable for debugging. */
  if (UNLIKELY(mempool_debug_memset)) {
    memset(addr, 255, pool->esize - POISON_REDZONE_SIZE);
  }
#endif

  if (pool->flag & BLI_MEMPOOL_ALLOW_ITER) {
#ifndef NDEBUG
    /* This will detect double free's. */
    BLI_assert(newhead->freeword != FREEWORD);
#endif
    newhead->freeword = FREEWORD;
  }

  newhead->next = pool->free;
  pool->free = newhead;

  BLI_asan_poison(newhead, pool->esize);

  pool->totused--;

  if (pool->totused < 0) {
    fprintf(stderr, "Corrupted mempool\n");
    fflush(stderr);
    abort();
  }

#ifdef WITH_MEM_VALGRIND
  VALGRIND_MEMPOOL_FREE(pool, addr);
#endif

  /* Nothing is in use; free all the chunks except the first. */
  if (UNLIKELY(pool->totused == 0) && (pool->chunks->next)) {
    const uint esize = pool->esize;
    BLI_freenode *curnode;
    uint j;
    BLI_mempool_chunk *first;

    first = pool->chunks;
    mempool_chunk_free_all(first->next, pool);
    first->next = NULL;
    pool->chunk_tail = first;

    mempool_update_chunktable(pool);

#ifdef USE_TOTALLOC
    pool->totalloc = pool->pchunk;
#endif

    /* Temp alloc so valgrind doesn't complain when setting free'd blocks 'next'. */
#ifdef WITH_MEM_VALGRIND
    VALGRIND_MEMPOOL_ALLOC(pool, CHUNK_DATA(first), pool->csize);
#endif

    curnode = CHUNK_DATA(first);
    pool->free = curnode;

    j = pool->pchunk;
    while (j--) {
      BLI_asan_unpoison(curnode, pool->esize - POISON_REDZONE_SIZE);
      BLI_freenode *next = curnode->next = NODE_STEP_NEXT(curnode);
      BLI_asan_poison(curnode, pool->esize);
      curnode = next;
    }

    BLI_asan_unpoison(curnode, pool->esize - POISON_REDZONE_SIZE);
    BLI_freenode *prev = NODE_STEP_PREV(curnode);
    BLI_asan_poison(curnode, pool->esize);

    curnode = prev;

    BLI_asan_unpoison(curnode, pool->esize - POISON_REDZONE_SIZE);
    curnode->next = NULL; /* terminate the list */
    BLI_asan_poison(curnode, pool->esize);

#ifdef WITH_MEM_VALGRIND
    VALGRIND_MEMPOOL_FREE(pool, CHUNK_DATA(first));
#endif
  }

  mempool_poison(pool);
}

int BLI_mempool_len(const BLI_mempool *pool)
{
  mempool_unpoison((BLI_mempool *)pool);
  int ret = pool->totused;
  mempool_poison((BLI_mempool *)pool);

  return ret;
}

void *BLI_mempool_findelem(BLI_mempool *pool, uint index)
{
  mempool_unpoison(pool);

  BLI_assert(pool->flag & BLI_MEMPOOL_ALLOW_ITER);

  if (index < (uint)pool->totused) {
    /* We could have some faster mem chunk stepping code inline. */
    BLI_mempool_iter iter;
    void *elem;
    BLI_mempool_iternew(pool, &iter);
    for (elem = BLI_mempool_iterstep(&iter); index-- != 0; elem = BLI_mempool_iterstep(&iter)) {
      /* pass */
    }

    mempool_poison(pool);
    return elem;
  }

  mempool_poison(pool);
  return NULL;
}

void BLI_mempool_as_table(BLI_mempool *pool, void **data)
{
  BLI_mempool_iter iter;
  void *elem;
  void **p = data;

  mempool_unpoison(pool);
  BLI_assert(pool->flag & BLI_MEMPOOL_ALLOW_ITER);
  mempool_poison(pool);

  BLI_mempool_iternew(pool, &iter);

  while ((elem = BLI_mempool_iterstep(&iter))) {
    *p++ = elem;
  }

  mempool_unpoison(pool);
  BLI_assert((int)(p - data) == pool->totused);
  mempool_poison(pool);
}

void **BLI_mempool_as_tableN(BLI_mempool *pool, const char *allocstr)
{
  void **data = MEM_mallocN((size_t)pool->totused * sizeof(void *), allocstr);
  BLI_mempool_as_table(pool, data);
  return data;
}

void BLI_mempool_as_array(BLI_mempool *pool, void *data)
{
  mempool_unpoison(pool);

  const uint esize = pool->esize - (uint)POISON_REDZONE_SIZE;
  BLI_mempool_iter iter;
  char *elem, *p = data;

  BLI_assert(pool->flag & BLI_MEMPOOL_ALLOW_ITER);

  mempool_poison(pool);

  BLI_mempool_iternew(pool, &iter);
  while ((elem = BLI_mempool_iterstep(&iter))) {
    memcpy(p, elem, (size_t)esize);
    p = NODE_STEP_NEXT(p);
  }
}

void *BLI_mempool_as_arrayN(BLI_mempool *pool, const char *allocstr)
{
  mempool_unpoison(pool);
  char *data = MEM_malloc_arrayN((size_t)pool->totused, pool->esize, allocstr);
  mempool_poison(pool);

  BLI_mempool_as_array(pool, data);
  return data;
}

void BLI_mempool_iternew(BLI_mempool *pool, BLI_mempool_iter *iter)
{
  BLI_assert(pool->flag & BLI_MEMPOOL_ALLOW_ITER);

  iter->pool = pool;
  mempool_unpoison(pool);
  iter->curchunk = pool->chunks;
  mempool_poison(pool);
  iter->curindex = 0;
}

static void mempool_threadsafe_iternew(BLI_mempool *pool, BLI_mempool_threadsafe_iter *ts_iter)
{
  BLI_mempool_iternew(pool, &ts_iter->iter);
  ts_iter->curchunk_threaded_shared = NULL;
}

ParallelMempoolTaskData *mempool_iter_threadsafe_create(BLI_mempool *pool, const size_t iter_num)
{
  BLI_assert(pool->flag & BLI_MEMPOOL_ALLOW_ITER);

  ParallelMempoolTaskData *iter_arr = MEM_mallocN(sizeof(*iter_arr) * iter_num, __func__);
  BLI_mempool_chunk **curchunk_threaded_shared = MEM_mallocN(sizeof(void *), __func__);

  mempool_threadsafe_iternew(pool, &iter_arr->ts_iter);

  *curchunk_threaded_shared = iter_arr->ts_iter.iter.curchunk;
  iter_arr->ts_iter.curchunk_threaded_shared = curchunk_threaded_shared;
  for (size_t i = 1; i < iter_num; i++) {
    iter_arr[i].ts_iter = iter_arr[0].ts_iter;
    *curchunk_threaded_shared = iter_arr[i].ts_iter.iter.curchunk =
        ((*curchunk_threaded_shared) ? (*curchunk_threaded_shared)->next : NULL);
  }

  return iter_arr;
}

void mempool_iter_threadsafe_destroy(ParallelMempoolTaskData *iter_arr)
{
  BLI_assert(iter_arr->ts_iter.curchunk_threaded_shared != NULL);

  MEM_freeN(iter_arr->ts_iter.curchunk_threaded_shared);
  MEM_freeN(iter_arr);
}

#if 0
/* unoptimized, more readable */

static void *bli_mempool_iternext(BLI_mempool_iter *iter)
{
  void *ret = NULL;

  if (iter->curchunk == NULL || !iter->pool->totused) {
    return ret;
  }

  ret = ((char *)CHUNK_DATA(iter->curchunk)) + (iter->pool->esize * iter->curindex);

  iter->curindex++;

  if (iter->curindex == iter->pool->pchunk) {
    iter->curindex = 0;
    iter->curchunk = iter->curchunk->next;
  }

  return ret;
}

void *BLI_mempool_iterstep(BLI_mempool_iter *iter)
{
  BLI_freenode *ret;

  do {
    ret = bli_mempool_iternext(iter);
  } while (ret && ret->freeword == FREEWORD);

  return ret;
}

#else /* Optimized version of code above. */

void *BLI_mempool_iterstep(BLI_mempool_iter *iter)
{
  if (UNLIKELY(iter->curchunk == NULL)) {
    return NULL;
  }

  mempool_unpoison(iter->pool);

  const uint esize = iter->pool->esize;
  BLI_freenode *curnode = POINTER_OFFSET(CHUNK_DATA(iter->curchunk), (esize * iter->curindex));
  BLI_freenode *ret;
  do {
    ret = curnode;

    BLI_asan_unpoison(ret, iter->pool->esize - POISON_REDZONE_SIZE);

    if (++iter->curindex != iter->pool->pchunk) {
      curnode = POINTER_OFFSET(curnode, esize);
    }
    else {
      iter->curindex = 0;
      iter->curchunk = iter->curchunk->next;
      if (UNLIKELY(iter->curchunk == NULL)) {
        BLI_asan_unpoison(ret, iter->pool->esize - POISON_REDZONE_SIZE);
        void *ret2 = (ret->freeword == FREEWORD) ? NULL : ret;

        if (ret->freeword == FREEWORD) {
          BLI_asan_poison(ret, iter->pool->esize);
        }

        mempool_poison(iter->pool);
        return ret2;
      }
      curnode = CHUNK_DATA(iter->curchunk);
    }
  } while (ret->freeword == FREEWORD);

  mempool_poison(iter->pool);
  return ret;
}

void *mempool_iter_threadsafe_step(BLI_mempool_threadsafe_iter *ts_iter)
{
  BLI_mempool_iter *iter = &ts_iter->iter;
  if (UNLIKELY(iter->curchunk == NULL)) {
    return NULL;
  }

  mempool_unpoison(iter->pool);

  const uint esize = iter->pool->esize;
  BLI_freenode *curnode = POINTER_OFFSET(CHUNK_DATA(iter->curchunk), (esize * iter->curindex));
  BLI_freenode *ret;
  do {
    ret = curnode;

    BLI_asan_unpoison(ret, esize - POISON_REDZONE_SIZE);

    if (++iter->curindex != iter->pool->pchunk) {
      curnode = POINTER_OFFSET(curnode, esize);
    }
    else {
      iter->curindex = 0;

      /* Begin unique to the `threadsafe` version of this function. */
      for (iter->curchunk = *ts_iter->curchunk_threaded_shared;
           (iter->curchunk != NULL) && (atomic_cas_ptr((void **)ts_iter->curchunk_threaded_shared,
                                                       iter->curchunk,
                                                       iter->curchunk->next) != iter->curchunk);
           iter->curchunk = *ts_iter->curchunk_threaded_shared) {
        /* pass. */
      }
      if (UNLIKELY(iter->curchunk == NULL)) {
        if (ret->freeword == FREEWORD) {
          BLI_asan_poison(ret, esize);
          mempool_poison(iter->pool);
          return NULL;
        }
        else {
          mempool_poison(iter->pool);
          return ret;
        }
      }
      /* End `threadsafe` exception. */

      iter->curchunk = iter->curchunk->next;
      if (UNLIKELY(iter->curchunk == NULL)) {
        if (ret->freeword == FREEWORD) {
          BLI_asan_poison(ret, iter->pool->esize);
          mempool_poison(iter->pool);
          return NULL;
        }
        else {
          mempool_poison(iter->pool);
          return ret;
        }
      }

      curnode = CHUNK_DATA(iter->curchunk);
    }

    if (ret->freeword == FREEWORD) {
      BLI_asan_poison(ret, iter->pool->esize);
    }
    else {
      break;
    }
  } while (true);

  mempool_poison(iter->pool);
  return ret;
}

#endif

void BLI_mempool_clear_ex(BLI_mempool *pool, const int elem_num_reserve)
{
  mempool_unpoison(pool);

  BLI_mempool_chunk *mpchunk;
  BLI_mempool_chunk *mpchunk_next;
  uint maxchunks;

  BLI_mempool_chunk *chunks_temp;
  BLI_freenode *last_tail = NULL;

#ifdef WITH_MEM_VALGRIND
  VALGRIND_DESTROY_MEMPOOL(pool);
  VALGRIND_CREATE_MEMPOOL(pool, 0, false);
#endif

  if (elem_num_reserve == -1) {
    maxchunks = pool->maxchunks;
  }
  else {
    maxchunks = mempool_maxchunks((uint)elem_num_reserve, pool->pchunk);
  }

  /* Free all after 'pool->maxchunks'. */
  mpchunk = mempool_chunk_find(pool->chunks, maxchunks - 1);
  if (mpchunk && mpchunk->next) {
    /* terminate */
    mpchunk_next = mpchunk->next;
    mpchunk->next = NULL;
    mpchunk = mpchunk_next;

    do {
      mpchunk_next = mpchunk->next;
      mempool_chunk_free(mpchunk, pool);
    } while ((mpchunk = mpchunk_next));
  }

  /* re-initialize */
  pool->free = NULL;
  pool->totused = 0;
#ifdef USE_TOTALLOC
  pool->totalloc = 0;
#endif

  chunks_temp = pool->chunks;
  pool->chunks = NULL;
  pool->chunk_tail = NULL;

  while ((mpchunk = chunks_temp)) {
    chunks_temp = mpchunk->next;
    last_tail = mempool_chunk_add(pool, mpchunk, last_tail);
  }

  mempool_poison(pool);
}

void BLI_mempool_clear(BLI_mempool *pool)
{
  BLI_mempool_clear_ex(pool, -1);
}

void BLI_mempool_destroy(BLI_mempool *pool)
{
  mempool_unpoison(pool);
  mempool_chunk_free_all(pool->chunks, pool);

  if (pool->memtag) {
    free(pool->memtag);
  }
  if (pool->memtag_chunk) {
    free(pool->memtag_chunk);
  }

  MEM_SAFE_FREE(pool->chunktable);

#ifdef WITH_MEM_VALGRIND
  VALGRIND_DESTROY_MEMPOOL(pool);
#endif

#ifdef HAVE_MEMPOOL_ASAN
  BLI_asan_unpoison(pool->poisoned, sizeof(pool->poisoned));
#endif

  MEM_freeN(pool);
}

#ifndef NDEBUG
void BLI_mempool_set_memory_debug(void)
{
  mempool_debug_memset = true;
}
#endif
