/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2020 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#include <assert.h>
#include <omp.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "grid_common.h"
#include "grid_constants.h"
#include "grid_library.h"

static grid_library_globals **per_thread_globals = NULL;
static bool library_initialized = false;
static grid_library_config config = {.backend = GRID_BACKEND_AUTO,
                                     .device_id = 0,
                                     .validate = false,
                                     .apply_cutoff = false,
                                     .queue_length = 8192};

#if !defined(_OPENMP)
#error "OpenMP is required. Please add -fopenmp to your C compiler flags."
#endif

/*******************************************************************************
 * \brief Initializes the grid library.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_init() {
  if (library_initialized) {
    printf("Error: Grid library was already initialized.\n");
    abort();
  }

  const int n = omp_get_max_threads();
  per_thread_globals = malloc(n * sizeof(grid_library_globals *));

// Using parallel regions to ensure memory is allocated near a thread's core.
#pragma omp parallel default(none) shared(per_thread_globals) num_threads(n)
  {
    const int ithread = omp_get_thread_num();
    per_thread_globals[ithread] = malloc(sizeof(grid_library_globals));
    memset(per_thread_globals[ithread], 0, sizeof(grid_library_globals));
  }

  library_initialized = true;
}

/*******************************************************************************
 * \brief Finalizes the grid library.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_finalize() {
  if (!library_initialized) {
    printf("Error: Grid library is not initialized.\n");
    abort();
  }

  for (int i = 0; i < omp_get_max_threads(); i++) {
    grid_sphere_cache_free(&per_thread_globals[i]->sphere_cache);
    free(per_thread_globals[i]);
  }
  free(per_thread_globals);
  per_thread_globals = NULL;
  library_initialized = false;
}

/*******************************************************************************
 * \brief Returns a pointer to the thread local sphere cache.
 * \author Ole Schuett
 ******************************************************************************/
grid_sphere_cache *grid_library_get_sphere_cache() {
  return &per_thread_globals[omp_get_thread_num()]->sphere_cache;
}

/*******************************************************************************
 * \brief Configures the grid library.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_set_config(const enum grid_backend backend,
                             const int device_id, const bool validate,
                             const bool apply_cutoff, const int queue_length) {
  config.backend = backend;
  config.device_id = device_id;
  config.validate = validate;
  config.apply_cutoff = apply_cutoff;
  config.queue_length = queue_length;
}

/*******************************************************************************
 * \brief Returns the library config.
 * \author Ole Schuett
 ******************************************************************************/
grid_library_config grid_library_get_config() { return config; }

/*******************************************************************************
 * \brief Adds given increment to counter specified by lp, backend, and kernel.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_counter_add(const int lp, const enum grid_backend backend,
                              const enum grid_library_kernel kernel,
                              const int increment) {
  const int back = backend - GRID_BACKEND_REF;
  const int idx = back * 4 * 20 + kernel * 20 + imin(lp, 19);
  per_thread_globals[omp_get_thread_num()]->counters[idx] += increment;
}

/*******************************************************************************
 * \brief Comperator passed to qsort to compare two counters.
 * \author Ole Schuett
 ******************************************************************************/
static int compare_counters(const void *a, const void *b) {
  return *(long *)b - *(long *)a;
}

/*******************************************************************************
 * \brief Prints statistics gathered by the grid library.
 * \author Ole Schuett
 ******************************************************************************/
void grid_library_print_stats(void (*mpi_sum_func)(long *, int),
                              const int mpi_comm,
                              void (*print_func)(char *, int),
                              const int output_unit) {
  if (!library_initialized) {
    printf("Error: Grid library is not initialized.\n");
    abort();
  }

  // Sum all counters across threads and mpi ranks.
  long counters[320][2] = {0};
  double total = 0.0;
  for (int i = 0; i < 320; i++) {
    counters[i][1] = i;
    for (int j = 0; j < omp_get_max_threads(); j++) {
      counters[i][0] += per_thread_globals[j]->counters[i];
    }
    mpi_sum_func(&counters[i][0], mpi_comm);
    total += counters[i][0];
  }

  // Sort counters.
  qsort(counters, 320, 2 * sizeof(long), &compare_counters);

  // Print counters.
  print_func("\n", output_unit);
  print_func(" ----------------------------------------------------------------"
             "---------------\n",
             output_unit);
  print_func(" -                                                               "
             "              -\n",
             output_unit);
  print_func(" -                                GRID STATISTICS                "
             "              -\n",
             output_unit);
  print_func(" -                                                               "
             "              -\n",
             output_unit);
  print_func(" ----------------------------------------------------------------"
             "---------------\n",
             output_unit);
  print_func(" LP    KERNEL             BACKEND                              "
             "COUNT     PERCENT\n",
             output_unit);

  const char *kernel_names[] = {"collocate ortho", "integrate ortho",
                                "collocate general", "integrate general"};
  const char *backend_names[] = {"REF", "CPU", "GPU", "HYBRID"};

  for (int i = 0; i < 320; i++) {
    if (counters[i][0] == 0)
      continue; // skip empty counters
    const double percent = 100.0 * counters[i][0] / total;
    const int idx = counters[i][1];
    const int back = idx / 80;
    const int kern = (idx % 80) / 20;
    const int lp = (idx % 80) % 20;
    char buffer[100];
    snprintf(buffer, sizeof(buffer), " %-5i %-17s  %-6s  %34li %10.2f%%\n", lp,
             kernel_names[kern], backend_names[back], counters[i][0], percent);
    print_func(buffer, output_unit);
  }

  print_func(" ----------------------------------------------------------------"
             "---------------\n",
             output_unit);
}

// EOF
