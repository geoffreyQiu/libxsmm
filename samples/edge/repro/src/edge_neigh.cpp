#include <iostream>
#include <string>
#include <cstdlib>

#ifdef PP_USE_OMP
#include <omp.h>
#else
#define omp_get_max_threads() 1
#define omp_get_num_threads() 1
#define omp_get_thread_num()  0
#endif

#include "parallel_global.hpp"
#include "AderDg.hpp"
#include "edge_setup.hpp"

#if defined(PP_T_KERNELS_XSMM) || defined(PP_T_KERNELS_XSMM_DENSE_SINGLE)
#include "libxsmm.h"
#endif


int main(int i_argc, char *i_argv[]) {
  std::cout << std::endl;
  std::cout << "EDGE Neighboring Update Reproducer" << std::endl;

  double                                  l_dT = 0.01;
  unsigned int                            l_nSteps;
  unsigned int                            l_nElements;
  unsigned int                         (* l_elFa)[ C_ENT[T_SDISC.ELEMENT].N_FACES ];
  t_faceChars                             l_faChars;
  unsigned int                         (* l_elFaEl)[ C_ENT[T_SDISC.ELEMENT].N_FACES ];
  unsigned short                       (* l_fIdElFaEl)[ C_ENT[T_SDISC.ELEMENT].N_FACES ];
  unsigned short                       (* l_vIdElFaEl)[ C_ENT[T_SDISC.ELEMENT].N_FACES ];
  t_dg                                    l_dg;
  t_matStar                            (* l_starM)[N_DIM];
  t_fluxSolver                         (* l_fluxSolvers)[ C_ENT[T_SDISC.ELEMENT].N_FACES ];
  real_base                            (* l_dofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS];
  real_base                            (* l_tInt)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS];
  edge::data::MmXsmmFused< real_base >    l_mm;
  unsigned int                            l_dummyUInt;

  // 1. Parse cmd arguments for hyper-parameters
  if ( i_argc == 3 ) {
    l_nSteps    = (unsigned int)atoi(i_argv[1]);
    l_nElements = (unsigned int)atoi(i_argv[2]);    
  } else {
    std::cout << "Usage: ./local {NUM_STEPS} {NUM_ELEMENTS} [-h|--help]\n" << std::endl;
    std::exit(1);
  }
  std::cout << "#Steps: " << l_nSteps << ", #Elements: " << l_nElements << std::endl;
  std::cout << std::endl;

  // 2. Set up structures
  setupDg( l_dg );
  setupStarM( l_nElements, l_starM );
  setupFluxSolv( l_nElements, l_fluxSolvers );
  setupKernel( l_mm );

  setupTensor( l_nElements, l_dofs, l_tInt );
#ifdef PP_USE_OMP
  #pragma omp parallel
  #pragma omp critical
#endif
  setupScratchMem( edge::parallel::g_scratchMem );

  // use one pseudo face for all elements
  l_elFa = (unsigned int (*)[ C_ENT[T_SDISC.ELEMENT].N_FACES ]) new unsigned int[ l_nElements * C_ENT[T_SDISC.ELEMENT].N_FACES ];
  for ( unsigned int l_el = 0; l_el < l_nElements; l_el++ ) {
    for ( unsigned int l_fa = 0; l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fa++ ) {
      l_elFa[l_el][l_fa] = 0;
    }
  }
  // zero init - force derive neighboring elmt from pseudo mesh
  l_faChars.spType = 0;

  /* setup pseudo mesh - l_elFaEl    : neighboring element id
   *                     l_fIdElFaEl : neighboring face id
   *                     l_vIdElFaEl : neighboring face orientation
   */
  l_elFaEl = (unsigned int (*)[ C_ENT[T_SDISC.ELEMENT].N_FACES ]) new unsigned int[ l_nElements * C_ENT[T_SDISC.ELEMENT].N_FACES ];
  l_fIdElFaEl = (unsigned short (*)[ C_ENT[T_SDISC.ELEMENT].N_FACES ]) new unsigned short[ l_nElements * C_ENT[T_SDISC.ELEMENT].N_FACES ];
  l_vIdElFaEl = (unsigned short (*)[ C_ENT[T_SDISC.ELEMENT].N_FACES ]) new unsigned short[ l_nElements * C_ENT[T_SDISC.ELEMENT].N_FACES ];
  for ( unsigned int l_el = 0; l_el < l_nElements; l_el++ ) {
    for ( unsigned int l_fa = 0; l_fa < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fa++ ) {
      l_elFaEl[l_el][l_fa] = (unsigned int)((unsigned int)rand() % l_nElements);
      l_fIdElFaEl[l_el][l_fa] = (unsigned short)(unsigned short)rand() % C_ENT[T_SDISC.ELEMENT].N_FACES;
      l_vIdElFaEl[l_el][l_fa] = 0;
    }
  }


  // 3. Run solvers
  unsigned long long l_start = libxsmm_timer_tick();
#ifdef PP_USE_OMP
  #pragma omp parallel firstprivate( l_nSteps, l_nElements, l_dT )  \
                       firstprivate( l_elFa, l_faChars )            \
                       firstprivate( l_dg, l_starM, l_fluxSolvers ) \
                       firstprivate( l_dofs, l_tInt )               \
                       firstprivate( l_mm )                         \
                       private( l_dummyUInt )
#endif
  {
    const unsigned int l_nThreads = omp_get_num_threads();
    const unsigned int l_tid = omp_get_thread_num();
    unsigned int l_firstEl = (unsigned int)((l_nElements + l_nThreads - 1) / l_nThreads) * l_tid;
    unsigned int l_lastEl = (unsigned int)((l_nElements + l_nThreads - 1) / l_nThreads) * (l_tid + 1);
    l_lastEl = std::min(l_lastEl, l_nElements);
    unsigned int l_numEl = l_lastEl - l_firstEl;

    for ( unsigned int l_step = 0; l_step < l_nSteps; l_step++ )
      edge::elastic::solvers::AderDg::neigh< unsigned int,
                                             real_base,
                                             edge::data::MmXsmmFused< real_base > > 
                                           ( l_firstEl,
                                             l_numEl,
                                             l_dummyUInt,
                                             l_dg,
                                            &l_faChars,
                                             l_fluxSolvers,
                                             l_elFa,
                                             l_elFaEl,
                                             nullptr,
                                             nullptr,
                                             l_fIdElFaEl,
                                             l_vIdElFaEl,
                                             l_tInt,
                                             nullptr,
                                             l_dofs,                                             
                                             l_mm           );
  }
  unsigned long long l_end = libxsmm_timer_tick();

  // 4. Print statistics
  double l_time = libxsmm_timer_duration(l_start, l_end);
  unsigned int l_neigh_flops[] =
  {
    792, 3078, 8532, 20502, 43848, 86832, 159516
  };
  unsigned long long l_flops = (unsigned long long)l_neigh_flops[ORDER-1] * PP_N_CRUNS * \
                               l_nElements * l_nSteps;
  double l_gflops = (double)l_flops / (l_time * 1000000000);
  std::cout << "Elapsed time: " << l_time << "s" << std::endl;
  std::cout << "Performance:  " << l_gflops << "GFLOPS" << std::endl;
  std::cout << std::endl;


  // 5. Clean up
  for ( unsigned int l_st = 0; l_st < (ORDER-1)*N_DIM; l_st++ )
    delete[] l_dg.mat.stiffT[l_st];
  for ( unsigned int l_sv = 0; l_sv < N_DIM; l_sv++ )
    delete[] l_dg.mat.stiff[l_sv];
  for ( unsigned int l_fl = 0; l_fl < C_ENT[T_SDISC.ELEMENT].N_FACES; l_fl++ )
    delete[] l_dg.mat.fluxL[l_fl];
  for ( unsigned int l_fn = 0; l_fn < N_FLUXN_MATRICES; l_fn++ )
    delete[] l_dg.mat.fluxN[l_fn];
  for ( unsigned int l_ft = 0; l_ft < C_ENT[T_SDISC.ELEMENT].N_FACES; l_ft++ )
    delete[] l_dg.mat.fluxT[l_ft];
  delete[] (t_matStar *)l_starM;
  delete[] (t_fluxSolver *)l_fluxSolvers;

  free(l_dofs);
  free(l_tInt);

#ifdef PP_USE_OMP
  #pragma omp parallel
  #pragma omp critical
#endif
  free( edge::parallel::g_scratchMem );

  delete[] (unsigned int *)l_elFa;
  delete[] (unsigned int *)l_elFaEl;
  delete[] (unsigned short *)l_fIdElFaEl;
  delete[] (unsigned short *)l_vIdElFaEl;
  
  return 0;
}
