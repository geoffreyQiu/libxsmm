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
  std::cout << "EDGE Local Update Reproducer" << std::endl;

  double                                  l_dT = 0.01;
  unsigned int                            l_nSteps;
  unsigned int                            l_nElements;
  t_elementChars                        * l_elChars; /* zero initialization */
  t_dg                                    l_dg;
  t_matStar                            (* l_starM)[N_DIM];
  t_fluxSolver                         (* l_fluxSolvers)[ C_ENT[T_SDISC.ELEMENT].N_FACES ];
  real_base                            (* l_dofs)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS];
  real_base                            (* l_tInt)[N_QUANTITIES][N_ELEMENT_MODES][N_CRUNS];
  edge::io::Receivers                     l_recvs;
  edge::data::MmXsmmFused< real_base >    l_mm;
  unsigned int                            l_dummyUInt;
  double                                  l_dummyDouble;

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

  // zero init - disable read/write recvs
  l_elChars = new t_elementChars[l_nElements];
  for ( unsigned int l_el = 0; l_el < l_nElements; l_el++ ) l_elChars[l_el].spType = 0;


  // 3. Run solvers
  unsigned long long l_start = libxsmm_timer_tick();
#ifdef PP_USE_OMP
  #pragma omp parallel firstprivate( l_nSteps, l_nElements, l_dT )  \
                       firstprivate( l_elChars )                    \
                       firstprivate( l_dg, l_starM, l_fluxSolvers ) \
                       firstprivate( l_dofs, l_tInt  )              \
                       firstprivate( l_mm )                         \
                       private( l_recvs, l_dummyUInt, l_dummyDouble )
#endif
  {
    const unsigned int l_nThreads = omp_get_num_threads();
    const unsigned int l_tid = omp_get_thread_num();
    unsigned int l_firstEl = (unsigned int)((l_nElements + l_nThreads - 1) / l_nThreads) * l_tid;
    unsigned int l_lastEl = (unsigned int)((l_nElements + l_nThreads - 1) / l_nThreads) * (l_tid + 1);
    l_lastEl = std::min(l_lastEl, l_nElements);
    unsigned int l_numEl = l_lastEl - l_firstEl;

    for ( unsigned int l_step = 0; l_step < l_nSteps; l_step++ )
      edge::elastic::solvers::AderDg::local< unsigned int,
                                             real_base,
                                             edge::data::MmXsmmFused< real_base > > 
                                           ( l_firstEl,
                                             l_numEl,
                                             l_dummyDouble,
                                             l_dT,
                                             l_dummyUInt,
                                             l_dummyUInt,
                                             nullptr,
                                             nullptr,
                                             l_elChars,
                                             l_dg,
                                             l_starM,
                                             l_fluxSolvers,
                                             l_dofs,
                                             l_tInt,
                                             nullptr,
                                             l_recvs,
                                             l_mm           );
  }
  unsigned long long l_end = libxsmm_timer_tick();

  // 4. Print statistics
  double l_time = libxsmm_timer_duration(l_start, l_end);
  unsigned int l_local_flops[] =
  {
    792, 3564, 11412, 31500, 77184, 173538, 360522
  };
  unsigned long long l_flops = (unsigned long long)l_local_flops[ORDER-1] * PP_N_CRUNS * \
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

  delete[] l_elChars;
  
  
  return 0;
}
