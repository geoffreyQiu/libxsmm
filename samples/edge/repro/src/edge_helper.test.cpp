#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include "constants.hpp"
#include "edge_helper.hpp"

int main( int i_argc, char *i_argv[] ) {
  {
    std::string l_cscFileName[] = 
    { 
      "../mats/tet4_4_stiffT_2_csc.mtx",
      "../mats/tet4_4_fluxL_2_csc.mtx",
      "../mats/tet4_4_fluxT_2_csc.mtx"
    };

    unsigned int l_cscMatSize[][3] = 
    { 
      {35, 35, 287},
      {35, 15, 363},
      {15, 35, 363}
    };

    std::vector< real_base >       l_matVal;
    std::vector< unsigned int > l_matColPtr;
    std::vector< unsigned int > l_matRowIdx;

    for ( unsigned int l_i = 0; l_i < 3; l_i++ ) {
      readSparseMatrixCsc( l_cscFileName[l_i],
                           l_matVal,
                           l_matColPtr,
                           l_matRowIdx );
      assert( l_matVal.size() == l_cscMatSize[l_i][2] );
      assert( l_matColPtr.size() == l_cscMatSize[l_i][1]+1 );
    }
  }

  {
    std::string l_csrFileName[] = 
    { 
      "../mats/tet4_4_stiffT_2_csr.mtx",
      "../mats/tet4_4_fluxL_2_csr.mtx",
      "../mats/tet4_4_fluxT_2_csr.mtx"
    };

    unsigned int l_csrMatSize[][3] = 
    { 
      {35, 35, 287},
      {35, 15, 363},
      {15, 35, 363}
    };

    std::vector< real_base >       l_matVal;
    std::vector< unsigned int > l_matRowPtr;
    std::vector< unsigned int > l_matColIdx;

    for ( unsigned int l_i = 0; l_i < 3; l_i++ ) {
      readSparseMatrixCsr( l_csrFileName[l_i],
                           l_matVal,
                           l_matRowPtr,
                           l_matColIdx );
      assert( l_matVal.size() == l_csrMatSize[l_i][2] );
      assert( l_matRowPtr.size() == l_csrMatSize[l_i][0]+1 );
    }
  }


  return 0;
}
