#ifndef EDGE_REPRODUCER_HELPER_HPP
#define EDGE_REPRODUCER_HELPER_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "constants.hpp"

int readSparseMatrixCsc ( std::string  const     i_fileName,
                          real_base          * & io_matNZ,
                          unsigned int       * & io_matColPtr,
                          unsigned int       * & io_matRowIdx ) {
  abort();
  return 0;
};

int readSparseMatrixCsc ( std::string            const   i_fileName,
                          std::vector< real_base >     & io_matVal,
                          std::vector< unsigned int >  & io_matColPtr,
                          std::vector< unsigned int >  & io_matRowIdx ) {
  std::ifstream l_fp( i_fileName );
  std::string l_lineBuf;

  unsigned int l_header = 0;
  unsigned int l_nEntries;  
  unsigned int l_nCols;  
  unsigned int l_nRows;  
  unsigned int l_row;
  unsigned int l_col; 
  real_base    l_entry;
  unsigned int l_nzCounter;
  unsigned int l_colCounter;
  int          l_errCheck;
  
  while (l_fp) {
    std::getline(l_fp, l_lineBuf);
    if ( l_lineBuf.length() == 0 || l_lineBuf[0] == '%' ) continue;
    if (l_header == 0) {
      l_errCheck= sscanf(l_lineBuf.c_str(), "%u %u %u", &l_nRows, &l_nCols, &l_nEntries);
      if (l_errCheck != 3) abort();

      io_matVal.resize(l_nEntries);
      io_matColPtr.resize(l_nCols+1);
      io_matRowIdx.resize(l_nEntries);

      l_nzCounter = 0;
      l_colCounter = 0;

      l_header = 1;
    } else {
#if PP_PRECISION==32
      l_errCheck= sscanf(l_lineBuf.c_str(), "%u %u %f", &l_row, &l_col, &l_entry);
#elif PP_PRECISION==64
      l_errCheck= sscanf(l_lineBuf.c_str(), "%u %u %lf", &l_row, &l_col, &l_entry);
#endif
      if (l_errCheck != 3) abort();

      io_matVal[l_nzCounter] = l_entry;
      io_matRowIdx[l_nzCounter] = l_row - 1;
      for ( unsigned int l_cc = l_colCounter; l_cc < l_col; l_cc++ ) {
        io_matColPtr[l_cc] = l_nzCounter;
      }
      l_nzCounter++;
      l_colCounter = l_col;
    }
  }
  assert ( l_nzCounter == l_nEntries );

  for ( unsigned int l_cc = l_colCounter; l_cc < l_nCols+1; l_cc++ ) {
        io_matColPtr[l_cc] = l_nzCounter;
  }
  assert ( io_matColPtr.back() == l_nEntries );

  return 0;
};

int readSparseMatrixCsr ( std::string            const   i_fileName,
                          std::vector< real_base >     & io_matVal,
                          std::vector< unsigned int >  & io_matRowPtr,
                          std::vector< unsigned int >  & io_matColIdx ) {
  std::ifstream l_fp( i_fileName );
  std::string l_lineBuf;

  unsigned int l_header = 0;
  unsigned int l_nEntries;  
  unsigned int l_nCols;  
  unsigned int l_nRows;  
  unsigned int l_row;
  unsigned int l_col; 
  real_base    l_entry;
  unsigned int l_nzCounter;
  unsigned int l_rowCounter;
  int          l_errCheck;
  
  while (l_fp) {
    getline(l_fp, l_lineBuf);
    if ( l_lineBuf.length() == 0 || l_lineBuf[0] == '%' ) continue;
    if (l_header == 0) {
      l_errCheck= sscanf(l_lineBuf.c_str(), "%u %u %u", &l_nRows, &l_nCols, &l_nEntries);
      if (l_errCheck != 3) abort();

      io_matVal.resize(l_nEntries);
      io_matRowPtr.resize(l_nRows+1);
      io_matColIdx.resize(l_nEntries);

      l_nzCounter = 0;
      l_rowCounter = 0;

      l_header = 1;
    } else {
#if PP_PRECISION==32
      l_errCheck= sscanf(l_lineBuf.c_str(), "%u %u %f", &l_row, &l_col, &l_entry);
#elif PP_PRECISION==64
      l_errCheck= sscanf(l_lineBuf.c_str(), "%u %u %lf", &l_row, &l_col, &l_entry);
#endif
      if (l_errCheck != 3) abort();

      io_matVal[l_nzCounter] = l_entry;
      io_matColIdx[l_nzCounter] = l_col - 1;
      for ( unsigned int l_rr = l_rowCounter; l_rr < l_row; l_rr++ ) {
        io_matRowPtr[l_rr] = l_nzCounter;
      }
      l_nzCounter++;
      l_rowCounter = l_row;
    }
  }
  assert ( l_nzCounter == l_nEntries );

  for ( unsigned int l_rr = l_rowCounter; l_rr < l_nRows+1; l_rr++ ) {
        io_matRowPtr[l_rr] = l_nzCounter;
  }
  assert ( io_matRowPtr.back() == l_nEntries );

  return 0;
};

int selectSubSparseMatrixCsc ( std::vector< real_base >    const & i_matVal,
                               std::vector< unsigned int > const & i_matColPtr,
                               std::vector< unsigned int > const & i_matRowIdx,
                               unsigned int                const   i_nSubRows,
                               unsigned int                const   i_nSubCols,
                               std::vector< real_base >          & o_subMatVal,
                               std::vector< unsigned int >       & o_subMatColPtr,
                               std::vector< unsigned int >       & o_subMatRowIdx
                             ) {
  unsigned int l_nCols = i_matColPtr.size() - 1;
  unsigned int l_nEntries = i_matVal.size();
  unsigned int l_subNzCounter = 0;

  o_subMatVal.resize(l_nEntries); 
  o_subMatColPtr.resize(i_nSubCols+1);
  o_subMatRowIdx.resize(l_nEntries);

  for ( unsigned int l_cc = 0; l_cc < l_nCols; l_cc++ ) {
    if ( l_cc >= i_nSubCols ) break;

    o_subMatColPtr[l_cc] = l_subNzCounter;
    for ( unsigned int l_nz = i_matColPtr[l_cc]; l_nz < i_matColPtr[l_cc+1]; l_nz++ ) {
      if ( i_matRowIdx[l_nz] >= i_nSubRows ) break;

      o_subMatVal[l_subNzCounter] = i_matVal[l_nz];
      o_subMatRowIdx[l_subNzCounter] = i_matRowIdx[l_nz];
      l_subNzCounter++;
    }
  }

  for ( unsigned int l_cc = l_nCols; l_cc < i_nSubCols + 1; l_cc) {
    o_subMatColPtr[l_cc] = l_subNzCounter;
  }

  o_subMatVal.resize(l_subNzCounter);
  o_subMatRowIdx.resize(l_subNzCounter);
                          
  return 0;
};

#endif
