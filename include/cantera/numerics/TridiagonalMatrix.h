/**
 *  @file TridiagonalMatrix.h
 *   Declarations for the class TridiagonalMatrix
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_TRIDIAGONALMATRIX_H
#define CT_TRIDIAGONALMATRIX_H

#include "cantera/base/ctexceptions.h"
#include "cantera/base/global.h"

namespace Cantera
{

//! A class for the general tridiagonal matrices, involving linear solving.
//! The class is based upon the LAPACK banded storage matrix format.
class TridiagonalMatrix
{
public:
    //! Base Constructor
    /*!
     * Create an \c 0 by \c 0 matrix, and initialize all elements to \c 0.
     */
    TridiagonalMatrix();

    //! Creates a banded matrix and sets all elements to zero
    /*!
     * Create an \c n by \c n  banded matrix, and initialize all elements to \c v.
     *
     * @param n   size of the square matrix
     * @param v   initial value of all matrix components.
     */
    TridiagonalMatrix(size_t n, double v = 0.0);

    ~TridiagonalMatrix() {};

    //! Resize the matrix problem
    /*!
     * All data is lost
     *
     * @param n   size of the square matrix
     * @param v   initial value of all matrix components.
     */
    void resize(size_t n, double v = 0.0);

    //! Fill or zero the matrix
    /*!
     *  @param v  Fill value, defaults to zero.
     */
    void bfill(double v = 0.0);

    //! Solve the matrix problem Ax = b
    /*!
     * @param b     INPUT RHS of the problem
     *              OUTPUT solution to the problem
     * @param nrhs  Number of right hand sides to solve
     * @param ldb   Leading dimension of `b`. Default is nColumns()
     * @returns a success flag. 0 indicates a success; ~0 indicates some error
     *     occurred, see the LAPACK documentation
     */
    int solve(vector_fp& b, size_t nrhs=1, size_t ldb=0);

    //! Return a changeable reference to element of the diagonal
    /*!
     * @param i  row
     * @returns a reference to the value of the matrix entry
     */
    inline double& valueDiagonal(size_t i) {
        return m_d[i];
    }

    //! Return the value of the diagonal
    /*!
     * This method does not alter the array.
     * @param i  row
     * @returns the value of the matrix entry
     */
    inline double valueDiagonal(size_t i) const {
        return m_d[i];
    }

    //! Return a changeable reference to element of the diagonal
    /*!
     * @param i  row
     * @returns a reference to the value of the matrix entry
     */
    inline double& valueSubdiagonal(size_t i) {
        return m_dl[i];
    }

    //! Return the value of the diagonal
    /*!
     * This method does not alter the array.
     * @param i  row
     * @returns the value of the matrix entry
     */
    inline double valueSubdiagonal(size_t i) const {
        return m_dl[i];
    }

    //! Return a changeable reference to element of the diagonal
    /*!
     * @param i  row
     * @returns a reference to the value of the matrix entry
     */
    inline double& valueSuperdiagonal(size_t i) {
        return m_du[i];
    }

    //! Return the value of the diagonal
    /*!
     * This method does not alter the array.
     * @param i  row
     * @returns the value of the matrix entry
     */
    inline double valueSuperDiagonal(size_t i) const {
        return m_du[i];
    }

    inline size_t size() const {
        return m_n;
    }

    inline size_t nColumns() const {
        return m_n;
    }

    inline size_t nRows() const {
        return m_n;
    }

protected:

    //! Number of rows and columns of the matrix
    size_t m_n;

    //! Matrix data

    vector_fp m_d; //!< diagonal n elements
    vector_fp m_dl;//!< sub-diagonal n-1 elements
    vector_fp m_du;//!< super-diagonal n-1 elements

};

}

#endif
