/*
 * matrix.h
 *
 *  Created on: 03-Jul-2014
 *      Author: amit
 */

#ifndef MATRIX_H_
#define MATRIX_H_

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "expm.h"
namespace math {

template<typename scalar_type> class matrix: public boost::numeric::ublas::matrix<
		scalar_type> {
public:
	typedef boost::numeric::ublas::matrix<scalar_type> ublas_matrix_impl;
	typedef typename boost::numeric::ublas::matrix<scalar_type>::size_type size_type;
	typedef boost::numeric::ublas::vector<scalar_type> ublas_vector_impl;
	typedef typename boost::numeric::ublas::matrix<scalar_type>::array_type array_type;

	matrix();
	matrix(size_type r, size_type c);
	//matrix(int row, int col);			will try later
	matrix(size_type r, size_type c, array_type data);
	void matrix_exponentiation(math::matrix<scalar_type>& res, double time_tau) const;
	void matrix_exponentiation(math::matrix<scalar_type>& res) const;
	void multiply(matrix& A, matrix& res);

	//calling matrix as minuend, passed as subtrahend and result as difference (minuend − subtrahend =	difference)
	void minus(matrix& B, matrix& res);
	void mult_vector(std::vector<scalar_type> v, std::vector<scalar_type> &res);
	//void transpose_assign();
	void transpose(matrix& res);
	/*
	 * Appends a column vector to the end of the calling matrix and returns the new resized matrix
	 */
//	void addColumn(std::vector <scalar_type> columnVector, math::matrix<scalar_type>& resized_matrix);
	scalar_type norm_inf();
	void matrix_copy(math::matrix<scalar_type>& destination);
	void matrix_join(const math::matrix<scalar_type> mat2,
			math::matrix<scalar_type>& joined_matrix);
	//void matrix_Identity(math::matrix<scalar_type>& newIdentityMatrix);
	/*
	 * inverse of a matrix : Returns True if Inverse Exists otherwise returns False
	 */
	bool inverse(math::matrix<scalar_type>& res);

	/*private:
	 ublas_matrix_impl my_matrix;*/
};

#include "matrix.cpp"

} // end of math namespace

#endif /* MATRIX_H_ */
