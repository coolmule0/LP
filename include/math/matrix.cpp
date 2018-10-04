/**
 * Matrix class with operations
 */

#ifndef MATRIX_CPP_
#define MATRIX_CPP_

#include <iostream>
#include "matrix.h"
#include "basic_functions.h"

using namespace std;

/** constructors */

template<typename scalar_type> math::matrix<scalar_type>::matrix() {
	ublas_matrix_impl();
}
template<typename scalar_type> math::matrix<scalar_type>::matrix(size_type r,
		size_type c) :
		ublas_matrix_impl(r, c) {
}

template<typename scalar_type> math::matrix<scalar_type>::matrix(size_type r,
		size_type c, array_type data) :
		ublas_matrix_impl(r, c, data) {
}

/**
 * returns the ublas matrix

 const void get_matrix_impl(ublas_matrix_impl &x) const {
 ublas_matrix_impl m(this->size1(), this->size2(), this->data());
 return m;
 } */
/**
 * returns the exponentiated matrix U = exp(A*time_tau)
 */

template<typename scalar_type> void math::matrix<scalar_type>::matrix_exponentiation(
		math::matrix<scalar_type>& res, double time_tau) const {

	/*ublas_matrix_impl emat(size1(),size2());
	 for(unsinged int i=0;i<size1()i++){
	 for(unsigned int j=0;j<size2();j++){
	 emat(i,j) = this
	 }
	 }*/
	ublas_matrix_impl m(this->size1(), this->size2(), this->data());
	m = expm_pad(m, time_tau);
	res = math::matrix<scalar_type>(m.size1(), m.size2(), m.data());

}
/**
 * returns the exponentiated matrix U = exp(A)
 */

template<typename scalar_type>
void math::matrix<scalar_type>::matrix_exponentiation(
		math::matrix<scalar_type>& res) const {

	/*ublas_matrix_impl emat(size1(),size2());
	 for(unsinged int i=0;i<size1()i++){
	 for(unsigned int j=0;j<size2();j++){
	 emat(i,j) = this
	 }
	 }*/
	ublas_matrix_impl m(this->size1(), this->size2(), this->data());
	m = expm_pad(m);
	res = math::matrix<scalar_type>(m.size1(), m.size2(), m.data());

}

/**
 * Matrix Multiply routine with ublas matrix. The matrix passed as parameter is multiplied
 * with the caller matrix and the result is assigned to the caller matrix.
 */

template<typename scalar_type> void math::matrix<scalar_type>::multiply(
		matrix& A, matrix& res) {

	ublas_matrix_impl m1(this->size1(), this->size2(), this->data());
	ublas_matrix_impl m2(A.size1(), A.size2(), A.data());

	m1 = boost::numeric::ublas::prod(m1, m2);
	res = math::matrix<scalar_type>(m1.size1(), m1.size2(), m1.data());
}

template<typename scalar_type> void math::matrix<scalar_type>::minus(
		matrix& B, matrix& res) {

	ublas_matrix_impl m1(this->size1(), this->size2(), this->data());
	ublas_matrix_impl m2(B.size1(), B.size2(), B.data());

	m1 = m1 - m2;	//matrix subtraction
	res = math::matrix<scalar_type>(m1.size1(), m1.size2(), m1.data());
}



/**
 * Matrix vector multiplication. The result is stored in the passed argument res
 */
template<typename scalar_type> void math::matrix<scalar_type>::mult_vector(
		std::vector<scalar_type> v, std::vector<scalar_type> &res) {
	ublas_matrix_impl m(this->size1(), this->size2(), this->data());
	ublas_vector_impl uv(v.size());
	for (unsigned int i = 0; i < v.size(); i++)
		uv(i) = v[i];
	uv = prod(m, uv);
	res.resize(uv.size());
	//cout<<"\n"<<v.size()<<endl;
	for (unsigned int i = 0; i < m.size1(); i++)	//m.size1()  or uv.size()
		res[i] = uv(i);
	//cout<<"\n"<<uv.size()<<endl;
}

/**
 * Implements the transpose and assings the result to the caller matrix.

 template<typename scalar_type> void math::matrix<scalar_type>::transpose_assign() {

 size_type r = this->size2();
 size_type c = this->size1();
 ublas_matrix_impl m(r, c);
 for (size_type i = 0; i < r; i++) {
 for (size_type j = 0; j < c; j++)
 m(i, j) = this(j, i);
 }
 this(r, c, m.data());
 }
 */

/*template<typename scalar_type> void math::matrix<scalar_type>::addColumn(
 std::vector<scalar_type> columnVector,
 math::matrix<scalar_type>& resized_matrix) {

 size_type r = this->size2();
 size_type c = this->size1();
 ublas_matrix_impl m(r, c);
 for (size_type i = 0; i < r; i++) {
 for (size_type j = 0; j < c; j++)
 m(i, j) = this->at_element(j, i);
 }
 res = math::matrix<scalar_type>(m.size1(), m.size2(), m.data());
 }*/

/**
 * Implements the transpose and assings the result to the caller matrix.
 */
template<typename scalar_type> void math::matrix<scalar_type>::transpose(
		matrix& res) {

	size_type r = this->size2();
	size_type c = this->size1();
	ublas_matrix_impl m(r, c);
	for (size_type i = 0; i < r; i++) {
		for (size_type j = 0; j < c; j++)
			m(i, j) = this->at_element(j, i);
	}
	res = math::matrix<scalar_type>(m.size1(), m.size2(), m.data());
}
template<typename scalar_type> void math::matrix<scalar_type>::matrix_copy(
		math::matrix<scalar_type>& destination) {
	ublas_matrix_impl m(this->size1(), this->size2(), this->data());
	destination = math::matrix<scalar_type>(m.size1(), m.size2(), m.data());
}

/*template<typename scalar_type> void math::matrix<scalar_type>::matrix_Identity(
		math::matrix<scalar_type>& newIdentityMatrix) {
	;
}*/

template<typename scalar_type> void math::matrix<scalar_type>::matrix_join(
		matrix mat2, matrix& joined_matrix) {
	size_type row, col, index_i;
	row = this->size1();
	col = this->size2();
	if (col == mat2.size2()) {
		ublas_matrix_impl m(this->size1(), this->size2(), this->data());
		matrix mat1;
		mat1 = math::matrix<scalar_type>(m.size1(), m.size2(), m.data());
		row = row + mat2.size1(); //only row will increase as the col=system_dimension so will not change
		matrix mat2_temp(row, col);
		//cout << "This is the new Rows = " << mat2_temp.size1() << "AmitJI" << endl;
		mat1.matrix_copy(mat2_temp);
		mat2_temp.resize(row, col, true);
		//cout << "This is the new Rows after matrix_copy and resize = " << mat2_temp.size1() << "AmitJI" << endl;
		//joined_matrix.resize(row, col);
		//mat1.matrix_copy(joined_matrix);		//not working here
		for (size_type i = mat1.size1(), index_i = 0; i < mat2_temp.size1();
				i++, index_i++) {
			for (size_type j = 0; j < mat2_temp.size2(); j++) {
				mat2_temp(i, j) = mat2(index_i, j);
			}
			//cout<<endl;
		}
		/*
		 cout << "\n\n\nBEFORE JOINING"<<endl;
		 for (unsigned int i=0;i<mat2_temp.size1();i++){
		 for (unsigned int j=0;j<mat2_temp.size2();j++)
		 std::cout << mat2_temp(i,j)<<"\t";
		 cout<<endl;
		 cout<<"Check it out Ra\t";
		 }
		 */
		joined_matrix = math::matrix<scalar_type>(mat2_temp.size1(),
				mat2_temp.size2(), mat2_temp.data());

		/*
		 cout << "\n\n\nAFTER JOINING"<<endl;
		 for (unsigned int i=0;i<joined_matrix.size1();i++){
		 for (unsigned int j=0;j<joined_matrix.size2();j++)
		 std::cout << joined_matrix(i,j)<<"\t";
		 cout<<endl;
		 }
		 */

	} else {
		cout
				<< "Matices are not of the same dimension:: number of columns do not match!!!";
	}

}

/**
 * returns the infinity norm of the matrix. inf norm of a matrix m is defined as
 * inf_norm = max(a_{i,j}), 0<=i<r, 0<=j<c
 *
 */
template<typename scalar_type> scalar_type math::matrix<scalar_type>::norm_inf() {
	scalar_type norm = 0;
	for (size_type i = 0; i < this->size1(); i++) {
		for (size_type j = 0; j < this->size2(); j++) {
			if (abs(this->at_element(size_type(i), size_type(j))) > norm)
				norm = abs(this->at_element(size_type(i), size_type(j)));
			//	if (abs(this(i, j)) > norm)
			//		norm = abs(this(i, j));
		}
	}
//	std::cout<<"math norm = "<<norm<<std::endl;
	return norm;
}

/**
 * Computes the inverse of the matrix.
 *
 * \return The inverse matrix if it exists. If the matrix is singular, then the parameter singular is set to true.
 * singular is set to false otherwise. The result is put in the argument res.
 * @notice This function works doesn't work for Rational type
 */
/* Matrix inversion routine.
 Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<typename scalar_type>
bool math::matrix<scalar_type>::inverse(math::matrix<scalar_type>& inverse) {
	//bool NonSingular=false;
	assert(this->size1() == this->size2());
	using namespace boost::numeric::ublas;
	typedef permutation_matrix<std::size_t> pmatrix;
	// create a working copy of the input
	ublas_matrix_impl A(this->size1(), this->size2(), this->data());
	// make A same as the current matrix.
	/*
	 cout << "\nAmit here\n";
	 for (int i = 0; i < A.size1(); i++) {
	 for (int j = 0; j < A.size2(); j++)
	 cout << A(i, j) << "\t";
	 cout << endl;
	 }
	 cout << "\nAmit there\n";
	 */
	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());
	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0) {
		//NonSingular=false;
		return false;
	}
	// create identity matrix of "inverse"
	inverse.assign(identity_matrix<scalar_type>(A.size1()));
	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);
	return true;	//	NonSingular=true;
}

#endif
