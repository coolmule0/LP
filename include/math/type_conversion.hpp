/*
 * type_conversion.hpp
 *
 *  Created on: Nov 12, 2009
 *      Author: frehse
 */

#ifndef TYPE_CONVERSION_HPP_
#define TYPE_CONVERSION_HPP_

#include <stdexcept>
#include <typeinfo>

#include "type_conversion.h"

template<typename result_type, typename scalar_type>
result_type converter<result_type, scalar_type>::convert(const scalar_type& x) {
	std::string name1 = typeid(x).name();
	result_type rtmp;
	std::string name2 = typeid(rtmp).name();
	throw std::runtime_error(
			"convert_element from " + name1 + " to " + name2
					+ " not implemented.");
	return result_type(); // this is a dummy return to avoid compiler warnings
}

/** Specialization for __float128 because typeinfo is missing (g++ bug).
 *
 * @see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=43622
 */
template<>
class converter<__float128, __float128> {
public:
	static __float128 convert(const __float128& x) {
		return x;
	}
	;
};

/** Specialization for __float128 because typeinfo is missing (g++ bug).
 *
 * @see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=43622
 */
template<typename result_type>
class converter<result_type, __float128> {
public:
	static result_type convert(const __float128& x) {
		result_type rtmp;
		std::string name2 = typeid(rtmp).name();
		throw std::runtime_error(
				"convert_element from __float128 to " + name2
						+ " not implemented.");
		return result_type(); // this is a dummy return to avoid compiler warnings
	}
};

/** Specialization for __float128 because typeinfo is missing (g++ bug).
 *
 * @see http://gcc.gnu.org/bugzilla/show_bug.cgi?id=43622
 */
template<typename scalar_type>
class converter<__float128, scalar_type> {
public:
	static __float128 convert(const scalar_type& x) {
		std::string name1 = typeid(x).name();
		throw std::runtime_error(
				"convert_element from " + name1
						+ " to __float128 not implemented.");
	return __float128(); // this is a dummy return to avoid compiler warnings
}
};

/** Convert double to long double */
template<>
class converter<long double, double> {
public:
	static long double convert(const double& x) {
		long double y(x);
		return y;
	}
	;
};

/** Convert long double to double */
template<>
class converter<double, long double> {
public:
	static double convert(const long double& x) {
		double y(x);
		return y;
	}
	;
};

/** Convert double to __float128 */
template<>
class converter<__float128, double> {
public:
	static __float128 convert(const double& x) {
		__float128 y(x);
		return y;
	}
	;
};

/** Convert __float128 to double */
template<>
class converter<double, __float128> {
public:
	static double convert(const __float128& x) {
		double y(x);
		return y;
	}
	;
};

/** Convert double to float */
template<>
class converter<float, double> {
public:
	static double convert(const double& x) {
		double y(x);
		return y;
	}
	;
};

#endif /* TYPE_CONVERSION_HPP_ */
