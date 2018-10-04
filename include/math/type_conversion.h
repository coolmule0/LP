#ifndef TYPE_CONVERSION_H_
#define TYPE_CONVERSION_H_

//#include "math/scalar_types/rational.h"
//#include "../utility/calc_string.h"

/** Helper class to allow partial specialization.
 * The default conversion is to throw (otherwise we'd have to implement conversion
 * for all kinds of combinations that are never used.) */
template<typename result_type, typename scalar_type>
class converter {
public:
	static result_type convert(const scalar_type& x);
};

/** This is an interface that calls the converter class.
 */
template<typename result_type, typename scalar_type> result_type convert_element(
		const scalar_type& x) {
	return converter<result_type, scalar_type>::convert(x);
}
;

/** Don't convert if scalar_type is the same as result_type. */
template<typename result_type>
class converter<result_type, result_type> {
public:
	static result_type convert(const result_type& x) {
		return x;
	}
	;
};

#include "type_conversion.hpp"

#endif /*TYPE_CONVERSION_H_*/
