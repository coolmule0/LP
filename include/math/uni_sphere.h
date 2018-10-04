#ifndef _UNI_SPHERE
#define _UNI_SPHERE
#include <vector>

namespace math {

/**
 *
 * Compute N points uniformly distributed on a hyper-sphere of dimension dim
 *
 * The algorithm works as follows:
 * firstly it generates N points with a uniform random distribution on the sphere
 * secondly it improves the repartition by iteratively picking new points,
 * exchanging them with their nearest neighbors if this increases the minimum distances
 * between the points.
 *
 * @param N number of points
 * @param dim number of dimensions
 * @param max_iter max number of random picking in the optimization phase
 * @param tol stopping value on the increase of the minimum distance between any two points
 *        (averaged on three iterations)
 *
 * @author Alexandre Donze
 */
std::vector<std::vector<double> > uni_sphere(int N, int dim, int max_iter,
		double tol);

}

#endif
