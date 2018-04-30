#include <iostream>
#include <limits>
#include "time.h"
#include "boost/random.hpp"
#include "uni_sphere.h"
#include<cmath>

namespace math {

using namespace std;
#define INF (std::numeric_limits<double>::infinity())

long double dist_pt_pt(const vector<double> &v1, const vector<double> &v2) {

	long double res = 0.;

	for (int i = 0; i < v1.size(); i++)
		res += (v1[i] - v2[i]) * (v1[i] - v2[i]);

	return std::sqrt(res);
}

double dist_pt_vec(const vector<double> &v,
		const std::vector<vector<double> > &vv, int &idx) {

	int i = 0;
	idx = 0;
	double d;
	double dmin = INF;

	for (i = 0; i < vv.size(); i++) {
		d = dist_pt_pt(v, vv[i]);
		if ((d < dmin) && (d != 0)) {
			dmin = d;
			idx = i;
		}
	}
	return dmin;
}

// look for the two closest points, return their indices and the distance to the second closest
double dist_pt_vec2(const vector<double> &v,
		const std::vector<vector<double> > &vv, int &index1, int & index2) {

	int i = 0; //, imin=0;
	double d;
	double dmin1 = INF;
	double dmin2 = INF;

	for (i = 0; i < vv.size(); i++) {
		d = dist_pt_pt(v, vv[i]);
		if ((d < dmin2)) {
			if (d < dmin1) {
				dmin2 = dmin1;
				index2 = index1;
				dmin1 = d;
				index1 = i;
			} else {
				dmin2 = d;
				index2 = i;
			}
		}
	}
	return dmin2;
}

std::vector<vector<double> > uni_sphere(int N, int dim, int max_iter,
		double tol) {

	std::vector<vector<double> > vout(N);
	vector<double> distv(N);
	int j = 0;

	/* initialize random generator : */

	boost::mt19937 rng;                 // produces randomness out of thin air
										// see pseudo-random number generators
	// default constructor of random generator always uses the same seed, see
	// http://www.boost.org/doc/libs/1_43_0/doc/html/boost/random/mersenne_twister.html

	boost::uniform_on_sphere<> rsphere(dim);      // distribution 

	boost::variate_generator<boost::mt19937&, boost::uniform_on_sphere<> > give_sphere_pt(
			rng, rsphere);             // glues randomness with mapping

	/* generate N points on the sphere */

	// cout << "Generate N points on the sphere"<< endl; // DEBUG__
	for (j = 0; j < N; j++) {

		vector<double> new_point = vector<double>(give_sphere_pt());
		vout.at(j) = new_point;

	}

	/* initialize the minimum distances between the points */

	vector<int> idistv(N); // indices of nearest neighbors
	int idx;
	double cur_min = INF;
	double prev_min = INF;
	double cur_max = 0;

	for (int i = 0; i < N; i++) {
		distv[i] = dist_pt_vec(vout[i], vout, idx);
		idistv[i] = idx;

		if (distv[i] < cur_min)
			cur_min = distv[i];
		if (distv[i] > cur_max)
			cur_max = distv[i];

	}

	/* start maximizing minimum distances */

	//std::cout << "\n\n Starting optimization \n----------------\n";
	//std::cout << "tol: "<< tol << endl; // DEBUG__
	int idx1, idx2;
	double d2;
	int nb_iter = 0;
	double err1 = 1., err2 = 1., curr_err = 1., simp = 1;

	while (1) {

		vector<double> candidate = vector<double>(give_sphere_pt());
		d2 = dist_pt_vec2(candidate, vout, idx1, idx2);

		if (d2 > distv[idx1]) { // we found a better candidate

			/* Replace the point and update the  minimum distances vector */

			// insert new point 
			vout[idx1] = candidate;
			distv[idx1] = d2;
			idistv[idx1] = idx2;

			// update distances for points that need to be updated 

			prev_min = cur_min;
			cur_min = INF;
			cur_max = 0;
			for (int i = 0; i < N; i++) {
				if (idistv[i] == idx1) {
					distv[i] = dist_pt_vec(vout[i], vout, idx);
					idistv[i] = idx;
				} else {
					d2 = dist_pt_pt(vout[i], candidate);
					if ((d2 < distv[i]) && (d2 != 0)) {
						distv[i] = d2;
						idistv[i] = idx1;
					}
				}

				if (distv[i] < cur_min)
					cur_min = distv[i];

				if (distv[i] > cur_max)
					cur_max = distv[i];
			}

			if (cur_min > prev_min) {
				err1 = err2;
				err2 = curr_err;
				curr_err = cur_min - prev_min;
				simp = (err1 + err2 + curr_err) / (3. * cur_min);
				//	cout << "cur_min: "<< cur_min << " cur_max: " << cur_max << " cur_max-cur_min: " << cur_max - cur_min << "smoothed improvt:" << simp << endl; // DEBUG__
			}

		}

		if ((nb_iter++ > max_iter) || (simp < tol)) {
//      std::cout << "cur_min: "<< cur_min << " cur_max: " << cur_max << " cur_max-cur_min: " << cur_max - cur_min << " smoothed improvt:" << simp << std::endl; // DEBUG__
			break;
		}

	} // end while(1)

	//std::cout << "\ndone\n"<< endl; // DEBUG__
	//std::cout << "distv: "<< distv<< endl; // DEBUG__

	return vout;

}

}

