#include <boost/filesystem.hpp>
#include <casm/CASM_global_definitions.hh>
#include "casmutils/exceptions.hpp"
#include "casmutils/structure.hpp"
#include <fstream>
#include <set>
namespace Viewpoint {
struct NeighborInfo {
	int index;
	double distance;
};

Eigen::Vector3d find_layer_direction(const Rewrap::Structure& struc,
				     int num_samples, double nn_bubble_radius,
				     std::set<std::string>& layering_species) {
	Eigen::Vector3d total_hkl(0.0, 0.0, 0.0);
	// repeat finding normal direction from random sample points
	for (size_t i = 0; i < num_samples; i += 1) {
		int rand_idx = rand() % (struc.basis.size() - 1);
		auto& pivot = struc.basis[rand_idx];
		// first find all neighbors in range
		std::vector<NeighborInfo> neighbors;
		for (size_t j = 0; j < struc.basis.size(); j++) {
			const auto& neighbor = struc.basis[j];
			double dist =
			    (neighbor.const_cart() - pivot.const_cart()).norm();
			auto end = layering_species.end();
			if (0 < dist && dist < nn_bubble_radius &&
			    (!layering_species.size() ||
			     (layering_species.find(pivot.occ_name()) != end && layering_species.find(neighbor.occ_name()) != end))) {
				neighbors.push_back(NeighborInfo{j, dist});
			}
		}
		// sort the neighbors
		std::sort(neighbors.begin(), neighbors.end(),
			  [](const NeighborInfo& a, const NeighborInfo& b) {
				  return a.distance < b.distance;
			  });
		// cut the number of neighbors
		if (neighbors.size() < 2) {
			continue;
		}
		// get normals from the sample point and two neighbors
		std::vector<Eigen::Vector3d> normals;
		for (size_t j = 0; j < neighbors.size() - 1; j++) {
			const auto& neighbor0 = struc.basis[neighbors[j].index];
			const auto& neighbor1 =
			    struc.basis[neighbors[j + 1].index];
			Eigen::Vector3d vec0 =
			    neighbor0.const_frac() - pivot.const_frac();
			Eigen::Vector3d vec1 =
			    neighbor1.const_frac() - pivot.const_frac();
			auto n = (vec0.cross(vec1)).normalized();
			if (n(2) < 0.0) {
				n = -1.0 * n;
			} else if (n(2) == 0.0) {
				if (n(1) < 0.0) {
					n = -1.0 * n;
				}
			}
			Eigen::Vector3d vec2 =
			    neighbor0.const_cart() - pivot.const_cart();
			Eigen::Vector3d vec3 =
			    neighbor1.const_cart() - pivot.const_cart();
			if (std::abs(vec2.dot(vec3)/vec2.norm()/vec3.norm()) >0.99){
				//This sample's test vectors are essentially parallel. Not using because it will ruin average.
				continue;
			}
			normals.push_back(n);
		
		}
		// average the normals
		Eigen::Vector3d averaged(0.0, 0.0, 0.0);
		for (const auto& n : normals) {
			averaged += n;
		}
		averaged /= double(normals.size());
		averaged.normalize();

		total_hkl += averaged;
	}

	// average the normals from each sample points
	total_hkl /= double(num_samples);
	total_hkl.normalize();
	// This is the layering plane in terms of h k l
	Eigen::Vector3d ret_val= total_hkl;
	// This is the layering plane in terms of u v w
	ret_val= struc.lattice().inv_lat_column_mat()*struc.lattice().reciprocal().lat_column_mat()*ret_val;
	return ret_val;
}
}
