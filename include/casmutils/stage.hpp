#ifndef UTILS_STAGE_HH
#define UTILS_STAGE_HH

#include "casmutils/definitions.hpp"
#include "casmutils/structure.hpp"
#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>
#include <iostream>

namespace Viewpoint
{

/// This function finds the most likely layering direction in the given
/// Structure. Layering direction is given in terms of crystal axes.
Eigen::Vector3d find_layer_direction(const Rewrap::Structure& struc, int n_samples, double nn_bubble_radius,
                                     std::set<std::string>& layering_species);

/// This function reorients a structure along the plane perpendicular to the miller direction given. 
/// the new a and b vectors lie in the plane.
Rewrap::Structure reoriented_struc(const Rewrap::Structure& struc, Eigen::Vector3i millers);

} // namespace Viewpoint

#endif
