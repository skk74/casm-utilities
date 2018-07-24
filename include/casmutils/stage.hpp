#ifndef STAGE_HH
#define STAGE_HH

#include "casmutils/definitions.hpp"
#include "casmutils/structure.hpp"
#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>
#include <iostream>

/// This function takes two structures and returns a structure that represents
/// the second structure in the frame of reference of the first structure.
/// This only considers the symmetry breaking distortions required to reach structure 2
/// Feeding the first structure (in the proper supercell) and the return value to an interpolator
/// would result in a visualization of the deformation pathway from structure 1 to structure 2.
Rewrap::Structure minimally_distorted_structure(Rewrap::Structure& ref_struc, Rewrap::Structure& deformed_struc);

/// This function takes two structures and returns a structure that represents
/// the second structure in the frame of reference of the first structure.
/// Feeding the first structure (in the proper supercell) and the return value to an interpolator
/// would result in a visualization of the deformation pathway from structure 1 to structure 2.
Rewrap::Structure distorted_structure(Rewrap::Structure& ref_struc, Rewrap::Structure& deformed_struc);

/// This function takes two aligned structures and returns a vector of structures
/// representing a linearly interpolated path between the input.
std::vector<Rewrap::Structure> interpolate(Rewrap::Structure& init_struc, Rewrap::Structure& final_struc, int n_images);

/// This function finds the deformation pathway between init_struc and final_struc
/// and then interpolates n_images along the path.
std::vector<Rewrap::Structure> deformation_pathway(Rewrap::Structure& init_struc, Rewrap::Structure& final_struc,
                                                   int n_images, bool minimal);

#endif
