#ifndef STAGE_HH
#define STAGE_HH
#include "casmutils/definitions.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/frankenstein.hpp"
#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/Structure.hh>
#include <iostream>

/// This function takes two structures and returns a structure that represents
/// the second structure in the frame of reference of the first structure.
/// This only considers the symmetry breaking distortions required to reach structure 2
/// Feeding the first structure (in the proper supercell) and the return value to an interpolator
/// would result in a visualization of the deformation pathway from structure 1 to structure 2.
Rewrap::Structure minimally_distorted_structure(const Rewrap::Structure& ref_struc,
                                                const Rewrap::Structure& deformed_struc, bool allow_va, double lattice_weight);

/// This function takes two structures and returns a structure that represents
/// the second structure in the frame of reference of the first structure.
/// Feeding the first structure (in the proper supercell) and the return value to an interpolator
/// would result in a visualization of the deformation pathway from structure 1 to structure 2.
Rewrap::Structure distorted_structure(const Rewrap::Structure& ref_struc, const Rewrap::Structure& deformed_struc, bool allow_va, double lattice_weight);

/// This function takes two aligned structures and returns a vector of structures
/// representing a linearly interpolated path between the input.
std::vector<Rewrap::Structure> interpolate(const Rewrap::Structure& init_struc, const Rewrap::Structure& final_struc,
                                           int n_images, double lattice_weight);

/// This function finds the deformation pathway between init_struc and final_struc
/// and then interpolates n_images along the path.
std::vector<Rewrap::Structure> deformation_pathway(const Rewrap::Structure& init_struc,
                                                   const Rewrap::Structure& final_struc, int n_images, bool minimal, bool allow_va, double lattice_weight);

/// Creates the matrix entry information for a single entry in a GUS-style matrix from host_struc and test_struc
/// the argument sym_break_only tells whether to remove the symmetry preserving part of the mapping score or not
/// Returns a tuple of mapping score, lattice score,  basis score, host point group, test point group and bool representing group-subgroup
std::tuple<double, double, double,std::string, std::string, bool> gus_entry(const Rewrap::Structure& host_struc, const Rewrap::Structure& test_struc,
                                  bool sym_break_only, double lattice_weight);

///This function searches all files in struc_folder and attempts to read them into a structure from a prim.json format
///The structure title becomes the file name with .json and any instances of POSCAR removed
std::vector<Rewrap::Structure> read_and_rename_json(const Rewrap::fs::path& struc_folder);

///This function searches all files in struc_folder and attempts to read them into a structure from a POSCAR format
///The structure title becomes the file name with any instances of POSCAR removed
std::vector<Rewrap::Structure> read_and_rename_poscar(const Rewrap::fs::path& struc_folder);

///This function replaces all occupants of a structure with the given vector of occupant names
Rewrap::Structure reassign_all_occs(const Rewrap::Structure &original, const std::set<std::string> &occ_list);

///This function returns a symmetrized version of the given structure with a relaxed tolerance
Rewrap::Structure symmetrize(const Rewrap::Structure &struc, double tol);

///This function takes a 2D hexagonal layer and returns all the equivalents (may have some repeats)
///due to operation of the lattice point group + basis shifts [1/3, 2/3, 0] and [2/3, 1/3, 0]
std::vector<Rewrap::Structure> enumerate_layer_equivalents(const Rewrap::Structure &original);

///This function takes a fcc or bcc structure and returns a reoriented structure in a cartesian 
/// frame of reference where <100> directions are along x,y, z along with all the variants of the
/// bain path distortion in the cartesian space.
std::vector<Rewrap::Structure> bainify(const Rewrap::Structure &original);

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
