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
Rewrap::Structure minimally_distorted_structure(const Rewrap::Structure& ref_struc,
                                                const Rewrap::Structure& deformed_struc);

/// This function takes two structures and returns a structure that represents
/// the second structure in the frame of reference of the first structure.
/// Feeding the first structure (in the proper supercell) and the return value to an interpolator
/// would result in a visualization of the deformation pathway from structure 1 to structure 2.
Rewrap::Structure distorted_structure(const Rewrap::Structure& ref_struc, const Rewrap::Structure& deformed_struc);

/// This function takes two aligned structures and returns a vector of structures
/// representing a linearly interpolated path between the input.
std::vector<Rewrap::Structure> interpolate(const Rewrap::Structure& init_struc, const Rewrap::Structure& final_struc,
                                           int n_images);

/// This function finds the deformation pathway between init_struc and final_struc
/// and then interpolates n_images along the path.
std::vector<Rewrap::Structure> deformation_pathway(const Rewrap::Structure& init_struc,
                                                   const Rewrap::Structure& final_struc, int n_images, bool minimal);

/// Creates the matrix entry information for a single entry in a GUS-style matrix from host_struc and test_struc
/// the argument sym_break_only tells whether to remove the symmetry preserving part of the mapping score or not
/// Returns a pair of mapping score and bool representing group-subgroup
std::pair<double, bool> gus_entry(const Rewrap::Structure& host_struc, const Rewrap::Structure& test_struc,
                                  bool sym_break_only);

///This function searches all files in struc_folder and attempts to read them into a structure from a prim.json format
///The structure title becomes the file name with .json and any instances of POSCAR removed
std::vector<Rewrap::Structure> read_and_rename_json(const Rewrap::fs::path& struc_folder);

///This function searches all files in struc_folder and attempts to read them into a structure from a POSCAR format
///The structure title becomes the file name with any instances of POSCAR removed
std::vector<Rewrap::Structure> read_and_rename_poscar(const Rewrap::fs::path& struc_folder);

///This function replaces all occupants of a structure with the given vector of occupant names
Rewrap::Structure reassign_all_occs(const Rewrap::Structure &original, const std::set<std::string> &occ_list);

#endif
