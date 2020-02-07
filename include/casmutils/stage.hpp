#include <tuple>
#include <unordered_set>
#include <casm/external/Eigen/Core>
namespace CASM {
	namespace xtal{
		class SimpleStructure;
	}

}

///\brief This function takes a test structure and maps it onto a parent
///structure \param _parent \parblock 	The parent crystal structure used as a
///reference \endparblock \param allowed_occupants \parblock
///  The allowed occupants on every site of the parent crystal
///\endparblock
///\param _test
///\parblock
///  The structure being mapped
///\endparblock
///\param max_va_concentration
///\parblock
///  The maximum allowed vacancy concentration during the assignment of basis
///  sites
///\endparblock
///\param lattice_weight
///\parblock
///  The tradeoff between the importance of lattice score vs. basis score in
///  mapping optimization 0 is all basis 1 is all lattice.
///\endparblock
/// returns A tuple of 
//		The unstrained superlattice of _parent as a 3x3 matrix of column vectors
//		The Left stretch tensor that takes this superlattice to test's lattice
//		The displacement that takes the supercell's basis to the basis of test.
std::tuple<Eigen::Matrix3d, Eigen::Matrix3d, Eigen::MatrixXd> structure_map(
    CASM::xtal::SimpleStructure &_parent,
    std::unordered_set<std::string> &allowed_occupants, CASM::xtal::SimpleStructure &_test,
    double max_va_concentration, double lattice_weight);
