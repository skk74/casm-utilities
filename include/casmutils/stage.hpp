#ifndef UTILS_STAGE_HH
#define UTILS_STAGE_HH

namespace Rewrap {

class Structure;

}

/// Symmetrizes the given structure, loosening the tolerance of symmetric equivalence
/// to tol. Returned structure is a copy with potentially altered lattice vectors
/// and basis positions. 
Rewrap::Structure symmetrize(const Rewrap::Structure& struc, double tol);

#endif
