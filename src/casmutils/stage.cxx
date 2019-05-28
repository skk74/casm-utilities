#include "include/casmutils/structure.hpp"
#include "include/casmutils/stage.hpp"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Structure.hh"



Rewrap::Structure symmetrize(const Rewrap::Structure& struc, double tol)
{
    std::set<int> point_group_sizes {1,2,3,4,6,8,12,16,24,48};
    Rewrap::Structure sym_struc=struc;
    int biggest = struc.factor_group().size();
    CASM::Structure tmp = struc;
    // a) symmetrize the lattice vectors
    CASM::Lattice lat = tmp.lattice();
    lat.symmetrize(tol);
    //lat.set_tol(tol);
    CASM::SymGroup pg;
    lat.generate_point_group(pg);
    if (point_group_sizes.find(pg.size()) == point_group_sizes.end())
    {
        return struc;
    }
    tmp.set_lattice(lat, CASM::FRAC);

    tmp.factor_group();
    // b) find factor group with same tolerance
    tmp.fg_converge(tol);
    if (point_group_sizes.find(tmp.point_group().size()) == point_group_sizes.end())
    {
        return struc;
    }

    // c) symmetrize the basis sites
    tmp.symmetrize(tmp.factor_group());
    if (tmp.factor_group().is_group(tol) && (tmp.factor_group().size() > biggest))
    {
        sym_struc=tmp;
        biggest = tmp.factor_group().size();
    }
    return sym_struc;
}

