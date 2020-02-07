#include <casm/CASM_global_definitions.hh>
#include <casm/crystallography/io/VaspIO.hh>
#include <casm/crystallography/Niggli.hh>
#include <casm/crystallography/Structure.hh>
#include <casm/crystallography/SimpleStructureTools.hh>
#include <boost/filesystem.hpp>
#include "casmutils/structure.hpp"
#include <fstream>

namespace Rewrap
{
Structure::Structure(CASM::xtal::Structure init_struc) : CASM::xtal::Structure(init_struc) {}

bool Structure::is_primitive() const { return CASM::xtal::Structure::is_primitive(); }

Structure Structure::primitive() const { return Simplicity::make_primitive(*this); }
}

namespace Simplicity
{
Rewrap::Structure make_primitive(const Rewrap::Structure& input)
{
    const CASM::xtal::Structure& casted_input(input);
    CASM::xtal::Structure true_prim;
    bool is_prim = casted_input.is_primitive(true_prim);
    return true_prim;
}

Rewrap::Structure make_niggli(const Rewrap::Structure& non_niggli)
{
    CASM::xtal::Structure niggli = non_niggli;
    CASM::xtal::Lattice lat_niggli = CASM::xtal::niggli(non_niggli.lattice(), CASM::TOL);
    niggli.set_lattice(lat_niggli, CASM::CART);
    return niggli;
}

void make_niggli(Rewrap::Structure* non_niggli)
{
    CASM::xtal::Lattice lat_niggli = CASM::xtal::niggli(non_niggli->lattice(), CASM::TOL);
    non_niggli->set_lattice(lat_niggli, CASM::CART);
    return;
}

void print_poscar(const Rewrap::Structure& printable, std::ostream& outstream)
{
    CASM::VaspIO::PrintPOSCAR p(CASM::xtal::make_simple_structure(printable));
    p.sort();
    p.print(outstream);
    return;
}

void write_poscar(const Rewrap::Structure& printable, const Rewrap::fs::path& filename)
{
    std::ofstream file_out(filename.string());
    print_poscar(printable, file_out);
    file_out.close();
    return;
}
}
