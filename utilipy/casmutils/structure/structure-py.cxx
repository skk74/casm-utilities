#include <boost/filesystem.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"
#include <casm/crystallography/Structure.hh>
#include <string>

//******************************************************************************************************//
//******************************************************************************************************//

namespace WrapPy
{
Rewrap::Structure from_poscar(const std::string& filename);
void to_poscar(const Rewrap::Structure& writeable, const std::string& filename);

//******************************************************************************************************//

std::string __str__(const Rewrap::Structure& printable)
{
    std::ostringstream sstream;
    Simplicity::print_poscar(printable, sstream);
    return sstream.str();
}

Rewrap::Structure from_poscar(const std::string& filename) { return Rewrap::Structure::from_poscar(filename); }

void to_poscar(const Rewrap::Structure& writeable, const std::string& filename)
{
    Simplicity::write_poscar(writeable, filename);
    return;
}

int num_atoms(const Rewrap::Structure& struc)
{
	return struc.basis.size();
}

double c_vec_length(const Rewrap::Structure& struc)
{
	return struc.lattice().length(2);
}

PYBIND11_MODULE(_structure, m)
{
    using namespace pybind11;
    using namespace WrapPy;

    m.doc() = "Raw python bindings for a re-wrapped CASM::Structure class.";

    class_<Rewrap::Structure>(m, "Structure")
        .def("is_primitive", &Rewrap::Structure::is_primitive)
        .def("primitive", &Rewrap::Structure::primitive)
        /* .def("__repr__", &Structure::to_string) */
        .def("__str__", &__str__);

    m.def("num_atoms", num_atoms);
    m.def("c_vec_length",c_vec_length);
    m.def("from_poscar", from_poscar);
    m.def("to_poscar", to_poscar);
    m.def("structure_map", gus_entry);
    m.def("symmetrize", symmetrize);
    /* m.def("make_niggli", (Rewrap::Structure(*)(const Rewrap::Structure&)) Simplicity::make_niggli); */
}
} // namespace WrapPy
