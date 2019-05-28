#include "casmutils/definitions.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"
#include <boost/program_options.hpp>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include <iostream>

namespace Utilities
{

void symmetrize_initializer(po::options_description& symmetrize_desc)
{
    UtilityProgramOptions::add_help_suboption(symmetrize_desc);
    UtilityProgramOptions::add_output_suboption(symmetrize_desc);

    symmetrize_desc.add_options()("structure,s", po::value<fs::path>()->required(),
                               "POS.vasp like file that you want to get the symmetrized structure for.");
    symmetrize_desc.add_options()("tolerance,t", po::value<double>()->default_value(1e-5),
                               "tolerance to symmetrize to.");

    return;
}
} // namespace Utilities

using namespace Utilities;

int main(int argc, char* argv[])
{
    Handler symmetrize_launch(argc, argv, symmetrize_initializer);

    /* if(symmetrize_launch.count("help") || symmetrize_launch.argc()<2) */
    if (symmetrize_launch.count("help"))
    {
        std::cout << symmetrize_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        symmetrize_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto struc_path = symmetrize_launch.fetch<fs::path>("structure");
    auto tol = symmetrize_launch.fetch<double>("tolerance");

    // Should all CASM calls be wrapped up?
    auto struc = Rewrap::Structure(struc_path);
    auto sym_struc = symmetrize(struc,tol);

    if (symmetrize_launch.vm().count("output"))
    {
        auto out_path = symmetrize_launch.fetch<fs::path>("output");
        Simplicity::write_poscar(sym_struc, out_path);
    }

    else
    {
        Simplicity::print_poscar(sym_struc, std::cout);
    }

    return 0;
}
