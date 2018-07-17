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

void superstructure_initializer(po::options_description& superstructure_desc)
{
    UtilityProgramOptions::add_help_suboption(superstructure_desc);
    UtilityProgramOptions::add_output_suboption(superstructure_desc);

    superstructure_desc.add_options()("structure,s", po::value<fs::path>()->required(),
                               "POS.vasp like file that you want to get the super structure for.");
    superstructure_desc.add_options()("transformation matrix,t", po::value<fs::path>()->required(),
                                      "path to the file with transformation matrix.");

    return;
}
}

using namespace Utilities;

int main(int argc, char* argv[])
{
    Handler superstructure_launch(argc, argv, superstructure_initializer);

    /* if(primify_launch.count("help") || primify_launch.argc()<2) */
    if (superstructure_launch.count("help"))
    {
        std::cout << superstructure_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        superstructure_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto structure_path = superstructure_launch.fetch<fs::path>("structure");
    auto transf_file_path = superstructure_launch.fetch<fs::path>("transformation matrix");

    //change this to Rewrap structure after john pushes a branch with constructor that takes path
    auto struc = CASM::Structure(structure_path);
    auto super_struc = Simplicity::make_super_struc(struc,transf_file_path);

    if (superstructure_launch.vm().count("output"))
    {
        auto out_path = superstructure_launch.fetch<fs::path>("output");
        Simplicity::write_poscar(super_struc, out_path);
    }

    else
    {
        Simplicity::print_poscar(super_struc, std::cout);
    }

    return 0;
}
