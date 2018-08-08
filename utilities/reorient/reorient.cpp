#include "casmutils/definitions.hpp"
#include "casmutils/frankenstein.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/structure.hpp"
#include "casmutils/stage.hpp"
#include <boost/program_options.hpp>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include <iostream>

namespace Utilities
{

void reorient_initializer(po::options_description& reorient_desc)
{
    UtilityProgramOptions::add_help_suboption(reorient_desc);
    UtilityProgramOptions::add_desc_suboption(reorient_desc);
    UtilityProgramOptions::add_output_suboption(reorient_desc);
    reorient_desc.add_options()("structure,s", po::value<fs::path>()->required(),
                                   "POS.vasp like file you want to "
                                   "reorient");
    reorient_desc.add_options()("miller-indices,m", po::value<std::vector<int>>()->multitoken()->required(),
                                   "The miller indices of the direction that is perpendicular to the "
                                   "planes you wish to reveal.");
    return;
}
} // namespace Utilities

using namespace Utilities;

int main(int argc, char* argv[])
{
    Handler reorient_launch(argc, argv, reorient_initializer);

    if (reorient_launch.count("help"))
    {
        std::cout << reorient_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        reorient_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto struc_path = reorient_launch.fetch<fs::path>("structure");
    auto struc = Rewrap::Structure(struc_path);
    auto vec = reorient_launch.fetch<std::vector<int>>("miller-indices");
    auto out_struc = Viewpoint::reoriented_struc(struc, Eigen::Map<Eigen::Vector3i>(&vec[0]));
    if (reorient_launch.vm().count("output"))
    {
        Simplicity::write_poscar(out_struc, reorient_launch.fetch<fs::path>("output"));
        return 0;
    }
    Simplicity::print_poscar(out_struc, std::cout);
    return 0;
}
