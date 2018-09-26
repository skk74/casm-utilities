#include "casmutils/definitions.hpp"
#include "casmutils/stage.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/structure.hpp"
#include <boost/program_options.hpp>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include <iostream>

namespace Utilities
{

void deformation_pathway_initializer(po::options_description& deformation_pathway_desc)
{
    UtilityProgramOptions::add_help_suboption(deformation_pathway_desc);
    UtilityProgramOptions::add_desc_suboption(deformation_pathway_desc);
    UtilityProgramOptions::add_output_suboption(deformation_pathway_desc);
    deformation_pathway_desc.add_options()("structures,s", po::value<std::vector<fs::path>>()->multitoken()->required(),
                                           "POS.vasp like files you want to "
                                           "get deformation pathway for");
    deformation_pathway_desc.add_options()("number-images,n", po::value<int>()->required(),
                                           "number of structures you want along the path");
    deformation_pathway_desc.add_options()(
        "minimal,m", "Give this flag if you only care about the symmetry breaking portion of the deformation path");
    deformation_pathway_desc.add_options()(
        "allow-vacancies,v", "Give this flag if you want vacancies to be considered during the map");
    deformation_pathway_desc.add_options()("weight,w", po::value<double>()->required(),
                                           "lattice-basis tradeoff that determines the shortest path for the best mapping");
    return;
}
} // namespace Utilities

using namespace Utilities;

int main(int argc, char* argv[])
{
    Handler deformation_pathway_launch(argc, argv, deformation_pathway_initializer);

    if (deformation_pathway_launch.count("help"))
    {
        std::cout << deformation_pathway_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        deformation_pathway_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto structure_paths = deformation_pathway_launch.fetch<std::vector<fs::path>>("structures");
    auto lattice_weight = deformation_pathway_launch.fetch<double>("weight");

    if (structure_paths.size() != 2)
    {
        std::cerr << "You must give exactly 2 --structures for deformation pathway, otherwise reference is ambiguous"
                  << std::endl;
        return 3;
    }

    auto init_struc = Rewrap::Structure(structure_paths[0]);
    auto final_struc = Rewrap::Structure(structure_paths[1]);
    auto strucs = deformation_pathway(init_struc, final_struc, deformation_pathway_launch.fetch<int>("number-images"),
                                      deformation_pathway_launch.count("minimal"),deformation_pathway_launch.count("allow-vacancies"),lattice_weight);
    if (deformation_pathway_launch.vm().count("output"))
    {
	auto out_path = deformation_pathway_launch.fetch<fs::path>("output");
        int count = 0;
        for (auto& item : strucs)
        {
            // TODO: what if directory doesn't exist?
            std::ostringstream ostr;
            ostr << std::setfill('0') << std::setw(2) << count;
            Simplicity::write_poscar(item, out_path / Rewrap::fs::path("image" + ostr.str() + "POSCAR"));
            count++;
        }
        return 0;
    }
    int count = 0;
        for (auto& item : strucs)
        {
            std::cout << "image " << count << std::endl;
            Simplicity::print_poscar(item, std::cout);
            count++;
        }
    return 0;
}
