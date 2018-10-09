#include "casmutils/definitions.hpp"
#include "casmutils/handlers.hpp"
#include "casmutils/stage.hpp"
#include "casmutils/structure.hpp"
#include <boost/program_options.hpp>
#include <casm/crystallography/Structure.hh>
#include <fstream>
#include <iostream>

namespace Utilities
{

void find_layer_initializer(po::options_description& find_layer_desc)
{
    UtilityProgramOptions::add_help_suboption(find_layer_desc);

    find_layer_desc.add_options()("structure,s", po::value<fs::path>()->required(),
                                  "POS.vasp like file that you want to get the layering direction for.");
    find_layer_desc.add_options()("trials,t", po::value<int>()->required(), "number of trial pivots for averaging");
    find_layer_desc.add_options()("radius,r", po::value<double>()->required(),
                                  "radius to consider as neighboring atoms");
    find_layer_desc.add_options()("layering-species,l", po::value<std::vector<std::string>>()->multitoken(),
                                  "if given, specifies which atom(s) construct layers");
    return;
}
} // namespace Utilities

using namespace Utilities;

int main(int argc, char* argv[])
{
    Handler find_layer_launch(argc, argv, find_layer_initializer);

    if (find_layer_launch.count("help"))
    {
        std::cout << find_layer_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        find_layer_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto path = find_layer_launch.fetch<fs::path>("structure");
    auto n_samples = find_layer_launch.fetch<int>("trials");
    auto nn_radius = find_layer_launch.fetch<double>("radius");
    std::vector<std::string> layer_species;
    if (find_layer_launch.count("layering-species"))
    {
        layer_species = find_layer_launch.fetch<std::vector<std::string>>("layering-species");
    }
    std::set<std::string> allowed_species(layer_species.begin(), layer_species.end());
    const auto& struc = Rewrap::Structure(path);
    std::cout << Viewpoint::find_layer_direction(struc, n_samples, nn_radius, allowed_species) << std::endl;

    return 0;
}
