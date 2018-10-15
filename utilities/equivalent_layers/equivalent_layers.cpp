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

void equivalent_layers_initializer(po::options_description& equivalent_layers_desc)
{
    UtilityProgramOptions::add_help_suboption(equivalent_layers_desc);
    UtilityProgramOptions::add_desc_suboption(equivalent_layers_desc);
    UtilityProgramOptions::add_output_suboption(equivalent_layers_desc);
    equivalent_layers_desc.add_options()("structure,s", po::value<fs::path>()->required(),
                                   "POS.vasp like file you want to find "
                                   "equivalent layers for");
    return;
}
} // namespace Utilities

using namespace Utilities;

int main(int argc, char* argv[])
{
    Handler equivalent_layers_launch(argc, argv, equivalent_layers_initializer);

    if (equivalent_layers_launch.count("help"))
    {
        std::cout << equivalent_layers_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        equivalent_layers_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto struc_path = equivalent_layers_launch.fetch<fs::path>("structure");
    auto struc = Rewrap::Structure(struc_path);
	std::vector<Rewrap::Structure> equivs = enumerate_layer_equivalents(struc); 
		int count = 0;
		for (auto &out_struc : equivs){
    if (equivalent_layers_launch.vm().count("output"))
    {   
        Simplicity::write_poscar(out_struc, equivalent_layers_launch.fetch<fs::path>("output") / fs::path( "equiv_" + std::to_string(count)));
		}
	else{
		std::cout << count << std::endl;
    Simplicity::print_poscar(out_struc, std::cout);
	}
		count++;
    }
    return 0;
}
