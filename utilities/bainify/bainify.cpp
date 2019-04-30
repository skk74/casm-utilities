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

void bainify_initializer(po::options_description& bainify_desc)
{
    UtilityProgramOptions::add_help_suboption(bainify_desc);
    UtilityProgramOptions::add_desc_suboption(bainify_desc);
    UtilityProgramOptions::add_output_suboption(bainify_desc);
    bainify_desc.add_options()("structure,s", po::value<fs::path>()->required(),
                                    "POS.vasp like file that you want to "
                                    "get the bain equivalents for.");

    return;
}
} // namespace Utilities

using namespace Utilities;

int main(int argc, char* argv[])
{
    double tol = 1e-5;
    Handler bainify_launch(argc, argv, bainify_initializer);

    if (bainify_launch.count("help"))
    {
        std::cout << bainify_launch.desc() << std::endl;
        return 1;
    }

    try
    {
        bainify_launch.notify();
    }

    catch (po::required_option& e)
    {
        std::cerr << e.what() << std::endl;
        return 2;
    }

    auto path = bainify_launch.fetch<fs::path>("structure");

    auto struc = Rewrap::Structure(path);
    std::vector<Rewrap::Structure> equivs;
	equivs = bainify(struc);
    
    if (bainify_launch.vm().count("output"))
    {
        auto out_path = bainify_launch.fetch<fs::path>("output");
        int count = 0;
        for (auto& item : equivs)
        {
			std::string title = "";
			switch (count)
			{
				case 0:
					title="OriginalPOSCAR";
					break;
				case 1:
					title="XBainPOSCAR";
					break;
				case 2:
					title="YBainPOSCAR";
					break;
				case 3:
					title="ZBainPOSCAR";
					break;
				default:
					std::cerr<< "too many equivs" << std::endl;	
					return 3;
			
			}
            // TODO: what if directory doesn't exist?
            Simplicity::write_poscar(item, out_path / Rewrap::fs::path(title));
            count++;
        }
    }

    else
    {
        int count = 0;
        for (auto& item : equivs)
        {
            std::string title = "";
			switch (count)
			{
				case 0:
					title="OriginalPOSCAR";
					break;
				case 1:
					title="XBainPOSCAR";
					break;
				case 2:
					title="YBainPOSCAR";
					break;
				case 3:
					title="ZBainPOSCAR";
					break;
				default:
					std::cerr<< "too many equivs" << std::endl;	
					return 3;
			
			}
			std::cout << title << std::endl;
			Simplicity::print_poscar(item, std::cout);
            count++;
        }
    }

    return 0;
}
