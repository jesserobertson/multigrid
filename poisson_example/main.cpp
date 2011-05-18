/*
    main.cpp (poisson)
    Jess Robertson, 2011-01-28
*/ 

#include <multigrid/multigrid.hpp>
#include <boost/program_options.hpp>     
#include "poisson.hpp"                  

using namespace mgrid;   
using namespace std;
namespace bpo = boost::program_options; 

int main (int argc, char *argv[]) {
    try{
        // Declare some option variables
        vector<double> aspectRatios;
        
        // Set up command-line options
        bpo::options_description 
            visibleOptions("Usage ./poisson <aspect ratios>\n\nOptions:");
        visibleOptions.add_options() \
            ("help", "prints this help message");    
        bpo::options_description hiddenOptions("Hidden options");
        hiddenOptions.add_options()("aspects", \
            bpo::value< vector<double> >(&aspectRatios), "aspect ratios"); 
        bpo::positional_options_description positional;
        positional.add("aspects", -1);  
        
        // Compose options for command line 
        bpo::options_description cmdline_options;
        cmdline_options.add(visibleOptions).add(hiddenOptions);
        
        // Parse command line variables, pass to variable map  
        bpo::variables_map varMap;
        bpo::store(
            bpo::command_line_parser(argc, argv).\
                options(cmdline_options).positional(positional).run(),
            varMap);
        bpo::notify(varMap);
        
        // Check for help flag
        if (varMap.count("help")) {
            std::cout << visibleOptions << std::endl;
            return 1;
        }
        
        // Run through given aspect ratios
        std::cout << "Running..." << std::endl;

        // Solve problems   
        foreach(double aspect, aspectRatios) {  
            Settings settings;
            settings.aspectRatio = aspect;
            std::auto_ptr<Poisson> problem(new Poisson(settings));
            problem->solve();
            problem->write(1, "poisson_");
        }    
        std::cout << "Finished!" << std::endl;
        return 0;
        
    } catch (std::exception& e) {
		std::cout << e.what() << std::endl;
        return 1;
    }      
}