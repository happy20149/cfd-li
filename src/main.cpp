#include <iostream>
#include <string>
#include "Communication.hpp"
#include "Simulation.hpp"

int main(int argn, char **args) {
    Params params;
    Communication::init_mpi(&params);
    std::string file_name;
    
    if (argn > 1) {
        file_name = args[1];  // Use the provided file name
    } else {
        // Default file path when no argument is provided
        file_name = "D:\\cfd-li\\example_cases\\StepFlowTurb\\StepFlowTurb.dat";
        std::cout << "No input file provided. Using default file: " << file_name << std::endl;
    }

    Simulation problem(file_name, argn, args, params);
    problem.simulate(params);
    Communication::finalize();
}
