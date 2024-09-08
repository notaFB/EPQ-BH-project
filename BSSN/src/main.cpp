#include <iostream>
#include "grid.hpp"
#include "fieldData.hpp"
#include "initial_conditions.hpp"
#include "evolution.hpp"
#include "boundary_conditions.hpp"
#include "utilities.hpp"
#include "config.hpp"

#include <iostream>
#include <chrono>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>  // Include this header for std::istringstream

std::size_t getMemoryUsage() {
    std::ifstream file("/proc/self/status");
    std::string line;
    std::size_t memory = 0;
    while (std::getline(file, line)) {
        if (line.substr(0, 6) == "VmRSS:") {  // Resident Set Size, memory in RAM (in KB)
            std::istringstream iss(line);
            std::string key;
            iss >> key >> memory;  // Skip the "VmRSS:" part, read the number
            memory *= 1024;  // Convert from KB to bytes
            break;
        }
    }
    return memory;
}

int main() {
    // Measure memory usage before allocation
    std::size_t memoryBefore = getMemoryUsage();
    std::cout << "Memory usage before: " << memoryBefore << " bytes" << std::endl;

    // Allocate the FieldData instance (singleton or normal allocation)
    FieldData& fieldData = FieldData::getInstance();  // If using a singleton pattern
    // FieldData fieldData;  // If not using singleton

    // Measure memory usage after allocation
    std::size_t memoryAfter = getMemoryUsage();
    std::cout << "Memory usage after: " << memoryAfter << " bytes" << std::endl;

    std::cout << "Total memory increase: " << (memoryAfter - memoryBefore) << " bytes" << std::endl;
    //Now in GB
    std::cout << "Total memory increase: " << (memoryAfter - memoryBefore) / (1024.0 * 1024.0 * 1024.0) << " GB" << std::endl;

    //return 0;
    //Now test Evolution::F_K, measure time
    SpatialSlice& slice = fieldData.getV1Slice();
    auto start = std::chrono::high_resolution_clock::now();

    InitialConditions::setInitialConditions(slice);

    //Save phi data to file
    std::ofstream file_pre("phi_pre_data.txt");
    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                file_pre << slice.phi(nx, ny, nz) << std::endl;
            }
        }
    }


    for(int i=0;i<50;i++) {
        Evolution::RK4_step();
        std::cout << "Step " << i << std::endl;
        //print phi(0,0,0) to check if its not nan
        std::cout << "phi(0,0,0) = " << slice.phi(0,0,0) << std::endl;
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time for Evolution::F_K: " << elapsed.count() << " s" << std::endl;
    
    //Save phi data to file
    std::ofstream file("phi_data.txt");
    for(int nx = 0; nx < N; nx++) {
        for(int ny = 0; ny < N; ny++) {
            for(int nz = 0; nz < N; nz++) {
                file << slice.phi(nx, ny, nz) << std::endl;
            }
        }
    }

    return 0;
}