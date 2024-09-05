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

    //Now test Evolution::F_K, measure time
    SpatialSlice& slice = fieldData.getCurrentSlice();
    SpatialSlice& sliceOut = fieldData.getFutureSlice();
    auto start = std::chrono::high_resolution_clock::now();
    Evolution::F_K(slice, sliceOut);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Elapsed time for Evolution::F_K: " << elapsed.count() << " s" << std::endl;
    

    return 0;
}