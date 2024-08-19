#include <iostream>
#include <math.h>

// include for fopen
#include <stdio.h>


// Define constants for the grid size and time steps
#define NX 100  // Number of grid points in the x-direction
#define NY 100  // Number of grid points in the y-direction
#define NT 1000  // Number of time steps

// Create a 3D array to store field data for each point in space and time
float fieldData[NX][NY][NT];

// Define the time step
float dt = 0.001;

// Define the spatial resolution
float dx = 1.0 / NX;
float dy = 1.0 / NY;

// Define the physical dimensions of the domain
float LX = 1.0, LY = 1.0;

float c = 1.0;

// Function to initialize the field data
void initializeField()
{
    // Loop over each grid point in the x and y directions
    for(int i = 0; i < NX; i++)
    {
        for(int j = 0; j < NY; j++)
        {
            // Calculate the physical coordinates
            float x = i * LX / NX;
            float y = j * LY / NY;

            // Initialize the field with a sine wave pattern at t=0 and t=1
            fieldData[i][j][0] = sin(x * 2.0 * M_PI) * sin(y * 2.0 * M_PI);
            

            // Use the wave equation to estimate the field at t=1 using a central difference in time
            // Assume that the initial velocity is zero, i.e., ∂f/∂t = 0 at t = 0.
            fieldData[i][j][1] = fieldData[i][j][0] + 0.5 * c * c * dt * dt * (
                (sin((x + dx) * 2.0 * M_PI) + sin((x - dx) * 2.0 * M_PI) - 2.0 * sin(x * 2.0 * M_PI)) / (dx * dx) +
                (sin((y + dy) * 2.0 * M_PI) + sin((y - dy) * 2.0 * M_PI) - 2.0 * sin(y * 2.0 * M_PI)) / (dy * dy)
            );
        }
    }
}

// Function to perform one leapfrog step of the simulation
void LeapfrogStep(int nt)
{
    // Loop over each interior grid point in the x and y directions
    for(int i = 1; i < NX - 1; i++)
    {
        for(int j = 1; j < NY - 1; j++)
        {
            // Calculate periodic boundary conditions
            int nxp1 = (i + 1) % NX;  // Next point in x, wrapping around
            int nxm1 = (i - 1 + NX) % NX;  // Previous point in x, wrapping around
            int nyp1 = (j + 1) % NY;  // Next point in y, wrapping around
            int nym1 = (j - 1 + NY) % NY;  // Previous point in y, wrapping around

            // Calculate the second spatial derivatives using central differences
            float d2fdx2 = (fieldData[nxp1][j][nt] - 2.0 * fieldData[i][j][nt] + fieldData[nxm1][j][nt]) / (dx * dx);
            float d2fdy2 = (fieldData[i][nyp1][nt] - 2.0 * fieldData[i][j][nt] + fieldData[i][nym1][nt]) / (dy * dy);

            // Update the field using the leapfrog method (second-order in time)
            fieldData[i][j][nt + 1] = 2.0 * fieldData[i][j][nt] - fieldData[i][j][nt - 1] + c * c * dt * dt * (d2fdx2 + d2fdy2);
        }
    }
}

int main()
{
    // Print a greeting message to the console
    
    initializeField();

    //print first couple of rows
    for(int i = 0; i < 10; i++)
    {
        for(int j = 0; j < 10; j++)
        {
            //std::cout << fieldData[i][j][0] << " ";
        }
        //std::cout << std::endl;
    }

    // Loop over each time step
    for(int nt = 1; nt < NT-1; nt++)
    {
        // Perform one Euler step of the simulation
        LeapfrogStep(nt);
    }



     //print first couple of rows
    for(int i = 0; i < 10; i++)
    {
        for(int j = 0; j < 10; j++)
        {
            std::cout << fieldData[i][j][99] << " ";
        }
        std::cout << std::endl;// << std::endl;
    }


    //output the data into a file to be parsed by python
    //the format should be:
    //first nt=0, then nt=1, etc... separated by "___"
    //each nt consists of a matrix of csv formatted values

    //open file
    FILE *f = fopen("output.txt", "w");
    if (f == NULL)
    {
        printf("Error opening file!\n");
        exit(1);
    }

    //write the data
    int NT_calc = 0;
    for(int nt = 0; nt < NT; nt+=10)
    {
        for(int i = 0; i < NX; i++)
        {
            for(int j = 0; j < NY; j++)
            {
                fprintf(f, "%f", fieldData[i][j][nt]);
                if(j < NY-1)
                {
                    fprintf(f, ",");
                }
            }
            fprintf(f, "\n");
        }
        if(nt < NT-11)
        {
            fprintf(f, "___\n");
        }

        NT_calc++;
    }
    std::cout << "NT_calc: " << NT_calc << std::endl;

    return 0;
}
