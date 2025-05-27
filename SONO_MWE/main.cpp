#include "../SONO.hpp"

#include "../curves/trefoil.hpp" //  provides the routine get_trefoil_15.

#include <chrono>

// Load name space.

// C++ way for typedefs
using Int    = SONO::Stepper::Int;
using Real   = SONO::Stepper::Real;

// Prepare a clock for timing.
using Clock  = std::chrono::high_resolution_clock;
using TimePt = std::chrono::time_point<Clock>;


constexpr Int AmbDim = SONO::Stepper::AmbDim;

int main( int argc, char** argv )
{
    // At the moment we just ignore the input arguments from the command line.
    (void) argc;
    (void) argv;
    
    std::cout << "Hello, there! This is a minimal working example for the class SONO:Stepper." << std::endl;
        
    std::cout << "Loading a simple trefoil curve..." << std::endl;
    SONO::Stepper sono( get_trefoil_15() );
    std::cout << "Loading done." << std::endl;
    
    sono.SetShiftIterationCount(1);
    sono.SetShift(0.37);
    sono.SetLengthIterationCount(3);
    sono.SetCollisionIterationCount(1);
    sono.SetStepSize(0.25);
    
    std::cout << "Number of vertices = " << sono.VertexCount() << "." << std::endl;
    
    {
        std::cout << "\nCoordinates before stepping." << std::endl;
        
        Real * x = sono.data();

        for( Int i = 0; i < sono.VertexCount(); ++i )
        {
            std::cout << "{ " << x[AmbDim * i + 0]
            << ", " << x[AmbDim * i + 1]
            << ", " << x[AmbDim * i + 2] << " }\n";
        }
    }
    
    const Int step_count = 1000000;
    
    {
        std::cout << "\nNow we step the curve for " << step_count << " steps..." << std::endl;
        
        const TimePt start_time = Clock::now();
        sono.Step(step_count);
        const TimePt stop_time  = Clock::now();
        double duration = std::chrono::duration<double>(stop_time - start_time).count();
        
        std::cout << "Stepping done. Time elapsed : " << duration << "s." << std::endl;
    }
    
    {
        std::cout << "\nCoordinates after stepping." << std::endl;
        
        Real * x = sono.data();
        
        for( Int i = 0; i < sono.VertexCount(); ++i )
        {
            std::cout << "{ " << x[AmbDim * i + 0]
            << ", " << x[AmbDim * i + 1]
            << ", " << x[AmbDim * i + 2] << " }\n";
        }
    }
    
    std::cout << "\nEnd of program." << std::endl;
}
