#pragma once

#include <stdio.h>
// Load the bread-and-butter container class std::vector
#include <vector>

#include <algorithm>

// Load very slow random number generator std::random_device.
#include <random>

// Load the fast pseudo random number generator pcg64.
#include "submodules/pcg-cpp/include/pcg_random.hpp"

namespace SONO
{
    // Very rough port of Keenan Crane's experimental implementation.
    // Many magic numbers for which I do not know how to adjust them.
    
    class Stepper
    {
        
    public:
        
        using Real = double;
        using Int  = std::size_t;

        using VertexContainer_T = std::vector<Real>;
        using PRNG_T = pcg64;
        
        static constexpr Int AmbDim = 3;
        
        Stepper() = default;
        
        ~Stepper() = default;
        
        explicit Stepper( const Real * vertex_coords_, const Int vertex_count_ )
        :   vertex_count  ( vertex_count_ )
        ,   vertex_coords ( vertex_count * AmbDim )
        ,   vertex_buffer ( vertex_count * AmbDim )
        {
            std::copy_n( vertex_coords_, vertex_count * AmbDim, &vertex_coords[0]);
            
            InitializeRandomEngine();
        }
        
        explicit Stepper( const std::vector<Real> & vertex_coords_ )
        :   Stepper( &vertex_coords_[0], static_cast<Int>(vertex_coords_.size())/3 )
        {}
        
    private:
        
        Int vertex_count = 0;
        
        VertexContainer_T vertex_coords;
        VertexContainer_T vertex_buffer;
        
        PRNG_T random_engine;
        
        Real average_overlap       = 0.0;
        Real shift                 = 0.37; // No idea what a good value would be.
        
        // `collision_damping` should be in the interval [0,1).
        // `collision_damping = 0` means no damping (as in the original SONO).
        // `collision_damping = 1` means no progress at all.
        Real collision_damping     = 0.0;
        Real target_overlap_factor = 0.005;
        Real step_size             = 0.25; // No idea what a good value would be.
        Real shrink_factor         = 0.999;
        
        Int  shift_iter_count      = 30;
        Int  length_iter_count     = 3;
        Int  collision_iter_count  = 1;
        
        static constexpr Real node_diam = 2.;

        
    private:
        
        void InitializeRandomEngine()
        {
            using RNG_T = std::random_device;
            
            static constexpr std::size_t PRNG_T_state_size = 4;
        
            static constexpr std::size_t seed_size  = (PRNG_T_state_size * sizeof(PRNG_T::result_type)) / sizeof(RNG_T::result_type);
            
            std::array<RNG_T::result_type,seed_size> seed;
            
            std::generate( seed.begin(), seed.end(), RNG_T() );
            
            std::seed_seq seed_sequence ( seed.begin(), seed.end() );
            
            random_engine = PRNG_T( seed_sequence );
        }
        
    public:
        
        Int VertexCount() const
        {
            return vertex_count;
        }
        
        const VertexContainer_T & VertexCoordinates() const
        {
            return vertex_coords;
        }

        VertexContainer_T & VertexCoordinates()
        {
            return vertex_coords;
        }
        
        const Real * data() const
        {
            return &vertex_coords[0];
        }
        
        Real * data()
        {
            return &vertex_coords[0];
        }
        
        Real CurveLength() const
        {
            return Stepper::CurveLength( &vertex_coords[0], vertex_count );
        }
        
        Real StepSize() const
        {
            return step_size;
        }
        
        void SetStepSize( const Real value )
        {
            step_size = value;
        }
        
        
        Int ShiftIterationCount() const
        {
            return shift_iter_count;
        }
        
        void SetShiftIterationCount( const Int value )
        {
            shift_iter_count = value;
        }
        
        Real Shift() const
        {
            return shift;
        }
        
        void SetShift( const Real value )
        {
            shift = value;
        }
        
        
        Int LengthIterationCount() const
        {
            return length_iter_count;
        }
        
        void SetLengthIterationCount( const Int value )
        {
            length_iter_count = value;
        }
        
        
        Int CollisionIterationCount() const
        {
            return collision_iter_count;
        }
        
        void SetCollisionIterationCount( const Int value )
        {
            collision_iter_count = value;
        }
        
        Real CollisionDamping() const
        {
            return collision_damping;
        }
        
        void SetCollisionDamping( const Real value )
        {
            collision_damping = value;
        }
        
        
        Real ShrinkFactor() const
        {
            return shrink_factor;
        }
        
        void SetShrinkFactor( const Real value )
        {
            shrink_factor = value;
        }
        
        
    public:
        
        // Take `step_count` steps of SONO algorithm to tighten curve.
        void Step( const Int step_count )
        {
            if( vertex_count <= Int(0) )
            {
                return;
            }
            
            try
            {
                for( Int step = 0; step < step_count; ++step )
                {
                    std::copy_n(
                        &vertex_coords[0], AmbDim * vertex_count, &vertex_buffer[0]
                    );
                    
                    // TODO: What would best order of these operations?
                    for( Int iter = 0; iter < collision_iter_count; ++iter )
                    {
                        UpdateCollisions( &vertex_buffer[0], vertex_count );
                    }
                    for( Int iter = 0; iter < length_iter_count; ++iter )
                    {
                        UpdateEdgeLengths( &vertex_buffer[0], vertex_count );
                    }
                    for( Int iter = 0; iter < shift_iter_count; ++iter )
                    {
                        Shift( shift, &vertex_buffer[0], vertex_count );
                    }
                    
                    // if there's not much overlap remaining, shrink the curve
                    // and target curve length (but not the diameter)
                    if( average_overlap < target_overlap_factor * node_diam )
                    {
                        Scale( shrink_factor, vertex_buffer.data(), vertex_count );
                    }
                    
                    // TODO: Remove LinearCombine.
                    
                    const double alpha = step_size;
                    const double beta  = Real(1) - step_size;
                    
                    double * x = &vertex_buffer[0];
                    double * y = &vertex_coords[0];
                    
                    for( Int i = 0; i < AmbDim * vertex_count; ++i )
                    {
                        y[i] = alpha * x[i] + beta  * y[i];
                    }
                }
            }
            catch( const std::exception & e )
            {
                std::cout << "Error: " << e.what() << std::endl;
                throw;
            }
            catch(...)
            {
                std::cout << "Error: Unknown exception." << std::endl;
                throw;
            }
        }
        
        
        static Real CurveLength( const Real * x, const Int n )
        {
            Real L  = 0;

            for( Int j = 1; j < n; ++j )
            {
                const Int i = j - 1;
                
                // get the current distance
                    const Real uij [3] = {
                    x[3 * j + 0] - x[3 * i + 0],
                    x[3 * j + 1] - x[3 * i + 1],
                    x[3 * j + 2] - x[3 * i + 2]
                };

                L += std::sqrt(uij[0] * uij[0] + uij[1] * uij[1] + uij[2] * uij[2]);
            }
            
            // Handle wrap-around.
            {
                Int i = n - 1;
                Int j = 0;
                
                // get the current distance
                const Real uij [3] = {
                    x[3 * j + 0] - x[3 * i + 0],
                    x[3 * j + 1] - x[3 * i + 1],
                    x[3 * j + 2] - x[3 * i + 2]
                };

                L += std::sqrt(uij[0] * uij[0] + uij[1] * uij[1] + uij[2] * uij[2]);
            }
            
            return L;
        }
        
        
        // Compute x[i] += alpha * (x[(i+1) % n] - x[i]).
        static void Shift( const Real alpha, Real * x, const Int n )
        {
            // Remember the position of the first node so that it is not overwritten
            const Real x_0 [3] = {x[0],x[1],x[2]};
            
            for( Int j = 1; j < n; ++j )
            {
                const Int i = j - 1;
                
                x[3 * i + 0] += alpha * (x[3 * j + 0] - x[3 * i + 0]);
                x[3 * i + 1] += alpha * (x[3 * j + 1] - x[3 * i + 1]);
                x[3 * i + 2] += alpha * (x[3 * j + 2] - x[3 * i + 2]);
            }
            
            // Handle wrap-around.
            {
                //const Int j = 0;
                const Int i = n - 1;
                
                x[3 * i + 0] += alpha * (x_0[0] - x[3 * i + 0]);
                x[3 * i + 1] += alpha * (x_0[1] - x[3 * i + 1]);
                x[3 * i + 2] += alpha * (x_0[2] - x[3 * i + 2]);
            }
        }
        
        // Compute x = alpha * x.
        static void Scale( const Real alpha, Real * x, const Int n )
        {
            for( Int i = 0; i < AmbDim * n; ++i )
            {
                x[i] *= alpha;
            }
        }
        
        // Classical axpby:
        // Compute y = alpha * x + beta * y.
        static void LinearCombine(
            const Real alpha, const Real * x, const Real beta, Real * y, const Int n
        )
        {
            for( Int i = 0; i < AmbDim * n; ++i )
            {
                y[i] = alpha * x[i] + beta * y[i];
            }
        }
        
        void UpdateEdgeLengths( Real * x, const Int n )
        {
            std::uniform_int_distribution<Int> uniform_dist( Int(0), n - Int(1) );
            std::uniform_int_distribution<Int> coin_flip   ( Int(0), Int(1) );
            
            Int i = uniform_dist( random_engine ); // starting node
            const Int dir = (coin_flip( random_engine ) ? Int(1) : n - Int(1) ); // direction of traversal
            
            const Real curve_length = CurveLength(x,n);

            const Real target_length = curve_length / static_cast<Real>(n);

            // iterate over curve, projecting each edge onto the correct length
            // (essentially via nonlinear Gauss-Seidel)
            for( Int k = 0; k < n; ++k )
            {
                // Get the next neighbor.
//                Int j = (i + 1) % n;
                Int j = i + dir;
                if( j >= n ) { j -= n; }

                // move the vertices to the target length

                // midpoint
                const Real m [3] = {
                    Real(0.5) * x[3 * i + 0] + Real(0.5) * x[3 * j + 0],
                    Real(0.5) * x[3 * i + 1] + Real(0.5) * x[3 * j + 1],
                    Real(0.5) * x[3 * i + 2] + Real(0.5) * x[3 * j + 2]
                };

                // get the current distance
                const Real u [3] = {
                    x[3 * j + 0] - x[3 * i + 0],
                    x[3 * j + 1] - x[3 * i + 1],
                    x[3 * j + 2] - x[3 * i + 2]
                };

                const Real d = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);

                const Real scale = (Real(0.5) * target_length) / d;

                x[3 * i + 0] = m[0] - scale * u[0];
                x[3 * i + 1] = m[1] - scale * u[1];
                x[3 * i + 2] = m[2] - scale * u[2];

                x[3 * j + 0] = m[0] + scale * u[0];
                x[3 * j + 1] = m[1] + scale * u[1];
                x[3 * j + 2] = m[2] + scale * u[2];

                // go to the next node around the curve
                i = j;
            }
        }
        
        void UpdateCollisions( Real * x, const Int n )
        {
            Real alpha = Real(1.) - collision_damping;
            Real beta  = collision_damping;
            
            constexpr Real delta   = 0.001; // "wiggle factor" between hard spheres
            constexpr Real epsilon = 0.0;
            
            constexpr Real D_eps = node_diam + epsilon;
            
            
            const Real curve_length = CurveLength(x,n);
            
            // Determine how many neighboring nodes to skip.
            const Real l = curve_length / static_cast<Real>(n);
            const Real D = node_diam;
            
            // TODO: Using the average edge length is a bit fishy as the polygon is never really equilateral.
            // TODO: One should rather employ the true arc distance here, which can be done efficiently, at least if we "freeze" the edge lengths before the loop.
            const Int skip = static_cast<Int>( std::ceil( M_PI * D / (Real(2) * l) ) );
            
            if( skip >= n )
            {
                // This is how an error message is thrown.
                throw std::runtime_error("skip >= n. There is no point in continuing.");
            }
            
            std::uniform_int_distribution<Int> uniform_dist( Int(0), n - Int(1) );
            Int i = uniform_dist( random_engine ); // starting node

            // Keep track of the average amount by which nodes overlap;
            // if no nodes overlap, this amount will default to zero.
            
            average_overlap = 0.;
            Int overlap_count = 0;
            
            // Iterate over curve, projecting each pair of nodes onto a non-intersecting configuration (essentially via nonlinear Gauss-Seidel).
            for( Int k = 0; k < n; ++k )
            {
                // TODO: skipping too much? too little? off by 1 on either end?
                for( Int m = skip + 1; m < n - skip; ++m )
                {

//                    Int j = (i+m) % n;
                    Int j = i + m;
                    if( j >= n ) { j-= n; };
                    
                    if( j == i ) continue;
                    
                    // distance vector u = x[j] - x[i]
                    const Real u [3] = {
                        x[3 * j + 0] - x[3 * i + 0],
                        x[3 * j + 1] - x[3 * i + 1],
                        x[3 * j + 2] - x[3 * i + 2]
                    };
                    
                    // current distance
                    const Real d = std::sqrt( u[0] * u[0] + u[1] * u[1] + u[2] * u[2] );
                    
                    if( d < D_eps )
                    {
                        average_overlap += std::max(Real(0),D - d);
                        ++overlap_count;
                        
                        // move the vertices to the target length
                        
                        // midpoint
                        const Real m [3] = {
                            Real(0.5) * x[3 * i + 0] + Real(0.5) * x[3 * j + 0],
                            Real(0.5) * x[3 * i + 1] + Real(0.5) * x[3 * j + 1],
                            Real(0.5) * x[3 * i + 2] + Real(0.5) * x[3 * j + 2]
                        };
                        
                        const Real scale = (0.5*(D+delta)) / d;
                        
                        const Real c_i [3] = {
                            m[0] - scale * u[0],
                            m[1] - scale * u[1],
                            m[2] - scale * u[2]
                        };
                        
                        // Applying some damping here (if damping > 0).
                        // The original SONO algorithm does not do it.
                        x[3 * i + 0] = alpha * c_i[0] + beta * x[3 * i + 0];
                        x[3 * i + 1] = alpha * c_i[1] + beta * x[3 * i + 1];
                        x[3 * i + 2] = alpha * c_i[2] + beta * x[3 * i + 2];

                        const Real c_j [3] = {
                            m[0] + scale * u[0],
                            m[1] + scale * u[1],
                            m[2] + scale * u[2]
                        };
                        
                        // Applying some damping here (if damping > 0).
                        // The original SONO algorithm does not do it.
                        x[3 * j + 0] = alpha * c_j[0] + beta * x[3 * j + 0];
                        x[3 * j + 1] = alpha * c_j[1] + beta * x[3 * j + 1];
                        x[3 * j + 2] = alpha * c_j[2] + beta * x[3 * j + 2];
                    }
                }
                
                // Go to the next node around the curve.
                // i = (i+1) % n;
                ++i;
                if( i >= n ) { i -= n; }
            }
            
            if( overlap_count > 0 )
            {
                average_overlap /= static_cast<Real>(overlap_count);
            }
        }
        
    }; // class Stepper
    
} // namespace SONO
