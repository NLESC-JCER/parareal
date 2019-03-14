// ------ language="C++" file="src/test-rk4.cc"
#include "methods.hh"
#include "types.hh"

#include <argagg/argagg.hpp>

#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>

using namespace pint;

// ------ begin <<harmonic-oscillator>>[0]
template <typename real_t, typename vector_t>
ODE<real_t, vector_t> harmonic_oscillator
    ( real_t omega_0
    , real_t zeta )
{
    return [=] (real_t t, vector_t const &y) {
        return vector_t
            ( y[1] 
            , -2 * zeta * omega_0 * y[1] - omega_0*omega_0 * y[0] );
    };
}
// ------ end

template <typename real_t, typename vector_t>
std::vector<vector_t> solve_iterative
    ( IterationStep<real_t, vector_t> step
    , std::vector<vector_t> const &y_0
    , std::vector<real_t> const &t
    , real_t abs_err
    , unsigned max_iter )
{
    std::vector<vector_t> y = y_0;

    for (unsigned i = 0; i < t.size(); ++i) {
        std::cout << t[i] << " " << y[i][0] << " " << y[i][1] << std::endl;
    }
    std::cout << "\n\n";

    for (unsigned i = 0; i < max_iter; ++i) {
        auto y_next = step(y, t);
        real_t max_err = 0.0;
        for (unsigned i = 0; i < t.size(); ++i) {
            real_t err = (y_next[i] - y[i]).norm();
            max_err = (err > max_err ? err : max_err);
        }
        y = y_next;

        std::cout << "# iteration=" << i + 1
                  << " max_abs_err=" << max_err << "\n";
        for (unsigned i = 0; i < t.size(); ++i) {
            std::cout << t[i] << " " << y[i][0] << " " << y[i][1] << std::endl;
        }
        std::cout << "\n\n";
    }
    return y;
}

int main(int argc, char **argv)
{
    argagg::parser argparser
        {{ { "help",   {"-h", "--help"}
           , "shows this help message", 0 }
         , { "omega0", {"--omega0"}
           , "undamped angular frequency (default 1.0)", 1 }
         , { "zeta",   {"--zeta"}
           , "damping ratio (default 0.5)", 1 }
         , { "n",      {"--n"}
           , "number of time slices (default 9)", 1 }
         , { "h",      {"--h"}
           , "size of time step in fine integrator (default 0.01)", 1 }
        }};
    
    argagg::parser_results args;
    try {
        args = argparser.parse(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    if (args["help"]) {
        std::cerr << "Parareal test case: harmonic oscillator\n";
        std::cerr << argparser;
        return EXIT_SUCCESS;
    }

    using real_t = double;
    using vector_t = Eigen::Vector2d;

    unsigned n    = args["n"].as<unsigned>(9);
    real_t h      = args["h"].as<real_t>(0.01);
    real_t omega0 = args["omega0"].as<real_t>(1.0);
    real_t zeta   = args["zeta"].as<real_t>(0.5);

    auto ts = linspace<real_t>(0, 15.0, n);
    auto ode = harmonic_oscillator<real_t, vector_t>(omega0, zeta);
    auto coarse = runge_kutta_4<real_t, vector_t>(ode);
    auto fine = iterate_step<real_t, vector_t>(coarse, h); 
    auto y_0 = solve(coarse, vector_t(1.0, 0.0), ts);

    auto t_ref = linspace<real_t>(0, 15.0, 100);
    auto y_ref = solve(fine, vector_t(1.0, 0.0), t_ref);
    for (unsigned i = 0; i < t_ref.size(); ++i) {
        std::cout << t_ref[i] << " " << y_ref[i][0] << " " << y_ref[i][1] << std::endl;
    }
    std::cout << "\n\n";

    auto y = solve_iterative
        ( parareal(coarse, fine)
        , y_0
        , ts
        , 1e-6
        , n );

    return EXIT_SUCCESS;
}
// ------ end
