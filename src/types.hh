// ------ language="C++" file="src/types.hh"
#pragma once
#include <functional>
#include <vector>
#include <cassert>
#include <iostream>

namespace pint
{
    // ------ begin <<ode-type>>[0]
    template <typename real_t, typename vector_t>
    using Function = std::function
        < vector_t
          ( real_t
          , vector_t const & ) >;
    // ------ end
    // ------ begin <<ode-type>>[1]
    template <typename real_t, typename vector_t>
    using ODE = std::function
        < vector_t
          ( real_t
          , vector_t const & ) >;
    // ------ end
    // ------ begin <<ode-type>>[2]
    template <typename real_t, typename vector_t, typename matrix_t>
    using Jacobian = std::function
        < std::tuple<vector_t, matrix_t>
          ( real_t
          , vector_t const & ) >;
    // ------ end
    // ------ begin <<method-type>>[0]
    template <typename real_t, typename vector_t>
    using Integral = std::function
        < vector_t
          ( vector_t const &
          , real_t
          , real_t ) >;
    
    // ------ begin <<solve-function>>[0]
    template <typename real_t, typename vector_t>
    std::vector<vector_t> solve
        ( Integral<real_t, vector_t> step
        , vector_t const &y_0
        , std::vector<real_t> const &t )
    {
        assert(t.size() > 0);
        std::vector<vector_t> y(t.size());
        y[0] = y_0;
        for (unsigned i = 1; i < t.size(); ++i) {
            y[i] = step(y[i-1], t[i-1], t[i]);
        }
        return y;
    }
    
    template <typename real_t>
    std::vector<real_t> linspace(real_t start, real_t end, unsigned n)
    {
        std::vector<real_t> x(n);
        for (unsigned i = 0; i < n; ++i) {
            x[i] = start + (end * i - start * i) / (n - 1);
        }
        return x;
    }
    // ------ end
    
    template <typename real_t, typename vector_t>
    using StepMethod = std::function
        < Integral<real_t, vector_t> 
          ( ODE<real_t, vector_t> ) >;
    
    template <typename real_t, typename vector_t>
    Integral<real_t, vector_t> iterate_step
        ( Integral<real_t, vector_t> step
        , real_t h )
    {
        return [=]
            ( vector_t const &y_0
            , real_t t_init
            , real_t t_end ) -> vector_t
        {
            real_t t = t_init;
            vector_t y = y_0;
            while (t + h < t_end) {
                y = step(y, t, t + h);
                t += h;
            }
            return step(y, t, t_end);
        };
    }
    // ------ end
}
// ------ end
