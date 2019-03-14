// ------ language="C++" file="src/methods.hh"
#pragma once
#include "types.hh"

namespace pint
{
    // ------ begin <<forward-euler-method>>[0]
    template <typename real_t, typename vector_t>
    Integral<real_t, vector_t> forward_euler
        ( ODE<real_t, vector_t> f )
    {
        return [=]
            ( vector_t const &y
            , real_t t_init
            , real_t t_end ) -> vector_t
        {
            return y + (t_end - t_init) * f(t_init, y);
        };
    }
    // ------ end
    // ------ begin <<runge-kutta-4-method>>[0]
    template <typename real_t, typename vector_t>
    Integral<real_t, vector_t> runge_kutta_4
        ( ODE<real_t, vector_t> f )
    {
        return [=]
            ( vector_t const &y
            , real_t t_init
            , real_t t_end ) -> vector_t
        {
            real_t   t  = t_init,
                     h  = t_end - t_init;
            vector_t k1 = h * f(t, y),
                     k2 = h * f(t + h/2, y + k1/2),
                     k3 = h * f(t + h/2, y + k2/2),
                     k4 = h * f(t + h, y + k3);
            return y + (k1 + 2*k2 + 2*k3 + k4) / 6;
        };
    }
    // ------ end
    // ------ begin <<parareal-method>>[0]
    template <typename real_t, typename vector_t>
    using IterationStep = std::function
        < std::vector<vector_t>
              ( std::vector<vector_t> const &
              , std::vector<real_t> const & ) >;
    
    template <typename real_t, typename vector_t>
    IterationStep<real_t, vector_t> parareal
        ( Integral<real_t, vector_t> coarse
        , Integral<real_t, vector_t> fine )
    {
        return [=]
            ( std::vector<vector_t> const &y
            , std::vector<real_t> const &t ) -> std::vector<vector_t>
        {
            unsigned m = t.size();
            std::vector<vector_t> y_next(m);
    
            y_next[0] = y[0];
            for (unsigned i = 1; i < m; ++i) {
                y_next[i] = coarse(y_next[i-1], t[i-1], t[i])
                          + fine(y[i-1], t[i-1], t[i])
                          - coarse(y[i-1], t[i-1], t[i]);
            }
            return y_next;
        };
    }
    // ------ end
}
// ------ end
