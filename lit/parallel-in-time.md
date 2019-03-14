---
title: Parallel-in-time integration in OpenFOAM
author: Johan Hidding
tangle:
    prefix: pint
---

# Introduction

We'll implement Parallel-in-time methods in OpenFOAM. OpenFOAM is an extensive cpp code base with its own peculiar style of cpp, lots of template libraries, different executables, and using a build system called `wmake`.

First we should have a naive implementation of a PinT method in plain cpp to make sure we understand what we're doing. We will implement Parareal on top of integrators available in the GNU Scientific Library.

# Prerequisites

I will use the following libraries:

* [GNU Scientific Library](http://gnu.org/s/gsl) generic scientific algorithms.
* [Eigen3](http://eigen.tuxfamily.org) templated vector and matrix types and linear algebra.
* [Threading Building Blocks](https://www.threadingbuildingblocks.org/) task based parallel cpp framework.

# Parareal

``` {.cpp file=src/types.hh}
#pragma once
#include <functional>
#include <vector>
#include <cassert>
#include <iostream>

namespace pint
{
    <<ode-type>>
    <<method-type>>
}
```

From Wikipedia:

> Parareal solves an initial value problem of the form
> 
> $$\dot{y}(t) = f(y(t), t), \quad y(t_0) = y_0 \quad \text{with} \quad t_0 \leq t \leq T.$$
>
> Here, the right hand side $f$ can correspond to the spatial discretization of a partial differential equation in a method of lines approach.

We can define the type of problem as follows

$$y_i: \mathbb{R} \to \mathbb{R},$$

``` {.cpp #ode-type}
template <typename real_t, typename vector_t>
using Function = std::function
    < vector_t
      ( real_t
      , vector_t const & ) >;
```

and the ODE in the form

$$\frac{dy_i(t)}{dt} = f_i(t, y_1(t), \dots, y_n(t))$$

``` {.cpp #ode-type}
template <typename real_t, typename vector_t>
using ODE = std::function
    < vector_t
      ( real_t
      , vector_t const & ) >;
```

and the Jacobian (when its needed by the method)

$$J_{ij} = \frac{\partial f_i(t, y(t))}{\partial y_j},$$

``` {.cpp #ode-type}
template <typename real_t, typename vector_t, typename matrix_t>
using Jacobian = std::function
    < std::tuple<vector_t, matrix_t>
      ( real_t
      , vector_t const & ) >;
```

> Parareal now requires a decomposition of the time interval $[t_0, T]$ into $P$ so-called time slices $[t_j, t_{j+1}]$ such that
>
> $$[t_0, T] = [t_0, t_1] \cup [t_1, t_2] \cup \ldots \cup [t_{P-1}, t_{P} ].$$
>
> Each time slice is assigned to one processing unit when parallelizing the algorithm, so that $P$ is equal to the number of processing units used for Parareal.
>
> Parareal is based on the iterative application of two methods for integration of ordinary differential equations. One, commonly labelled ${\mathcal {F}}$, should be of high accuracy and computational cost while the other, typically labelled ${\mathcal {G}}$, must be computationally cheap but can be much less accurate. Typically, some form of Runge-Kutta method is chosen for both coarse and fine integrator, where ${\mathcal {G}}$ might be of lower order and use a larger time step than ${\mathcal {F}}$. If the initial value problem stems from the discretization of a PDE, ${\mathcal {G}}$ can also use a coarser spatial discretization, but this can negatively impact convergence unless high order interpolation is used. The result of numerical integration with one of these methods over a time slice $[t_{j}, t_{j+1}]$ for some starting value $y_{j}$ given at $t_{j}$ is then written as
>
> $$y = \mathcal{F}(y_j, t_j, t_{j+1})\ {\rm or}\ y = \mathcal{G}(y_j, t_j, t_{j+1}).$$

The function $\mathcal{G}$ in this instance is created from the ODE, so we say $\mathcal{G}$ has the type `Integral`, while the solver (returning an integral) will have type `StepMethod`.

``` {.cpp #method-type}
template <typename real_t, typename vector_t>
using Integral = std::function
    < vector_t
      ( vector_t const &
      , real_t
      , real_t ) >;

<<solve-function>>

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
```

The easiest example is the forward Euler method.

$$y_{i + 1} = y_i + h f(t_i, y_i)$$

``` {.cpp #forward-euler-method}
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
```

In general we will want to solve the ODE for a range of times. The following convenience function takes an integral, a $y_0$ and a range of times $t_i$, and returns $y_i$.

``` {.cpp #solve-function}
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
```

> Serial time integration with the fine method would then correspond to a step-by-step computation of
>
> $$y_{j+1} = \mathcal{F}(y_j, t_j, t_{j+1}), \quad j=0, \ldots, P-1.$$
>
> Parareal instead uses the following iteration
>
> $$y_{j+1}^{k+1} = \mathcal{G}(y^{k+1}_j, t_j, t_{j+1}) + \mathcal{F}(y^k_j, t_j, t_{j+1}) - \mathcal{G}(y^k_j, t_j, t_{j+1}), \quad j=0, \ldots, P-1, \quad k=0, \ldots, K-1,$$
>
> where $k$ is the iteration counter. As the iteration converges and $y^{k+1}_j - y^k_j \to 0$, the terms from the coarse method cancel out and Parareal reproduces the solution that is obtained by the serial execution of the fine method only. It can be shown that Parareal converges after a maximum of $P$ iterations. For Parareal to provide speedup, however, it has to converge in a number of iterations significantly smaller than the number of time slices, that is $K \ll P$.
>
> In the Parareal iteration, the computationally expensive evaluation of $\mathcal{F}(y^k_j, t_j, t_{j+1})$ can be performed in parallel on $P$ processing units. By contrast, the dependency of $y^{k+1}_{j+1}$ on $\mathcal{G}(y^{k+1}_j, t_j, t_{j+1})$ means that the coarse correction has to be computed in serial order.

``` {.cpp #parareal-method}
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
```

# Runge-Kutta Method

The fourth order Runge-Kutta method is given as follows

$$\begin{aligned}
k_1 &= h f(t_i, y_i)\\
k_2 &= h f(t_i + h/2, y_i + k_1 / 2)\\
k_3 &= h f(t_i + h/2, y_i + k_2 / 2)\\
k_4 &= h f(t_i + h, y_i + k_3)\\
y_{i+1} &= y_i + (k_1 + 2k_2 + 2k_3 + k_4) / 6
\end{aligned}$$

``` {.cpp #runge-kutta-4-method}
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
```

# Harmonic Oscillator

Let's test the integrator on a damped oscillator

$$y'' + 2\zeta \omega_0 y' + \omega_0^2 y = 0,$$

where $\omega_0 = \sqrt{k/m}$ and $\zeta = c / 2\sqrt{mk}$, $k$ being the spring constant, $m$ the test mass and $c$ the friction constant.

To solve this second order ODE we need to introduce a second variable to solve for. Say $q = y$ and $p = y'$.

$$\begin{aligned}
    q' &= p\\
    p' &= -2\zeta \omega_0 p + \omega_0^2 q
\end{aligned}$$

``` {.cpp #harmonic-oscillator}
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
```

## Synthesis

``` {.cpp file=src/methods.hh}
#pragma once
#include "types.hh"

namespace pint
{
    <<forward-euler-method>>
    <<runge-kutta-4-method>>
    <<parareal-method>>
}
```

``` {.cpp file=src/test-rk4.cc}
#include "methods.hh"
#include "types.hh"

#include <argagg/argagg.hpp>

#include <vector>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <Eigen/Dense>

using namespace pint;

<<harmonic-oscillator>>

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
```

## Test results

``` {.gnuplot file=images/plot1.gp}
set term svg font "Bitstream Charter"
set output 'images/plot1.svg'

plot '< ./parareal' \
       i 0 w l lc 'black' lw 2 t 'runge-kutta-4, h=0.01', \
    '' i 1 w lp t'parareal n=0', \
    '' i 2 w lp t'parareal n=1', \
    '' i 3 w lp t'parareal n=2', \
    '' i 4 w lp ls 6 t'parareal n=3'
```

# Building

``` {.makefile file=Makefile}
# allow spaces for indenting
.RECIPEPREFIX +=

# find sources
build_dir = ./build
cc_files = $(shell find ./src -name *.cc)
obj_files = $(cc_files:%.cc=$(build_dir)/%.o)
dep_files = $(obj_files:%.o=%.d)

# set compiler
compile = g++
link = g++

# libfmt
fmtlib_lflags = -lfmt

# eigen3
eigen_lflags = $(shell pkg-config --libs eigen3)
eigen_cflags = $(shell pkg-config --cflags eigen3)

# compile and link flags
compile_flags = -O3 -std=c++17 -Wall -Werror $(eigen_cflags)
link_flags = $(fmt_lflags) $(eigen_lflags)

# rules
.PHONY: clean build

build: parareal

clean:
    rm -rf $(build_dir)

-include $(dep_files)

$(build_dir)/%.o: %.cc Makefile
    @mkdir -p $(@D)
    $(compile) $(compile_flags) -MMD -c $< -o $@

parareal: $(obj_files)
    @mkdir -p $(@D)
    $(link) $^ $(link_flags) -o $@
```
