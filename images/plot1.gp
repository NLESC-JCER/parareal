## ------ language="Gnuplot" file="images/plot1.gp"
set term svg font "Bitstream Charter"
set output 'images/plot1.svg'

plot '< ./parareal' \
       i 0 w l lc 'black' lw 2 t 'runge-kutta-4, h=0.01', \
    '' i 1 w lp t'parareal n=0', \
    '' i 2 w lp t'parareal n=1', \
    '' i 3 w lp t'parareal n=2', \
    '' i 4 w lp ls 6 t'parareal n=3'
## ------ end
