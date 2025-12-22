# Plot Averaged g(r) from Production Runs
set terminal pngcairo size 1000,600 enhanced font 'Verdana,10'
set output 'averaged_msd_comparison.png'

set title "Averaged Mean self-Intermediate (ScatterinHard Spheres)"
set xlabel "Time (tau)"
set ylabel "MSD"
set grid
set key top right

# Define colors for different Phi
set style line 1 lc rgb '#E6194B' lw 2 # Red
set style line 2 lc rgb '#3CB44B' lw 2 # Green
set style line 3 lc rgb '#FFE119' lw 2 # Yellow
set style line 4 lc rgb '#4363D8' lw 2 # Blue
set style line 5 lc rgb '#F58231' lw 2 # Orange
set style line 6 lc rgb '#cf21b8ff' lw 2 # Orange
set style line 7 lc rgb '#5b31f5b7' lw 2 # Orange

set logscale
set format "10^{%L}"

set key left top
set xrange [0.001:100]

# Plot command
# Assuming directory structure: ../HardSphere/phi_0.X/average_RadialDist.dat
# We plot specific known Phis for clarity
plot "../HardSphere/phi_0.1/average_MSD.dat" using 1:2 with lines ls 1 title "Phi = 0.1", \
     "../HardSphere/phi_0.2/average_MSD.dat" using 1:2 with lines ls 2 title "Phi = 0.2", \
     "../HardSphere/phi_0.3/average_MSD.dat" using 1:2 with lines ls 3 title "Phi = 0.3", \
     "../HardSphere/phi_0.4/average_MSD.dat" using 1:2 with lines ls 4 title "Phi = 0.4", \
     "../HardSphere/phi_0.5/average_MSD.dat" using 1:2 with lines ls 5 title "Phi = 0.5", \
     "../HardSphere/phi_0.55/average_MSD.dat" using 1:2 with lines ls 6 title "Phi = 0.55", \
     "../HardSphere/phi_0.56/average_MSD.dat" using 1:2 with lines ls 7 title "Phi = 0.56"
