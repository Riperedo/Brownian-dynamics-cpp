# Plot all g(r) files
set terminal pngcairo size 1000,600 enhanced font 'Verdana,10'
set output 'msd_comparison.png'

set title "Mean Square Displacement Comparison"
set xlabel "Time (tau)"
set ylabel "MSD"
set grid
set key outside right

set logscale
set format "10^{\%T}"

# Plot all RadialDist_*.dat files found
# "noenhanced" helps with filenames containing underscores
# we filter out 'RadialDist.dat' if it's identical or just plot everything matching
files = system("ls MSD_*.dat 2>/dev/null")

if (strlen(files) == 0) {
    print "No MSD_*.dat files found."
} else {
    plot for [file in files] file using 1:2 with lines lw 2 title file noenhanced
}
