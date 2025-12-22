# Plot all g(r) files
set terminal pngcairo size 1000,600 enhanced font 'Verdana,10'
set output 'g_r_comparison.png'

set title "Radial Distribution Function Comparison"
set xlabel "Distance r (sigma)"
set ylabel "g(r)"
set grid
set key outside right

# Plot all RadialDist_*.dat files found
# "noenhanced" helps with filenames containing underscores
# we filter out 'RadialDist.dat' if it's identical or just plot everything matching
files = system("ls RadialDist_*.dat 2>/dev/null")

if (strlen(files) == 0) {
    print "No RadialDist_*.dat files found."
} else {
    plot for [file in files] file using 1:2 with lines lw 2 title file noenhanced
}
