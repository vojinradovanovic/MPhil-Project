# --- Constants & Configuration ---
    x_min = 0.0
    x_max = 0.1

    y_min = 0.0
    y_max = 0.1

    X     = 1   # Column index for X centroid
    Y     = 2   # Column index for Y centroid
    VAL = 8   # Column index for the solution variable

    file = "diffusion"
    data = "../../data/" . file . ".dat"

# --- Data Statistical Analysis ---
    stats data using VAL prefix "val" nooutput

# --- Terminal & Output Settings ---
    set term pngcairo size 1000, 700 enhanced font "Helvetica Bold,24"
    set output "../../visualisations/" . file . "_dirichlet.png"

# --- Plot Styles, Borders & Tics ---
    set border 15 lw 1.5 lc rgb "#333333"  
    set tics out nomirror scale 0.75       
    set mxtics 2                           
    set mytics 2
    set size square

# --- Tick Label Font Reductions ---
    set xtics font "Helvetica Bold,32"
    set ytics font "Helvetica Bold,32"
    set cbtics font "Helvetica Bold,32"

    cb_step = (val_max-val_min)/5.0
    cb_step = real(sprintf("%.1g", cb_step))
    set cbtics cb_step
    set xtics 0.05
    set ytics 0.05

# --- Layout Margin Adjustments ---
    set rmargin 5
    set lmargin 8
    set bmargin 4

# --- Axis Ranges ---
    set xrange [x_min:x_max]
    set yrange [y_min:y_max]
    set cbrange [val_min:val_max]

# --- Color Palette ---
    set palette viridis

# --- Labels, Titles & Colorbar ---
    set xlabel "X" font "Helvetica Bold,32" tc rgb "#111111" offset 0,-0.5
    set ylabel "Y" font "Helvetica Bold,32" tc rgb "#111111" offset -0.5,0
    set cblabel "ϕ value" font "Helvetica Bold,32" tc rgb "#111111" offset 3,0
    set format cb "%.3g"

# --- Legend / Key ---
    set key off

# --- Plot Command ---
    plot data using X:Y:VAL with image