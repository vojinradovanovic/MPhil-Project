# --- Constants & Configuration ---
    x_min = 0
    x_max = 0.1

    y_min = 0
    y_max = 0.1

    # columns
    X = 1           # X centroid
    Y = 2           # Y centroid
    P = 3
    VX_START = 4    # vector origin
    VY_START = 5
    VX_MOD = 6      # vector length coordinates
    VY_MOD = 7
    V_MOD = 8

    array names[8] = ["x centroid", "y centroid", "pressure", "vx start coordinate", "vy start coordinate", \
                        "vx value", "vy value", "v modulus"]

    frame = 99
    steps = 1

    file = "simple_end"
    data = "../../data/" . file . ".dat"

# --- Set plotting column ---
    COL = P

# --- Data Statistical Analysis ---
    stats data using COL name "col" nooutput

# --- Terminal & Output Settings ---
    set term pngcairo enhanced font "Helvetica Bold,24" size 1100, 600
    set output "../../visualisations/" . file . ".png"

# --- Plot Styles, Borders & Tics ---
    set border 15 lw 1.5 lc rgb "#333333"  
    set tics out nomirror scale 0.75       
    set mxtics 2                           
    set mytics 2

# --- Tick Label Font & Format Configuration ---
    set xtics font "Helvetica Bold,32"
    set ytics font "Helvetica Bold,32"
    set cbtics font "Helvetica Bold,32"

    cb_step = (col_max - col_min) / 5.0
    set cbtics cb_step
    set xtics 0.05 offset 0,-0.5
    set ytics 0.05 offset -0.5,0

# --- Layout Margin Adjustments ---
    set rmargin 14
    set lmargin 14
    set bmargin 4

# --- Axis Ranges ---
    set xrange [x_min:x_max]
    set yrange [y_min:y_max]
    set cbrange [col_min:col_max]

# --- Color Palette ---
    set palette defined (col_min "blue", 1 "white", col_max "red")

# --- Labels, Titles & Colorbar ---
    set xlabel "X" font "Helvetica Bold,32" tc rgb "#111111" offset 0,-0.5
    set ylabel "Y" font "Helvetica Bold,32" tc rgb "#111111" offset -3.5,0
    set cblabel names[COL] font "Helvetica Bold,32" tc rgb "#111111" offset 5,0

# --- Legend / Key ---
    set key off

# --- Plot Command ---
    plot data using X:Y:COL index frame with image, \
         data using VX_START:VY_START:VX_MOD:VY_MOD index frame every 7 \
             with vectors lw 2 lc rgb "black" size screen 0.02,10