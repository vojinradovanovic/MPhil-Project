# --- Constants & Configuration ---
    x_min = 0
    x_max = 0.4

    y_min = 0
    y_max = 0.1

    # Columns
    X        = 1    # X centroid
    Y        = 2    # Y centroid
    P        = 3
    VX_START = 4    # vector origin
    VY_START = 5
    VX_MOD   = 6    # vector length coordinates
    VY_MOD   = 7
    V_MOD    = 8    # Velocity Modulus (Solution Variable)

    file = "advection_diffusion"
    data = "../../data/" . file . ".dat"

# --- Data Statistical Analysis ---
    stats data using V_MOD index 0 prefix "stats" nooutput
    v_mod_max = stats_max
    v_mod_min = 0

# --- Terminal & Output Settings ---
    set term pngcairo enhanced font "Helvetica Bold,24" size 1200, 300
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

    cb_step = (v_mod_max - v_mod_min) / 2
    cb_step = real(sprintf("%.1g", cb_step))
    set cbtics cb_step

    set xtics 0.1 offset 0,-0.5
    set ytics 0.05 offset -0.5,0
    set format cb "%.3g"

# --- Layout Margin Adjustments ---
    set rmargin 5
    set lmargin 10
    set bmargin 4

# --- Axis Ranges ---
    set xrange [x_min:x_max]
    set yrange [y_min:y_max]
    set cbrange [0:v_mod_max]

# --- Color Palette ---
    set palette defined ( \
        0  "#000004", \
        1  "#160b39", \
        2  "#420a68", \
        3  "#6a1c6e", \
        4  "#932667", \
        5  "#bc3754", \
        6  "#dd513a", \
        7  "#f37c1b", \
        8  "#f9a620", \
        9  "#f6d142", \
        10 "#fcffa4"  \
    )

# --- Labels, Titles & Colorbar ---
    set xlabel "X" font "Helvetica Bold,32" tc rgb "#111111" offset 0,-0.5
    set ylabel "Y" font "Helvetica Bold,32" tc rgb "#111111" offset -0.5,0
    set cblabel "|{/Helvetica-BoldItalic u}| value" font "Helvetica Bold,32" tc rgb "#111111" offset 2.5,0

# --- Legend / Key ---
    set key off

# --- Plot Command ---
    plot data using X:Y:V_MOD index 0 with image, \
         data using VX_START:VY_START:VX_MOD:VY_MOD index 0 \
             with vectors lw 2 lc rgb "black" size screen 0.02,10