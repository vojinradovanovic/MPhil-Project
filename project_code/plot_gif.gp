# --- Constants & Configuration ---
    x_min = 0
    x_max = 0.5

    y_min = 0
    y_max = 0.2

    # columns
    X = 1           # X centroid
    Y = 2           # Y centroid
    P = 3
    VX_START = 4    # vector origin
    VY_START = 5
    VX_LEN = 6      # vector length coordinates
    VY_LEN = 7
    V_MOD = 8
    K = 9
    EPS = 10
    NU_TURB = 11

    array names[11] = ["x centroid", "y centroid", "pressure", "vx start coordinate", "vy start coordinate", \
                        "vx value", "vy value", "|u| value", "k", "{/symbol e}", "{/symbol n}"]

    frames = 10
    steps = 1

    file = "turbulence_end"
    data = "../data/" . file . ".dat"

# --- Set plotting column ---
    COL = V_MOD

# --- Data Statistical Analysis ---
    stats data using COL name "col" nooutput

# --- Terminal & Output Settings ---
    set term gif animate delay 20 font "Helvetica Bold,24" size 1200, 300
    set output "../visualisations/" . file . ".gif"

# --- Plot Styles, Borders & Tics ---
    set border 15 lw 1.5 lc rgb "#333333"  
    set tics out nomirror scale 0.75       
    set mxtics 2                           
    set mytics 2

# --- Tick Label Font & Format Configuration ---
    set xtics font "Helvetica Bold,32"
    set ytics font "Helvetica Bold,32"
    set cbtics font "Helvetica Bold,32"

    cb_step = (col_max - col_min) / 2
    set cbtics cb_step

    set xtics 0.1 offset 0,-0.5
    set ytics 0.05 offset -0.5,0
    set format cb "%.5f"

# --- Layout Margin Adjustments ---
    set rmargin 10
    set lmargin 10
    set bmargin 4

# --- Axis Ranges ---
    set xrange [x_min:x_max]
    set yrange [y_min:y_max]
    set cbrange [col_min:col_max]

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
    set cblabel names[COL] font "Helvetica Bold,32" tc rgb "#111111" offset 3.5,0

# --- Legend / Key ---
    set key off

# --- Plot ---
    do for [t_i=0:frames:steps] {
        plot    data using X:Y:COL index t_i with points pt 5 palette, \
                data using VX_START:VY_START:VX_LEN:VY_LEN index t_i every 35  \
                    with vectors lw 2 lc rgb "black" size screen 0.02,10
    }