# --- Constants & Configuration ---
    x_min = -11
    x_max = 220

    y_min = -25.5
    y_max = 25.5

    # columns
    X = 1           # X centroid
    Y = 2           # Y centroid
    P = 3
    UX = 4
    UY = 5
    U_MOD = 6
    K = 7
    EPS = 8
    NU_TURB = 9

    array names[9] = ["x centroid", "y centroid", "Δp", "Δu_x", "Δu_y", "|Δu|", "Δk", "Δε", "Δν_t"]

    array file_names[9] = ["x_centroid", "y_centroid", "pressure", "velocity_x", "velocity_y", "velocity", "k", "epsilon", "nu_t"]

    file = "comparison"
    data = "../data/" . file . ".dat"

do for [COL=P:NU_TURB] {
    # --- Data Statistical Analysis ---
        stats data using COL prefix "col" nooutput

    # --- Terminal & Output Settings ---
        set term pngcairo enhanced font "Helvetica Bold,24" size 1200, 300
        set output "../visualisations/" . file . "_" . file_names[COL] . ".png"

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

        set xtics 40 offset 0,0
        set ytics 25 offset 0.5,0
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
        max_abs = (abs(col_min) > abs(col_max)) ? abs(col_min) : abs(col_max)

        v(gray) = col_min + gray * (col_max - col_min)  # colour function

        set palette model RGB functions \
            (v(gray) > 0 ? 1.0 : 1.0 - abs(v(gray))/max_abs), \
            (1.0 - abs(v(gray))/max_abs), \
            (v(gray) > 0 ? 1.0 - abs(v(gray))/max_abs : 1.0)

    # --- Labels, Titles & Colorbar ---
        set xlabel "X (mm)" font "Helvetica Bold,32" tc rgb "#111111" offset 0,0
        set ylabel "Y (mm)" font "Helvetica Bold,32" tc rgb "#111111" offset 0,0
        set cblabel names[COL] font "Helvetica Bold,32" tc rgb "#111111" offset 3.5,0

        # add black bw_step
        set object 1 rectangle from x_min, y_min to 0, 0 fillstyle solid noborder fillcolor rgb "black" front

    # --- Legend / Key ---
        set key off
        # unset xtics
        # unset ytics
        # unset border

    # --- Plot Command ---
        plot data using X:Y:COL with points pt 5 palette
}