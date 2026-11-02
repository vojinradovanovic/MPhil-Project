# --- Constants & Configuration ---
    x_min = -11
    x_max = 220

    y_min = -25.5
    y_max = 25.5

    # columns
    X = 1           # X centroid
    Y = 2           # Y centroid
    P = 3
    UX_START = 4    # vector origin
    UY_START = 5
    UX_LEN = 6      # vector arrow lengths
    UY_LEN = 7
    UX = 8          # velocity component values
    UY = 9
    U_MOD = 10
    K = 11
    EPS = 12
    NU_TURB = 13

    # other
    HEAT_PALETTE = 1
    RB_PALETTE = 2

    NO_ARROWS = 1
    WITH_ARROWS = 2

    array names[13] = ["x centroid", "y centroid", "pressure", "ux start coordinate", "uy start coordinate", \
                        "ux arrow length", "uy arrow length", "u_x", "u_y", "|u| value", "k", "ε", "ν_t"]

    array file_names[13] = ["x_centroid", "y_centroid", "pressure", "ux_start_coordinate", "uy_start_coordinate", \
                            "ux_length", "uy_length", "ux", "uy", "u_mod", "k", "epsilon", "nu_t"]

    array plotting_columns[7] = [P, UX, UY, U_MOD, K, EPS, NU_TURB]
    array plotting_palettes[7] = [HEAT_PALETTE, RB_PALETTE, RB_PALETTE, HEAT_PALETTE, HEAT_PALETTE, HEAT_PALETTE, HEAT_PALETTE]
    array plotting_arrows[7] = [NO_ARROWS, NO_ARROWS, NO_ARROWS, WITH_ARROWS, NO_ARROWS, NO_ARROWS, NO_ARROWS]

    file = "turbulence_end"
    data = "../data/" . file . ".dat"

do for [i=1:7] {
    # --- Set plotting column, plot and colour palette ---
        COL = plotting_columns[i]
        PALETTE = plotting_palettes[i]
        ARROWS = plotting_arrows[i]

    # --- Data Statistical Analysis ---
        set xrange [*:*]
        set yrange [*:*]
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
        set format cb "%.3f"

    # --- Layout Margin Adjustments ---
        set rmargin 10
        set lmargin 10
        set bmargin 4

    # --- Axis Ranges ---
        set xrange [x_min:x_max]
        set yrange [y_min:y_max]
        set cbrange [col_min:col_max]

    # --- Color Palette ---
        if(PALETTE == HEAT_PALETTE) {

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

            unset object 1

        } else if(PALETTE == RB_PALETTE) {

            max_abs = (abs(col_min) > abs(col_max)) ? abs(col_min) : abs(col_max)

            v(gray) = col_min + gray * (col_max - col_min)  # colour function

            set palette model RGB functions \
                (v(gray) > 0 ? 1.0 : 1.0 - abs(v(gray))/max_abs), \
                (1.0 - abs(v(gray))/max_abs), \
                (v(gray) > 0 ? 1.0 - abs(v(gray))/max_abs : 1.0)

            # add black bw_step
                set object 1 rectangle from x_min, y_min to 0, 0 fillstyle solid noborder fillcolor rgb "black" front

        }

    # --- Labels, Titles & Colorbar ---
        set xlabel "X (mm)" font "Helvetica Bold,32" tc rgb "#111111" offset 0,0
        set ylabel "Y (mm)" font "Helvetica Bold,32" tc rgb "#111111" offset 0,0
        set cblabel names[COL] font "Helvetica Bold,32" tc rgb "#111111" offset 3.5,0

    # --- Legend / Key ---
        set key off
        # unset xtics
        # unset ytics
        # unset border

    # --- Plot Command ---
        if(ARROWS == WITH_ARROWS) {

            plot data using X:Y:COL with points pt 5 palette, \
                data using UX_START:UY_START:UX_LEN:UY_LEN every 31 \
                    with vectors lw 2 lc rgb "blue" size screen 0.02,10

        } else if(ARROWS == NO_ARROWS) {

            plot data using X:Y:COL with points pt 5 palette

        }
}