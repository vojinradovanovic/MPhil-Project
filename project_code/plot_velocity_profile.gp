# --- Constants & Configuration ---
    x_min = 0
    x_max = 220

    y_min = -25.5
    y_max = 25.5

    H = 25.5
    u_inlet = 13.3
    cell_size_x = 11./5
    plot_tolerance = cell_size_x/2.1
    scale = 18
    n_sampling_points = 9
    array sampling_points[n_sampling_points] = [0.5, 1, 2, 3, 4, 5, 6, 7, 8]

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

    array names[13] = ["x centroid", "y centroid", "pressure", "ux start coordinate", "uy start coordinate", \
                        "ux arrow length", "uy arrow length", "u_x", "u_y", "|u| value", "k", "ε", "ν_t"]

    array file_names[13] = ["x_centroid", "y_centroid", "pressure", "ux_start_coordinate", "uy_start_coordinate", \
                            "ux_length", "uy_length", "ux", "uy", "u_mod", "k", "epsilon", "nu_t"]

    file = "turbulence_end"
    data = "../data/" . file . ".dat"

# --- Set plotting column ---
    COL = UX

    # display plotting details
    print("X profile for column " . file_names[COL])
    print sprintf("Inlet velocity (m/s): %.1f",  u_inlet)
    print sprintf("Step height (mm): %.1f",  H)
    print sprintf("Cell length (mm): %.1f", cell_size_x)

# --- Data Statistical Analysis ---
    stats data using COL prefix "col" nooutput

# --- Terminal & Output Settings ---
    set term pngcairo enhanced font "Helvetica Bold,24" size 1400, 400
    set output "../visualisations/" . "x_profile" . "_" . file_names[COL] . ".png"

# ==============================================================================
# MULTIPLOT CONFIGURATION
# ==============================================================================
    set multiplot

    # --- 1. Define Common Vertical Alignment & Right Panel Screen Position ---
    top_margin    = 0.90
    bottom_margin = 0.28

    gap = 0.035            # Gap between left and right plot boxes (3.5% of canvas width)

    right_lmargin = 0.16 + gap  # Left edge of right plot box (screen fraction 0.0 to 1.0)
    right_rmargin = 0.95        # Right edge of right plot box
    right_plot_width = right_rmargin - right_lmargin

    # --- 2. Calculate Exact Left Panel Width for 1:1 Scale Matching ---
    left_u_max = 1.2       # Max U/U_0 on left plot range
    right_domain_mm = x_max - x_min

    # Left box width needed so 1.0 unit of U/U_0 equals 'scale' mm on the right
    left_plot_width = right_plot_width * (scale * left_u_max) / right_domain_mm

    left_rmargin = 0.16
    left_lmargin = left_rmargin - left_plot_width

# ------------------------------------------------------------------------------
# PANEL 1: INLET PROFILE (LEFT)
# ------------------------------------------------------------------------------
    # Tics
    unset xtics
    unset ytics
    set xtics ("0" 0, "1" 1) font "Helvetica Bold,32"
    set ytics ("-1" -H, "0" 0, "1" H) font "Helvetica Bold,32"
    set grid xtics lt 1 lc rgb "black" lw 1

    # --- Layout Margin Adjustments ---
    # Lock inner plot box to exact screen coordinates
    set lmargin at screen left_lmargin
    set rmargin at screen left_rmargin
    set bmargin at screen bottom_margin
    set tmargin at screen top_margin

    # --- Axis Ranges for Normalized Velocity U/U_0 ---
    set xrange [0:left_u_max]
    set yrange [y_min:y_max]

    # --- Labels & Titles ---
    set xlabel "U/U_{inlet}" font "Helvetica Bold,32" tc rgb "#111111" offset 1.5,0
    set ylabel "Y/H" font "Helvetica Bold,32" tc rgb "#111111" offset 0,0

    # Plot inlet profile
    inlet_x = -1

    plot sprintf("< awk '$%d >= %f && $%d <= %f' %s | sort -k%d,%dn", \
                X, inlet_x - plot_tolerance, X, inlet_x + plot_tolerance, data, Y, Y) \
        using (column(COL)/u_inlet):Y \
        with lines lw 2.5 lc rgb "red" notitle


# ------------------------------------------------------------------------------
# PANEL 2: DOWNSTREAM PROFILES (RIGHT)
# ------------------------------------------------------------------------------
    # Tics
    unset xtics
    unset ytics
    set xtics ("0" 0, "" 0.5*H) font "Helvetica Bold,32"
    do for [i=1:n_sampling_points] {
        set xtics add (sprintf("%d", sampling_points[i]) sampling_points[i]*H)
    }
    set ytics ("" y_min+H)  # Draws tic marks without repeating numbers

    set grid xtics ytics lt 1 lc rgb "black" lw 1

    # --- Layout Margin Adjustments ---
    # Lock inner plot box to exact screen coordinates
    set lmargin at screen right_lmargin
    set rmargin at screen right_rmargin
    set bmargin at screen bottom_margin
    set tmargin at screen top_margin

    # --- Axis Ranges ---
    set xrange [x_min:x_max]
    set yrange [y_min:y_max]

    # --- Labels & Titles ---
    set xlabel "X/H" font "Helvetica Bold,32" tc rgb "#111111" offset 0,0
    unset ylabel                 # Y-axis label is already on the left panel

    # # add black bw_step
    # set object 1 rectangle from x_min, y_min to 0, 0 fillstyle solid noborder fillcolor rgb "black" front

    # --- Legend / Key ---
    set key off
    # unset border

    # --- Plot Command ---
    plot for [i=1:n_sampling_points] \
            sprintf("< awk '$%d >= %f && $%d <= %f' %s | sort -k%d,%dn", \
                    X, \
                    sampling_points[i]*H - plot_tolerance, \
                    X, \
                    sampling_points[i]*H + plot_tolerance, \
                    data, \
                    Y, \
                    Y) \
            using ((sampling_points[i]*H) + scale * (column(COL)/u_inlet)):Y \
                with lines lw 2.5 lc rgb "red" notitle

                # with linespoints lw 1 lc rgb "red" ps .5 pt 7 notitle

# ==============================================================================
# END MULTIPLOT
# ==============================================================================
    unset multiplot