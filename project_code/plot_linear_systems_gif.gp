# ---------- helper functions ----------
function $make_range(start, stop) << EOD
    array arr[stop-start]
    do for [i=start:(stop-1)] { arr[i-start+1] = i }
    return arr
EOD

function $sub(str, prev_expr, expr) << EOD
    array temp = split(str, prev_expr)
    return join(temp, expr)
EOD

function $plot(file, col_vars, prefix, output_suffixes, legend_suffixes, flip_plot, plot_bounds) << EOD
    data = "../data/" . file . ".dat"

    # initialisations
    legend_prefix = prefix
    output_prefix = $sub(legend_prefix, "_", "")

    # read columns
    do for [i = 1 : |col_vars|] {

        # statistics
        stats data using col_vars[i] nooutput

        if(STATS_max == 0) {
            STATS_max = 1
        }

        param = 0
        if(STATS_min>0) {
            param = STATS_min
        }

        # set up plot
        set term gif animate delay 200
        set output "../visualisations/linear_systems/" . output_prefix . "_" . output_suffixes[i] . ".gif"
        set xrange [plot_bounds[1]:plot_bounds[2]]
        set yrange [plot_bounds[3]:plot_bounds[4]]
        set key tmargin right

        if(flip_plot == 1) {
            flip = -1
        }
        else {
            flip = 1
        }

        do for [t_n=0:(frames-1):steps] {
            plot data using X:( flip * column(Y) ):col_vars[i] \
                    index t_n with image title (legend_prefix . " " . legend_suffixes[i])
        }

        unset xrange
        unset yrange
    }

    return 0
EOD

# ---------- global constants ----------
    num_vars = 3
    num_plot_bounds = 4

    array prefixes[num_vars] = ["v_x", "v_y", "p"]

    # global columns
    X = 1
    Y = 2

    frames = 1
    steps = 1

# ---------- plot velocity and pressure matrices ----------
    # constants
    x_id_min = 0
    x_id_max = 99

    y_id_min = -99
    y_id_max = 0

    num_plots = 1

    # columns
    array VX_COLS[num_plots] = [3]
    array VY_COLS[num_plots] = [4]
    array P_COLS[num_plots] = [3]

    array files[num_vars] = ["v_matrix_data", "v_matrix_data", "p_matrix_data"]

    array output_file_suffixes[num_plots] = ["matrix"]
    array legend_title_suffixes[num_plots] = ["matrix"]

    flip_plot = 1
    array plot_bounds[num_plot_bounds] = [x_id_min, x_id_max, y_id_min, y_id_max]

    # plot matrix data
    evaluate $plot(files[1], VX_COLS, prefixes[1], output_file_suffixes, legend_title_suffixes, flip_plot, plot_bounds)
    evaluate $plot(files[2], VY_COLS, prefixes[2], output_file_suffixes, legend_title_suffixes, flip_plot, plot_bounds)
    evaluate $plot(files[3], P_COLS, prefixes[3], output_file_suffixes, legend_title_suffixes, flip_plot, plot_bounds)

# ---------- plot velocity and pressure vectors ----------
    # constants
    x_min = 0
    x_max = 0.1

    y_min = 0
    y_max = 0.1

    num_plots = 4

    # columns
    VX_COLS = $make_range(3+0   *num_plots, 3+1 *num_plots)
    VY_COLS = $make_range(3+1   *num_plots, 3+2 *num_plots)
    P_COLS = $make_range(3, 3+num_plots)

    array files[num_vars] = ["v_vector_data", "v_vector_data", "p_vector_data"]

    array legend_title_suffixes[num_plots] = ["diagonal coefficients", "rhs vector", "residuals before", "residuals after"]
    array output_file_suffixes[num_plots]
    do for [i=1:|output_file_suffixes|] { output_file_suffixes[i] = $sub(legend_title_suffixes[i], " ", "_") }

    flip_plot = 0
    array plot_bounds[num_plot_bounds] = [x_min, x_max, y_min, y_max]

    # plot vector data
    evaluate $plot(files[1], VX_COLS, prefixes[1], output_file_suffixes, legend_title_suffixes, flip_plot, plot_bounds)
    evaluate $plot(files[2], VY_COLS, prefixes[2], output_file_suffixes, legend_title_suffixes, flip_plot, plot_bounds)
    evaluate $plot(files[3], P_COLS, prefixes[3], output_file_suffixes, legend_title_suffixes, flip_plot, plot_bounds)