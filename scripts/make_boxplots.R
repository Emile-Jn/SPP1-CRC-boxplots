# Author: Emile Johnston
# University: Technische Universität Wien
# Date: 23 July 2025

# This R script makes boxplot figures from the json files in data/json/ for the paper.

# Use pacman to load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(ggplot2, jsonlite, tools)

#' Convert a p-value to significance stars
#'
#' Returns conventional asterisk notation for statistical significance.
#'
#' @param p A numeric p-value.
#' @return A character string: "", ".", "*", "**", or "***"
p_to_stars <- function(p) {
    if (p < 0.001) {
        return("***")
    } else if (p < 0.01) {
        return("**")
    } else if (p < 0.05) {
        return("*")
    } else if (p < 0.1) {
        return(".")
    } else {
        return("")  # not significant
    }
}

# diverse colours for the boxplots
base_colours <- c(
    "#1f77b4",  # Blue
    "#ff7f0e",  # Orange
    "#2ca02c",  # Green
    "#d62728",  # Red
    "#9467bd",  # Purple
    "#8c564b",  # Brown
    "#17becf"   # Cyan
)

#' Lighten a colour
#'
#' Function to lighten a colour, to visually differentiate case and control
#'
#' @param colour A character string representing a colour (e.g., "#RRGGBB").
#' @param factor A numeric factor by which to lighten the colour (default is 1.5).
#' @return A character string representing the lightened colour.
lighten_colour <- function(colour, factor = 1.5) {
    col_rgb <- col2rgb(colour)
    col_light <- pmin(255, col_rgb * factor)
    rgb(t(col_light), maxColorValue = 255)
}

#' Insert a skipline
#'
#' make names of clinical trial data (P4HA1 and SPP1) into 2 lines instead of 1 for better
#' plotting
#'
#' @param name A string, e.g. "P4HA1  I vs. IV_case"
#' @return the name with a skipline, e.g. "P4HA1\nI vs. IV_case"
insert_skipline <- function(name) {
    if (grepl("vs", name)) { # only clinical trials contain "vs"
        # replace the first space with a newline
        return(sub("1  ", "1\n", name, fixed = TRUE))
    } else {
        return(name)
    }
}


#' Make and save boxplot figures from a list of arrays.
#'
#' This function produces a figure in which each array in the given list is displayed as a
#' boxplot. The figures are not displayed but simply saved in plots/paper_plots/
#'
#' @param list A named list of numeric vectors, where each vector represents a case or
#' control group of data for one gene.
#' @param file_name The name under which to save the file, if `unified` is `TRUE`. If false,
#' it's the name of the folder where individual plots will be saved.
#' @param unified Whether the different genes in list are plotted in one figure (TRUE) or in
#' separate figures (FALSE)
#' @return nothing, creates directories if needed and saves plots.
boxplots <- function(list, file_name, unified) {

    n <- length(list)
    if (n < 2) {
        warning("List must contain at least two data elements.")
        return(NULL)
    }
    if (is.null(names(list))) {
        names(list) <- paste0("Group", seq_along(list))
    }
    num_pairs <- floor(n / 2)
    if (num_pairs < 1) {
        warning("The list needs at least one case-control pair.")
        return(NULL)
    }
    if (!unified) {
        dir.create(file.path("plots", "individual plots", file_name),
                   showWarnings = FALSE,
                   recursive = TRUE)
    } else {
        dir.create(file.path("plots", "paper_plots"), showWarnings = FALSE)
    }

    # Define spacing variables to control layout
    space_within_pair <- 0.6 # Controls distance between case/control boxes
    space_between_pairs <- 0.9 # Controls distance between different pairs
    box_width <- 0.5

    # Data preparation
    plot_data_list <- list()
    annotation_data_list <- list()
    x_axis_breaks <- c()
    all_colours <- c()
    current_x <- 1 # Starting position for the first box

    # For each gene
    for (i in 1:num_pairs) {
        case_idx <- (i * 2) - 1
        ctrl_idx <- i * 2
        pair_name <- sub("_case$", "", names(list)[case_idx]) # gene name

        group1 <- list[[case_idx]] # case group data
        group2 <- list[[ctrl_idx]] # control group data

        if (length(group1) == 0 || length(group2) == 0) {
            warning("Skipping pair '", pair_name, "' as one group is empty.")
            next
        }

        # 1. Calculate numeric x-positions for custom spacing
        x_case <- current_x
        x_ctrl <- current_x + space_within_pair
        x_center <- current_x + (space_within_pair / 2)
        current_x <- x_ctrl + space_between_pairs # Set start for the next pair

        # 2. Create data frame for the current pair
        df <- data.frame(
            value = c(group1, group2),
            group = factor(rep(c("case", "control"), c(length(group1), length(group2)))),
            pair = pair_name,
            x_pos = rep(c(x_case, x_ctrl), c(length(group1), length(group2)))
        )
        plot_data_list[[pair_name]] <- df

        # 3. Prepare data for annotations (stars and pair names)
        p_val <- t.test(group1, group2)$p.value
        annotation_data_list[[pair_name]] <- data.frame(
            x = x_center,
            stars = p_to_stars(p_val),
            name = insert_skipline(pair_name)
        )

        # 4. Store x-axis tick positions
        x_axis_breaks <- c(x_axis_breaks, c(x_case, x_ctrl))

        # 5. Define and store colours for this pair
        colour_case <- base_colours[i]
        colour_ctrl <- grDevices::adjustcolor(colour_case, alpha.f = 0.6)
        all_colours[paste(pair_name, "case", sep = ".")] <- colour_case
        all_colours[paste(pair_name, "control", sep = ".")] <- colour_ctrl
    }

    if (length(plot_data_list) == 0) {
        warning("No valid data to plot.")
        return(NULL)
    }

    # Plotting
    if (unified) {
        # Combine all data into single data frames
        plot_data <- do.call(rbind, plot_data_list)
        annot_data <- do.call(rbind, annotation_data_list)
        annot_data$vjust_name <- ifelse(grepl("\n", annot_data$name), 1.7, 3.5)

        # Create a unique fill identifier for each box
        plot_data$fill_id <- interaction(plot_data$pair, plot_data$group)

        # Define integer breaks for the y-axis gridlines
        y_data_range <- range(plot_data$value, na.rm = TRUE)
        y_grid_breaks <- seq(floor(y_data_range[1]), ceiling(y_data_range[2]), by = 1)

        # Create the plot
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_pos, y = value, group = x_pos)) +
            ggplot2::geom_boxplot(ggplot2::aes(fill = fill_id), width = box_width) +
            ggplot2::scale_fill_manual(values = all_colours) +

            ggplot2::geom_text(data = annot_data, ggplot2::aes(x = x, y = Inf, label = stars), vjust = 1.5, size = 8, inherit.aes = FALSE) +
            ggplot2::geom_text(data = annot_data, ggplot2::aes(x = x, y = Inf, label = name, vjust = vjust_name), fontface = "bold", size = 6, inherit.aes = FALSE) +

            # Use the calculated positions for the x-axis ticks and labels
            ggplot2::scale_x_continuous(breaks = x_axis_breaks, labels = rep(c("case", "control"), num_pairs)) +

            # Add extra space at the top of the plot for annotations
            ggplot2::scale_y_continuous(
                expand = ggplot2::expansion(mult = c(0.05, 0.32)),
                breaks = y_grid_breaks
            ) +

            ggplot2::labs(y = "Log2FC") +

            ggplot2::theme_classic(base_size = 14) +
            ggplot2::theme(
                legend.position = "none",
                axis.title.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(), # Remove x-axis ticks for a cleaner look
                # Add solid, thin, grey gridlines for the major y-axis breaks
                panel.grid.major.y = ggplot2::element_line(color = "grey85", linewidth = 0.3, linetype = "solid")
            )

        # Save the plot
        plot_file <- file.path("plots", "paper_plots", paste0(file_name, ".png"))
        plot_width_inches <- 1.3 * num_pairs + 1.1 # Dynamically adjust width
        ggplot2::ggsave(plot_file, p, width = plot_width_inches, height = 5, units = "in", dpi = 600)

    } else {
        for (name in names(plot_data_list)) {
            df_sub <- plot_data_list[[name]]
            annot_sub <- annotation_data_list[[name]]

            colour_vals <- c(
                "case" = all_colours[[paste(name, "case", sep = ".")]],
                "control" = all_colours[[paste(name, "control", sep = ".")]]
            )

            p_indiv <- ggplot2::ggplot(df_sub, ggplot2::aes(x = group, y = value, fill = group)) +
                ggplot2::geom_boxplot(width = 0.6) +
                ggplot2::scale_fill_manual(values = colour_vals) +
                ggplot2::ggtitle(annot_sub$name) +
                ggplot2::annotate("text", x = 1.5, y = max(df_sub$value) * 1.05, label = annot_sub$stars, size = 6) +
                ggplot2::theme_minimal(base_size = 14) +
                ggplot2::theme(
                    legend.position = "none",
                    axis.title.x = ggplot2::element_blank(),
                    panel.grid.major.x = ggplot2::element_blank()
                ) +
                ggplot2::ylab("Log2FC")

            plot_file <- file.path("plots", "individual plots", file_name, paste0(name, ".png"))
            ggplot2::ggsave(plot_file, p_indiv, width = 500, height = 600, units = "px", dpi = 150)
        }
    }
}

#' Make boxplots for all data
#'
#' Iterate over all json files, for each file, extract the list of arrays and produce the
#' boxplots by calling the function boxplots() .
#'
#' @param unified Whether to create a unified plot for all pairs (TRUE) or individual
#' plots for each pair (FALSE).
#' @return nothing, creates directories if needed and saves plots.
make_all_boxplots <- function(unified = TRUE) {
    json_files <- list.files(path = "data/json", pattern = "\\.json$", full.names = TRUE)

    if (length(json_files) == 0) {
        message("No JSON files found in the current directory.")
        return()
    }

    for (file in json_files) {
        message("Processing: ", file)

        # Read JSON and convert to list
        data <- fromJSON(file)

        if (!is.list(data)) {
            warning("Skipping file ", file, ": not a list.")
            next
        }

        # Ensure data is a named list of vectors (expected format for boxplots)
        if (is.null(names(data))) {
            names(data) <- paste0("Group", seq_along(data))
        }

        file_name <- file_path_sans_ext(basename(file))

        tryCatch({
            # boxplots(data, title = tools::file_path_sans_ext(file))
            boxplots(data, file_name, unified = unified)
        }, error = function(e) {
            warning("Error while processing ", file, ": ", e$message)
        })
    }
}

# Run this to create all boxplots for the paper.
make_all_boxplots()

# For reproducibility:
writeLines(capture.output(sessionInfo()), "session_info.txt")
