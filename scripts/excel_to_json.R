# Author: Emile Johnston
# University: Technische Universität Wien
# Date: 23 July 2025

# This R script converts the case vs control data from an Excel file into JSON format for
# easier data manipulation.

# Use pacman to load required packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(readxl, jsonlite, tools)

# This function reads a specific sheet from an Excel file, makes a case array and a control
# array for each gene, and returns the list of arrays.
# case values are above the 'control' row and control values are below the 'control' row.
read_sheet <- function(file_path, sheet = 1) {
    # Read the full sheet without setting column names
    data <- read_excel(file_path, sheet = sheet, col_names = FALSE)

    # Find the first row that contains "control"
    control_row <- which(apply(data, 1, function(row) any(grepl("control", row, ignore.case = TRUE))))[1]

    if (is.na(control_row)) {
        stop("No 'control' row found in the sheet.")
    }

    # Check if row 2 contains any strings
    second_row <- unlist(data[2, ])
    non_na_values <- second_row[!is.na(second_row)]
    cat("non_na_values in row 2:", non_na_values, "\n")
    numeric_check <- suppressWarnings(as.numeric(non_na_values))
    cat("numeric_check for row 2:", numeric_check, "\n")
    row2_has_string <- !isTRUE(all(!is.na(numeric_check)))
    cat("Row 2 has string values:", row2_has_string, "\n")

    # Initialize column_names
    column_names <- c()

    if (row2_has_string) { # only clinical cases
        # Construct composite column names from row 1 and row 2
        first_row <- as.character(unlist(data[1, ]))
        second_row <- as.character(unlist(data[2, ]))

        current_group <- NA # the gene, P4HA1 or SPP1

        # For each phase
        for (i in seq_along(first_row)) {
            top <- first_row[i]
            bottom <- second_row[i]

            # Update current group if top cell is not empty
            if (!is.na(top) && top != "") {
                current_group <- top
            }

            # Only add if current group is defined and bottom is not empty
            if (!is.na(current_group) && !is.na(bottom) && bottom != "") {
                column_names <- c(column_names, paste0(current_group, "  ", bottom))
            }
        }

        # Remove first two rows (they were headers)
        data <- data[-c(1, 2), ]

        # Remove columns that are completely empty (all NA or blank "")
        cat("Shape before removing empty columns:", dim(data), "\n")
        data <- data[, colSums(!is.na(data) & data != "") > 0]
        cat("Shape after removing empty columns:", dim(data), "\n")
    } else { # all 5 other sheets
        # Use first row as column names
        column_names <- as.character(unlist(data[1, ]))
        data <- data[-1, ]
    }

    # Recalculate control row in updated data
    control_row <- which(apply(data, 1, function(row) any(grepl("control", row, ignore.case = TRUE))))[1]

    if (is.na(control_row)) {
        stop("No 'control' row found in the sheet after removing headers.")
    }

    # Reset row indices
    row_count <- nrow(data)

    # Initialize result list
    result <- list()

    col_index <- 1

    # For each gene
    for (col_name in column_names) {
        # Skip columns that no longer exist due to filtering
        if (col_index > ncol(data)) break

        col_data <- data[[col_index]]

        # Extract and clean case and control values
        case_values <- col_data[1:(control_row - 1)]
        case_values <- as.numeric(case_values[!is.na(case_values)])

        control_values <- col_data[(control_row + 1):row_count]
        control_values <- as.numeric(control_values[!is.na(control_values)])

        # Store in result list
        result[[paste0(col_name, "_case")]] <- case_values
        result[[paste0(col_name, "_control")]] <- control_values

        col_index <- col_index + 1
    }
    return(result)
}

# This function iterates over all sheets in the excel file, for each it uses the
# read_sheet function to create a list of arrays, and saves a json file that contains this
# list.
process_all_sheets_to_json <- function(file_path, output_dir = "data/json") {
    # Get all sheet names
    sheets <- excel_sheets(file_path)

    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }

    # Loop through each sheet
    for (sheet_name in sheets) {
        cat("Processing sheet:", sheet_name, "\n")

        # Get the list of arrays from the sheet
        result <- tryCatch({
            read_sheet(file_path, sheet = sheet_name)
        }, error = function(e) {
            warning(paste("Skipping sheet", sheet_name, "due to error:", e$message))
            return(NULL)
        })

        if (!is.null(result)) {

            # Build JSON file path
            json_file <- file.path(output_dir, paste0(sheet_name, ".json"))

            # Write to JSON
            write_json(result, path = json_file, pretty = TRUE, auto_unbox = TRUE)
        }
    }
}

process_all_sheets_to_json('data/case_vs_control.xlsx')