# =============================================================================
# R/00_io.R
#
# Helpers for consistent figure and table saving, dual format exports
# (vector PDF + 300 DPI PNG), and automatic manifest generation.
# =============================================================================

.manifest_env <- new.env(parent = emptyenv())
.manifest_env$entries <- list()

#' Save a ggplot figure in both PDF and PNG formats
#'
#' @param plot_obj A ggplot or patchwork object.
#' @param base_name Primary filename (without extension), e.g. "Fig1_perturbation_mechanics".
#' @param alt_names Optional character vector of additional alias filenames to save
#'   (e.g., "fig_perturbation_mechanics" to match LaTeX \includegraphics directly).
#' @param caption One-line descriptive caption for the manifest.
#' @param width Width in inches.
#' @param height Height in inches.
#' @param supplementary Logical; if TRUE, saves to output/figures/supplementary.
save_figure <- function(plot_obj, base_name, alt_names = NULL, caption = "",
                        width = 8, height = 5, supplementary = FALSE) {
  target_dir <- if (supplementary) DIR_FIG_SUPP else DIR_FIG_MAIN

  all_names <- unique(c(base_name, alt_names))

  for (nm in all_names) {
    pdf_path <- file.path(target_dir, paste0(nm, ".pdf"))
    png_path <- file.path(target_dir, paste0(nm, ".png"))

    ggplot2::ggsave(pdf_path, plot = plot_obj, width = width, height = height,
                    device = "pdf", bg = "white")
    ggplot2::ggsave(png_path, plot = plot_obj, width = width, height = height,
                    dpi = 300, bg = "white")
  }

  .manifest_env$entries[[length(.manifest_env$entries) + 1]] <- list(
    type        = if (supplementary) "Supplementary Figure" else "Main Figure",
    file        = paste0(base_name, ".pdf"),
    path        = file.path(target_dir, paste0(base_name, ".pdf")),
    caption     = caption
  )

  invisible(plot_obj)
}

#' Save a data frame as a CSV table
#'
#' @param df Data frame or tibble to save.
#' @param base_name Filename without extension (e.g. "Table1_baseline_results").
#' @param caption One-line descriptive caption for the manifest.
#' @param supplementary Logical; if TRUE, saves to output/tables/supplementary.
save_table <- function(df, base_name, caption = "", supplementary = FALSE) {
  target_dir <- if (supplementary) DIR_TAB_SUPP else DIR_TAB_MAIN
  csv_path   <- file.path(target_dir, paste0(base_name, ".csv"))

  readr::write_csv(df, csv_path)

  .manifest_env$entries[[length(.manifest_env$entries) + 1]] <- list(
    type        = if (supplementary) "Supplementary Table" else "Main Table",
    file        = paste0(base_name, ".csv"),
    path        = csv_path,
    caption     = caption
  )

  invisible(df)
}

#' Write manifest Markdown file summarizing all generated artifacts
write_manifest <- function(output_file = file.path(DIR_OUT, "manifest.md")) {
  entries <- .manifest_env$entries
  if (length(entries) == 0) return(invisible(NULL))

  lines <- c(
    "# Reproduction Artifacts Manifest",
    "",
    "Generated automatically by `run_all.R`.",
    "",
    "| Type | File | Description |",
    "| :--- | :--- | :--- |"
  )

  for (e in entries) {
    rel_path <- gsub("\\\\", "/", e$path)
    lines <- c(lines, sprintf("| %s | [`%s`](%s) | %s |", e$type, e$file, rel_path, e$caption))
  }

  lines <- c(lines, "", paste0("*Total artifacts generated: ", length(entries), "*"))
  writeLines(lines, output_file)
  invisible(output_file)
}
