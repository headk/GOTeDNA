#' Displays heatmap of scaled monthly detection probabilities for taxon
#'
#' @description This function displays scaled species detection probabilities
#' of the chosen taxonomic group with probabilities produced with
#' [scale_newprob()]. All primers that are available in the data are
#' displayed.
#' NOTE: interpolated `(i.e., missing)` data are not used in this
#' representation.
#'
#' @param scaledprobs (required, data.frame) Normalized detection
#' probabilities as returned by [scale_newprob()].
#'
#' @author Anais Lacoursiere-Roussel \email{Anais.Lacoursiere@@dfo-mpo.gc.ca}
#' @rdname hm_fig
#' @export
#' @examples
#' \dontrun{
#' newprob <- calc_det_prob(D_mb)
#' scaledprobs <- scale_newprob(D_mb, newprob)
#' hm_fig(scaledprobs)
#' }
hm_fig <- function(scaledprobs, selected_taxon_level = "species") {
  # Preprocess data
  df <- scaledprobs %>%
    filter(!is.na(year)) %>%
    group_by(year, month, .data[[selected_taxon_level]]) %>%   # note: use .data here
    summarise(scaleP = mean(scaleP, na.rm = TRUE), .groups = "drop") %>%
    rename("Detection rate" = "scaleP") %>%
    mutate(
      taxon = as.character(.data[[selected_taxon_level]]),     # works now
      Month = factor(month, levels = 1:12, labels = month.abb)
    )

  taxa_list <- unique(df[[selected_taxon_level]])
  n_taxa <- length(taxa_list)
  plots <- lapply(seq_along(taxa_list), function(i) {
    sp <- taxa_list[i]
    df_sp <- df %>% filter(taxon == sp)
    df_sp$zero_layer <- ifelse(df_sp$`Detection rate` == 0, 0, NA)

    # Only show x-axis labels on the bottom subplot if <= 3 subplots
    show_xticks <- if (n_taxa > 3) TRUE else (i == n_taxa)

    plot_ly(
      df_sp,
      x = ~Month,
      y = ~factor(year, levels = rev(unique(year))),
      z = ~`Detection rate`,
      type = "heatmap",
      colors = viridis::viridis(100, direction = -1),
      zmin = 0.00001,
      zmax = 1,
      showscale = FALSE,
      na.color = "white"
    ) %>%
      add_trace(
        z = ~zero_layer,
        type = "heatmap",
        colorscale = list(c(0, "lightgrey"), c(1, "lightgrey")),
        showscale = FALSE,
        na.color = "transparent"
      ) %>%
      layout(
        xaxis = list(
          title = "",
          tickangle = -45,
          showgrid = FALSE,
          showticklabels = show_xticks
        ),
        yaxis = list(title = "", autorange = "reversed", showgrid = FALSE),
        margin = list(l = 50, r = 20, t = 50, b = 50),
        annotations = list(
          list(
            x = 0,
            y = 1,
            text = sp,
            xref = "paper",
            yref = "paper",
            xanchor = "left",
            yanchor = "bottom",
            showarrow = FALSE,
            font = list(size = 14, color = "black")
          )
        )
      )
  })


  n_taxa <- length(taxa_list)

  # Set base height per subplot
  height_per_subplot <- 275
  extra_spacing <- 30  # extra height per subplot for spacing

  # Total height scales with number of taxa
  total_height <- height_per_subplot * n_taxa + extra_spacing * max(0, n_taxa - 1)

  subplot(
    plots,
    nrows = n_taxa,
    shareX = FALSE,  # <- important
    shareY = FALSE,
    titleY = FALSE,
    margin = .03
  ) %>%
    layout(height = total_height)
  }

