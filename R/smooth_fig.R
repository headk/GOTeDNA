#' Display species monthly detection probability by primer and region.
#'
#' @description This function displays the species optimal detection period as
#' monthly detection probabilities with LOESS smoothing over years.
#'
#' @param data (required, data.frame): Data.frame read in with [read_data()].
#'
#' @author Anais Lacoursiere-Roussel \email{Anais.Lacoursiere@@dfo-mpo.gc.ca}
#' @rdname smooth_fig
#' @export
#' @examples
#' \dontrun{
#' data <- gotedna_data$metabarcoding |> dplyr::filter(
#'      species == "Acartia longiremis",
#'      primer == "COI1"
#'  )
#' smooth_fig(data = data)
#' }
smooth_fig <- function(scaledprobs) {
  # Use the month-level scaled probability directly
  data <- scaledprobs %>%
    dplyr::filter(!is.na(scaleP), !is.na(year)) %>%
    dplyr::mutate(month = as.numeric(month),
                  year = as.factor(year))  # ensure year is factor for coloring

  Dsummary24 <- Dsummary12 <- data
  Dsummary12$month <- data$month + 12
  Dsummary24$month <- data$month + 24

  Dsummary_comb <- rbind(data, Dsummary12, Dsummary24)

  loessmod <- loess(scaleP ~ month, Dsummary_comb, span = 3 / 12)

  NEW <- data.frame(month = seq(1, 35, 0.1))
  NEW$PRED <- predict(loessmod, newdata = NEW$month)

  NEW2 <- NEW[NEW$month > 12 & NEW$month <= 24, ]
  NEW2$month <- NEW2$month - 12

  ggplot2::ggplot() +
    ggplot2::geom_hline(ggplot2::aes(yintercept = y), data.frame(y = c(0:4)/4), color = "lightgrey") +
    ggplot2::geom_vline(ggplot2::aes(xintercept = x), data.frame(x = 0:12), color = "lightgrey") +
    ggplot2::geom_path(data = NEW2, ggplot2::aes(x = month, y = PRED), colour = "blue") +
    ggplot2::geom_point(data = data, ggplot2::aes(x = month, y = scaleP, col = year), alpha = 0.9, size = 5) +
    ggplot2::coord_polar(clip = "off") +
    ggplot2::labs(col = "Year", x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::scale_colour_manual(values = palette("Alphabet")) +
    ggplot2::scale_x_continuous(limits = c(0.5, 12.5), breaks = 1:12, labels = month.abb) +
    ggplot2::guides(colour = ggplot2::guide_legend(label.position = "left", label.hjust = 1)) +
    ggplot2::scale_y_continuous(limits = c(-0.1, 1.01), breaks = c(0, 0.25, 0.50, 0.75, 1)) +
    theme_circle
}



