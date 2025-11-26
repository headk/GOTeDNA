#' Normalize detection probabilities with min-max method (range 011) by month
#' for each year separately.
#'
#' @description This function normalizes (i.e., scales) monthly detection
#' probabilities for each selected taxon (e.g. species), primer, and year that were calculated with
#' the previous function, [calc_det_prob()]. Outputs fed into figure and window
#' calculation functions.
#' * NOTE: Currently this function only works for metabarcoding data.
#'
#' @param data (required, data.frame) Data.frame imported with [read_data()]. Required
#' to join taxonomic information.
#' @param newprob (required, list) detection probabilities aggregated per month and year
#' [calc_det_prob()].
#'
#' @return Grouped data.frame with 15 columns (though lower taxonomic levels are not included beyond the selected taxonomic level):
#' * `id`: unique taxon;primer;year identifier
#' * `month`:
#' * `detect`: number of detections
#' * `nondetect`: number of non-detections
#' * `scaleP`: detection probability scaled to range 0-1
#' * `protocol_ID`:
#' * `species`:
#' * `primer`:
#' * `year`:
#' * `phylum`:
#' * `class`:
#' * `order`:
#' * `family`:
#' * `genus`:
#' * `fill`:
#'
#' @author Anais Lacoursiere-Roussel \email{Anais.Lacoursiere@@dfo-mpo.gc.ca}
#' @rdname scale_newprob
#' @export
#' @examples
#' \dontrun{
#' newprob <- calc_det_prob(gotedna_data$metabarcoding, "Scotian Shelf")
#' scale_newprob(data = gotedna_data$metabarcoding, newprob)
#' }
scale_newprob <- function(data, newprob, selected_taxon_level = "species") {
  CPscaled <- lapply(newprob, function(species_list) {

  # 1. Extract the maximum p across all elements for this species
  all_p <- unlist(lapply(species_list, function(df) df$p))
  max_p <- max(all_p, na.rm = TRUE)

  # 2. Apply minimal change: scale by max_p instead of scale_prop
  lapply(species_list, function(y) {
    data.frame(y) |>
      dplyr::mutate(
        scaleP = dplyr::case_when(
          p == 1 ~ 1,
          p == 0 ~ 0,
          TRUE   ~ p / max_p
        )
      )
    })
  })
  DFmo <- lapply(CPscaled$newP_agg, function(x) {
    out <- data.frame(
      month = 1:12,
      detect = NA_integer_,
      nondetect = NA_integer_,
      scaleP = NA_real_
    )
    out$detect[x$month] <- x$nd
    out$nondetect[x$month] <- x$n - x$nd
    out$scaleP[x$month] <- x$scaleP
    out
  }) |>
    do.call(what = rbind) |>
    dplyr::mutate(
      id = rep(names(CPscaled$newP_agg), each = 12)
    ) |>
    dplyr::select(id, month, detect, nondetect, scaleP) |>
    dplyr::tibble()
  row.names(DFmo) <- NULL

  DFmo[c("protocol_ID", selected_taxon_level, "primer")] <- stringr::str_split_fixed(DFmo$id, ";", 3)

  taxa_columns <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  cols_to_keep <- taxa_columns[1:which(taxa_columns == selected_taxon_level)]

  DFmo <- DFmo |>
    dplyr::left_join(
      unique(data[, cols_to_keep]),
      by = selected_taxon_level,
      multiple = "first"
    )

  # Interpolate missing months
  DFmo$det_int <- NA
  DFmo$nd_int <- NA
  DFmo$fill <- NA # add column

  for (taxon in unique(DFmo$id)) {
    DF1 <- DFmo[DFmo$id == taxon, ]

    # then add code for interpolation that starts with DF2 = .....
    # add dataframe above and below to help will fills for jan and dec. Needed to have 4 copyies because of max function used below
    DF2 <- rbind(
      cbind(DF1, data.frame(G = 1)),
      cbind(DF1, data.frame(G = 2)),
      cbind(DF1, data.frame(G = 3)),
      cbind(DF1, data.frame(G = 4))
    )

    DF2$det_int <- DF2$detect
    DF2$nd_int <- DF2$nondetect
    DF2$fill <- DF2$scaleP

    # which months are NA and define groups with sequential NAs
    month_na_id <- which(is.na(DF2$detect))
    nagroups <- cumsum(c(1, abs(month_na_id[-length(month_na_id)] - month_na_id[-1]) > 1))

    # identify which NA groups are in G = 2 or 3 (ignore 1 and 2)
    nagroupsG <- list()
    for (i in unique(nagroups)) {
      nagroupsG[[i]] <- max(DF2$G[month_na_id[which(nagroups == i)]])
    }
    nagroupsGv <- unlist(nagroupsG)
    nagroupsf <- which(nagroupsGv %in% 2:3)

    # loop over final NA groups and fill in using average
    for (i in unique(nagroupsf)) {
      DF2$det_int[month_na_id[which(nagroups == i)]] <- (DF2$detect[min(
        month_na_id[which(nagroups == i)]) - 1] + DF2$detect[max(
          month_na_id[which(nagroups == i)]) + 1])/2

      DF2$nd_int[month_na_id[which(nagroups == i)]] <- (DF2$nondetect[min(
        month_na_id[which(nagroups == i)]) - 1] + DF2$nondetect[max(
          month_na_id[which(nagroups == i)]) + 1])/2

      DF2$fill[month_na_id[which(nagroups == i)]] <- (DF2$scaleP[min(
        month_na_id[which(nagroups == i)]) - 1] + DF2$scaleP[max(
          month_na_id[which(nagroups == i)]) + 1])/2
    }
    # then put values from DF3 back into DF$sp.pr. This assumes that the months are all in the correct order (jan to dec) in DF3 and test_interp
    # DF3 is final DF with fills
    DF3 <- DF2[DF2$G == 2, ]

    DFmo$det_int[DFmo$id == taxon] <- DF3$det_int
    DFmo$nd_int[DFmo$id == taxon] <- DF3$nd_int
    DFmo$fill[DFmo$id == taxon] <- DF3$fill
  }

  # Pscaled_month <- DFmo %>%
  #  dplyr::ungroup()

  # scale and interpolate each month separately
  DFyr <- lapply(CPscaled$newP_yr, function(x) {
    out <- data.frame(
      month = 1:12,
      detect = NA_integer_,
      nondetect = NA_integer_,
      scaleP = NA_real_
    )
    out$detect[x$month] <- x$nd
    out$nondetect[x$month] <- x$n - x$nd
    out$scaleP[x$month] <- x$scaleP
    out
  }) |>
    do.call(what = rbind) |>
    dplyr::mutate(
      id = rep(names(CPscaled$newP_yr), each = 12)
    ) |>
    dplyr::select(id, month, detect, nondetect, scaleP) |>
    dplyr::tibble()
  row.names(DFyr) <- NULL

  DFyr[c("protocol_ID", selected_taxon_level, "primer", "year")] <- stringr::str_split_fixed(DFyr$id, ";", 4)

  taxa_columns <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus", "species")
  cols_to_keep <- taxa_columns[1:which(taxa_columns == selected_taxon_level)]

  DFyr <- DFyr |>
    dplyr::left_join(
      unique(data[, cols_to_keep]),
      by = selected_taxon_level,
      multiple = "first"
    )

  # Interpolate missing months
  DFyr$fill <- NA # add column

  for (taxon in unique(DFyr$id)) {
    DF1 <- DFyr[DFyr$id == taxon, ]

    # then add code for interpolation that starts with DF2 = .....
    # add dataframe above and below to help will fills for jan and dec. Needed to have 4 copies because of max function used below
    DF2 <- rbind(
      cbind(DF1, data.frame(G = 1)),
      cbind(DF1, data.frame(G = 2)),
      cbind(DF1, data.frame(G = 3)),
      cbind(DF1, data.frame(G = 4))
    )

    DF2$fill <- DF2$scaleP

    # which months are NA and define groups with sequential NAs
    month_na_id <- which(is.na(DF2$scaleP))
    nagroups <- cumsum(c(1, abs(month_na_id[-length(month_na_id)] - month_na_id[-1]) > 1))

    # identify which NA groups are in G = 2 or 3 (ignore 1 and 2)
    nagroupsG <- list()
    for (i in unique(nagroups)) {
      nagroupsG[[i]] <- max(DF2$G[month_na_id[which(nagroups == i)]])
    }
    nagroupsGv <- unlist(nagroupsG)
    nagroupsf <- which(nagroupsGv %in% 2:3)

    # loop over final NA groups and fill in using average
    for (i in unique(nagroupsf)) {
      DF2$fill[month_na_id[which(nagroups == i)]] <- (DF2$scaleP[min(month_na_id[which(nagroups == i)]) - 1] + DF2$scaleP[max(month_na_id[which(nagroups == i)]) + 1]) / 2
    }
    # then put values from DF3 back into DF$id.sp.pr.yr. This assumes that the months are all in the correct order (jan to dec) in DF3 and test_interp
    # DF3 is final DF with fills
    DF3 <- DF2[DF2$G == 2, ]
    DFyr$fill[DFyr$id == taxon] <- DF3$fill
  }

  DFmo$year = NA
  scaledprobs = dplyr::bind_rows(DFmo, DFyr)
#  list(Pscaled_month = DFmo, Pscaled_year = DFyr)
  return(scaledprobs)
}

