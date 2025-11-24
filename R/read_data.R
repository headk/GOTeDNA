#' *New proposal of this page*
#' Read and format metabarcoding metadata and data from OBIS 
#'
#' @description Function that reads for any data in OBIS that has a DNA Derived Data 
#' extension. Will compile occurrence (core) and DNA Derived Data (extension) files 
#' into a single dataframe. Data frames are then formatted appropriately for use in 
#' subsequence analysis and visualizations (e.g., date reformatting, merging metadata 
#' and data).
#' * Note that template column names are in the format of Darwin Core Archive
#' (DwC-A) using Darwin Core (DwC) data standards where possible.
#'
#' @return A tibble with 25 columns:
#' * `protocol_ID`
#' * `protocolVersion`
#' * `samp_name`
#' * `eventID`
#' * `primer`
#' * `species`
#' * `domain`
#' * `kingdom`
#' * `phylum`
#' * `class`
#' * `order`
#' * `family`
#' * `genus`
#' * `concentration`: provided when choose.method = "qPCR"
#' * `pcr_primer_lod` : provided when choose.method = "qPCR"
#' * `organismQuantity`: provided when choose.method = "metabarcoding"
#' * `date`
#' * ecodistrict`
#' * `LClabel` : Local Contexts label to denote First Nations data sovereignty
#' * `decimalLatitude`
#' * `decimalLongitude`
#' * `station`
#' * `year`
#' * `month`
#' * `detected`
#' * `msct` :logical, where minimum sequence copy threshold = 10
#' * `ownerContact` : email of data owner/steward
#' * `bibliographicCitation` : DOI reference, if applicable
#'
#' @author Anais Lacoursiere-Roussel \email{Anais.Lacoursiere@@dfo-mpo.gc.ca}
#' @rdname read_data
#' @export
#' @examples
#' \dontrun{
#' D_mb <- read_data(
#'  choose.method = "metabarcoding", path.folder = "./inst/testdata"
#' )
#' }

read_data <- function(dataset_ids, join_by = c("auto", "occurrenceID", "id")) {
  library(robis)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(purrr)

  dataset_ids <- as.character(dataset_ids)
  join_by <- match.arg(join_by)

  obis_list <- purrr::map(dataset_ids, function(ds) {

    message("Pulling OBIS dataset: ", ds)

    # ---- 1. Check that dataset has DNADerivedData extension ----
    ds_meta <- robis::dataset(datasetid = ds)
    exts <- tolower(unlist(ds_meta$extensions))

    if (!"dnaderiveddata" %in% exts) {
      warning("Dataset ", ds, " has no DNADerivedData extension; skipping.")
      return(NULL)
    }

    # ---- 2. Pull occurrence records with DNADerivedData ----
    rec <- robis::occurrence(
      datasetid     = ds,
      extensions    = "DNADerivedData",
      hasextensions = "DNADerivedData"
    )

    if (nrow(rec) == 0L) {
      warning("No occurrence records returned for dataset ", ds)
      return(NULL)
    }

    # We'll keep both occurrenceID and id around
    core_occ <- rec %>%
      distinct(occurrenceID, .keep_all = TRUE)
    core_id <- rec %>%
      distinct(id, .keep_all = TRUE)

    # DNADerivedData extension (includes `id` by default; keep occurrenceID too if present)
    dna_only <- robis::unnest_extension(rec, "DNADerivedData")

    # ---- 3. Decide how to join core + extension ----
    join_choice <- join_by
    if (join_choice == "auto") {
      can_occ <- "occurrenceID" %in% names(core_occ) &&
        "occurrenceID" %in% names(dna_only) &&
        any(!is.na(core_occ$occurrenceID) & !is.na(dna_only$occurrenceID))

      can_id  <- "id" %in% names(core_id) &&
        "id" %in% names(dna_only) &&
        any(!is.na(core_id$id) & !is.na(dna_only$id))

      if (can_occ) {
        join_choice <- "occurrenceID"
      } else if (can_id) {
        join_choice <- "id"
      } else {
        stop("Neither occurrenceID nor id can be used to join core and DNADerivedData for dataset ",
             ds, ".")
      }
    }

    if (join_choice == "occurrenceID") {
      core_dna <- core_occ %>%
        left_join(dna_only, by = "occurrenceID")
    } else {  # "id"
      core_dna <- core_id %>%
        left_join(dna_only, by = "id")
    }

    core_dna <- core_dna %>%
      filter(!is.na(samp_name)) %>%  # keep only real samples
      mutate(
        datasetID_obis = ds,
        # eventDate usually ISO; strip time & parse
        eventDate_chr   = as.character(eventDate),
        eventDate_clean = suppressWarnings(
          ymd(substr(eventDate_chr, 1, 10))
        ),
        year  = year(eventDate_clean),
        month = month(eventDate_clean),
        decimalLatitude  = suppressWarnings(as.numeric(decimalLatitude)),
        decimalLongitude = suppressWarnings(as.numeric(decimalLongitude))
      )

    # Safe versions if fields are missing
    core_dna <- core_dna %>%
      mutate(
        station = if ("samplingStation" %in% names(.)) samplingStation else NA_character_,
        ownerContact = if ("ownerInstitutionCode" %in% names(.)) ownerInstitutionCode else NA_character_,
        bibliographicCitation = if ("bibliographicCitation" %in% names(.)) bibliographicCitation else NA_character_
      )

    # ---- 4. Metabarcoding detection + primer ----
    core_dna <- core_dna %>%
      mutate(
        organismQuantity = suppressWarnings(as.numeric(organismQuantity)),
        detected = case_when(
          !is.na(organismQuantity) & organismQuantity > 0 ~ 1L,
          TRUE ~ 0L
        ),
        primer = dplyr::coalesce(target_subfragment, target_gene)
      )

    # ---- 5. Return in GOTeDNA_df-like shape ----
    out <- core_dna %>%
      transmute(
        protocol_ID           = NA_character_,
        protocolVersion       = NA_real_,
        samp_name             = as.character(samp_name),
        primer                = paste (target_gene, pcr_primer_name_forward, pcr_primer_name_reverse, sep = "|"),
        scientificName        = scientificName,
        kingdom               = kingdom,
        phylum                = phylum,
        class                 = class,
        order                 = order,
        family                = family,
        genus                 = genus,
        date                  = eventDate_clean,
        LClabel               = NA_character_,
        decimalLatitude       = decimalLatitude,
        decimalLongitude      = decimalLongitude,
        station               = station,
        year                  = year,
        month                 = month,
        organismQuantity      = organismQuantity,
        concentration         = NA_real_,
        pcr_primer_lod        = NA_real_,
        detected              = detected,
        ownerContact          = ownerContact,
        bibliographicCitation = bibliographicCitation,
        datasetID_obis        = datasetID_obis
      )

    out
  })

  obis_list <- purrr::compact(obis_list)

  if (length(obis_list) == 0L) {
    warning("No OBIS datasets with DNADerivedData returned any records.")
    return(tibble::tibble())
  }

  GOTeDNA_df <- dplyr::bind_rows(obis_list) %>%
    dplyr::rename(species = scientificName)

  rownames(GOTeDNA_df) <- NULL
  GOTeDNA_df
}
