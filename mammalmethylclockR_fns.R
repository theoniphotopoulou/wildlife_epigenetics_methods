# MammalMethylClock functions

# MammalMethylClockR needs a bunch of other bioconductor packages to install, avoid
# this by sourcing these key functions here

alignDatToInfo <- function(info, dat, SAMPVAR.info, SAMPVAR.dat) {
  if (!(SAMPVAR.info %in% colnames(info))) {
    stop(paste0("SAMPVAR.info=", SAMPVAR.info, " not found"))
  }
  if (!(SAMPVAR.dat %in% colnames(dat))) {
    stop(paste0("SAMPVAR.dat=", SAMPVAR.dat, " not found"))
  }
  info_alignvar <- info[, SAMPVAR.info]
  dat_alignvar <- dat[, SAMPVAR.dat]
  matchid <- match(info_alignvar, dat_alignvar)
  matchid.new <- matchid[which(!is.na(matchid))]
  info.new <- info[which(!is.na(matchid)), ]
  dat.new <- dat[matchid.new, ]
  ys <- info.new
  xs <- as.matrix(dat.new[, which(colnames(dat.new) != SAMPVAR.dat)])
  rownames(xs) <- NULL
  xs <- xs[, which(grepl("cg", colnames(xs)))] # removes probes that don't start with 'cg'
  if (is.vector(xs)) {
    xs <- t(as.matrix(xs))
  }
  yxs.list <- list(ys, xs)
  return(yxs.list)
}

predictAge <- function(xs.new, ys.new, tissue.names, species.name) {
  if (is.null(ys.new)) {
    stop("ys.new cannot be NULL")
  }
  if (length(species.name) >= 2) {
    warning("'species.name' has length > 1, so only the first element will be used")
    species.name <- species.name[1]
  }
  # DELETE LATER (USED FOR CODE TESTING):
  # devtools::load_all()
  # tissue.names = c("Blood","Ear")
  # species.name = "Ovis aries"
  # i=7
  if ("ALL" %in% tissue.names) {
    rows_match <- unique(unlist(lapply(c(species.name, "ALL"), grep, clocks_metadatabase$SpeciesNames)))
    if ("ALL" == species.name) {
      stop("Please select one species or at least one tissue type")
    }
  } else if ("ALL" == species.name) {
    rows_match <- unique(unlist(lapply(c(tissue.names, "ALL"), grep, clocks_metadatabase$Tissues)))
  } else {
    rows_match <- intersect(
      unique(unlist(lapply(c(tissue.names, "ALL"), grep, clocks_metadatabase$Tissues))),
      unique(unlist(lapply(c(species.name, "ALL"), grep, clocks_metadatabase$SpeciesNames)))
    )
  }
  if (length(rows_match) == 0) {
    stop("No tissue or species matches found")
  }
  # DELETE LATER (USED FOR CODE TESTING):
  # clocks_metadatabase[rows_match,]
  # clocks_metadatabase[rows_match,"Tissues"]
  # clocks_metadatabase[rows_match,"SpeciesNames"]
  ys.output <- ys.new
  for (i in 1:length(rows_match)) {
    clock_coefficients <- clock_coefficients_list[[match(
      clocks_metadatabase$CoefficientsName[rows_match[i]],
      names(clock_coefficients_list)
    )]]
    in.valbeta <- na.omit(clock_coefficients[, 1 + c(0, clocks_metadatabase$CoefficientsColumnNum[rows_match[i]])])
    rownames(in.valbeta) <- NULL
    attr(in.valbeta, "na.action") <- NULL
    PREDVAR_i <- unlist(strsplit(colnames(in.valbeta)[2], "Coef."))[2]
    if (!grepl("Relative", PREDVAR_i)) {
      PREDVAR_i <- paste0("DNAmAge.", PREDVAR_i)
    } else {
      PREDVAR_i <- paste0("DNAmRelAge.", PREDVAR_i)
    }
    NAME_fun_inv <- clocks_metadatabase$InverseTransform[rows_match[i]]
    fun_inv <- get(NAME_fun_inv)
    
    clock_species_vec <- unlist(strsplit(clocks_metadatabase$SpeciesNames[rows_match[i]], ";"))
    if (any(grepl(species.name, clock_species_vec))) {
      index_species_match <- grep(species.name, clock_species_vec)
      fun_VAR1_val <- as.numeric(unlist(strsplit(as.character(clocks_metadatabase$fun_VAR1[rows_match[i]]), ";")))[index_species_match]
      fun_VAR2_val <- as.numeric(unlist(strsplit(as.character(clocks_metadatabase$fun_VAR2[rows_match[i]]), ";")))[index_species_match]
    } else {
      # captures edge case, where clock_species_vec == "ALL"
      anage_row_number <- which(species.name == anAge_table$SpeciesLatinName)
      anage_column_name.fun_VAR1 <- clocks_metadatabase$fun_VAR1[rows_match[i]]
      anage_column_name.fun_VAR2 <- clocks_metadatabase$fun_VAR2[rows_match[i]]
      # valid options for clocks_metadatabase$fun_VAR1[rows_match[i]], clocks_metadatabase$fun_VAR2[rows_match[i]]:
      # - "Gestation.Incubation..days." (will be converted to years)
      # - "averagedMaturity.yrs"
      # - "maxAgeCaesar"
      if (species.name %in% anAge_table$SpeciesLatinName) {
        # if species is IN anAge
        if (!is.na(anage_column_name.fun_VAR1) && anage_column_name.fun_VAR1 == "Gestation.Incubation..days.") {
          fun_VAR1_val <- anAge_table[anage_row_number, anage_column_name.fun_VAR1] / 365
        } else {
          fun_VAR1_val <- anAge_table[anage_row_number, anage_column_name.fun_VAR1]
        }
        if (!is.na(anage_column_name.fun_VAR2) && anage_column_name.fun_VAR2 == "Gestation.Incubation..days.") {
          fun_VAR2_val <- anAge_table[anage_row_number, anage_column_name.fun_VAR2] / 365
        } else {
          fun_VAR2_val <- anAge_table[anage_row_number, anage_column_name.fun_VAR2]
        }
      } else {
        # if species is NOT IN anAge
        fun_VAR1_val <- NA
        fun_VAR2_val <- NA
        warning("Species not in anAge database")
      }
    }
    
    ys.new.prediction_i <- predictClockSimple(in.valbeta, xs.new,
                                              fun_inv = fun_inv,
                                              fun_VAR1_val = fun_VAR1_val, fun_VAR2_val = fun_VAR2_val
    )
    if (length(ys.new.prediction_i) > 0) {
      ys.output[, PREDVAR_i] <- ys.new.prediction_i
      attr(ys.output[, PREDVAR_i], "dimnames") <- NULL
    }
  }
  
  return(ys.output)
}

fun_identity <- function(x, ...) x

predictClockSimple <- function(in.valbeta, xs.new, fun_inv = fun_identity, fun_VAR1_val = NA, fun_VAR2_val = NA) {
  if (ncol(in.valbeta) != 2) {
    stop("in.valbeta doesn't follow correct format")
  }
  
  ones.column <- data.frame(V1 = rep(1, nrow(xs.new)))
  colnames(ones.column) <- as.character(in.valbeta[1, 1])
  
  xs.new_complete <- t(na.omit(t(xs.new))) # remove data columns (probes) that have missing values
  if (ncol(xs.new) > ncol(xs.new_complete)) {
    warning(paste0(ncol(xs.new) - ncol(xs.new_complete), " columns of xs.new removed due to missing values"))
  }
  
  in.valbeta_trimmed <- in.valbeta[sort(c(1, which(as.character(in.valbeta[, 1]) %in% colnames(xs.new_complete)))), ]
  if (nrow(in.valbeta_trimmed) - 1 == 0) {
    warning("methylation data (xs argument) is missing all clock coefficients in columns, or has missing data in all corresponding columns")
    return(rep(NA, nrow(xs.new)))
  }
  num_beta.not_in_New <- nrow(in.valbeta) - nrow(in.valbeta_trimmed)
  if (num_beta.not_in_New > 0) {
    warning(paste0(nrow(in.valbeta_trimmed) - 1, " clock probes found in xs.new, ", num_beta.not_in_New, " clock probes not found or have missing values in xs.new"))
  }
  xs.new_trimmed <- cbind(ones.column, xs.new_complete)
  xs.new_trimmed <- xs.new_trimmed[, c(1, which(colnames(xs.new_trimmed) %in% as.character(in.valbeta_trimmed[, 1])))]
  xs.new_trimmed <- dplyr::select(xs.new_trimmed, as.character(in.valbeta_trimmed[, 1])) # sort xs.new to match in.valbeta
  ys.new.prediction.raw <- as.vector(as.matrix(xs.new_trimmed) %*% in.valbeta_trimmed[, 2])
  if (length(ys.new.prediction.raw) > 0) {
    ys.new.prediction <- fun_inv(ys.new.prediction.raw, fun_VAR1_val, fun_VAR2_val)
  } else {
    ys.new.prediction <- ys.new.prediction.raw
  }
  
  return(ys.new.prediction)
}

