#' Count functional groups of a compound
#'
#' Using a mol file for compound of interest, return number of
#' many functional groups in a dataframe
#'
#' @param compound_id character string that is 5 digits prepended with a "C"
#' @param pathway_id character string that is 5 digits prepended with "map"
#' @param path relative path to location to download data
#'
#' @return single row dataframe with columns for numbers of different functional groups and basic compound details
#' @export
get_fx_groups <- function(compound_id, pathway_id, path){
  rowname <- n <- phosphoric_acid <- phosphoric_ester <- rings_aromatic <- phenol <- hydroxyl_groups <- carbon_dbl_bonds <- NULL
  mol_path <- paste0(path, "/", pathway_id, "/", compound_id, ".mol")
  compound_sdf <- ChemmineR::read.SDFset(sdfstr = mol_path)
  kegg_data <- KEGGREST::keggGet(compound_id)
  groups <- data.frame(t(ChemmineR::groups(compound_sdf, groups = "fctgroup",
                                           type = "countMA")))
  rings <- data.frame(t(ChemmineR::rings(compound_sdf, type = "count", arom = TRUE)))
  atoms <- data.frame(t(unlist(ChemmineR::atomcount(compound_sdf))))
  carbon_bond_data <- data.frame(ChemmineR::conMA(compound_sdf)[[1]]) %>%
    dplyr::select(tidyselect::contains("C_")) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(stringr::str_detect(rowname, "C_")) %>%
    tibble::column_to_rownames(var = "rowname")
  carbon_dbl_count <- data.frame(all = unlist(carbon_bond_data)) %>%
    dplyr::count(all) %>%
    dplyr::filter(all == 2) %>%
    dplyr::select(n) %>%
    dplyr::mutate(n = n / 2)
  if(nrow(carbon_dbl_count) == 0) {
    carbon_dbl_count <- tibble::add_row(carbon_dbl_count, n = 0)
  }
  # *_pattern are SMARTS strings: https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
  phenol_pattern <- "[OX2H][cX3]:[c]"
  nitrate_pattern <- "[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]"
  amine_pattern <- "[NX3;H2,H1;!$(NC=O)]"
  amide_pattern <- "[NX3][CX3](=[OX1])[#6]"
  nitro_pattern <- "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]"
  phosphoric_acid_pattern <- "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]"
  phosphoric_ester_pattern <- "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]"
  sulfate_pattern <- "[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]"
  sulfonate_pattern <- "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]"
  thiol_pattern <- "[#16X2H]"
  carbothioester_pattern <- "S([#6])[CX3](=O)[#6]"
  fx_groups_df <- data.frame(pathway = pathway_id,
                             compound = compound_id,
                             formula = kegg_data[[1]]$FORMULA,
                             name = kegg_data[[1]]$NAME[1],
                             mass = ifelse(!is.null(kegg_data[[1]]$MOL_WEIGHT),
                                           as.numeric(kegg_data[[1]]$MOL_WEIGHT), NA),
                             carbons = ifelse("CMP1.C" %in% colnames(atoms),
                                              atoms$CMP1.C, 0),
                             ketones = groups$RCOR,
                             aldehydes = groups$RCHO,
                             hydroxyl_groups = groups$ROH,
                             carbox_acids = groups$RCOOH,
                             peroxide = NA,
                             hydroperoxide = NA,
                             nitrate = ChemmineR::smartsSearchOB(compound_sdf, nitrate_pattern),
                             nitro = ChemmineR::smartsSearchOB(compound_sdf, nitro_pattern),
                             carbon_dbl_bonds = carbon_dbl_count$n,
                             rings = rings$RINGS,
                             rings_aromatic = rings$AROMATIC,
                             phenol = ChemmineR::smartsSearchOB(compound_sdf, phenol_pattern, uniqueMatches = FALSE),
                             nitrophenol = NA,
                             nitroester = NA,
                             ester = groups$RCOOR,
                             ether_alicyclic = NA,
                             ether_aromatic = NA,
                             amine_primary = groups$RNH2,
                             amine_secondary = groups$R2NH,
                             amine_tertiary = groups$R3N,
                             amine_aromatic = NA,
                             amines = ChemmineR::smartsSearchOB(compound_sdf, amine_pattern),
                             amides = ChemmineR::smartsSearchOB(compound_sdf, amide_pattern),
                             phosphoric_acid = ChemmineR::smartsSearchOB(compound_sdf, phosphoric_acid_pattern),
                             phosphoric_ester = ChemmineR::smartsSearchOB(compound_sdf, phosphoric_ester_pattern),
                             sulfate = ChemmineR::smartsSearchOB(compound_sdf, sulfate_pattern),
                             sulfonate = ChemmineR::smartsSearchOB(compound_sdf, sulfonate_pattern),
                             thiol = ChemmineR::smartsSearchOB(compound_sdf, thiol_pattern),
                             carbothioester = ChemmineR::smartsSearchOB(compound_sdf, carbothioester_pattern),
                             oxygens = ifelse("CMP1.O" %in% colnames(atoms),
                                              atoms$CMP1.O, 0),
                             chlorines = ifelse("CMP1.Cl" %in% colnames(atoms),
                                                atoms$CMP1.Cl, 0),
                             nitrogens = ifelse("CMP1.N" %in% colnames(atoms),
                                                atoms$CMP1.N, 0),
                             sulfurs = ifelse("CMP1.S" %in% colnames(atoms),
                                              atoms$CMP1.S, 0),
                             phosphoruses = ifelse("CMP1.P" %in% colnames(atoms),
                                                   atoms$CMP1.P, 0),

                             bromines = ifelse("CMP1.Br" %in% colnames(atoms),
                                              atoms$CMP1.Br, 0),
                             iodines = ifelse("CMP1.I" %in% colnames(atoms),
                                              atoms$CMP1.I, 0),
                             fluorines = ifelse("CMP1.F" %in% colnames(atoms),
                                              atoms$CMP1.F, 0))
  fx_groups_df <- fx_groups_df %>%
    # to fix double counting of rings, aromatic rings, phenols, hydroxyls, carbon double bonds, and phosphoric acids/esters
    dplyr::mutate(phenol = ifelse(rings !=0 & rings_aromatic != 0 & phenol > 1, (phenol/2) - (hydroxyl_groups - 1), phenol),
                  rings = ifelse(rings !=0 & rings_aromatic != 0, rings - rings_aromatic, rings),
                  hydroxyl_groups = hydroxyl_groups - phenol,
                  carbon_dbl_bonds = ifelse(carbon_dbl_bonds != 0 & rings_aromatic != 0, carbon_dbl_bonds - (rings_aromatic * 3), carbon_dbl_bonds),
                  carbon_dbl_bonds = ifelse(carbon_dbl_bonds < 0, 0, carbon_dbl_bonds),
                  phosphoric_acid = ifelse(phosphoric_acid != 0 & phosphoric_ester != 0, phosphoric_acid - phosphoric_ester, phosphoric_acid))
  return(fx_groups_df)
}
