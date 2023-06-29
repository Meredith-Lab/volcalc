#' Count functional groups of a compound
#'
#' Return functional group counts relevant to volatility calculation for
#' specified compound
#' 
#' @note Users do not typically need to interact with this function, but it is
#'   used internally by `calc_vol()`
#'
#' @param compound_sdf object created by `ChemmineR::read.SDFset()`
#'
#' @return a tibble with columns of basic compound info and functional group
#'   counts.
#'
#' @examples
#' \dontrun{
#' get_df <- get_mol_kegg(compound_ids = "C16181", dir = tempdir())
#' mol <- ChemmineR::read.SDFset(get_df$mol_path)
#' get_fx_groups(mol)
#' }
#' 
#' @export
get_fx_groups <-
  function(compound_sdf) {
  
    #TODO: Should get_fx_groups be vectorized?  Vectorization to work on
    #multiple compounds could happen at three levels: 1) get_fx_groups() could
    #accept a vector of SDFset objects. 2) SDFset objects can contain multiple
    #compounds, so it could accept a single SDFset with multiple compounds. 3)
    #it could be restricted to a single SDFset with a single molecule and then
    #be applied to multiple compunds in `calc_vol()`
    
    if(length(compound_sdf) != 1) {
      stop("SDFset objects must contain a single molecule only")
      #this is partly because of type instability of groups() https://github.com/girke-lab/ChemmineR/issues/15
    }
              
    #assign variables to quiet devtools::check()
    rowname <- n <- phosphoric_acid <- phosphoric_ester <- rings_aromatic <- phenol <- hydroxyl_groups <- carbon_dbl_bonds <- NULL
    
  groups <- data.frame(t(ChemmineR::groups(compound_sdf,
    groups = "fctgroup",
    type = "countMA"
  )))
  rings <- data.frame(t(ChemmineR::rings(compound_sdf, type = "count", arom = TRUE)))
  atoms <- data.frame(t(unlist(ChemmineR::atomcount(compound_sdf))))
  carbon_bond_data <- data.frame(ChemmineR::conMA(compound_sdf)[[1]]) %>%
    dplyr::select(dplyr::contains("C_")) %>%
    tibble::rownames_to_column() %>%
    dplyr::filter(stringr::str_detect(rowname, "C_")) %>%
    tibble::column_to_rownames(var = "rowname")
  if (nrow(carbon_bond_data) == 0){
    carbon_dbl_count <- tibble::tibble(n = 0)
  } else {
    carbon_dbl_count <- 
      data.frame(all = unlist(carbon_bond_data)) %>%
      dplyr::count(all) %>%
      dplyr::filter(all == 2) %>%
      dplyr::select(n) %>%
      dplyr::mutate(n = n / 2)
  }
  if (nrow(carbon_dbl_count) == 0) {
    carbon_dbl_count <- tibble::add_row(carbon_dbl_count, n = 0)
  }
  # *_pattern are SMARTS strings: https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html
  peroxide_pattern <- "[OX2,OX1-][OX2,OX1-]"
  phenol_pattern <- "[OX2H][cX3]:[c]"
  nitrate_pattern <- "[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]"
  amine_pattern <- "[NX3;H2,H1;!$(NC=O)]"
  amide_pattern <- "[NX3][CX3](=[OX1])[#6]"
  nitro_pattern <- "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]"
  ether_pattern <- "[OD2]([#6])[#6]"
  phosphoric_acid_pattern <- "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]"
  phosphoric_ester_pattern <- "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]"
  sulfate_pattern <- "[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]"
  sulfonate_pattern <- "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]"
  thiol_pattern <- "[#16X2H]"
  carbothioester_pattern <- "S([#6])[CX3](=O)[#6]"
  fx_groups_df <- data.frame(
    formula = ChemmineR::propOB(compound_sdf)$formula,
    name = ChemmineR::propOB(compound_sdf)$title, #TODO replace empty string with NA
    mass = ChemmineR::propOB(compound_sdf)$MW, #TODO need to replace with NA if empty?
    carbons = ifelse("CMP1.C" %in% colnames(atoms),
      atoms$CMP1.C, 0
    ),
    ketones = groups$RCOR,
    aldehydes = groups$RCHO,
    hydroxyl_groups = groups$ROH,
    carbox_acids = groups$RCOOH,
    peroxide = ChemmineR::smartsSearchOB(compound_sdf, peroxide_pattern),
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
    ether = ChemmineR::smartsSearchOB(compound_sdf, ether_pattern),
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
      atoms$CMP1.O, 0
    ),
    chlorines = ifelse("CMP1.Cl" %in% colnames(atoms),
      atoms$CMP1.Cl, 0
    ),
    nitrogens = ifelse("CMP1.N" %in% colnames(atoms),
      atoms$CMP1.N, 0
    ),
    sulfurs = ifelse("CMP1.S" %in% colnames(atoms),
      atoms$CMP1.S, 0
    ),
    phosphoruses = ifelse("CMP1.P" %in% colnames(atoms),
      atoms$CMP1.P, 0
    ),
    bromines = ifelse("CMP1.Br" %in% colnames(atoms),
      atoms$CMP1.Br, 0
    ),
    iodines = ifelse("CMP1.I" %in% colnames(atoms),
      atoms$CMP1.I, 0
    ),
    fluorines = ifelse("CMP1.F" %in% colnames(atoms),
      atoms$CMP1.F, 0
    )
  )
  fx_groups_df <- 
    fx_groups_df %>%
    # to fix double counting of rings, aromatic rings, phenols, hydroxyls, carbon double bonds, and phosphoric acids/esters
    dplyr::mutate(
      phenol = ifelse(rings != 0 & rings_aromatic != 0 & phenol > 1, (phenol / 2) - (hydroxyl_groups - 1), phenol),
      rings = ifelse(rings != 0 & rings_aromatic != 0, rings - rings_aromatic, rings),
      hydroxyl_groups = hydroxyl_groups - phenol,
      carbon_dbl_bonds = ifelse(carbon_dbl_bonds != 0 & rings_aromatic != 0, carbon_dbl_bonds - (rings_aromatic * 3), carbon_dbl_bonds),
      carbon_dbl_bonds = ifelse(carbon_dbl_bonds < 0, 0, carbon_dbl_bonds),
      phosphoric_acid = ifelse(phosphoric_acid != 0 & phosphoric_ester != 0, phosphoric_acid - phosphoric_ester, phosphoric_acid)
    )
  tibble::as_tibble(fx_groups_df)
}
