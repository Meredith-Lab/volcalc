#' Count compound functional groups
#'
#' Returns functional group counts relevant to calculating estimated volatility
#' for specified compounds. Users will not typically interact with this function
#' directly, but rather by using [calc_vol()].
#' 
#' @note This function currently does not capture hydroperoxide, nitrophenol,
#'   nitroesther, alicyclic ether, aromatic ether, or armoatic amine groups.
#'   Contributions of SMARTS strings to capture these groups are welcome.
#'
#' @param compound_sdf a [ChemmineR::SDFset] object returned by
#'   [ChemmineR::read.SDFset()] or [ChemmineR::smiles2sdf()], for example.
#'
#' @return a tibble with columns of basic compound info and functional group
#'   counts.
#' @seealso [calc_vol()]
#' @examples
#' mol_path <- mol_example("C16181.mol")
#' sdf <- ChemmineR::read.SDFset(mol_path)
#' fx_groups <- get_fx_groups(sdf)
#' 
#' @export
get_fx_groups <- function(compound_sdf) {
  
  if(length(compound_sdf) != 1) {
    stop("SDFset objects must contain a single molecule only")
    # this is partly because of type instability of groups():
    # https://github.com/girke-lab/ChemmineR/issues/15
    #TODO the above bug has been fixed in the dev version, which may cause breaking changes here
  }
    
  #assign variables to quiet devtools::check()
  rowname <- n <- phosphoric_acid <- phosphoric_ester <- rings_aromatic <- hydroxyl_aromatic <- hydroxyl_groups <- carbon_dbl_bonds <- NULL
  
  groups <- 
    tibble::as_tibble_row(ChemmineR::groups(compound_sdf,
                                            groups = "fctgroup",
                                            type = "countMA")) %>%
    dplyr::mutate(dplyr::across(dplyr::everything(), as.integer))
  rings <- data.frame(t(ChemmineR::rings(compound_sdf, type = "count", arom = TRUE, inner = TRUE)))
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
  hydroxyl_aromatic_pattern <- "[OX2H]c"
  nitrate_pattern <- "[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]"
  amide_pattern <- "[NX3][CX3](=[OX1])[#6]"
  nitro_pattern <- "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]"
  ether_pattern <- "[OD2]([#6])[#6]"
  phosphoric_acid_pattern <- "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]"
  phosphoric_ester_pattern <- "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]"
  sulfate_pattern <- "[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]"
  sulfonate_pattern <- "[$([#16X4](=[OX1])(=[OX1])([#6])[OX2H0]),$([#16X4+2]([OX1-])([OX1-])([#6])[OX2H0])]"
  thiol_pattern <- "[#16X2H]"
  carbothioester_pattern <- "S([#6])[CX3](=O)[#6]"
  
  #TODO make these column names as specific as possible.  E.g. instead of "hydroxyl_groups" it should be "alkyl_hydroxyls" (what we want to capture) or "total_hydroxyls" (what is currently captured).  Instead of "phenol" it should be "aromatic_hydroxyls".
  fx_groups_df <- data.frame(
    formula = ChemmineR::propOB(compound_sdf)$formula,
    #TODO should name be moved to `calc_vol`? `formula` also?
    name = ChemmineR::propOB(compound_sdf)$title,
    mass = ChemmineR::propOB(compound_sdf)$MW, #TODO need to replace with NA if empty?
    #TODO these columns should all be integer
    carbons = ifelse("CMP1.C" %in% colnames(atoms),
                     atoms$CMP1.C, 0L
    ),
    carbons_asa = NA_integer_, #carbon number on the acid-side of amide
    rings_aromatic = rings$AROMATIC,
    rings = rings$RINGS, #TODO: call this rings_total?
    carbon_dbl_bonds = carbon_dbl_count$n, #TODO: this should be only non-aromatic double bonds
    CCCO_aliphatic_ring = NA_integer_, # C=C-C=O in a non-aromatic ring
    hydroxyl_groups = groups$ROH, #TODO: this is total, should be just aliphatic for SIMPOL.1
    aldehydes = groups$RCHO,
    ketones = groups$RCOR,
    carbox_acids = groups$RCOOH,
    ester = groups$RCOOR,
    ether = ChemmineR::smartsSearchOB(compound_sdf, ether_pattern),
    ether_alicyclic = NA_integer_,
    ether_aromatic = NA_integer_,
    nitrate = ChemmineR::smartsSearchOB(compound_sdf, nitrate_pattern),
    nitro = ChemmineR::smartsSearchOB(compound_sdf, nitro_pattern),
    hydroxyl_aromatic = ChemmineR::smartsSearchOB(compound_sdf, hydroxyl_aromatic_pattern, uniqueMatches = FALSE),
    amine_primary = groups$RNH2,
    amine_secondary = groups$R2NH,
    amine_tertiary = groups$R3N,
    amine_aromatic = NA_integer_,
    amides = ChemmineR::smartsSearchOB(compound_sdf, amide_pattern),
    #TODO should have primary, secondary, tertiary amides
    carbonylperoxynitrate = NA_integer_,
    peroxide = ChemmineR::smartsSearchOB(compound_sdf, peroxide_pattern),
    hydroperoxide = NA_integer_,
    carbonylperoxyacid = NA_integer_,
    nitrophenol = NA_integer_,
    nitroester = NA_integer_,

    phosphoric_acid = ChemmineR::smartsSearchOB(compound_sdf, phosphoric_acid_pattern),
    phosphoric_ester = ChemmineR::smartsSearchOB(compound_sdf, phosphoric_ester_pattern),
    sulfate = ChemmineR::smartsSearchOB(compound_sdf, sulfate_pattern),
    sulfonate = ChemmineR::smartsSearchOB(compound_sdf, sulfonate_pattern),
    thiol = ChemmineR::smartsSearchOB(compound_sdf, thiol_pattern),
    carbothioester = ChemmineR::smartsSearchOB(compound_sdf, carbothioester_pattern),
    oxygens = ifelse("CMP1.O" %in% colnames(atoms),
                     as.integer(atoms$CMP1.O), 0
    ),
    chlorines = ifelse("CMP1.Cl" %in% colnames(atoms),
                       as.integer(atoms$CMP1.Cl), 0L
    ),
    nitrogens = ifelse("CMP1.N" %in% colnames(atoms),
                       as.integer(atoms$CMP1.N), 0L
    ),
    sulfurs = ifelse("CMP1.S" %in% colnames(atoms),
                     as.integer(atoms$CMP1.S), 0L
    ),
    phosphoruses = ifelse("CMP1.P" %in% colnames(atoms),
                          as.integer(atoms$CMP1.P), 0L
    ),
    bromines = ifelse("CMP1.Br" %in% colnames(atoms),
                      as.integer(atoms$CMP1.Br), 0L
    ),
    iodines = ifelse("CMP1.I" %in% colnames(atoms),
                     as.integer(atoms$CMP1.I), 0L
    ),
    fluorines = ifelse("CMP1.F" %in% colnames(atoms),
                       as.integer(atoms$CMP1.F), 0L
    )
  ) %>% 
    #TODO should this be moved to `calc_vol?`. It's only relevant when from = "mol_path"
    dplyr::mutate(name = ifelse(.data$name == "", NA_character_, .data$name))
  
  fx_groups_df <- 
    fx_groups_df %>%
    # to fix double counting of rings, aromatic rings, hydroxyls, carbon double bonds, and phosphoric acids/esters
    #TODO clarify this in documentation.  E.g. "rings" doesn't include phenols and other aromatic rings, "peroxides" doesn't include hydroperoxides (eventually)
    dplyr::mutate(
      rings = ifelse(rings != 0 & rings_aromatic != 0, rings - rings_aromatic, rings),
      hydroxyl_groups = hydroxyl_groups - hydroxyl_aromatic,
      carbon_dbl_bonds = ifelse(carbon_dbl_bonds != 0 & rings_aromatic != 0, carbon_dbl_bonds - (rings_aromatic * 3), carbon_dbl_bonds),
      carbon_dbl_bonds = ifelse(carbon_dbl_bonds < 0, 0, carbon_dbl_bonds),
      phosphoric_acid = ifelse(phosphoric_acid != 0 & phosphoric_ester != 0, phosphoric_acid - phosphoric_ester, phosphoric_acid)
    )
  tibble::as_tibble(fx_groups_df)
}
