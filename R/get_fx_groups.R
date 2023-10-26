#' Count compound functional groups
#'
#' Returns functional group counts relevant to calculating estimated volatility
#' for specified compounds. Users will not typically interact with this function
#' directly, but rather by using [calc_vol()].
#' 
#' @note This function currently does **not** capture the following functional
#'   groups used in SIMPOL.1:
#' 
#' - carbon number on the acid-side of amide
#' - C=C-C=O in a non-aromatic ring
#' - alicyclic ether
#' - aromatic ether
#' - aromatic amine
#' - carbonylperoxynitrate
#' - hydroperoxide
#' - carbonylperoxyacid
#' - nitrophenol
#' - nitroesther
#' 
#' Contributions of SMARTS strings to capture these groups are welcome.
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

  # For now at least, this code only works with SDFset objects that contain single molecules.
  # TODO: make this function work with SDFset objects with multiple molecules?
  if (length(compound_sdf) != 1) {
    stop("SDFset objects must contain a single molecule only")
  }
  
  chem_groups <- ChemmineR::groups(compound_sdf,
                                   groups = "fctgroup",
                                   type = "countMA")
  
  # Handle different behavior of ChemmineR::groups() depending on version when
  # length(compound_sdf) == 1. See more here:
  # https://github.com/girke-lab/ChemmineR/issues/15
  if (utils::packageVersion("ChemmineR") < "3.53.1") {
    groups <- tibble::as_tibble_row(chem_groups)
  } else {
    groups <- tibble::as_tibble(chem_groups) 
  }
    
  #assign variables to quiet devtools::check()
  rowname <- n <- phosphoric_acid <- phosphoric_ester <- rings_aromatic <- hydroxyl_aromatic <- hydroxyl_total <- carbon_dbl_bonds <- NULL
  
  #convert counts to integer
  groups <- groups %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.integer))
  rings <- data.frame(t(ChemmineR::rings(compound_sdf, type = "count", arom = TRUE, inner = TRUE)))
  atoms <- atomcount2tibble(ChemmineR::atomcount(compound_sdf))
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
  carbon_dbl_bonds_pattern <- "C=C" #non-aromatic carbon double bonds
  CCCO_pattern <- "C(C=C[AR1])(=O)[AR1]" #C=C-C=O in a non-aromatic ring
  ether_pattern <- "[OD2]([C!R1!X1])[C!R1!X1]" #TODO disambiguate ether and esther.
  ether_alicyclic_pattern <- "[OD2]([C!R0])[C!R0]"
  ether_aromatic_pattern <- "[OD2]([cX2])[cX2]"
  nitro_pattern <- "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]"
  hydroxyl_aromatic_pattern <- "[OX2H]c"
  nitrate_pattern <- "[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]"
  amine_aromatic_pattern <-  "[NX3;!$(NO)]c" 
  # amide_total_pattern <- "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]"
  amide_primary_pattern <- "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H2]"
  amide_secondary_pattern <- "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H1][#6;!$(C=[O,N,S])]"
  amide_tertiary_pattern <- "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])]"
  peroxide_pattern <- "[OX2D2][OX2D2]" 
  hydroperoxide_pattern <- "[OX2][OX2H,OX1-]" #TODO this captures peroxyacids too
  carbonylperoxyacid_pattern <- "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][$([OX2H]),$([OX1-])]"
  nitroester_pattern <- "[OX2][N+]([O-])=O"
  nitrophenol_pattern <- "([OX2H][cr6]).([cr6$([NX3](=O)=O),$([NX3+](=O)[O-])])" #TODO Defined as # of hydroxyl on an aromatic ring when there's also a nitro on that ring.  I.e. nitro groups get counted, but count of aromatic hydroxyls goes into nitrophenol instead. Still not working
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
    exact_mass = ChemmineR::exactMassOB(compound_sdf),
    molecular_weight = ChemmineR::propOB(compound_sdf)$MW, #TODO need to replace with NA if empty?
    #TODO these columns should all be integer
    carbons = ifelse("C" %in% colnames(atoms),
                     atoms$C, 0L
    ),
    carbons_asa = NA_integer_, #carbon number on the acid-side of amide
    rings_aromatic = as.integer(rings$AROMATIC),
    rings_total = as.integer(rings$RINGS),
    carbon_dbl_bonds = ChemmineR::smartsSearchOB(compound_sdf, carbon_dbl_bonds_pattern),
    CCCO_aliphatic_ring = ChemmineR::smartsSearchOB(compound_sdf, CCCO_pattern), # C=C-C=O in a non-aromatic ring
    hydroxyl_total = groups$ROH, #this is total, need just aliphatic for SIMPOL.1, corrected below
    aldehydes = groups$RCHO,
    ketones = groups$RCOR,
    carbox_acids = groups$RCOOH,
    ester = groups$RCOOR,
    ether = ChemmineR::smartsSearchOB(compound_sdf, ether_pattern),
    ether_alicyclic = ChemmineR::smartsSearchOB(compound_sdf, ether_alicyclic_pattern),
    ether_aromatic = ChemmineR::smartsSearchOB(compound_sdf, ether_aromatic_pattern),
    nitrate = ChemmineR::smartsSearchOB(compound_sdf, nitrate_pattern),
    nitro = ChemmineR::smartsSearchOB(compound_sdf, nitro_pattern),
    hydroxyl_aromatic = ChemmineR::smartsSearchOB(compound_sdf, hydroxyl_aromatic_pattern, uniqueMatches = FALSE),
    amine_primary = groups$RNH2,
    amine_secondary = groups$R2NH,
    amine_tertiary = groups$R3N,
    amine_aromatic = ChemmineR::smartsSearchOB(compound_sdf, amine_aromatic_pattern),
    amide_primary = ChemmineR::smartsSearchOB(compound_sdf, amide_primary_pattern),
    amide_secondary = ChemmineR::smartsSearchOB(compound_sdf, amide_secondary_pattern),
    amide_tertiary = ChemmineR::smartsSearchOB(compound_sdf, amide_tertiary_pattern),
    carbonylperoxynitrate = NA_integer_,
    peroxide = ChemmineR::smartsSearchOB(compound_sdf, peroxide_pattern),
    hydroperoxide = ChemmineR::smartsSearchOB(compound_sdf, hydroperoxide_pattern),
    carbonylperoxyacid = ChemmineR::smartsSearchOB(compound_sdf, carbonylperoxyacid_pattern),
    nitrophenol = NA_integer_,
    # nitrophenol = ChemmineR::smartsSearchOB(compound_sdf, nitrophenol_pattern),
    nitroester = ChemmineR::smartsSearchOB(compound_sdf, nitroester_pattern),
  
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
      # rings_aliphatic = ifelse(rings != 0 & rings_aromatic != 0, rings - rings_aromatic, rings),
      rings_aliphatic = rings_total - rings_aromatic,
      hydroxyl_aliphatic = hydroxyl_total - hydroxyl_aromatic,
      phosphoric_acid = ifelse(phosphoric_acid != 0 & phosphoric_ester != 0, phosphoric_acid - phosphoric_ester, phosphoric_acid),
      #TODO probably should change name of `ether` to `ether_alkyl`
      ether_total = sum(ether, ether_alicyclic, ether_aromatic)
    )
  tibble::as_tibble(fx_groups_df)
}
