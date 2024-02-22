#' Count compound functional groups
#'
#' Returns functional group counts relevant to calculating estimated volatility
#' for specified compounds. Users will not typically interact with this function
#' directly, but rather by using [calc_vol()].
#' 
#' @note This function currently does **not** capture the carbon number on the
#'   acid-side of amide, one of the functional groups used in SIMPOL.1.
#'   Contributions of SMARTS strings or other methods to capture this
#'   "functional group" are welcome.
#'
#' @param compound_sdf a [ChemmineR::SDFset] object returned by
#'   [ChemmineR::read.SDFset()] or [ChemmineR::smiles2sdf()], for example.
#'
#' @returns A tibble with columns of basic compound info and functional group
#'   counts.
#' @seealso [calc_vol()]
#' @examples
#' if (rlang::is_installed("ChemmineOB")) {
#'   mol_path <- mol_example()[1]
#'   sdf <- ChemmineR::read.SDFset(mol_path)
#'   get_fx_groups(sdf)
#' }
#' 
#' @export
get_fx_groups <- function(compound_sdf) {
  # Check that ChemmineOB is installed
  rlang::check_installed("ChemmineOB", action = BiocManager::install)
  
  # For now at least, this code only works with SDFset objects that contain
  # single molecules. 
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
  rowname <- n <- NULL
  
  #convert counts to integer
  groups <-
    groups %>% 
    dplyr::mutate(dplyr::across(dplyr::everything(), as.integer))
  rings <- 
    data.frame(t(ChemmineR::rings(compound_sdf, type = "count", arom = TRUE, inner = TRUE)))
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
  # ether_alkyl_pattern <- "[OD2]([C!R1])[C!R1]" #currently unused--ether_alkly calculated as total - other ethers
  ether_alicyclic_pattern <- "[OD2]([C!R0])[C!R0]"
  ether_aromatic_pattern <- "O(c)[C,c]" #only one of the carbons has to be aromatic
  nitro_pattern <- "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]"
  hydroxyl_aromatic_pattern <- "[OX2H]c"
  nitrate_pattern <- "[$([NX3](=[OX1])(=[OX1])O),$([NX3+]([OX1-])(=[OX1])O)]"

  #TODO need patterns for amines that don't pick up amides
  amine_primary_pattern <- "[NX3;H2;!$(NC=[!#6]);!$(NC#[!#6])][#6X4]"
  amine_secondary_pattern <- "[NX3H1!$(NC=[!#6])!$(NC#[!#6])]([#6X4])[#6X4]"
  amine_tertiary_pattern <- "[NX3H0!$(NC=[!#6])!$(NC#[!#6])]([#6X4])([#6X4])[#6X4]"
  amine_aromatic_pattern <-  "[NX3;!$(NO)]c" 

  amide_primary_pattern <- "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H2]"
  amide_secondary_pattern <- "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H1][#6;!$(C=[O,N,S])]"
  amide_tertiary_pattern <- 
    "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3H0]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])]"
  
  # amide_total_pattern <- "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[#7X3;$([H2]),$([H1][#6;!$(C=[O,N,S])]),$([#7]([#6;!$(C=[O,N,S])])[#6;!$(C=[O,N,S])])]"
  
  carbonylperoxynitrate_pattern <- "*C(=O)OO[N+1](=O)[O-1]"
  peroxide_pattern <- "[OX2D2][OX2D2]"  #this captures carbonylperoxynitrates too
  hydroperoxide_pattern <- "[OX2][OX2H,OX1-]" #this captures peroxyacids too
  carbonylperoxyacid_pattern <- "[CX3;$([R0][#6]),$([H1R0])](=[OX1])[OX2][$([OX2H]),$([OX1-])]" 
  nitroester_pattern <- "C(=O)(OC)C~[NX3](-,=[OX1])-,=[OX1]" 
  # This captures OH groups on a ring that also has a nitro group (para, ortho, or meta).  Need to correct aromatic hydroxyl count later.
  nitrophenol_pattern <- 
    "[OX2H][$(c1ccccc1[$([NX3](=O)=O),$([NX3+](=O)[O-])]),$(c1cccc(c1)[$([NX3](=O)=O),$([NX3+](=O)[O-])]),$(c1ccc(cc1)[$([NX3](=O)=O),$([NX3+](=O)[O-])])]"
  phosphoric_acid_pattern <-
    "[$(P(=[OX1])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)]),$([P+]([OX1-])([$([OX2H]),$([OX1-]),$([OX2]P)])([$([OX2H]),$([OX1-]),$([OX2]P)])[$([OX2H]),$([OX1-]),$([OX2]P)])]"
  phosphoric_ester_pattern <-
    "[$(P(=[OX1])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)]),$([P+]([OX1-])([OX2][#6])([$([OX2H]),$([OX1-]),$([OX2][#6])])[$([OX2H]),$([OX1-]),$([OX2][#6]),$([OX2]P)])]"
  sulfate_pattern <-
    "[$([#16X4](=[OX1])(=[OX1])([OX2H,OX1H0-])[OX2][#6]),$([#16X4+2]([OX1-])([OX1-])([OX2H,OX1H0-])[OX2][#6])]"
  #sulfonate groups; sulfonate ions, and conjugate acid, sulfonic acids
  sulfonate_pattern <-
    "[#16X4](=[OX1])(=[OX1])([#6])[*$([O-1]),*$([OH1]),*$([OX2H0])]"
  thiol_pattern <- "[#16X2H]"
  carbothioester_pattern <- "S([#6])[CX3](=O)[#6]"
  
  fx_groups_df <- 
    dplyr::tibble(
      formula = ChemmineR::propOB(compound_sdf)$formula,
      #TODO should name be moved to `calc_vol`? `formula` also?
      name = ChemmineR::propOB(compound_sdf)$title,
      exact_mass = ChemmineR::exactMassOB(compound_sdf),
      molecular_weight = ChemmineR::propOB(compound_sdf)$MW
    ) %>% 
    dplyr::mutate(
      carbons = atoms[["C"]] %||% 0L,
      carbons_asa = NA_integer_, #carbon number on the acid-side of amide
      rings_aromatic = as.integer(rings$AROMATIC),
      rings_total = as.integer(rings$RINGS),
      rings_aliphatic = NA_integer_, #calculated below
      carbon_dbl_bonds_aliphatic = ChemmineR::smartsSearchOB(compound_sdf, carbon_dbl_bonds_pattern),
      CCCO_aliphatic_ring = ChemmineR::smartsSearchOB(compound_sdf, CCCO_pattern), # C=C-C=O in a non-aromatic ring
      hydroxyl_total = groups$ROH, 
      hydroxyl_aromatic = ChemmineR::smartsSearchOB(compound_sdf, hydroxyl_aromatic_pattern, uniqueMatches = FALSE),
      hydroxyl_aliphatic = NA_integer_, #calculated below
      aldehydes = groups$RCHO,
      ketones = groups$RCOR,
      carbox_acids = groups$RCOOH,
      ester = groups$RCOOR,
      ether_total = groups$ROR,
      # ether_alkyl = ChemmineR::smartsSearchOB(compound_sdf, ether_alkyl_pattern),
      ether_alkyl = NA_integer_,
      ether_alicyclic = ChemmineR::smartsSearchOB(compound_sdf, ether_alicyclic_pattern),
      ether_aromatic = ChemmineR::smartsSearchOB(compound_sdf, ether_aromatic_pattern),
      nitrate = ChemmineR::smartsSearchOB(compound_sdf, nitrate_pattern),
      nitro = ChemmineR::smartsSearchOB(compound_sdf, nitro_pattern),
      amine_primary   = ChemmineR::smartsSearchOB(compound_sdf, amine_primary_pattern),
      amine_secondary = ChemmineR::smartsSearchOB(compound_sdf, amine_secondary_pattern),
      amine_tertiary  = ChemmineR::smartsSearchOB(compound_sdf, amine_tertiary_pattern),
      amine_aromatic = ChemmineR::smartsSearchOB(compound_sdf, amine_aromatic_pattern),
      amide_primary = ChemmineR::smartsSearchOB(compound_sdf, amide_primary_pattern),
      amide_secondary = ChemmineR::smartsSearchOB(compound_sdf, amide_secondary_pattern),
      amide_tertiary = ChemmineR::smartsSearchOB(compound_sdf, amide_tertiary_pattern),
      carbonylperoxynitrate = ChemmineR::smartsSearchOB(compound_sdf, carbonylperoxynitrate_pattern),
      peroxide = ChemmineR::smartsSearchOB(compound_sdf, peroxide_pattern),
      hydroperoxide = ChemmineR::smartsSearchOB(compound_sdf, hydroperoxide_pattern),
      carbonylperoxyacid = ChemmineR::smartsSearchOB(compound_sdf, carbonylperoxyacid_pattern),
      nitrophenol = ChemmineR::smartsSearchOB(compound_sdf, nitrophenol_pattern),
      nitroester = ChemmineR::smartsSearchOB(compound_sdf, nitroester_pattern),
      
      # Additional groups from Meredith et al. 2023
      phosphoric_acids = ChemmineR::smartsSearchOB(compound_sdf, phosphoric_acid_pattern),
      phosphoric_esters = ChemmineR::smartsSearchOB(compound_sdf, phosphoric_ester_pattern),
      sulfates = ChemmineR::smartsSearchOB(compound_sdf, sulfate_pattern),
      sulfonates = ChemmineR::smartsSearchOB(compound_sdf, sulfonate_pattern),
      thiols = ChemmineR::smartsSearchOB(compound_sdf, thiol_pattern),
      carbothioesters = ChemmineR::smartsSearchOB(compound_sdf, carbothioester_pattern),
      oxygens   = atoms[["O"]] %||% 0L,
      chlorines = atoms[["Cl"]] %||% 0L,
      nitrogens = atoms[["N"]] %||% 0L,
      sulfurs   = atoms[["S"]] %||% 0L,
      phosphoruses = atoms[["P"]] %||% 0L,
      bromines = atoms[["Br"]] %||% 0L,
      iodines = atoms[["I"]] %||% 0L,
      fluorines = atoms[["F"]] %||% 0L
    ) %>% 
    
    # Corrections & Calculations
    # The order these happen in matters!
    dplyr::mutate(
      rings_aliphatic = .data$rings_total - .data$rings_aromatic,
      hydroxyl_aliphatic = .data$hydroxyl_total - .data$hydroxyl_aromatic, 
      ether_alkyl = .data$ether_total - .data$ether_alicyclic - .data$ether_aromatic,
      
      #hydroperoxide pattern also picks up peroxyacids
      hydroperoxide = .data$hydroperoxide - .data$carbonylperoxyacid,
      #peroxide, nitrate, and ester patterns also pick up carbonylperoxynitrate group
      ester = .data$ester - .data$carbonylperoxynitrate,
      nitrate = .data$nitrate - .data$carbonylperoxynitrate,
      peroxide = .data$peroxide - .data$carbonylperoxynitrate,
      #phosphoric ester also matches phosphoric acid
      phosphoric_acids = .data$phosphoric_acids - .data$phosphoric_esters,
      #according to SIMPOL.1 paper, nitrophenol shouldn't count aromatic hydroxyls that are part of the nitrophenol group separately.
      hydroxyl_aromatic = .data$hydroxyl_aromatic - .data$nitrophenol,
      #according to SIMPOL.1 paper, nitroester shouldn't count the esters that are part of the nitroester group separately.
      ester = .data$ester - .data$nitroester
    ) %>% 
    # some of the columns created by ChemmineR are named vectors sometimes,
    # strip names for consistency
    dplyr::mutate(
      dplyr::across(dplyr::everything(), function(x) stats::setNames(x, NULL))
    ) %>% 
    
    #TODO should this be moved to `calc_vol?`. It's only relevant when from = "mol_path"
    dplyr::mutate(name = ifelse(.data$name == "", NA_character_, .data$name))
  
  #return
  fx_groups_df
}
