# Data provenance and acknowledgements

This directory contains organism-specific annotation tables used by TermineR
to classify and annotate proteolytic processing sites. The tables are bundled
with the package so that normal analyses do not need to contact external
services at run time.

## UniProt-derived processing annotations

The `*_uniprot_processing` tables were prepared from protein-processing
features obtained from the UniProt database through its API. During the data
preparation workflow, the `drawProteins` R/Bioconductor package was used as a
helper for accessing and parsing UniProt API responses. No `drawProteins`
source code is bundled in TermineR.

Please acknowledge and cite:

- UniProt Consortium. UniProtKB and the UniProt databases. UniProt data are
  provided under the [Creative Commons Attribution 4.0 International (CC BY
  4.0) license](https://www.uniprot.org/help/license). The UniProt license
  applies to copyrightable parts of the database; particular records may also
  involve other rights.
- Brennan P. *drawProteins: a Bioconductor/R package for reproducible and
  programmatic generation of protein schematics*. F1000Research. 2018;7:1105.
  See the [drawProteins repository](https://github.com/brennanpincardiff/drawProteins)
  and its Bioconductor documentation. The Bioconductor package is licensed
  MIT + file LICENSE.

The exact retrieval and construction dates were not retained consistently for
each organism-specific annotation bundle. These tables should therefore be
treated as fixed snapshots distributed with the corresponding TermineR
release, rather than as current or automatically refreshed UniProt data.

## MEROPS-derived annotations

The `merops_sites`, `merops_pssm`, and `merops_accession_to_protease` objects
were prepared from MEROPS database files. The MEROPS data were processed
directly by TermineR scripts; `drawProteins` was not used for these objects.
The dataset-level sources and access information are documented in
the corresponding object documentation in `R/data.R`.

Please acknowledge and cite the MEROPS database and its maintainers. MEROPS
requests citation of the following publication together with the database URL.
The MEROPS database is provided under the GNU Library General Public License;
this applies to the database content and is separate from the GPL-3 license of
the TermineR source code.

> Rawlings ND, Barrett AJ, Thomas PD, Huang X, Bateman A, Finn RD. The MEROPS
> database of proteolytic enzymes, their substrates and inhibitors in 2017 and
> a comparison with peptidases in the PANTHER database. *Nucleic Acids
> Research*. 2018;46:D624-D632.

See the [MEROPS database](https://www.ebi.ac.uk/merops/), its [citation
guidance](https://www.ebi.ac.uk/merops/about/citing.shtml), and its
[availability and licence information](https://www.ebi.ac.uk/merops/about/availability.shtml).

## Rebuilding the data

The scripts used to prepare these objects are in `scr/`. Rebuilding them may
produce different results because UniProt and MEROPS are updated over time.
For a reproducible rebuild, record the upstream database release, query or
download URL, retrieval date, and the TermineR commit used for the build.
