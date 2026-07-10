# TermineR 1.4.1

## Minor improvements

- Added human CaspSites processing annotations as packaged data.
- Integrated human-only CaspSites exact-site annotation into `annotate_neo_termini()`.
- Added a reproducible CaspSites data preparation script and tests for the packaged data and annotation output.

# TermineR 1.4.0

## Minor improvements

- Added UniProt processing annotations for Synechocystis sp. PCC 6803.
- Added `organism = "synechocystis"` support in `annotate_neo_termini()`.
- Added N-terminal Formyl support across adapters.
- Added packaged TargetP processing annotations for human, mouse, Arabidopsis thaliana, rat, and yeast.
- Added packaged TargetP processing annotations for Medicago truncatula, Rhizobium melitoli, pig, human with isoforms, Escherichia coli, Caenorhabditis elegans, and Synechocystis sp. PCC 6803.
- Expanded the TargetP remote registry to all UniProt-supported organisms that do not yet have packaged TargetP annotations.
