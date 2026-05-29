# TermineR 1.4.0

## Minor improvements

- Added UniProt processing annotations for Synechocystis sp. PCC 6803.
- Added `organism = "synechocystis"` support in `annotate_neo_termini()`.
- Added N-terminal Formyl support across adapters.
- Added packaged TargetP processing annotations for human, mouse, Arabidopsis thaliana, rat, and yeast.
- Added packaged TargetP processing annotations for Medicago truncatula, Rhizobium melitoli, pig, human with isoforms, Escherichia coli, Caenorhabditis elegans, and Synechocystis sp. PCC 6803.
- Added automatic TargetP-aware `processing_type` annotation in `annotate_neo_termini()` when TargetP data is available.
- Added SSH-based TargetP automation helper for downloading UniProt FASTA files locally, running TargetP on Linux, and retrieving outputs.
- Expanded the TargetP remote registry to all UniProt-supported organisms that do not yet have packaged TargetP annotations.
- Added a Linux-native TargetP batch runner for downloading UniProt FASTA files and running TargetP directly on the Linux host.
