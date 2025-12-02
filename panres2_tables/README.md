# PanRes2 Ontology Export Tables

This directory contains simplified TSV files exported from the PanRes2 ontology.

`PanGenes.tsv`  
Contains one row per PanGene, including gene length, PanGeneCluster membership, source databases, original gene names, resistance annotations, and a binary flag indicating whether the gene was discarded.

`PanProteins.tsv`  
Lists all PanProtein sequences with their protein lengths, originating PanGene, and assigned PanProteinCluster membership.

`PanStructures.tsv`  
Contains all predicted PanStructure identifiers, the PanProteinCluster from which each structure originates, and the PanStructureCluster(s) to which they belong.
