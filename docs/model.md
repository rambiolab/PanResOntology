### Classes
#### Class: `Resistance`
- **Description**: Overall category for resistance in the PanRes database.

#### Class: `Gene`
- **Description**: The major component of a reference sequence database: the genes. 
- **Parent Class**: `Resistance`
- **Children**:
  - Class: `PanGene`
    - **Description**: This is a gene identifier that follows the pan_ naming scheme. 
      - **Children**:
      - Class: `AntimicrobialResistanceGene`
      - Class: `BiocideResistanceGene`
      - Class: `DiscardedResistanceGene`
      - Class: `MetalResistanceGene`
  - Class: `OriginalGene`
    - **Description**: This is the original name extracted from the fasta header for each individual gene.
  - Class: `PanGeneCluster`
    - **Description**: Represents a cluster of PanGenes.

#### Class: `Database`
- **Description**: Represents a database from which genes are originally sourced.
- **Parent Class**: `Resistance`
- **Children**: 
  - Class: `AMRFinderPlus`
  - Class: `ARGANNOT`
  - Class: `CARD`
  - Class: `CsabaPal`
  - Class: `MegaRes`
  - Class: `MetalRes`
  - Class: `ResFinder`
  - Class: `ResFinderFG`
  - Class: `BacMet`

#### Class: `ResistanceType`
- **Description**: Represents different types of resistance.
- **Parent Class**: `Resistance`
- **Children**:
  - Class: `AntibioticResistance`
    - **Description**: Represents antibiotic resistance.
    - **Children**:
      - Class: `AntibioticResistanceClass`
        - **Description**: Represents classes of antibiotic resistance.
      - Class: `AntibioticResistancePhenotype`
        - **Description**: Represents phenotypes of antibiotic resistance.
      - Class `AntibioticResistanceMechanism`
        - **Description**: Represents mechanisms of antibiotic resistance.
  - Class: `BiocideResistance`
    - **Description**: Represents biocide resistance.
    - **Children**:
      - Class: `BiocideClass`
        - **Description**: Represents classes of biocide resistance.
      - Class: `Biocide`
        - **Description**: Represents biocides involved in resistance.
  - Class: `MetalResistance`
    - **Description**: Represents metal resistance.
    - **Children**:
    - Class: `MetalClass`
      -  **Description**: Represents classes of metal resistance.
    - Class: `Metal`
      - **Description**: Represents metals involved in resistance.
  - Class: `UnclassifiedResistanceClass`
    - **Description**: Represents unclassified resistance classes.
    - **Children**:
      - Class: `UnclassifiedResistance`
        - **Description**: Represents unclassified resistance.

#### Class: `Protein`
- **Description**: Genes are translated into proteins, which was added as a new component of the PanRes database in version 2.0.
- **Parent Class**: `Resistance`
- **Children**:
  - Class: `PanProtein`
    - **Description**: Represents a PanProtein.
  - Class: `PanProteinCluster`
    - **Description**: Represents a cluster of PanProteins.
  - Class: `Structure`
    - **Description**: Overall category for protein structures in the PanRes database.
    - **Children**:
    - Class: `PanStructure`
      - **Description**: Represents a 3D protein structure of PanProtein.
    - Class: `PanStructureCluster`
      - **Description**: Represents a cluster of 3D predicted PanProtein clusters.

### Object Properties 
#### Property: `folds_to`
- **Description**: Links a PanProtein to its predicted 3D structure.
- **Connects**: `PanProtein` -> `PanStructure`

#### Property: `has_mechanism_of_resistance` 
- **Description**: Links a gene to its annotated mechanism of resistance.
- **Connects**: `PanGene` or `OriginalGene` -> `ResistanceMechaninsm`

#### Property: `same_as` 
- **Description**: Indicates that a PanGene corresponds to a specific OriginalGene from a source database.
- **Connects**: `PanGene`-> `OriginalGene`

#### Property: `translates_to` 
- **Description**: Links a PanGene to its translated protein product.
- **Connects**: `PanGene`-> `PanProtein`

#### Property: `has_pan_name` 
- **Description**: Inverse of same_as property. Links an OriginalGene to its standardized PanGene identifier.
- **Connects**: `OriginalGene`-> `PanGene`

#### Property: `has_predicted_phenotype` 
- **Description**: Links a gene to its predicted resistance phenotype.
- **Connects**: `PanGene` or `OriginalGene` -> `AntibioticResistancePhenotype` or `Biocide` or `Metal` or `UnclassifiedResistance`

#### Property: `has_resistance_class` 
- **Description**: Links a gene or phenotype to its broader resistance class.
- **Connects**: `PanGene` or `OriginalGene` or `AntibioticResistancePhenotype`  -> `AntibioticResistanceClass` or `BiocideClass` or `MetalClass` or `UnclassifiedResistanceClass`

#### Property: `is_discarded` 
- **Description**: Indicates that a PanGene was removed during database curation.
- **Connects**: `PanGene` -> `DiscardedPanGene`

#### Property: `is_from_database` 
- **Description**: Indicates the source database from which a gene originates.
- **Connects**: `PanGene` or `OriginalGene`-> `DiscardedPanGene`

#### Property: `found_in` 
- **Description**: Indicates which database contains a given resistance class.
- **Connects**: `AntibioticResistanceClass` or `BiocideClass` or `MetalClass` or `UnclassifiedResistanceClass` -> `Database`

#### Property: `member_of` 
- **Description**: Links a gene, protein, or structure to its corresponding cluster.
- **Connects**: `PanGene` or `PanProtein` or `PanStructure` -> `PanGeneCluster` or `PanProteinCluster` or `PanStructureCluster`

### Annotation Properties
#### Property: `has_length` [int]
- **Description**: Annotation property to specify the length of PanGenes and PanProteins.
- **Domain**: `PanGene`, `PanProtein`

#### Property: `accession` [str]
- **Description**: Annotation property to specify the accession number of genes and proteins.
- **Domain**: `PanGene`, `OriginalGene`, `Protein`

#### Property: `card_link` [str]
- **Description**: Specifies the corresponding CARD identifier or URL for a gene entry, if available.
- **Domain**: `PanGene`, `OriginalGene`

#### Property: `has_members` [int]
- **Description**: Specifies the number of members contained within a given cluster.
- **Domain**: `PanGeneCluster`, `PanProteinCluster`, `PanStructureCluster`

#### Property: `is_ecoli_homolog` [bool]
- **Description**: Indicates whether a PanProtein has a homolog in E. coli.
- **Domain**: `PanProtein`

#### Property: `original_fasta_header` [int]
- **Description**: Stores the original FASTA header associated with a PanGene entry.
- **Domain**: `PanGene`