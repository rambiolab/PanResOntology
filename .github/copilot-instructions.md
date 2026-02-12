# PanRes Ontology - AI Coding Agent Instructions

## Project Overview
PanRes Ontology is an OWL-based knowledge base that consolidates antimicrobial resistance gene annotations from 10+ external databases (ResFinder, CARD, MegaRes, AMRFinderPlus, etc.) into a unified semantic structure. The ontology represents genes, proteins, resistance phenotypes/classes, and databases as interconnected entities.

## Architecture & Data Flow

### Core Components
- **`code/model.py`**: Defines the OWL ontology schema using `owlready2` - classes for Genes (PanGene/OriginalGene), Databases, ResistanceTypes (Antibiotic/Metal/Biocide), and functional properties linking them
- **`code/PanResOntology.py`**: Main build script orchestrating ontology construction in sequence:
  1. Load schema from model.py
  2. Load resistance targets (drugs/metals/biocides) from targets.xlsx via `targets.py`
  3. Ingest data from 10 external databases sequentially via `code/databases/*.py` modules
  4. Clean up unused classes and reclassify genes based on phenotypes
  5. Export to `ontology/panres_v2.owl`

- **`code/databases/`**: One module per external database (e.g., `panres.py`, `card.py`, `megares.py`)
  - Each implements `add_*_annotations()` function following a common pattern:
    - Parse TSV/CSV data file
    - Map external gene names to PanGene identifiers
    - Extract resistance targets (phenotype/class/mechanism)
    - Call `gene_target()` in `targets.py` to link genes to targets
    - Log operations via `loguru` logger

- **`code/functions.py`**: Reusable utilities:
  - `get_or_create_instance/subclass()`: Idempotent OWL entity creation with name normalization (spaces/hyphens → underscores)
  - `clean_gene_name()`: Database-specific header parsing (handles `|` delimiters, parentheses, regex patterns)
  - Query helpers (`get_genes_from_database()`, `class_to_genes()`, `summarise_classes()`)
  - Visualization utilities using Graphviz

- **`code/export.py`**: CLI tool extracting gene metadata to CSV with configurable columns (phenotype, class, accession, etc.)

### Key Data Structures
```
PanGene (unified identifier: "pan_123") 
├─ same_as → [OriginalGene] (database-specific names)
├─ is_from_database → [Database]
├─ has_predicted_phenotype → [Phenotype]
├─ has_resistance_class → [ResistanceClass]
├─ has_mechanism_of_resistance → [Mechanism]
├─ accession → [str]
└─ pubmed → [str]

ResistanceType hierarchy:
├─ AntibioticResistancePhenotype (e.g., "Ampicillin")
│  └─ is_a → AntibioticResistanceClass (e.g., "Beta_Lactam")
├─ Metal (e.g., "Mercury")
│  └─ metal_symbol, metal_comment
└─ Biocide (e.g., "Triclosan")
   └─ is_a → BiocideClass
```

## Developer Workflows

### Building the Ontology from Scratch
```bash
# Install dependencies
pip install owlready2==0.44 loguru==0.7.3 pandas==1.5.3 numpy==1.24.3

# Ensure all data files are present in data/
python code/PanResOntology.py
# Output: ontology/panres_v2.owl
# Logs: panres_messages.log
```

### Querying the Ontology
See `notebooks/Queries.ipynb` for examples. Basic pattern:
```python
from owlready2 import get_ontology
import model
from functions import get_genes_from_database

onto = get_ontology("ontology/panres_v2.owl").load()
genes_df = get_genes_from_database(onto, "CARD")
```

### Adding a New Database
1. Create `code/databases/newdb.py` with `add_newdb_annotations(file, onto, logger)` function
2. Map external gene names to existing PanGenes using `get_instance()` or create new ones
3. Extract targets and call `gene_target(gene, og, target, onto, db_name)` from `targets.py`
4. Import in `PanResOntology.py` and call function in build sequence
5. Update `data/README.md` with data source location

### Exporting Gene Data
```bash
python code/export.py -f ontology/panres_v2.owl -o output.csv --columns name accession has_predicted_phenotype is_from_database
```

## Critical Patterns & Conventions

### Naming & IRI Resolution
- All entity names use underscore naming: spaces/hyphens converted via `.replace(" ", "_").replace("-", "_")`
- IRIs constructed as: `onto.base_iri + "EntityName"`
- Search via: `onto.search_one(iri="*name_pattern")` or `onto.search(type=Class)`
- Use `get_instance(onto, name)` and `get_or_create_instance()` to abstract this

### Gene Target Assignment Pattern
Database modules call `gene_target(gene, og, target, onto, db_name)` which:
1. Normalizes target name to underscore format
2. Distinguishes phenotypes (single drug/biocide/metal) from classes and mechanisms
3. Handles drug combinations (e.g., "Ampicillin+Sulbactam") by splitting on `+`
4. Establishes bidirectional links: gene ↔ phenotype, phenotype ↔ class
5. Records source database via `target_instance.found_in`

### Database-Specific Gene Name Parsing
Each database has distinct FASTA header format. Use `clean_gene_name(gene_name, db)` which applies regex for:
- **AMRFinderPlus**: Extract field [5] from pipe-delimited header
- **CARD**: Field [5] before `[` bracket
- **MegaRes**: Field [4]
- **ARG-ANNOT**: Everything after first `)` 
- **MetalRes**: First space-delimited field
- Add new patterns as `elif db == 'newdb': return ...` branches

### Ontology Class Hierarchy Pruning
After loading all databases, `targets.py:remove_unused_subclasses_with_property()` deletes ontology classes (e.g., drug phenotypes) that have no gene instances linked via specified property, reducing bloat.

### Logging
Use `loguru` logger (configured in `PanResOntology.py` to write to `panres_messages.log`):
```python
logger.success("Message")  # Green, for successful operations
logger.info("Message")      # Blue, for general info
logger.error("Message")     # Red, for errors
```

## File Organization & Responsibilities
- **`data/*.tsv/.csv`**: Input data files from external sources (see `data/README.md` for URLs)
- **`data/discarded/`**: Genes removed due to sequence validation failures
- **`data/proteins/`**: FASTA sequences and clustering output
- **`ontology/`**: Output OWL files
- **`notebooks/`**: Jupyter examples for querying ontology
- **`docs/model.md`**: Extended documentation of classes/properties (incomplete—marked as TODO)

## Important Notes

### Data Integrity
- PanRes v1.0.0 is the seed dataset; all new genes mapped against this
- DiscardedPanGenes tracked separately (not removed) for traceability
- Gene clusters (50-90% identity) stored in protein FASTA header metadata for protein-level analysis

### Known Limitations / TODOs
- Protein translations not yet fully linked to genes (`translates_to` property exists but incomplete)
- PubMed IDs extraction via `accession_to_pubmed()` not yet integrated into pipeline
- Protein structures (`folds_to` property) schema defined but unused
- Gene alt names from some databases not fully populated

### Reasoner Configuration
Uses `sync_reasoner(debug=0, infer_property_values=True)` at build time (not Pellet) to materialize inferred relationships. Commented-out Pellet option available but slower.

## External Dependencies & References
- **owlready2**: OWL manipulation library (note: version 0.44 pinned)
- **loguru**: Structured logging with context
- **pandas**: Data frame operations for TSV/CSV parsing
- **Graphviz**: Visualization (optional, only used in notebooks)
- **NCBI E-utilities**: `accession_to_pubmed()` uses `esearch|elink|efetch` CLI tools (not integrated)
