# Brownaming v2.0.0

Taxonomy‑aware protein naming by fast homology search with DIAMOND. Brownaming assigns concise names to predicted proteins by prioritizing homologs from the closest available taxa and only expanding outward in the taxonomy when needed.

## About
Goal: give each query protein the most specific, biologically meaningful name supported by homology.

## Features
* Iterative taxonomy expansion: start at the target species; if no hit, move up (genus → family → order ...)
* Prioritize candidates by: (1) taxonomy distance (only varies across iterations) then (2) DIAMOND bitscore, then (3) % identity
* If nothing found up to configured rank → "Uncharacterized protein"
* Fully local: UniProt Swiss‑Prot + TrEMBL DIAMOND database with taxonomy metadata
* Optional exclusion of unwanted taxa

## Requirements
* Linux environment for database build (tested on Ubuntu; search also works under WSL)
* DIAMOND ≥ 2.1.x compiled with taxonomy support
* curl, awk, gzip, coreutils (for the provided build script)
* Disk space: ≈ 120–150 GB (TrEMBL + Swiss‑Prot FASTA + DIAMOND index + temporary files)
* RAM: ≥ 32 GB recommended for faster makedb sorting (less will still work, slower)
* Python ≥ 3.9
* Conda or Mamba (recommended for dependency management)

### Python Dependencies
See `environment.yml` for the complete list of dependencies. Main packages include:
* BioPython
* NumPy, Pandas
* scikit-learn
* openpyxl
* matplotlib
* requests

## Installation

### Database Setup (Required for Both Conda and Docker)

Before using Brownaming, you must create the local DIAMOND database:

1. Clone repository:
```bash
git clone <repo_url>
cd Brownaming
```

2. Create or edit `config.json`:
```json
{
    "local_db_path": "/path/to/brownaming_db"
}
```

3. Run the database build script:
```bash
./create_local_db.sh
```

What it does:
* Downloads UniProt Swiss‑Prot + TrEMBL (current release)
* Extracts TaxIDs from FASTA headers (OX=)
* Generates `taxonmap.tsv`, taxonomy JSON caches (parent/rank/children)
* Builds two DIAMOND databases:
  - full (Swiss‑Prot + TrEMBL)
  - swissprot (Swiss‑Prot only)

Duration: ~8 h

---

### Option 1: Conda/Mamba

1. Create conda environment:
```bash
conda env create -f environment.yml
conda activate brownaming
```

#### Running Brownaming

```bash
# Basic run
python main.py -p /path/to/query.fasta -s 83333 --threads 16

# SwissProt only
python main.py -p /path/to/query.fasta -s 83333 --swissprot-only

# Resume run
python main.py --resume 2026-02-24-14-30-83333
```

Results saved in: `runs/YYYY-MM-DD-HH-MM-TAXID/`

#### Improving Runtime Predictions

After running several analyses, update the time prediction model:

```bash
# Activate environment (if not already active)
conda activate brownaming

# Generate updated dataset from all run logs
cd time_prediction_model/
python create_data.py

# Retrain the model
python train_model.py
```

The more analyses you run, the more accurate the time estimates become.
---

### Option 2: Docker

Docker provides an isolated, reproducible environment.

#### Build Docker image
```bash
docker build -t brownaming .
```

#### Configuration

The wrapper script `brownaming-compose` automatically:
- Reads database path from `config.json`
- Detects and mounts the directory containing your input file

Make the script executable:
```bash
chmod +x brownaming-compose
```

#### Running Brownaming

Use **absolute paths** for your input files:

```bash
# Basic run
./brownaming-compose run --rm brownaming \
  python main.py -p /absolute/path/to/query.fasta -s 83333 --threads 16

# SwissProt only
./brownaming-compose run --rm brownaming \
  python main.py -p /absolute/path/to/query.fasta -s 83333 --swissprot-only

# Resume run
./brownaming-compose run --rm brownaming \
  python main.py --resume 2026-02-24-14-30-83333
```

Results saved in: `runs/YYYY-MM-DD-HH-MM-TAXID/`

#### Improving Runtime Predictions

After running several analyses, update the time prediction model:

```bash
# Generate updated dataset from all run logs
./brownaming-compose run --rm brownaming \
  python time_prediction_model/create_data.py

# Retrain the model
./brownaming-compose run --rm brownaming \
  python time_prediction_model/train_model.py
```

The updated model is immediately available - no rebuild needed! The more analyses you run, the more accurate the time estimates become.


## Command‑Line Arguments
Required:
* -p / --proteins <file> : Query protein FASTA
* -s / --species <taxid> : NCBI TaxID of target species (root of initial search)

Optional:
* --threads <N> : DIAMOND threads (default: all)
* --last-tax <taxid> : Stop expanding after this specific TaxID is reached.
* --ex-tax <taxid> : TaxID to exclude. For multiple exclusions, use this flag multiple times; each instance excludes the specified taxon and its subtree.
* --swissprot-only: Run DIAMOND searches only on the SwissProt database.

## Outputs
* **_query_file_name_**_brownamed.fasta : FASTA file with updated headers containing the assigned names.  
* **_query_file_name_**_diamond_results.xlsx : Excel table listing, for each query protein, the match used for naming, including homology scores (identity, evalue, bitscore, ...), and the rank and name of the last common ancestor.
* **_query_file_name_**_brownaming_stats.png : Statistics figure showing the progression through taxonomic ranks.
* **YYYY-MM-DD-HH-MM-TAXID.log** : Complete log file of the run (in the run directory).

## Updating the Databases
Re-run:
./create_local_db.sh --refresh

