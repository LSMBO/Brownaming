#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/config.json"

# Command-line argument handling
REFRESH_MODE=false
for arg in "$@"; do
  case $arg in
    --refresh)
      REFRESH_MODE=true
      echo "[INFO] Refresh mode enabled - reinstalling database"
      shift
      ;;
  esac
done

for bin in jq curl diamond awk gunzip; do
  command -v "$bin" >/dev/null || { echo "Missing dependency: $bin" >&2; exit 1; }
done

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Config file not found: $CONFIG_FILE" >&2
  exit 1
fi

LOCAL_DB_PATH=$(jq -r '.local_db_path' "$CONFIG_FILE")
if [[ -z "$LOCAL_DB_PATH" || "$LOCAL_DB_PATH" == "null" ]]; then
    echo "local_db_path not found in config file" >&2
    exit 1
fi

# Remove old database if refresh mode is enabled
if [[ "$REFRESH_MODE" == true && -d "$LOCAL_DB_PATH" ]]; then
    echo "[INFO] Removing existing database at $LOCAL_DB_PATH"
    rm -rf "$LOCAL_DB_PATH"
fi

mkdir -p "${LOCAL_DB_PATH}/uniprot"/{fasta,taxonomy,mapping,diamond}

echo "[INFO] Download Swiss-Prot"
uniprot_sprot_file="${LOCAL_DB_PATH}/uniprot/fasta/uniprot_sprot.fasta"
curl -L --fail --retry 3 https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz \
  | gunzip -c > $uniprot_sprot_file

echo "[INFO] Download TrEMBL (~5h)"
uniprot_all_file="${LOCAL_DB_PATH}/uniprot/fasta/uniprot_all.fasta"
curl -L --fail --retry 3 https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz \
  | gunzip -c > $uniprot_all_file

echo "[INFO] Concatenate"
cat "${LOCAL_DB_PATH}/uniprot/fasta/uniprot_sprot.fasta" >> "${LOCAL_DB_PATH}/uniprot/fasta/uniprot_all.fasta"

echo "[INFO] Build taxonmap.tsv"
taxonmap_file="${LOCAL_DB_PATH}/uniprot/mapping/taxonmap.tsv"
printf "accession\taccession.version\ttaxid\tgi\n" > $taxonmap_file
awk -v OFS='\t' '
  /^>/{
    acc=""; tax="";
    split($1,a,"|"); if (length(a)>=2) acc=a[2];
    if (match($0,/OX=([0-9]+)/,m)) tax=m[1];
    if (acc!="" && tax!="") printf "%s\t%s\t%s\t0\n", acc, acc, tax;
  }
' "${LOCAL_DB_PATH}/uniprot/fasta/uniprot_all.fasta" >> $taxonmap_file

lines=$(($(wc -l < $taxonmap_file) - 1))
if (( lines == 0 )); then
  echo "Empty taxonmap.tsv (no headers parsed)" >&2
  exit 1
fi

echo "[INFO] Download taxonomy dump"
nodes="${LOCAL_DB_PATH}/uniprot/taxonomy/nodes.dmp"
names="${LOCAL_DB_PATH}/uniprot/taxonomy/names.dmp"
curl -L --fail --retry 3 -o "${LOCAL_DB_PATH}/uniprot/taxonomy/taxdump.tar.gz" \
  https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
tar -C "${LOCAL_DB_PATH}/uniprot/taxonomy" -xzf "${LOCAL_DB_PATH}/uniprot/taxonomy/taxdump.tar.gz" $names $nodes

if [[ ! -s "$nodes" || ! -s "$names" ]]; then
  echo "Taxonomy files missing" >&2
  exit 1
fi

echo "[INFO] Make diamond db (SwissProt only)"
uniprot_sprot_building="${LOCAL_DB_PATH}/uniprot/fasta/uniprot_sprot.building"
diamond makedb -p "$(nproc)" \
  -d $uniprot_sprot_building \
  --in $uniprot_sprot \
  --taxonmap $taxonmap_file \
  --taxonnodes $nodes \
  --taxonnames $names
diamond dbinfo -d $uniprot_sprot_building > /dev/null \
  && mv -f "$uniprot_sprot_building.dmnd" "$LOCAL_DB_PATH/uniprot/diamond/uniprot_sprot.dmnd"


echo "[INFO] Make diamond db (SwissProt and TrEMBL) (~3h)"
uniprot_all_building="${LOCAL_DB_PATH}/uniprot/fasta/uniprot_all.building"
diamond makedb -p "$(nproc)" \
  -d $uniprot_all_building \
  --in $uniprot_all \
  --taxonmap $taxonmap_file \
  --taxonnodes $nodes \
  --taxonnames $names
diamond dbinfo -d $uniprot_all_building > /dev/null \
  && mv -f "$uniprot_all_building.dmnd" "$LOCAL_DB_PATH/uniprot/diamond/uniprot_all.dmnd"

echo "[INFO] Generate taxonomy JSON helpers"
python "${SCRIPT_DIR}/create_taxonomy_json.py"

echo "[DONE]"
echo
if [[ "$REFRESH_MODE" == true ]]; then
  echo "Database successfully refreshed in $LOCAL_DB_PATH"
else
  echo "Database successfully created in $LOCAL_DB_PATH"
  echo "To refresh the database in the future, use: $0 --refresh"
fi
