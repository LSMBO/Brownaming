import argparse
import os, re
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
import copy
from datetime import datetime
import utils, homology, excel, stats

parser = argparse.ArgumentParser(description="Brownaming: Propagating Sequence Names for Similar Organisms")
parser.add_argument('-p', '--proteins', help='FASTA file of query proteins')
parser.add_argument('-s', '--species', type=int, help='Taxonomy ID of the target species')
parser.add_argument('--threads', type=int, default=None, help='Number of threads (default: all available)')
parser.add_argument('--last-tax', type=int, default=None, help="(Taxonomy ID) Last taxonomic group for which homology searches will be performed")
parser.add_argument('--ex-tax', type=int, action='append', help='Taxonomy ID exclude from the research')
parser.add_argument('--swissprot-only', action='store_true', help='Use only SwissProt database for homology searches')
parser.add_argument('--local-db', help='Path to local database (optional if defined in LOCAL_DB_PATH env var)')
parser.add_argument('--working-dir', help='Final output directory (optional, run still executes in runs/YYYY-MM-DD-HH-MM-TAXID)')
parser.add_argument('--resume', help='Resume a previous run using the run ID')
args = parser.parse_args()


def error_exit(message, run_id=None):
    if run_id:
        print(f"[ERROR] [run_id={run_id}] {message}")
    else:
        print(f"[ERROR] {message}")
    exit(1)


RUN_ID = None
state = None
final_output_dir = None

if args.resume:
    RUN_ID = str(args.resume)
    run_working_dir = utils.working_dir(RUN_ID)

    if not os.path.isdir(run_working_dir):
        error_exit("Run directory not found in Brownaming/runs.", RUN_ID)

    state_args, state = utils.load_state(RUN_ID)
    print(state['step'])
    if not state_args:
        error_exit("Could not load resume state from state_args.json.", RUN_ID)

    query_fasta = state_args.get('proteins')
    target_taxid = state_args.get('species')
    args.ex_tax = state_args.get('ex_tax')
    args.last_tax = state_args.get('last_tax')
    args.swissprot_only = state_args.get('swissprot_only', False)
    args.local_db = state_args.get('local_db')
    args.threads = state_args.get('threads')
    final_output_dir = state_args.get('working_dir')

    logger = utils.setup_logger(RUN_ID)
    logger.info(f"Resuming Brownaming with run ID: {RUN_ID}")

    if state:
        assigned = state['assigned']
        pending = state['pending']
        curr_tax = state['curr_tax']
        prev_group = state['prev_group']
        step = state['step']
        stats_data = state['stats_data']
        elapsed = state['elapsed']
        timer_start = state['timer_start']
        query_ids = state['query_ids']
        estimated_runtime_list = state['estimated_runtime_list']
        dbsizes = state['dbsizes']
        saved_args = state.get('args')
        if saved_args:
            args.ex_tax = getattr(saved_args, 'ex_tax', args.ex_tax)
            args.last_tax = getattr(saved_args, 'last_tax', args.last_tax)
            args.swissprot_only = getattr(saved_args, 'swissprot_only', args.swissprot_only)
            args.threads = getattr(saved_args, 'threads', args.threads)
        estimated_runtime = sum(estimated_runtime_list)
        estimated_hours = int(estimated_runtime // 60)
        estimated_minutes = int(estimated_runtime % 60)
        logger.info(f"Resuming from step {step} with {len(pending)} pending sequences")
        logger.info(f"Estimated remaining runtime: {estimated_hours:02d}:{estimated_minutes:02d} (hh:mm)")

else:
    query_fasta = args.proteins
    target_taxid = args.species

    if not os.path.isfile(query_fasta):
        error_exit(f"File not found: {query_fasta}")
    if target_taxid is None:
        error_exit("Target species taxonomy ID is required.")

    timestamp = datetime.now().strftime("%Y-%m-%d-%H-%M")
    RUN_ID = f"{timestamp}-{target_taxid}"

    if args.working_dir:
        final_output_dir = os.path.abspath(args.working_dir)
        args.working_dir = final_output_dir
        print(f"[INFO] Final output directory requested for run_id {RUN_ID}: {final_output_dir}")

    utils.create_run(RUN_ID)
    utils.save_state_args(args, RUN_ID)

    logger = utils.setup_logger(RUN_ID)
    logger.info(f"Starting the Brownaming process with run ID: {RUN_ID}")

working_directory = utils.working_dir(RUN_ID)
output_fasta_file = working_directory + '/' + os.path.basename(query_fasta).replace('.fasta', '_brownamed.fasta').replace('.faa', '_brownamed.fasta')
output_stats_file = working_directory + '/' + os.path.basename(query_fasta).replace('.fasta', '_brownaming_stats.png').replace('.faa', '_brownaming_stats.png')
output_excel_file = working_directory + '/' + os.path.basename(query_fasta).replace('.fasta', '_diamond_results.xlsx').replace('.faa', '_diamond_results.xlsx')
state_file = os.path.join(working_directory, f"state.pkl")
# save_interval = 15 * 60
save_interval = 5
next_save = save_interval

if args.local_db:
    utils.LOCAL_DB_PATH = args.local_db
    print(f"[INFO] Using local database path from command line argument: {args.local_db}")
else:
    utils.LOCAL_DB_PATH = utils.set_local_db_path()
    args.local_db = utils.LOCAL_DB_PATH
if not args.local_db:
    error_exit("Local database path must be provided either through --local-db argument or set in config.json.", RUN_ID)
                
utils.PARENT = utils.set_parent_dict()
parent = utils.get_parent_dict()

utils.RANK = utils.set_rank_dict()
rank = utils.get_rank_dict()

utils.CHILDREN = utils.set_children_dict()
children = utils.get_children_dict()

utils.TAXID_TO_NAME = utils.set_taxid_to_scientificname()
taxid2name = utils.get_taxid_to_scientificname()

utils.TAXID_TO_DBSIZE = utils.set_taxid_to_dbsize()

excluded_tax = []
if args.ex_tax:
    for tax in args.ex_tax:
        excluded_tax += utils.get_children(tax)

if not args.resume or not state:
    query_ids = [rec.id for rec in SeqIO.parse(query_fasta, 'fasta')]
    estimated_runtime, estimated_runtime_list, dbsizes = utils.estimate_runtime(len(query_ids), target_taxid, last_tax=args.last_tax, swissprot_only=args.swissprot_only)
    estimated_hours = int(estimated_runtime // 60)
    estimated_minutes = int(estimated_runtime % 60)
    logger.info(f"Estimated total runtime: {estimated_hours:02d}:{estimated_minutes:02d} (hh:mm)")

    assigned = {}
    pending = set(query_ids)
    curr_tax = target_taxid
    prev_group = None
    step = 0
    stats_data = {}
    timer_start = time.time()

while curr_tax is not None and pending:
    step += 1
    tmp_fasta = os.path.join(working_directory, f".pending_{os.getpid()}_{step}.fasta")
    curr_tax_name = taxid2name.get(str(curr_tax), "unknown")
    curr_tax_rank = rank.get(str(curr_tax), 'unknown')
    
    if len(pending) == len(query_ids):
        tmp_fasta = query_fasta
        n_written = len(pending)
    else:
        n_written = utils.write_pending_fasta(query_fasta, pending, tmp_fasta)

    if n_written > 0:
        logger.info(
            f"Step {step}: Searching among {dbsizes[step-1]} sequences of {curr_tax_name} "
            f"({curr_tax} ; {curr_tax_rank}) with {n_written} pending sequences "
            f"(estimated runtime={estimated_runtime_list[step-1]:.2f} minutes)..."
        )
        stats_data[f"Step {step}"] = {
            'dbsize': dbsizes[step-1],
            'taxon_name': curr_tax_name,
            'taxon_id': curr_tax,
            'rank': curr_tax_rank,
            'nb_query': n_written,
            'estimated_runtime': f"{estimated_runtime_list[step-1]:.2f}"
        }
        input_taxon_list = homology.build_taxon_list(curr_tax, prev_group)
        if not input_taxon_list:
            logger.info(f"Step {step}: Subject database empty, continue to upper taxon")
            stats_data[f"Step {step}"]['prots_with_hit'] = stats_data.get(f"Step {step-1}", {}).get('prots_with_hit', 0)
        else:
            hits = homology.run_diamond(
                RUN_ID,
                tmp_fasta,
                input_taxon_list,
                (curr_tax, curr_tax_name, curr_tax_rank),
                threads=args.threads,
                max_targets=50,
                mode="more-sensitive",
                excluded_tax=excluded_tax,
                swissprot_only=args.swissprot_only
            )
            best = homology.select_best_by_priority(hits, target_taxid, step)
            assigned.update(best)
            logger.info(f"Step {step}: Found a satisfying hit for {len(assigned)} proteins")
            stats_data[f"Step {step}"]['prots_with_hit'] = len(assigned)
            pending -= {key for key, value in best.items() if len(value) < 3}

        if tmp_fasta != query_fasta:
            try:
                os.remove(tmp_fasta)
            except OSError:
                pass

    prev_group = curr_tax
    if curr_tax == args.last_tax or curr_tax == 131567:
        curr_tax = None
    else:
        curr_tax = parent.get(str(curr_tax))

    elapsed = time.time() - timer_start
    logger.info(f"Elapsed time: {elapsed/60:.2f} minutes")
    stats_data[f"Step {step}"]['elapsed_time'] = f"{elapsed/60:.2f}"
    print('elapsed :', elapsed, '>= next_save :', next_save)
    if elapsed >= next_save:
        utils.save_state(
            state_file,
            assigned,
            pending,
            curr_tax,
            prev_group,
            step,
            stats_data,
            elapsed,
            query_fasta,
            target_taxid,
            query_ids,
            estimated_runtime_list,
            dbsizes,
            args
        )
        next_save = ((elapsed // save_interval) + 1) * save_interval
        print('next save :', next_save)
stats.generate_combined_figure(stats_data, output_file=output_stats_file)

output_data = {
    "Query accession": [],
    "Subject accession": [],
    "Subject description": [],
    "Subject species (taxid)": [],
    "Subject species (name)": [],
    "Gene Name": [],
    "Bitscore": [],
    "Evalue": [],
    "Identity (%)": [],
    "Similarity (%)": [],
    "Query coverage (%)": [],
    "Subject coverage (%)": [],
    "Common ancestor (rank)": [],
    "Common ancestor (taxID)": [],
    "Common ancestor (name)": [],
    "Hit found": []
}
output_top3 = copy.deepcopy(output_data)

for qid in query_ids:
    if qid in assigned:
        output_data = excel.add_hit(output_data, assigned[qid][0])
        for hit in assigned[qid]:
            output_top3 = excel.add_hit(output_top3, hit)
    else:
        output_data = excel.add_no_hit(output_data, qid)
        output_top3 = excel.add_no_hit(output_top3, qid)

excel.write_excel(output_data, output_excel_file)
excel.add_sheet(output_top3, output_excel_file, "Top3 hits")

output_records = []
for record in SeqIO.parse(query_fasta, "fasta"):
    new_description = "Uncharacterized protein"
    if record.id in assigned:
        re_description_search = re.findall(r" .* OS=", assigned[record.id][0].get("stitle", ""))
        if len(re_description_search) != 0:
            record_description = re_description_search[0][1:-4]
            new_description = f"{record_description} FROM {taxid2name.get(str(assigned[record.id][0].get('staxid')), '')}"

    rec = SeqRecord(
        Seq(str(record.seq).upper()),
        id=record.id,
        description=new_description
    )
    output_records.append(rec)

with open(output_fasta_file, "w") as f:
    SeqIO.write(output_records, f, "fasta")

if final_output_dir:
    internal_run_dir = utils.working_dir(RUN_ID)
    destination_dir = os.path.abspath(final_output_dir)
    if os.path.abspath(internal_run_dir) != destination_dir:
        if os.path.exists(destination_dir):
            error_exit(f"Cannot move completed run to '{destination_dir}' because destination already exists.", RUN_ID)
        destination_parent = os.path.dirname(destination_dir)
        if destination_parent:
            os.makedirs(destination_parent, exist_ok=True)
        shutil.move(internal_run_dir, destination_dir)
        logger.info(f"Run {RUN_ID} moved to custom output directory: {destination_dir}")