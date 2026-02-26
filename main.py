import argparse
import os, re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
import copy
import logging
from datetime import datetime

parser = argparse.ArgumentParser(description="Brownaming: Propagating Sequence Names for Similar Organisms")
parser.add_argument('-p', '--proteins', help='FASTA file of query proteins')
parser.add_argument('-s', '--species', type=int, help='Taxonomy ID of the target species')
parser.add_argument('--threads', type=int, default=None, help='Number of threads (default: all available)')
parser.add_argument('--last-tax', type=int, default=None, help="(Taxonomy ID) Last taxonomic group for which homology searches will be performed")
parser.add_argument('--ex-tax', type=int, action='append', help='Taxonomy ID exclude from the research')
parser.add_argument('--swissprot-only', action='store_true', help='Use only SwissProt database for homology searches')
parser.add_argument('--local-db', help='Path to local database (optional if defined in LOCAL_DB_PATH env var)')
parser.add_argument('--resume', help='Resume a previous run using the run ID')
args = parser.parse_args()

if args.local_db:
    print(f"[INFO] Using local database path from command line argument: {args.local_db}")
    os.environ['LOCAL_DB_PATH'] = args.local_db

# Import utils after setting LOCAL_DB_PATH so it can read the environment variable
import utils, homology, excel, stats


if args.resume:
    RUN_ID = args.resume
    if not os.path.isdir(utils.working_dir(RUN_ID)):
        print(f"[ERROR] Run ID not found: {RUN_ID}")
        exit()
    # Setup logger for resume
    logger = utils.setup_logger(RUN_ID)
    logger.info(f"Resuming Brownaming with run ID: {RUN_ID}")
    state_args, state = utils.load_state(RUN_ID)
    if state_args:
        query_fasta = state_args.get('proteins')
        target_taxid = state_args.get('species')
        args.ex_tax = state_args.get('ex_tax')
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
        args = state['args']
        estimated_runtime = sum(estimated_runtime_list)
        estimated_hours = int(estimated_runtime // 60)
        estimated_minutes = int(estimated_runtime % 60)
        logger.info(f"Resuming from step {step} with {len(pending)} pending sequences")
        logger.info(f"Estimated remaining runtime: {estimated_hours:02d}:{estimated_minutes:02d} (hh:mm)")
else:
    query_fasta = args.proteins
    target_taxid = args.species
    if not os.path.isfile(query_fasta):
        print(f"[ERROR] File not found: {query_fasta}")
        exit()
    if target_taxid is None:
        print(f"[ERROR] Target species taxonomy ID is required.")
        exit()
        
    # Create run ID with format: yyyy-mm-dd-hh-mm-taxid
    timestamp = datetime.now().strftime("%Y-%m-%d-%H-%M")
    RUN_ID = f"{timestamp}-{target_taxid}"
    utils.create_run(RUN_ID)
    utils.save_state_args(args, RUN_ID)
    # Setup logger
    logger = utils.setup_logger(RUN_ID)
    logger.info(f"Starting the Brownaming process with run ID: {RUN_ID}")

working_directory = utils.working_dir(RUN_ID)
output_fasta_file = working_directory + '/' + os.path.basename(query_fasta).replace('.fasta', '_brownamed.fasta').replace('.faa', '_brownamed.fasta')
output_stats_file = working_directory + '/' + os.path.basename(query_fasta).replace('.fasta', '_brownaming_stats.png').replace('.faa', '_brownaming_stats.png')
output_excel_file = working_directory + '/' + os.path.basename(query_fasta).replace('.fasta', '_diamond_results.xlsx').replace('.faa', '_diamond_results.xlsx')
state_file = os.path.join(working_directory, f"state.pkl")
save_interval = 15 * 60
next_save = save_interval
parent = utils.get_parent_dict()
rank = utils.get_rank_dict()
children = utils.get_children_dict()
taxid2name = utils.get_taxid_to_scientificname()

excluded_tax = []
if args.ex_tax:
    for tax in args.ex_tax:
        excluded_tax += utils.get_children(tax, children)

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
    curr_tax_rank = rank.get(str(curr_tax),'unknown')
    
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
    if curr_tax == args.last_tax or curr_tax == 131567: # 131567: cellular organisms
        curr_tax = None
    else:
        curr_tax = parent.get(str(curr_tax))
        
    elapsed = time.time() - timer_start
    logger.info(f"Elapsed time: {elapsed/60:.2f} minutes")
    stats_data[f"Step {step}"]['elapsed_time'] = f"{elapsed/60:.2f}"

    if elapsed >= next_save:
        utils.save_state(state_file, assigned, pending, curr_tax, prev_group, step, stats_data, elapsed, 
                  query_fasta, target_taxid, query_ids, 
                  estimated_runtime_list, dbsizes, args)
            
        next_save = ((elapsed // save_interval) + 1) * save_interval

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
        re_description_search = re.findall(r" .* OS=", assigned[record.id][0].get("stitle",""))
        if (len(re_description_search)!=0):
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

os.remove(os.path.join('runs', RUN_ID, 'state_args.json'))