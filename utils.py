import json
import os
import re
from Bio import SeqIO
import numpy as np
import pickle
import requests
import time
import pandas as pd
import sklearn

# Global variable for custom working directory base path
# If None, defaults to script_dir()/runs/
WORKING_DIR_BASE = None

# Read local_db_path at module import time (global variable)
# Priority: LOCAL_DB_PATH env var > local_db_path in config.json
local_db_path = os.environ.get("LOCAL_DB_PATH")

if not local_db_path:
    config_file = os.path.join(os.path.dirname(__file__), "config.json")
    if os.path.exists(config_file):
        with open(config_file, "r") as f:
            local_db_path = json.load(f).get("local_db_path", None)

if not local_db_path:
    print("[ERROR] Database path not found. Set LOCAL_DB_PATH env var or local_db_path in config.json")
    exit()
    

parent_path = os.path.join(local_db_path, "taxonomy", "parent.json")
if os.path.exists(parent_path):
    with open(parent_path, 'r') as f:
        PARENT = json.load(f)
else:
    print(f"[ERROR] '{parent_path}' not found. Please build parent.json with load_nodes.py.")
    exit()
    
rank_path = os.path.join(local_db_path, "taxonomy", "rank.json")
if os.path.exists(rank_path):
    with open(rank_path, 'r') as f:
        RANK = json.load(f)
else:
    print(f"[ERROR] '{rank_path}' not found. Please build rank.json with load_nodes.py.")
    exit()
    
children_path = os.path.join(local_db_path, "taxonomy", "children.json")
if os.path.exists(children_path):
    with open(children_path, 'r') as f:
        CHILDREN = json.load(f)
else:
    print(f"[ERROR] '{children_path}' not found. Please build children.json with load_nodes.py.")
    exit()

taxid2name_path = os.path.join(local_db_path, "taxonomy", "taxid2scientific_name.json")
if os.path.exists(taxid2name_path):
    with open(taxid2name_path, 'r') as f:
        TAXID_TO_NAME = json.load(f)
else:
    print(f"[ERROR] '{taxid2name_path}' not found. Please build taxid2scientific_name.json with create_local_db.sh.")
    exit()

taxid2dbsize_path = os.path.join(local_db_path, "taxonomy", "taxid2dbsize.json")
if os.path.exists(taxid2dbsize_path):
    with open(taxid2dbsize_path, 'r') as f:
        TAXID_TO_DBSIZE = json.load(f)


def create_run(run_id):
    working_directory = working_dir(run_id)
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)
    return run_id


def setup_logger(run_id):
    """Setup logger that writes to both console and log file."""
    import logging
    
    logger = logging.getLogger('brownaming')
    logger.setLevel(logging.INFO)
    
    # Clear existing handlers
    logger.handlers.clear()
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_format = logging.Formatter('[%(levelname)s] %(message)s')
    console_handler.setFormatter(console_format)
    logger.addHandler(console_handler)
    
    # File handler
    log_file = os.path.join(working_dir(run_id), f'{run_id}.log')
    file_handler = logging.FileHandler(log_file, mode='a')
    file_handler.setLevel(logging.INFO)
    file_format = logging.Formatter('[%(levelname)s] %(message)s')
    file_handler.setFormatter(file_format)
    logger.addHandler(file_handler)
    
    return logger


def script_dir():
    return os.path.dirname(os.path.abspath(__file__))

def working_dir(run_id):
    if WORKING_DIR_BASE is not None:
        if os.path.basename(WORKING_DIR_BASE) == run_id:
            return WORKING_DIR_BASE
        return os.path.join(WORKING_DIR_BASE, run_id)
    return os.path.join(script_dir(), 'runs', run_id)

def get_parent_dict():
    return PARENT

def get_rank_dict():
    return RANK

def get_children_dict():
    return CHILDREN

def get_taxid_to_scientificname():    
    return TAXID_TO_NAME

def get_taxid_to_dbsize():
    return TAXID_TO_DBSIZE

def get_db_dmnd(swissprot_only):
    if swissprot_only:
        return os.path.join(local_db_path, "diamond", "uniprot_sprot.dmnd")
    return os.path.join(local_db_path, "diamond", "uniprot_all.dmnd")

def gene_name_from_stitle(stitle):
    # UniProt: look for " GN=gene_name "
    match = re.search(r" GN=([^ ]+)", stitle)
    if match:
        return match.group(1).strip()
    return ""

def get_children(taxid):
    all_children = {taxid}
    stack = [taxid]
    while stack:
        current = str(stack.pop())
        if current in CHILDREN:
            for child in CHILDREN[current]:
                if child not in all_children:
                    all_children.add(child)
                    stack.append(child)
    return all_children

def write_pending_fasta(src_faa, pending_ids, out_path):
    n = 0
    with open(out_path, 'w') as out_f:
        for rec in SeqIO.parse(src_faa, 'fasta'):
            if rec.id in pending_ids:
                SeqIO.write(rec, out_f, 'fasta')
                n += 1
    return n


def estimate_runtime(nb_query, target_taxid, last_tax=None, swissprot_only=False):
    predicted_times = []
    dbsizes = []
    curr_tax = target_taxid
    prev_group = None
    sum_previous_dbsize = 0
    while curr_tax:
        dbsize_count = count_sequence_from_taxid(curr_tax)
        current_dbsize = dbsize_count['swissprot' if swissprot_only else 'total']
        previous_dbsize = dbsizes[-1] if dbsizes else 0
        sum_previous_dbsize += previous_dbsize
        dbsizes.append(current_dbsize-sum_previous_dbsize)        
        
        predicted_time = predict_diamond_time(nb_query, dbsizes[-1])
        predicted_times.append(max(0.0, predicted_time))
        
        prev_group = curr_tax
        if curr_tax == last_tax or curr_tax == 131567: # 131567: cellular organisms
            curr_tax = None
        else:
            curr_tax = PARENT.get(str(curr_tax))
    
    return sum(predicted_times), predicted_times, dbsizes

def count_sequence_from_taxid(taxid):
    url = f"https://rest.uniprot.org/taxonomy/search?query=(tax_id:{taxid})&format=json&fields=statistics"

    response = requests.get(url)
    if response and response.json()["results"]:
        results = response.json()["results"]
        if results:
            swissprot_count = results[0].get("statistics", {}).get("reviewedProteinCount", 0)
            trembl_count = results[0].get("statistics", {}).get("unreviewedProteinCount", 0)
            total_count = swissprot_count + trembl_count
            return {
                "swissprot": swissprot_count,
                "total": total_count,
            }
    return {}

def predict_diamond_time(nb_query, dbsize):
    model_path = os.path.join(script_dir(), 'diamond_time_model.pkl')
    if not os.path.exists(model_path):
        # Model is optional - if not found, return a default estimate
        # Simple heuristic: ~0.1 second per query per 100k sequences in DB
        return (nb_query * dbsize / 100000) / 60  # Convert to minutes
    
    with open(model_path, 'rb') as f:
        model = pickle.load(f)
    
    features = pd.DataFrame([[nb_query, dbsize]], columns=['nb_query', 'dbsize'])    
    predicted_time = model.predict(features)[0]
    return predicted_time

def save_state_args(args, run_id):
    args_dict = vars(args)
    args_path = os.path.join(working_dir(run_id), 'state_args.json')
    with open(args_path, 'w') as f:
        json.dump(args_dict, f, indent=4)

def save_state(state_file, assigned, pending, curr_tax, prev_group, step, stats_data, elapsed, query_fasta, target_taxid, query_ids, estimated_runtime_list, dbsizes, args):
    for step, step_data in stats_data.items():
        if 'elapsed_time' not in step_data:
            stats_data.pop(step)
    state = {
        'assigned': assigned,
        'pending': pending,
        'curr_tax': curr_tax,
        'prev_group': prev_group,
        'step': step,
        'stats_data': stats_data,
        'elapsed': elapsed,
        'timer_start': time.time() - elapsed,
        'query_fasta': query_fasta,
        'target_taxid': target_taxid,
        'query_ids': query_ids,
        'estimated_runtime_list': estimated_runtime_list,
        'dbsizes': dbsizes,
        'args': args
    }
    with open(state_file, 'wb') as f:
        pickle.dump(state, f)
    print(f"[INFO] State saved at elapsed time: {elapsed/60:.2f} minutes", flush=True)

def load_state(run_id):
    try:
        state_args_file = os.path.join(working_dir(run_id), 'state_args.json')
        with open(state_args_file, 'r') as f:
            state_args = json.load(f)

        state_file = os.path.join(working_dir(run_id), 'state.pkl')
        if os.path.exists(state_file):
            with open(state_file, 'rb') as f:
                state = pickle.load(f)
            print(f"[INFO] Loaded state from elapsed time: {state['elapsed']/60:.2f} minutes", flush=True)
            
        else:
            state = None
        return state_args, state
    
    except (FileNotFoundError, pickle.UnpicklingError) as e:
        print(f"[ERROR] Could not load state: {e}", flush=True)
        return None
