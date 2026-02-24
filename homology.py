import utils
import os
import shutil
import subprocess

def which_or_die(bin_name):
    path = shutil.which(bin_name)
    if not path:
        print(f"[ERROR] {bin_name} not found in PATH.", flush=True)
        exit()
    return path

def build_taxon_list(curr_tax, excluded_tax):
    taxon_list = []
    if utils.CHILDREN.get(str(curr_tax)):
        for child in utils.CHILDREN[str(curr_tax)]:
            if child != excluded_tax:
                taxon_list.append(child)
    else:
        taxon_list.append(curr_tax)
    return taxon_list

def run_diamond(run_id, query_fasta, taxonlist, group, threads=None, max_targets=50, mode="more-sensitive", excluded_tax=[], swissprot_only=False):
    if group[0] in excluded_tax:
        return []
    
    diamond = which_or_die("diamond")
    out_path = os.path.join(utils.working_dir(run_id), f".diamond_tmp_{os.getpid()}.tsv")
    args = [
        diamond, "blastp",
        "-d", utils.get_db_dmnd(swissprot_only),
        "-q", query_fasta,
        "-k", str(max_targets),
        "-e", "1e-5",
        "-p", str(threads or os.cpu_count() or 1),
        "--" + mode,
        "-f", "6",
        "qseqid","sseqid","pident","ppos","length","evalue","bitscore","qlen","slen","staxids","stitle",
        "-o", out_path
    ]
    args.extend(["--taxonlist", ",".join(str(t) for t in taxonlist)])
   
    print("[INFO] Running DIAMOND:\n", " ".join(args), flush=True)
    proc = subprocess.Popen(args, stdout=None, stderr=subprocess.PIPE, text=True)
    _, stderr = proc.communicate()
    if proc.returncode != 0:
        msg = stderr.strip() or "Unknown error"
        print(f"[ERROR] DIAMOND failed: {msg}", flush=True)
        exit()

    hits = parse_diamond_tsv(out_path, group, excluded_tax)
    try:
        os.remove(out_path)
    except OSError:
        pass
    return hits


def parse_diamond_tsv(path, ancestor, excluded_tax):
    hits = []
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            qseqid = parts[0]
            sseqid = parts[1]
            pident = float(parts[2])
            ppos = float(parts[3])
            alen = int(parts[4])
            evalue = float(parts[5])
            bits = float(parts[6])
            qlen = int(parts[7])
            slen = int(parts[8])
            stax_raw = parts[9]
            staxid = None
            if stax_raw:
                staxid = int(stax_raw.split(";")[0])
            if staxid in excluded_tax:
                continue
            stitle = parts[10] if len(parts) > 10 else ""
            hits.append({
                "qseqid": qseqid, "sseqid": sseqid, "pident": pident,
                "ppos": ppos, "alen": alen, "evalue": evalue, "bits": bits,
                "qlen": qlen, "slen": slen, "staxid": staxid, "stitle": stitle,
                "common_ancestor_taxid": ancestor[0], "common_ancestor_name": ancestor[1],
                "common_ancestor_rank": ancestor[2]
            })
    return hits

def select_best_by_priority(hits, target_taxid, step,
                            min_pid=0, min_qcov=0, min_scov=0, min_bits=50.0):
    best_per_query = {}
    for h in hits:
        if h["qlen"] <= 0 or h["slen"] <= 0:
            continue
        qcov = h["alen"] / h["qlen"]
        scov = h["alen"] / h["slen"]
        if h["pident"] < min_pid or qcov < min_qcov or scov < min_scov or h["bits"] < min_bits:
            continue
        h["_qcov"] = qcov
        h["_scov"] = scov
        
        # Prioriry: Highest bitscore, then higher identity
        key = (step, -h["bits"], -h["pident"])
        h["_key"] = key
        q = h["qseqid"]

        if q not in best_per_query:
            best_per_query[q] = []
        if len(best_per_query[q]) < 3:
            best_per_query[q].append(h)
        elif key < best_per_query[q][-1]["_key"]:
            best_per_query[q][-1] = h
            best_per_query[q].sort(key=lambda x: x["_key"])
            
    return best_per_query