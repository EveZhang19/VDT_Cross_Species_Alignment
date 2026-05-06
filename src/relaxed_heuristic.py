import pandas as pd
import networkx as nx
import os, itertools, math

BASE_DIR = os.path.expanduser("")
GRN_DIR = os.path.join(BASE_DIR, "")
ANNO_DIR = os.path.join(BASE_DIR, "")
ORTHO_FILE = os.path.join(BASE_DIR, "")
SPECIES = ['maize', 'sorghum', 'ath']

BETA = 0.4
GAMMA = 0.6
EPSILON = 0.01 

GLOBAL_MAX_ST = 0.2478

discovered_modules = []
trajectory_data = []

grn_graphs = {sp: nx.DiGraph() for sp in SPECIES}
for sp in SPECIES:
    df = pd.read_csv(os.path.join(GRN_DIR, f"Ensemble_GRN_{sp}.csv"), sep='\t')
    for _, row in df.iterrows(): grn_graphs[sp].add_edge(row['Gene1'], row['Gene2'], weight=row['Ensemble_Pred'])

og_df = pd.read_csv(ORTHO_FILE, sep='\t')

col_mz = [c for c in og_df.columns if 'maize' in c.lower() or 'zea' in c.lower()][0]
col_sb = [c for c in og_df.columns if 'sorghum' in c.lower() or 'bicolor' in c.lower()][0]
col_at = [c for c in og_df.columns if 'ath' in c.lower() or 'arabidopsis' in c.lower()][0]

core_og_genes = {}
for _, row in og_df.iterrows():
    og_id = row['Orthogroup']
    mz_v = [g.strip() for g in str(row[col_mz]).split(',') if g.strip() in grn_graphs['maize']] if pd.notna(row[col_mz]) else []
    sb_v = [g.strip() for g in str(row[col_sb]).split(',') if g.strip() in grn_graphs['sorghum']] if pd.notna(row[col_sb]) else []
    at_v = [g.strip() for g in str(row[col_at]).split(',') if g.strip() in grn_graphs['ath']] if pd.notna(row[col_at]) else []
    if mz_v and sb_v and at_v: core_og_genes[og_id] = {'mz': mz_v, 'sb': sb_v, 'at': at_v}

gene2go = {'maize': {}, 'sorghum': {}, 'ath': {}}
def load_ensembl_go(sp, filename):
    df = pd.read_csv(os.path.join(ANNO_DIR, filename), sep='\t')
    g_col = [c for c in df.columns if 'Gene' in c][0]
    go_col = [c for c in df.columns if 'GO' in c and 'accession' in c.lower()][0]
    for _, row in df.dropna(subset=[go_col]).iterrows():
        g_id = str(row[g_col]).strip()
        if g_id not in gene2go[sp]: gene2go[sp][g_id] = set()
        gene2go[sp][g_id].add(str(row[go_col]).strip())

load_ensembl_go('maize', "")
load_ensembl_go('sorghum', "")
load_ensembl_go('ate', "")


def calc_SF(mz_g, sb_g, at_g):
    mz_go, sb_go, at_go = gene2go['maize'].get(mz_g, set()), gene2go['sorghum'].get(sb_g, set()), gene2go['ath'].get(at_g, set())
    if not mz_go and not sb_go and not at_go: return 0.0
    def jaccard(s1, s2): return len(s1.intersection(s2))/len(s1.union(s2)) if s1 and s2 else 0.0
    return (jaccard(mz_go, sb_go) + jaccard(mz_go, at_go) + jaccard(sb_go, at_go)) / 3.0

all_tuples = {}
for og_id, genes in core_og_genes.items():
    for combo in itertools.product(genes['mz'], genes['sb'], genes['at']):
        node_id = f"{combo[0]}|{combo[1]}|{combo[2]}"
        all_tuples[node_id] = {'og': og_id, 'mz': combo[0], 'sb': combo[1], 'at': combo[2], 'sf': calc_SF(combo[0], combo[1], combo[2])}


# -----------------  Seed  -----------------
seed = "Zm00001eb337450|Sobic.009G072100|AT2G25900" 
diverse_seeds = [seed]

used_mz = {seed.split('|')[0]}
used_sb = {seed.split('|')[1]}
used_at = {seed.split('|')[2]}

sorted_tuples = sorted(all_tuples.keys(), key=lambda k: all_tuples[k]['sf'], reverse=True)

for t in sorted_tuples:
    mz, sb, at = t.split('|')
    if mz not in used_mz and sb not in used_sb and at not in used_at:
        diverse_seeds.append(t)
        used_mz.add(mz)
        used_sb.add(sb)
        used_at.add(at)
        
    if len(diverse_seeds) == 10: 
        break

print("Diverse Seeds Selected:")
for i, s in enumerate(diverse_seeds):
    print(f"  Seed {i+1}: {s}")

seeds = diverse_seeds
# --------------------------------------------------------

for idx, seed in enumerate(seeds):
    current_module = [seed]
    
    def eval_mod(mod_nodes):
        sum_SF = sum(all_tuples[n]['sf'] for n in mod_nodes)
        sum_ST_raw = 0.0
        
        for u in mod_nodes:
            for v in mod_nodes:
                if u == v: continue
                u_mz, u_sb, u_at = all_tuples[u]['mz'], all_tuples[u]['sb'], all_tuples[u]['at']
                v_mz, v_sb, v_at = all_tuples[v]['mz'], all_tuples[v]['sb'], all_tuples[v]['at']
                w_mz = grn_graphs['maize'][u_mz][v_mz]['weight'] if grn_graphs['maize'].has_edge(u_mz, v_mz) else 0.0
                w_sb = grn_graphs['sorghum'][u_sb][v_sb]['weight'] if grn_graphs['sorghum'].has_edge(u_sb, v_sb) else 0.0
                w_at = grn_graphs['ath'][u_at][v_at]['weight'] if grn_graphs['ath'].has_edge(u_at, v_at) else 0.0
                
                nonzero_count = sum([1 for w in [w_mz, w_sb, w_at] if w > 0])
                if nonzero_count >= 2:
                    sum_ST_raw += ((w_mz + w_sb + w_at) / 3.0) * (nonzero_count ** 2)

        sum_ST_norm = sum_ST_raw / GLOBAL_MAX_ST if GLOBAL_MAX_ST > 0 else 0.0

        # return BETA * sum_ST_norm + GAMMA * sum_SF
        return (BETA * sum_ST_norm + GAMMA * sum_SF) * math.log2(len(mod_nodes) + 1)

    current_score = eval_mod(current_module)
    step = 1
    trajectory_data.append({'Module_ID': f"CumEps_{idx+1}", 'Step': step, 'Node_Added': seed, 'Local_Max_Score': current_score})
    
    while True:
        best_cand, best_score = None, current_score
        
        frontier = set()
        for node in current_module:
            mz_g, sb_g, at_g = all_tuples[node]['mz'], all_tuples[node]['sb'], all_tuples[node]['at']
            mz_neighbors = set(grn_graphs['maize'].successors(mz_g)).union(grn_graphs['maize'].predecessors(mz_g)) if mz_g in grn_graphs['maize'] else set()
            sb_neighbors = set(grn_graphs['sorghum'].successors(sb_g)).union(grn_graphs['sorghum'].predecessors(sb_g)) if sb_g in grn_graphs['sorghum'] else set()
            at_neighbors = set(grn_graphs['ath'].successors(at_g)).union(grn_graphs['ath'].predecessors(at_g)) if at_g in grn_graphs['ath'] else set()
            
            for cand, cand_data in all_tuples.items():
                if cand not in current_module:
                    neighbor_votes = 0
                    if cand_data['mz'] in mz_neighbors: neighbor_votes += 1
                    if cand_data['sb'] in sb_neighbors: neighbor_votes += 1
                    if cand_data['at'] in at_neighbors: neighbor_votes += 1
                    if neighbor_votes >= 2:
                        frontier.add(cand)
        # ---------------------------------------------
        
        for cand in frontier:
            score = eval_mod(current_module + [cand])
            if score > best_score:
                best_score = score
                best_cand = cand
                
        if best_cand:
            improvement_ratio = (best_score - current_score) / current_score if current_score > 0 else 1.0
            
            if improvement_ratio < EPSILON:
                print(f"  -> Module {idx+1} hit Diminishing Returns boundary at size {len(current_module)} (Improvement {improvement_ratio:.4f} < {EPSILON}). Stopping.")
                break 
                
            current_module.append(best_cand)
            current_score = best_score
            step += 1
            trajectory_data.append({'Module_ID': f"CumEps_{idx+1}", 'Step': step, 'Node_Added': best_cand, 'Local_Max_Score': current_score})
        else:
            print(f"  -> Module {idx+1} reached absolute topological end at size {len(current_module)}.")
            break
            
    discovered_modules.append({'Module_ID': f"CumEps_{idx+1}", 'Score': current_score, 'Size': len(current_module), 'Tuples': ", ".join(current_module)})


pd.DataFrame(discovered_modules).to_csv(os.path.join(GRN_DIR, "Ultimate_CumEps_Summary.tsv"), sep='\t', index=False)
pd.DataFrame(trajectory_data).to_csv(os.path.join(GRN_DIR, "Ultimate_CumEps_Trajectory.tsv"), sep='\t', index=False)
print("Alignment Complete!")
