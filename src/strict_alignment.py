import pandas as pd
import os
import networkx as nx
import itertools

BASE_DIR = os.path.expanduser("")
GRN_DIR = os.path.join(BASE_DIR, "")
ORTHO_FILE = os.path.join(BASE_DIR, "")
SPECIES = ['maize', 'sorghum', 'ath']

print(">>> Running Baseline: Strict Alignment")

grn_graphs = {sp: nx.DiGraph() for sp in SPECIES}
for sp in SPECIES:
    df = pd.read_csv(os.path.join(GRN_DIR, f"Ensemble_GRN_{sp}.csv"), sep='\t')
    for _, row in df.iterrows():
        grn_graphs[sp].add_edge(row['Gene1'], row['Gene2'], weight=row['Ensemble_Pred'])

og_df = pd.read_csv(ORTHO_FILE, sep='\t')
col_mz = [c for c in og_df.columns if 'maize' in c.lower() or 'zea' in c.lower()][0]
col_sb = [c for c in og_df.columns if 'sorghum' in c.lower() or 'bicolor' in c.lower()][0]
col_at = [c for c in og_df.columns if 'ath' in c.lower() or 'arabidopsis' in c.lower()][0]

core_og_genes = {}
for _, row in og_df.iterrows():
    og_id = row['Orthogroup']
    mz_g = [g.strip() for g in str(row[col_mz]).split(',')] if pd.notna(row[col_mz]) else []
    sb_g = [g.strip() for g in str(row[col_sb]).split(',')] if pd.notna(row[col_sb]) else []
    at_g = [g.strip() for g in str(row[col_at]).split(',')] if pd.notna(row[col_at]) else []
    
    mz_v = [g for g in mz_g if g in grn_graphs['maize']]
    sb_v = [g for g in sb_g if g in grn_graphs['sorghum']]
    at_v = [g for g in at_g if g in grn_graphs['ath']]
    
    if mz_v and sb_v and at_v:
        core_og_genes[og_id] = {'maize': mz_v, 'sorghum': sb_v, 'ath': at_v}


alignment_graph = nx.DiGraph()
tuple_nodes = []
for og_id, genes in core_og_genes.items():
    combos = list(itertools.product(genes['maize'], genes['sorghum'], genes['ath']))
    for combo in combos:
        node_id = f"{combo[0]}|{combo[1]}|{combo[2]}"
        alignment_graph.add_node(node_id, og=og_id, mz=combo[0], sb=combo[1], at=combo[2])
        tuple_nodes.append(node_id)

for u in tuple_nodes:
    for v in tuple_nodes:
        if u == v: continue
        u_mz, u_sb, u_at = alignment_graph.nodes[u]['mz'], alignment_graph.nodes[u]['sb'], alignment_graph.nodes[u]['at']
        v_mz, v_sb, v_at = alignment_graph.nodes[v]['mz'], alignment_graph.nodes[v]['sb'], alignment_graph.nodes[v]['at']
        
        if grn_graphs['maize'].has_edge(u_mz, v_mz) and \
           grn_graphs['sorghum'].has_edge(u_sb, v_sb) and \
           grn_graphs['ath'].has_edge(u_at, v_at):
            
            w_mz = grn_graphs['maize'][u_mz][v_mz]['weight']
            w_sb = grn_graphs['sorghum'][u_sb][v_sb]['weight']
            w_at = grn_graphs['ath'][u_at][v_at]['weight']
            alignment_graph.add_edge(u, v, weight=w_mz*w_sb*w_at, mz_w=w_mz, sb_w=w_sb, at_w=w_at)

modules = [alignment_graph.subgraph(c) for c in nx.weakly_connected_components(alignment_graph) if len(c) >= 2]
summary_data, edges_data = [], []

for idx, mod in enumerate(sorted(modules, key=lambda g: sum(d['weight'] for u,v,d in g.edges(data=True)), reverse=True)):
    mod_id = f"StrictMod_{idx+1}"
    summary_data.append({
        'Module_ID': mod_id, 'Num_Tuples': mod.number_of_nodes(), 'Num_Edges': mod.number_of_edges(),
        'Score': sum(d['weight'] for u,v,d in mod.edges(data=True)), 'Tuples': ", ".join(list(mod.nodes()))
    })
    for u, v, d in mod.edges(data=True):
        edges_data.append({
            'Module_ID': mod_id, 'Source': u, 'Target': v, 'Joint_Score': d['weight'],
            'Mz_w': d['mz_w'], 'Sb_w': d['sb_w'], 'At_w': d['at_w']
        })

pd.DataFrame(summary_data).to_csv(os.path.join(GRN_DIR, "Strict_Aligned_Summary.tsv"), sep='\t', index=False)
pd.DataFrame(edges_data).to_csv(os.path.join(GRN_DIR, "Strict_Aligned_Edges.tsv"), sep='\t', index=False)
print(f"Baseline Complete! Saved {len(modules)} strict modules.")
