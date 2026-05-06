# Relaxed Topological Alignment of Cross-Species GRNs

This repository contains the custom Python and R source code for the relaxed topological alignment algorithm. This framework is designed to extract functionally coherent, conserved regulatory modules from noisy, multi-species Gene Regulatory Networks (GRNs). 

Traditional strict alignment algorithms frequently fail under the many-to-many orthology mappings inherent in evolutionary datasets. Our algorithm resolves this by formulating subgraph discovery as a multi-objective optimization problem—balancing sequence homology, network topology, and functional coherence—navigated via a greedy seed-and-extend heuristic with a dynamic epsilon-stopping condition.

## 📁 Repository Structure

* `src/`: Core algorithm implementation.
    * `strict_alignment.py`: Baseline topological intersection algorithm.
    * `relaxed_heuristic.py`: Main algorithm featuring the epsilon-stopping expansion.
* `analysis/`: Scripts to reproduce the statistical analysis and figures presented in the manuscript.
    * `trajectory_visualization.py`: Generates the alignment trajectory and early-stopping curves.
    * `upset_convergence.py`: UpSet plot generation for heuristic convergence testing.
* `results/`: Output directories for generated subgraphs, statistical summaries, and high-resolution PDF figures.

## ⚙️ Prerequisites and Installation

The pipeline relies on several standard Python and R libraries for graph computation and statistical analysis. 

**Python Dependencies:**
* `pandas`, `numpy`, `scipy`, `statsmodels`
* `networkx` (for core topological operations)
* `matplotlib`, `seaborn`, `upsetplot`, `adjustText` (for publication-quality visualization)

## 🚀 Usage 

Detailed execution instructions and parameter definitions for the command-line interface are provided in the `src/` directory. A standard run of the relaxed alignment pipeline follows this sequence:

1.  **Initialize the Graph:** Load species-specific ensemble GRNs and the orthology mapping table.
2.  **Seed Selection:** Identify highly conserved strict core tuples.
3.  **Heuristic Expansion:** Execute the relaxed alignment to dynamically recruit functionally coherent downstream targets.
4.  **Validation:** Run downstream enrichment scripts to analyze the biological relevance of the extracted super-module.


## 🙏 Acknowledgements

This algorithm and the corresponding stress-response case studies were developed with support from the NSF project: *Unraveling the origin of vegetative desiccation tolerance in vascular plants*.
We would like to express our sincere gratitude to Dr. Luis Herrera-Estrella and Dr. Liqing Zhang  for providing valuable feedback and suggestions.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
