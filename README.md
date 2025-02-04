# Leveraging Comparative Phylogenetics for Evolutionary Medicine  

## Setup  
Before running any scripts, execute **`checkLibrary.R`** to ensure all required R packages are installed.  

---

## R Scripts for Visualizations
Plots from these scripts are store in plots directory

### `arcEvoModelClass.R`  
- Generates **Figure 1** by fitting groups in the phylogeny to various evolutionary models.  
- Highlights the best-fit model for each group on the tree.  

### `simVizCombo.R`  
- Creates bar graph error rate visualizations for **Figure 3**.  

### `contMapPagel.R`  
- Produces **Figure 2**, an **Ancestral State Reconstructed Phylogeny** with **Threshold Model Pie Charts**.  

---

## Experiment Scripts  

Each experiment is organized into three directories:  
- **50**, **150**, **350** → Indicates the number of species included in each simulation study.  
- Each directory contains similar files, adjusted for sample size.  

### Evolutionary Model Experiments  
Scripts for testing different evolutionary models:  
- **`expOU.R`**  
- **`expBrown.R`**  
- **`expPagel.R`**  
- **`expWhite.R`**  

### **Experiment Workflow (Repeated for 1000 Iterations):**  
1. **Simulate necropsy data** based on real-world distributions.  
2. **Generate trait data** following the selected evolutionary model.  
3. **Test relationships** using three statistical methods:  
   - `Compar.GEE` (binomial regression)  
   - `PGLS.Sey`  
   - `PGLS.Sey.Optim`  

---


## Addtional Files

###  min1LH.csv
Necropsy data for minimum of 1 species

###  min20-2022.05.16.csv
Necropsy data for minimum of 20 species

### min1.nwk
Newick tree for min1 species

### min20Fixed516.nwk
Newick tree for min20 species


