## Testing distance measures on Tabula muris data

S. Anders, 2018-06


*To properly display the math formulae, creat an HTML file with:* `pandoc README.md -o README.html -s --mathjax
`

### Project overview

Here, I try out some ideas on distance emasures for single-cell RNA-Seq data on the Dropseq data from the Tabula muris project.


### Files

- `README.md` -- The present file

- `Tm_droplet_make_h5ad.ipynd` -- PythonJupyter notebook, showing how to create an `.h5ad` file from the Tabula muris Dropseq data. An `.h5ad` file is an HDF5 file that holds an `AnnData` object as defined by the *Scanpy* project.

- `Tabula_muris_10x.h5ad` -- The Scanpy Anndata file with the Tabula muris Dropseq data. This file can be generated with the above script, or downloaded [here](http://www.zmbh.uni-heidelberg.de/Anders/div/Tabula_muris_10x.h5ad)

- `similarity_by_blocks.ipynb` -- A notebook to create a cell-to-cell similarity-score matrix and save it in a (huge) HDF5 file

- `endo_neighbours.ipynb` -- A notebook trying out ways to use this similarity matrix to call cell types.

- and various outdated files that I should delete

### The similarity score

For now, just a brief summary

- I write $K_{ij}$ for the UMI count for gene $i$ in cell $j$.

- First, I transform all count values to their square roots. (This is the variance-stabilizing transformation for Poisson-distributed data; sometimes called the Anscombe transformation.)

- Then I normalize each cell to unit L2 norm, i.e., with $K_{ij}$ the count for gene $i$ in cell $j$, I use $y_{ij} = \frac{\sqrt{K_{ij}} }{ \sqrt{\sum_i K_{ij}} }$ as expression measure. Coincidentally this normalization can also be understood as L1 normalization of the counts (i.e., simply dividing by the total UMI count), *followed* by the square root transformation.

- As a similarity score, I simply use the scalar product of the cells' expression vectors: $s_{jj'} = \sum_i y_{ij} y_{ij'}$

- As the expression vectors have unit L2 norm, the similarity score varies between 0 (orthogonal) and 1 (colinear).

- The similarity score $s_{jj'}$ has direct relation to two distance metrics, the Euclidean distance $d_{jj'}^\text{E}$ and the cosine distance (or "angular distance") $d_{jj'}^\text{A}$:

- For the Euclidean distance, we have: $d_{jj'}^\text{E} = \sqrt{1 - s_{jj'}/2}$. (To see, just expand the squares in $(d^\text{E}_{jj'})^2 = \sum_i(y_{ij}-y_{ij'})^2$ and remember that $\sum_i y_{ij}^2=1$.) 

- The cosine of the angle between the expression vectors of cells $j$ and $j'$ is equal to their scalar product (because they are unit vectors). We therefor define the angular distance as $d_{jj'}^\text{A}=\frac{2}{\pi} \arccos s_{jj'}$. Dividing the arc cosine by $\pi/2$ maps the possible angles (colinear to orthogonal) to the unit interval.

- Note that $d^\text{E}_{jj'} \approx d^\text{A}_{jj'}$. The two values are quite similar up to around 0.7.

- My hope it that this very simple distance, which has no tuning parameter, performs as well as the one used in Seurat, which requires selecting highly variable genes and deciding on how many principle components to use.

- The core idea is just the Anscombe transformation. This ensures that each gene contributes the same amount of Poisson noise. The problem with Seurat's approach is, in my opinion, that after log transformation, the lowly expressed genes contribute too much noise. The main effect of selecting highly variably genes is 
merely that one gets rid of the weakly expressed genes.


## Nearest neighbor classification

[...]