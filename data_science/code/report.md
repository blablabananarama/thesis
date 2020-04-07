

# Differential expression analysis 

## Data

The data used in this notebook comes from the supplementary material of this paper: blabla
It is resampled from counts that represent means and standard deviations.

## Why EdgeR?
- Uses count data 
- Normalization suited for count data, no between probe normalization as for 
 
## Results

![N up and downregulated](../results/figures/go_analysis/differential-expression_basic-proteomics-dataset_n-up-downregulated05.png "Up and down")

*Number of up- and downregulated genes decided by fishers combined method with an alpha cutoff of 0.05*

![N up and downregulated](../results/figures/go_analysis/differential-expression_basic-proteomics-dataset_n-up-downregulated01.png "Up and down")

*Number of up- and downregulated genes decided by fishers combined method with an alpha cutoff of 0.01*

![N up and downregulated](../results/figures/go_analysis/differential-expression_basic-proteomics-dataset_n-up-downregulated001.png "Up and down")

*Number of up- and downregulated genes decided by fishers combined method with an alpha cutoff of 0.001*


# Go term analysis with piano

## Why the given parameters?
- fishers combined test: includes the 

## Results

![downregulated](../results/figures/go_analysis/go-p-values-fisher-combined_basic-proteomics-dataset_upregulated.png "down") ![downregulated](../results/figures/go_analysis/go-p-values-fisher-combined_basic-proteomics-dataset_downregulated.png "up") ![upregulated](../results/figures/go_analysis/go-p-values-fisher-combined_basic-proteomics-dataset_up-and-downregulated.png "both")


# Go term analysis with Goatools
