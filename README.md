# Polygenic_Background_Rare_Variant_Axis [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/gamazonlab/Polygenic_Background_Rare_Variant_Axis/blob/master/LICENSE) 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4767933.svg)](https://doi.org/10.5281/zenodo.4767933)


## Reference

Contextualizing genetic risk score for disease screening and rare variant discovery. Dan Zhou, Dongmei Yu, Jeremiah M. Scharf, Carol A. Mathews, Lauren McGrath, Edwin Cook, S Hong Lee, Lea K. Davis<sup>§</sup>, Eric R. Gamazon<sup>§</sup>. Nature Communications.

<sup>§</sup>Co-corresponding authors:  Eric R. Gamazon and Lea K. Davis

## Authors of methodology and code

For questions, please feel free to contact

Dan Zhou zdangm@gmail.com;
Eric R. Gamazon ericgamazon@gmail.com  

## Description

We provide a set of scripts for simulating, estimating, and visualizing the relationship between common, small effect variant Polygenic Burden (PB) and Large-Effect Variant (LEV). We present a mathematical and summary statistics based framework that utilizes the relationship between PB and LEV among cases for some methodologically-and clinically-relevant applications. Details can be found in our manuscript ‘Contextualizing genetic risk score for disease screening and rare variant discovery’.

## Software dependency

R packages ggplot2 and optparse will be needed for the simulations and the reproduction of the figures. Both plink1.9 and plink2 are needed for PB-LEV simulation. Plink1.9 is needed for genotype simulation (plink2 does not support it anymore) and plink2 is needed for generating polygenic score on standardized genotype data (plink1.9 don’t have the option).


## PB-LEV-SCAN_1.0.r (the main script)

### Arguments

### #path

--plink_path path to plink1.9 [required, default = plink]

--plink2_path path to plink2 [required, default = plink2]

--maf_distribution_path Path to the maf distribution file (estimated from empirical data) [required]

--ldsc_path Path to the ldsc file (estimated from empirical data) [required]

--tmp_folder path to a folder for saving intermediate files [required, no default value]

--out_folder path for output folder [required]

--out_prefix prefix for result files. Result files will be named as ${out_folder}/${out_prefix}_*.* [required]
### #model

--genetic_architecture the assumed genetic architecture for common variants. Options include (1) “polygenic” which denotes a highly polygenic (almost infinitesimal) architecture, (2) “NegativeSelection” denotes an architecture consistent with negative selection, and (3) “LDAK” denotes LD-adjusted kinship. See Methods for details. [default = polygenic] 

--disease_model the assumed disease model. Options include: (1) “LTM” which is short for liability-threshold model, and (2) “logit” denotes the logit risk model. See Methods for details. [default = LTM]

### #parameter

--n_simu number of simulations [default = 500]

--sample_size sample size including cases and controls [default = 10000]

--h2_PB heritability of Polygenic Burden [default = 0.3]

--freq_LEV minor allele frequency of Large-Effect Variant(s). [default = 0.01]

--h2_LEV heritability of Large-Effect Variant(s). Please only provide the estimated heritability (h2_LEV) or the odds ratio (OR_LEV). If both of them are provided, OR_LEV will be ignored. [default = 0.01]

--OR_LEV odds ratio of Large-Effect Variant(s). Please only provide the estimated heritability (h2_LEV) or the odds ratio (OR_LEV). If both of them are provided, OR_LEV will be ignored. [optional, no default value]

--prevalence the estimated prevalence of disease/traits in the general population [default = 0.01]

--pi0 the proportion of non-causal variants [default = 0]

--seed seed for genotype simulation, sampling, and other stochastic procedures. [default = 0]
### #other options

--generate_a_figure generate a figure comparing the number of LEV carriers among patients with different polygenic burdens. [default = TRUE]

--clean_tmp remove the tmp folder and all temporary files. If FALSE, the intermediate datasets (including simulated genotype files ${tmp_folder}/tmp_geno.bed/.bim/.fam) will not be removed after simulations. [default = TRUE]


### Simulation results
${out_path}/${out_prefix}_simu.txt includes raw simulation results. 

n_simu: simulation id

wilcox_p: P-value of two-sample Wilcoxon test testing, among cases, whether the polygenic risk for LEV carriers is lower than non-carriers (one-side test).

p1-p10: the estimated number of LEV carriers per 1000 cases from group 1 to 10. The 10 groups are equally sized based on the polygenic burden (PB). p1 and p10 denote the group with the lowest and the highest PB, respectively.

${out_path}/${out_prefix}.txt includes summarized simulation results. It reports all the parameters used for simulations. The utility/power (when pi!=0) of observing negative correlation of PB-LEV is calculated as the proportion of times seeing significant (Wilcoxon test) difference polygenic risk between LEV carriers and non-carriers among cases. The file also provides mean and sd for the number of LEV carriers per 1000 cases for bin1 to bin10.

${out_path}/${out_prefix}_LEV_carriers_per_1000_cases_comparison_figure.pdf is a figure visualizing the distribution of polygenic risk and the number of LEV carriers per 1000 cases. The mean ± 1sd of the number of LEV carriers per 1000 cases is displayed in blue circles and bars.

### Notes
For common variant (MAF > 0.01) simulations, the distribution of MAF among all variants was estimated from a large-scale biobank. Users can modify the distribution from their specific studies. Please refer to the plink webpage for the format of the distribution file. 
The ld-score is supposed to be estimated from the same biobank. The script will automatically match the pattern between MAF and ld-score based on the provided empirical data.
The default parameters are all assigned. Empirical parameters from UK Biobank traits are summarized from three resources including (1) Ben Neale Lab [http://www.nealelab.is/uk-biobank], (2) [Wang, Quanli, et al. bioRxiv 2020](https://www.biorxiv.org/content/10.1101/2020.12.13.422582v1.abstract), and (3) [Lu, Tianyuan, et al. Genetics in Medicine 2020](https://www.nature.com/articles/s41436-020-01007-7). The summarized data are now available in the supplementary tables of our paper.

### Example
Rscript ${dir_script}/PB-LEV-SCAN_1.0.r \
--tmp_folder ${dir_tmp_folder} \
--OR_LEV 2 \
--freq_LEV 0.01 \
--h2_PB 0.3 \
--n_simu 50  \
--sample_size 1000 \
--genetic_architecture LDAK \
--maf_distribution_path ${dir_maf}/maf_distribution.txt \
--ldsc_path ${dir_ldsc}/ldsc.txt.gz  \
--out_folder ${dir_output} \
--out_prefix random_run


## Other scripts

### schematic_figure.r
The script is a schematic figure generator for the central hypothesis underlying the relationship between PB and LEVs. By varying the heritability attributable to the PB, the heritability due to the LEV, the number of LEVs, allele frequency of LEV, sample size, and disease prevalence, users can customize the simulation settings to make predictions or to conduct statistical analysis of empirical data. A figure with panel (a) showing the distribution of the liability, a latent trait resulting from the cumulative effects of PB, LEVs, and the residual component, and panel (b) showing the relationship between PB and LEVs (and therefore, as a specific example, the extent to which the PB can be used for pathogenic variant discovery) will be generated.

### simulation_utility_power_and_OR_comparison_LEV_cSEV.r
This script simulates the PB and LEV burden and estimates the utility, power, and OR of LEV and the PB under three classes of genetic architectures (polygenic model, model consistent with negative selection, and LD-adjusted kinship) and two disease models (liability-threshold model and logit model). The computational time largely depends on the sample size, the number of parameters varied, and the number of values assumed for each parameter.

### negativeSelection_GeneticArchitecture.r
This script generates Q-Q plots that facilitate comparison of the significance of the PB and LEV correlation among three classes of genetic architectures: the polygenic model, model consistent with negative selection, and LD-adjusted kinship.

### subtype.r
By assuming the presence of (two) disease subtypes, using a polygenic subtype heterogeneity parameter lambda, under a linear model of the underlying liability, the script generates a series of figures showing the robustness of the PB and LEV correlation. 


