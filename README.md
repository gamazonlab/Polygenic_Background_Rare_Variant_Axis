# Polygenic_Background_Rare_Variant_Axis [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/gamazonlab/Polygenic_Background_Rare_Variant_Axis/blob/master/LICENSE) 

## Reference

Contextualizing genetic risk score for disease screening and rare variant discovery. Dan Zhou, Dongmei Yu, Jeremiah M. Scharf, Carol A. Mathews, Lauren McGrath, Edwin Cook, S Hong Lee, Lea K. Davis<sup>§</sup>, Eric R. Gamazon<sup>§</sup>

<sup>§</sup>Co-corresponding authors:  Eric R. Gamazon and Lea K. Davis

## Authors of methodology and code

For questions, please feel free to contact

Dan Zhou zdangm@gmail.com;
Eric R. Gamazon ericgamazon@gmail.com  

## Description

We provide a set of scripts for simulating, estimating, and visualizing the relationship between common, small effect variant polygenic background (cSEV-PB) and large effect variant (LEV). Details can be found in our manuscript ‘Contextualizing genetic risk score for disease screening and rare variant discovery’.

Software dependency: R (>3.4.0) and ggplot2 will be needed for the simulations and the reproduction of the figures.

### schematic_figure.r
The script is a schematic figure generator for the central hypothesis underlying the relationship between cSEV-PB and LEVs. By varying the heritability attributable to the cSEV-PB, the heritability due to the LEV, the number of LEVs, allele frequency of LEV, sample size, and disease prevalence, users can customize the simulation settings to make predictions or to conduct statistical analysis of empirical data. A figure with panel (a) showing the distribution of the liability, a latent trait resulting from the cumulative effects of cSEV-PB, LEVs, and the residual component, and panel (b) showing the relationship between cSEV-PB and LEVs (and therefore, as a specific example, the extent to which the cSEV-PB can be used for pathogenic variant discovery) will be generated.

### simulation_utility_power_and_OR_comparison_LEV_cSEV.r
This script simulates the cSEV-PB and LEV burden and estimates the utility, power, and OR of LEV and of the cSEV-PB under three classes of genetic architectures (polygenic model, model consistent with negative selection, and LD-adjusted kinship) and two disease models (liability threshold model and logit model). The computational time largely depends on the sample size, number of parameters varied, and number of values assumed for each parameter.

### simulation_cSEV-OR_comparison_between_LEV_carriers_noncarriers.r
The script simulates the cSEV-PB and LEV burden and estimates the OR of the cSEV-PB in LEV carriers and non-carriers.

### negativeSelection_GeneticArchitecture.r
This script generates Q-Q plots that facilitate comparison of the significance of the cSEV-PB and LEV correlation among three classes of genetic architectures: the polygenic model, model consistent with negative selection, and LD-adjusted kinship.

### subtype.r
By assuming the presence of (two) disease subtypes, using a polygenic subtype heterogeneity parameter lambda, under a linear model of the underlying liability , the script generates a series of figures showing the robustness of the cSEV-PB and LEV correlation. 

### Figures.r
This script generates the figures in the manuscript. 

