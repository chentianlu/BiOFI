BiOFI
================
 an R package for metabolome and microbiome feature identification
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

<!-- badges: start -->

<!-- badges: end -->

Interactions between the metabolome and microbiome are pervasive within the human body, where microbial compositions and functions are intricately mirrored by metabolome profiles. Simultaneously, the metabolic processes and expressions are under the influence and regulation of the microbiome. Various methods and tools have been developed by us and others for the identification of association pairs. Nonetheless, existing tools typically rely on correlation analysis alone, often overlooking crucial factors such as differences in study groups, the complex topology of association networks among features, the biological context of these features, and their ranking based on abundance. These tools tend to separately examine functional and compositional features, and there is a notable absence of automated solutions to bridge the gap between microbes, functions, and metabolites. Furthermore, the large number of significant features and correlations presents formidable challenges for subsequent experimental inquiries. An imperative requirement exists for an algorithm capable of synthesizing results from diverse perspectives and unifying them into a comprehensive variable importance score. Here, we designed an R package, BiOFI, to calculate a comprehensive importance score (important feature score, IFS) for each feature (a microbe, a metabolite, or a function), to link microbe and metabolite by their shared function (based on embedded library derived from KEGG), and to rank correlated pairs within specific functions by the important pair score (IPS). Users can choose a small selection of the most crucial features and pairs for further biological research.


## Installation

You can install the development version of BiOFI from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("chentianlu/BiOFI",force = TRUE)
library("BiOFI")
```

## Example

# `MicPathmatch`

A funciton that is used to get pathways correspondent with the microbes
in your own 16S or Metagenome data.The result includes two parts that
one contains related detailed strain information and pathways of
microbes, and the other only contains names of microbes you want to
predict and related pathways.

An example of MicPathmatch:

``` r
data('MicPathmatch.eg')
MicPathmatch_res <- MicPathmatch(microbes = MicPathmatch.eg)
```

# `Meta2pathway`

A function that can transfer metabolite abundance into related pathway
data.you can get data containing metabolite abundance and pathway that
metabolites belong to.  
Note: Please use C number in KEGG database to replace the metabolite
names.For example, C00002 means ATP or Adenosine 5’-triphosphate in
KEGG. Try your best to match the C number of metabolites for predicting
pathway that metabolites belong to more accurately. if the metabolites
do not have C numbers, keep the original names.

An example of Meta2pathway:

``` r
data('Meta2pathway.eg')
Meta2pathway_res <- Meta2pathway(metadf = Meta2pathway.eg)
```

# `IFS`

A function that calculates and selects the key nodes,the formula of IFS
is: IFS= nDFS×DFS + nDS×DS+nES×ES+nAS×AS; nDFS,nDS,nES,nAS indict the
weight value of DFS,DS,ES and AS, and nDFS+nDS+nES+nAS=1; DFS means
Difference Score; DS means Degree Score; ES means Edge Score; AS means
Abundance Score.

An example of IFS:

``` r
data('micro.eg')
data('metabo.eg')
data('groupInfo.eg')
data('confounder.eg')
IFScore=IFS(microApath = micro.eg, metaApath = metabo.eg,
conf=confounder.eg,groupInfo=groupInfo.eg)
```

# `MMfunc`

A function that selects microbe and metabolites on the same pathway
based on the data of Nodescore.The result will be saved as HTML file.

An example of MMfunc:

``` r
data('micro.eg')
data('metabo.eg')
data('groupInfo.eg')
data('confounder.eg')
IFScore=IFS(microApath = micro.eg, metaApath = metabo.eg,
	groupInfo=groupInfo.eg)
MMfunc_res <- MMfunc(IFS = IFScore)
```

# `TarNet`

A function used to plot a network that consists of microbes and
metabolites in same pathway.This function offers a network of microbes
and metabolites in same pathways.

An example of TarNet:

``` r
data('micro.eg')
data('metabo.eg')
data('groupInfo.eg')
data('confounder.eg')
IFScore=IFS(microApath = micro.eg, metaApath = metabo.eg,
	groupInfo=groupInfo.eg)
MMfunc_res <- MMfunc(IFS = IFScore)
TarNet_res <- TarNet(IFS = IFScore, mm2path = MMfunc_res)
```

# `IPS`

A function used to calculate the importance pair score (IPS) of a
network based on the return of the TarNet function.

An example of IPS:

``` r
data('micro.eg')
data('metabo.eg')
data('groupInfo.eg')
data('confounder.eg')
IFScore=IFS(microApath = micro.eg, metaApath = metabo.eg,
	groupInfo=groupInfo.eg)
MMfunc_res <- MMfunc(IFS = IFScore)
TarNet_res <- TarNet(IFS = IFScore, mm2path = MMfunc_res)
IPScore <- IPS(TarnetPairs = TarNet_res,IFS = IFScore)
```
