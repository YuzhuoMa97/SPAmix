# SPAmix
A scalable, accurate, and universal analysis framework using individual-level allele frequency for large-scale genetic association studies in an admixed population

Please do not hesitate to contact me (yuzhuoma@stu.pku.edu.cn) if you are interested in SPAmix or other retrospective saddlepoint approximation methods desigend for GWAS.

# Introduction of Retrospective saddlepoint approximation (Retrospective SPA) approach in GWAS

Retrospective saddlepoint approximation (Retrospective SPA) is a method applied in GWAS. For a score statistic (S=G^TR), retrospective SPA strategy considers genotypes as random variables and approximate the distribution of score statistics conditional on phenotype and covariates.

Retrospective saddlepoint approximation (Retrospective SPA) methods were first proposed in the master's thesis ([马雨茁.经验鞍点近似方法及其在全基因组关联分析中的应用研究.2022.山东大学,MA thesis.doi:10.27272/d.cnki.gshdu.2022.002946.](https://kns.cnki.net/kcms2/article/abstract?v=jkwd3qsBIEKwkKkgMuimTLSEojAEBaWSJzCAd3uOCepX09aaYi1Vhn87HddxnsydAW9MGQHzgdF9Nw93IZ_DZCdJbGAX3C13DfGxpW58VBV273z1eVlg75Je1akPxIDc5iiSpz46iutS1tt9m3MJRg==&uniplatform=NZKPT&language=CHS) 
**DOI：	10.27272/d.cnki.gshdu.2022.002946**.)

Based on the idea in the above master's thesis, we have applied retrospective saddlepoint approximation to several methods including SPAmix (since 2020), SPAmix+ (based on SPAmix since 2024), SPAGxECCT (based on SPAmix since 2021), and SPAGxEmix+ (based on SPAmix since 2024). 



**If you utilized the retrospective saddlepoint approximation method in your proposed methods or tools, please acknowledge and respect the original ideas presented in the two pioneering works (SPAGxECCT and SPAmix). Additionally, kindly cite the original papers (SPAGxECCT and SPAmix) or the master's thesis ([马雨茁.经验鞍点近似方法及其在全基因组关联分析中的应用研究.2022.山东大学,MA thesis.doi:10.27272/d.cnki.gshdu.2022.002946.](https://kns.cnki.net/kcms2/article/abstract?v=jkwd3qsBIEKwkKkgMuimTLSEojAEBaWSJzCAd3uOCepX09aaYi1Vhn87HddxnsydAW9MGQHzgdF9Nw93IZ_DZCdJbGAX3C13DfGxpW58VBV273z1eVlg75Je1akPxIDc5iiSpz46iutS1tt9m3MJRg==&uniplatform=NZKPT&language=CHS) in accordance with academic standards.**

MLA format citation: [1]马雨茁.经验鞍点近似方法及其在全基因组关联分析中的应用研究.2022.山东大学,MA thesis.doi:10.27272/d.cnki.gshdu.2022.002946.

Given that the two pivotal papers  

(1) **SPAmix** (A scalable, accurate, and universal analysis framework using individual-level allele frequency for large-scale genetic association studies in an admixed population)

and

(2) **SPAGxECCT** (A scalable and accurate framework for large-scale genome-wide gene-environment interaction analysis and its application to time-to-event and ordinal categorical traits)

have not yet been published  (**despite more than three years since retrospective saddlepoint approximation method was first proposed in 2021 (in SPAmix) and published in 2022 [马雨茁.经验鞍点近似方法及其在全基因组关联分析中的应用研究.2022.山东大学,MA thesis.doi:10.27272/d.cnki.gshdu.2022.002946.](https://kns.cnki.net/kcms2/article/abstract?v=jkwd3qsBIEKwkKkgMuimTLSEojAEBaWSJzCAd3uOCepX09aaYi1Vhn87HddxnsydAW9MGQHzgdF9Nw93IZ_DZCdJbGAX3C13DfGxpW58VBV273z1eVlg75Je1akPxIDc5iiSpz46iutS1tt9m3MJRg==&uniplatform=NZKPT&language=CHS)), please consult the relevant authors and supervisors.**

Suggestions or comments on retrospective saddlepoint approximation methods are also welcome.

# SPAmix can control for population admixture


# SPAmix+ is an extension of SPAmix
In the SPAmix paper (**A scalable, accurate, and universal analysis framework using individual-level allele frequency for large-scale genetic association studies in an admixed population**), we proposed SPAmix framework that is applicable to admixed populations and can incorporate local ancestry information into analyses. However, SPAmix is only designed for unrelated individuals and cannot account for sample relatedness.

**We are preparing the SPAmix+ paper (A scalable and accurate analysis framework accounting for sample relatedness and population structure for large-scale genetic association studies in an admixed population).** In this paper, we extend SPAmix to SPAmix+, a universal analysis framework that accounts for both population structure and sample relatedness. SPAmix+ allows for the incorporation of admixed individuals into analyses. In addition, SPAmix can be extended to SPAmix+(local), which has the ability to incorporates local ancestry.

# The SPAmix paper and the SPAmix+ paper are two distinct papers

**Please note that the SPAmix paper and the SPAmix+ paper are two distinct papers. Unfortunately, although the SPAmix algorithm was completed in 2021, its submission was delayed for some reasons. Consequently, the SPAmix+ algorithm will be presented in a separate paper rather than within the SPAmix paper.** 


