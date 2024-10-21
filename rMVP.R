library(rMVP)
# Full-featured function (Recommended)
MVP.Data(fileVCF="../02Genotype/278Sweet_corn_V5.miss05.maf001.recode.vcf",
         #filePhe="Phenotype.txt",
         fileKin=F,
         filePC=F,
         out="mvp.vcf"
)

# Only convert genotypes
#MVP.Data.VCF2MVP("mvp.vcf", out='mvp') 

# calculate from mvp_geno_file
MVP.Data.PC(TRUE, mvp_prefix='hyb_mvp.plink', pcs.keep=5)
# calculate from mvp_geno_file
MVP.Data.Kin(TRUE, mvp_prefix='hyb_mvp.plink', out='hyb_mvp.plink')


Kinship <- attach.big.matrix("purned_MVP/hyb_mvp.plink.kin.desc")
genotype <- attach.big.matrix("purned_MVP/hyb_mvp.plink.geno.desc")
phenotype <- read.table("transformed_BGEM_alpha_hyb_LN.txt",head=TRUE,sep="\t")
map <- read.table("purned_MVP/hyb_mvp.plink.geno.map" , head = TRUE)

for(i in 2:ncol(phenotype)){
  imMVP <- MVP(
    phe=phenotype[, c(1, i)],
    geno=genotype,
    map=map,
    K=Kinship,
    #CV.GLM=Covariates,
    #CV.MLM=Covariates,
    #CV.FarmCPU=Covariates,
    nPC.GLM=3,
    nPC.MLM=3,
    nPC.FarmCPU=3,
    priority="speed",
    #ncpus=10,
    vc.method="BRENT",
    maxLoop=10,
    method.bin="static",
    #permutation.threshold=TRUE,
    #permutation.rep=100,
    threshold=0.05, ##0.05/marker size, a cutoff line on manhattan plot
    method=c("GLM", "MLM", "FarmCPU"),
    p.threshold=0.00001
    
  )
  gc()
}