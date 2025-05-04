rm(list = ls())
gc()
setwd("/home/kangziyi/comparsionDL/RealData/rainbow")
source("utils.r")
set.seed(42)
randomseed = sample(.Machine$integer.max,50)#.Machine$integer.max is 2147483647
percent = 0
#=================================
crap_geno = fread("myGD.txt",header = TRUE)
ID = crap_geno$taxa
crap_geno = as.matrix(crap_geno[,-1])
bool_na = crap_geno==5
crap_geno[bool_na] = 9

fwrite(as.data.table(cbind(ID,crap_geno)),col.names = FALSE,quote = FALSE,file = "genotype.txt",sep=" ")

system("AlphaImpute2 -genotypes genotype.txt -out imp -pop_only -maxthreads 16")

impgeno= fread("imp.genotypes",sep = " ")

idforimp = impgeno$V1

impgeno = impgeno[,-1]

impgeno = as.matrix(impgeno)

impgeno = apply(impgeno,2,function(x){
    
    x[which(!x%in%c(0,1,2))] = NA
    
    avrgeno = mean(x,na.rm = TRUE)

    x[is.na(x)] = avrgeno

    return(x)
})


rownames(impgeno) = idforimp

impgeno = impgeno[match(ID,rownames(impgeno)),]

genodf = as.data.table(cbind(ID,impgeno))

genodtput = makeped(z = impgeno)

mapdt= data.table(chr = rep(1,ncol(impgeno)),id = colnames(crap_geno),pos = as.integer(colnames(crap_geno)), gl = rep(0,ncol(impgeno)))

fwrite(genodtput,file="hib.ped",col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

fwrite(genodf,file="geno.csv",col.names = TRUE,row.names = FALSE, quote = FALSE,sep = ",")

fwrite(mapdt, file = "hib.map",col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")

system("plink --file hib --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib")

system("hiblup --make-xrm --threads 32 --bfile hib --add --out Gmat_b")

#phenotype
crap_pheno = fread("phenotype.csv",sep = ",")

colnames(crap_pheno)[1] = "ID1"

if(nrow(crap_pheno)!=nrow(impgeno)){
crap_pheno = crap_pheno[match(ID,crap_pheno$ID1),]
}
colnames(crap_pheno)[1] = "ID"

colnames(crap_pheno)[colnames(crap_pheno)=="Binary_survival"] = "phe1"
colnames(crap_pheno)[colnames(crap_pheno)=="Final_weight"] = "phe2"

if(percent>0){
if(percent == 0.1) set.seed(randomseed[50])
if(percent == 0.3) set.seed(randomseed[49])
if(percent == 0.5) set.seed(randomseed[48])
sample_num = round(length(crap_pheno$ID)*percent)
sample_ID = sample(crap_pheno$ID,sample_num )
crap_pheno[ID%in%sample_ID,phe2:=-999]
}

pheno1 = crap_pheno[,.(ID,phe1)]
colnames(pheno1) = c("ID","phenotype")

pheno2 = crap_pheno[,.(ID,phe2)]
colnames(pheno2) = c("ID","phenotype")

fwrite(pheno1,file = "phenotype_t.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

fwrite(pheno2,file = "phenotype_c.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

rep = 10
fold = 5

count = 0
for(r in 1:rep){

    count = count + 1

    set.seed(randomseed[count])

    fid = crap_pheno$ID
    n = length(fid)
    n = floor(n/fold)*fold
    assign_ids = sample(fid, n)
    fold_idx = sample(rep(c(1:fold), length.out = n))
    crap_pheno[ID%in%assign_ids,nfold:=fold_idx]

    for(k in 1:fold){
      colnm = paste("r",r,"k",k,sep = "_")

      crap_pheno[nfold==k,(colnm):="test"]

      val_fold = ifelse(k < fold, k + 1, 1)

      crap_pheno[nfold == val_fold,(colnm):="val"]

      crap_pheno[!nfold %in% c(k,val_fold),(colnm):="train"]
    }
}
crap_pheno[, nfold := NULL]

fwrite(crap_pheno,file = "phe_cv.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

#=================================


phe = fread("phe_cv.csv",sep = ",")
r = 10
k = 5
cl = clnm(r,k)
phe_test = phe[phe[[cl]] == "test",]
phe_val = phe[phe[[cl]] == "val",]
fwrite(phe_test,file = "phe_test.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
fwrite(phe_val,file = "phe_val.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
phe[phe2==-999,phe2 := NA]
phe[ID %in% c(phe_val$ID,phe_test$ID),phe1 := NA]
phe[ID %in% c(phe_val$ID,phe_test$ID),phe2 := NA]
fwrite(phe,file = "phe_hi.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")