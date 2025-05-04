rm(list = ls())
gc()

dir = "/home/kangziyi/poster/BaseMapHaplo/"
mapFile = "mergedMap.csv"
haploFile = "mergedHaplo.csv"
reldegree = 0.3 #0.3 0.5 0.7

set.seed(42)
randomseed = sample(.Machine$integer.max,50)#.Machine$integer.max is 2147483647

source("utils.r")

num <- data.frame(chr = c("NC_088853.1",
                          "NC_088854.1",
                          "NC_088855.1",
                          "NC_088856.1",
                          "NC_088857.1",
                          "NC_088858.1",
                          "NC_088859.1",
                          "NC_088860.1",
                          "NC_088861.1",
                          "NC_088862.1"),
                  len = c(76070991,
                          61469542,
                          61039741,
                          57946171,
                          57274926,
                          56905015,
                          53672946,
                          51133819,
                          50364239,
                          37310742)
)

genMap = fread(paste(dir,mapFile,sep=""),sep = ",")

Haplo = fread(paste(dir,haploFile,sep=""),sep = ",",header = TRUE)

colnames(Haplo) = paste("Site_",colnames(Haplo),sep = "")

genMap[,site:=colnames(Haplo)]

map_list = list()
haplo_list = list()
qtl_pos = list()
snp_pos = list()

for(i in 1:nrow(num)){
    sites = genMap[chr == num$chr[i],site]

    Qtl_sites = genMap[chr == num$chr[i]&QTL==TRUE,site]

    Snp_sites = genMap[chr == num$chr[i]&QTL!=TRUE,site]

    map_list[[i]] = genMap[chr == num$chr[i],pos]
    map_list[[i]] = map_list[[i]]/num$len[i]

    tempIndex = match(sites,colnames(Haplo))

    haplo_list[[i]] = Haplo[,.SD,.SDcols = colnames(Haplo)[tempIndex] ]

    qtl_pos[[i]] = which(sites %in% Qtl_sites)

    snp_pos[[i]] = which(sites %in% Snp_sites)

}

founderPop <- newMapPop(genMap=map_list, haplotypes=haplo_list)
SP <- SimParam$new(founderPop)

SP$invalidQtl <- snp_pos
SP$invalidSnp <- qtl_pos

SP$addTraitADE(nQtlPerChr = sapply(qtl_pos, length),
              mean=c(0,0),var = c(1,1),
              meanDD = c(1,1),varDD = c(0.5,0.5),
              relAA = 1,
              corA = matrix(c(1,reldegree,reldegree,1),nrow = 2))

SP$setVarE(h2 = c(0.17,0.25)) # 0.17 0.25 in the range of heritability for growth, meat yield, survival, etc
SP$addSnpChip(nSnpPerChr = sapply(snp_pos, length)) # all non-QTL SNPs saved from simulation
SP$setSexes("yes_sys") # at the time of breeding, all individuals will only be one sex

pop_founder = newPop(founderPop, simParam=SP)

varD(pop_founder)/varP(pop_founder)
varD(pop_founder)/varA(pop_founder)
varAA(pop_founder)/varP(pop_founder)
varAA(pop_founder)/varA(pop_founder)
#varAA(pop_founder)
#varD(pop_founder)
#varA(pop_founder)
calmaf(pop_founder)



nDam = 30
nSire = 15
nCrosses = 30
nProgenyPerCross = 100
nProgeny = nProgenyPerCross
pop <- selectCross(pop_founder,
                   nFemale = nDam,nMale = nSire,
                   nCrosses = nCrosses,nProgeny = nProgenyPerCross,
                   use = "rand",
                   simParam = SP)


#tv <- quantile(pop@pheno[,1], probs = 0.4)
varD(pop)/varP(pop)
varD(pop)/varA(pop)
varAA(pop)/varP(pop)
varAA(pop)/varA(pop)
#varAA(pop)
#varD(pop)
#varA(pop)
calmaf(pop)
#tv = pop@pheno[,1][order(pop@pheno[,1],decreasing =FALSE)][pop@nInd*0.4]


ped = data.table()
phe = data.table()

for(g in c(1:3)){
  
  ped = rbind(ped,data.table(ID = pop@id,father = pop@father,mother = pop@mother))
  phe = rbind(phe, data.table(ID = pop@id,env = rep(1,pop@nInd),
                              gv1 = pop@gv[,1],gv2 = pop@gv[,2],
                              phe1 = pop@pheno[,1],phe2 = pop@pheno[,2],
                              FamilyID = paste(pop@father,"_",pop@mother,sep = "")))
  
  #pop = pop[(pop@pheno[,1]>tv)]

  females = selectWithinFam(pop,nInd=4,use = "rand",sex = "F")
  males = selectWithinFam(pop,nInd=2,use = "rand",sex = "M")
  candidate = c(females,males)
  
  #pop = selectInd(pop,nInd = 1000,use = "rand")
  
  pop <- selectCross(candidate,
                     nFemale = nDam,nMale = nSire,
                     nCrosses = nCrosses,nProgeny = nProgenyPerCross,
                     use = "rand",trait = 1,
                     simParam = SP)
  
}

varD(pop)/varP(pop)
varD(pop)/varA(pop)
varAA(pop)/varP(pop)
varAA(pop)/varA(pop)
#varAA(pop)
#varD(pop)
#varA(pop)
calmaf(pop)

ped = rbind(ped,data.table(ID = pop@id,father = pop@father,mother = pop@mother))
phe = rbind(phe, data.table(ID = pop@id,env = rep(1,pop@nInd),
                            gv1 = pop@gv[,1],gv2 = pop@gv[,2],
                            phe1 = pop@pheno[,1],phe2 = pop@pheno[,2],
                            FamilyID = paste(pop@father,"_",pop@mother,sep = "")))
tv <- quantile(pop@pheno[,1], probs = 0.4)
phe$phe_liability = phe$phe1
phe[phe1>tv,phe1:=1]
phe[phe1!=1,phe1:=0]

pop@id = paste("ind_",pop@id,sep = "")
phe$ID = paste("ind_",phe$ID,sep = "")
ped$ID = paste("ind_",ped$ID,sep = "")
ped$father = paste("ind_",ped$father,sep = "")
ped$mother = paste("ind_",ped$mother,sep = "")

ped = visPedigree::tidyped(ped = ped, cand = as.character(pop@id))
phe = phe[ID %in% pop@id,]

fwrite(ped[,1:3],file = "ped.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

fwrite(phe,file = "phe.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

pheno1 = phe[,.(ID,phe1,FamilyID)]
colnames(pheno1) = c("ID","phenotype","FamilyID")

pheno2 = phe[,.(ID,phe2,FamilyID)]
colnames(pheno2) = c("ID","phenotype","FamilyID")

fwrite(pheno1,file = "phenotype_t.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

fwrite(pheno2,file = "phenotype_c.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

genodt = pullSnpGeno(pop)

needname = rownames(genodt)

genodf = as.data.table(cbind(needname,genodt))

mapdt = makemap()

genodtput = makeped(z = genodt)

fwrite(genodtput,file="hib.ped",col.names = FALSE,row.names = FALSE, quote = FALSE,sep = "\t")

fwrite(genodf,file="geno.csv",col.names = TRUE,row.names = FALSE, quote = FALSE,sep = ",")

fwrite(mapdt, file = "hib.map",col.names = FALSE,row.names = FALSE, quote = FALSE, sep = " ")



#system('hiblup --make-xrm --threads 32 --pedigree ped.csv --add --dom --out Amat_b')
system('hiblup --make-xrm --threads 32 --pedigree ped.csv --add --out Amat_b')

system("plink --file hib --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib")

#system("hiblup --make-xrm --threads 32 --bfile hib --add --dom --out Gmat_b")
system("hiblup --make-xrm --threads 32 --bfile hib --add --out Gmat_b")


phe = fread("phe.csv",sep = ",")

rep = 10
fold = 5

count = 0
for(r in 1:rep){

    count = count + 1

    set.seed(randomseed[count])

    for(f in unique(phe$FamilyID)){
      fid = phe[FamilyID==f,ID]
      n = length(fid)
      n = floor(n/fold)*fold
      assign_ids = sample(fid, n)
      fold_idx = sample(rep(c(1:fold), length.out = n))
      phe[FamilyID==f,nfold:=fold_idx]
      remaining_ids <- setdiff(fid, assign_ids)
      if(length(remaining_ids) > 0){
        phe[ID %in% remaining_ids, nfold := NA]
      }
    }

    for(k in 1:fold){
      colnm = paste("r",r,"k",k,sep = "_")

      phe[nfold==k,(colnm):="test"]

      val_fold = ifelse(k < fold, k + 1, 1)

      phe[nfold == val_fold,(colnm):="val"]

      phe[!nfold %in% c(k,val_fold),(colnm):="train"]
    }

}
phe[, nfold := NULL]

fwrite(phe,file = "phe_cv.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

# nIndiv = nrow(phe)
# nInd_per_fold = floor(nIndiv/fold)
# for(r in 1:rep){
#   for(k in 1:fold){

#     count = count + 1
#     set.seed(randomseed[count])
#     colnm = paste("r",r,"k",k,sep = "_")

#     test_ind = c()
#     for(f in unique(phe$FamilyID)){
#       fid = phe[FamilyID==f,ID]
#       test_ind = c(test_ind,sample(fid,nInd_per_fold/length(unique(phe$FamilyID))))
#     }
#     phe[phe$ID %in% test_ind,(colnm):="test"]

#     phe_trainval = phe[!phe$ID %in% test_ind,] 

#     val_ind = c()
#     for(f in unique(phe_trainval$FamilyID)){
#       fid = phe_trainval[FamilyID==f,ID]
#       val_ind = c(val_ind,sample(fid,nInd_per_fold/length(unique(phe_trainval$FamilyID))))
#     }

#     phe[phe$ID %in% val_ind,(colnm):="val"] 

#     phe[!phe$ID %in% c(test_ind,val_ind),(colnm):="train"] 
#   }
# }
# fwrite(phe,file = "phe_cv.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")


phe = fread("phe_cv.csv",sep = ",")
r = 10
k = 5
cl = clnm(r,k)
phe_test = phe[phe[[cl]] == "test",]
phe_val = phe[phe[[cl]] == "val",]
fwrite(phe_test,file = "phe_test.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
fwrite(phe_val,file = "phe_val.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

tv_2_8 <- quantile(phe[phe1==1,phe2], probs = 0.8)

tv_2_5 <- quantile(phe[phe1==1,phe2], probs = 0.5)

tv_2_2 <- quantile(phe[phe1==1,phe2], probs = 0.2)


phe[ID %in% c(phe_val$ID,phe_test$ID),phe1 := NA]
phe[ID %in% c(phe_val$ID,phe_test$ID),phe2 := NA]
fwrite(phe,file = "phe_hi_continue.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")


phe = fread("phe_cv.csv",sep = ",")
phe[phe1==1 & phe2>tv_2_8,phe2:=tv_2_8]
phe[phe1==0,phe2 := NA]
phe[ID %in% c(phe_val$ID,phe_test$ID),phe1 := NA]
phe[ID %in% c(phe_val$ID,phe_test$ID),phe2 := NA]
fwrite(phe,file = "phe_hi.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")


phe = fread("phe_cv.csv",sep = ",")
phe[phe1==1 & phe2<=tv_2_5,phe2:=0]
phe[phe1==1 & phe2>tv_2_5,phe2:=1]
phe[phe1==0,phe2 := 0]
phe[ID %in% c(phe_val$ID,phe_test$ID),phe1 := NA]
phe[ID %in% c(phe_val$ID,phe_test$ID),phe2 := NA]
fwrite(phe,file = "phe_hi_th.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")



phe = fread("phe_cv.csv",sep = ",")

phe[, phe2 := fcase(
  phe1 == 1 & phe2 <= tv_2_2, 1,
  phe1 == 1 & phe2 > tv_2_2 & phe2 <= tv_2_5, 2,
  phe1 == 1 & phe2> tv_2_5 & phe2 <= tv_2_8, 3,
  phe1 == 1 & phe2 > tv_2_8, 4,
  phe1 == 0, 0
)]
phe[ID %in% c(phe_val$ID,phe_test$ID),phe1 := NA]
phe[ID %in% c(phe_val$ID,phe_test$ID),phe2 := NA]
fwrite(phe,file = "phe_hi_soft.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

#======#======#======#======#======#======#======#======#======#======#======#======#======#======#======#======
phename = "phe_hi.csv"

domA = ""
domG = ""

#domG = ",Gmat_b.GD"
#domA = ",Amat_b.PD"

system(paste("hiblup",
    "--multi-trait",
       "--pheno",
       phename,
        "--pheno-pos 5 6" ,
         paste("--xrm Gmat_b.GA",domG,sep=""),
         "--vc-method AI",
        "--ai-maxit 30",
         "--threads 32",
         "--ignore-cove",
         "--out G_hib",sep = " "))

system(paste("hiblup",
    "--single-trait",
       "--pheno",
       phename,
        "--pheno-pos 5" ,
         paste("--xrm Gmat_b.GA",domG,sep=""),
         "--vc-method AI",
        "--ai-maxit 30",
         "--threads 32",
         "--out G_hib_ST",sep = " "))

system(paste("hiblup",
    "--single-trait",
       "--pheno",
       phename,
        "--pheno-pos 6" ,
        paste("--xrm Gmat_b.GA",domG,sep=""),
         "--vc-method AI",
        "--ai-maxit 30",
         "--threads 32",
         "--out G_hib_SC",sep = " "))

system(paste("hiblup",
    "--multi-trait",
       "--pheno",
       phename,
        "--pheno-pos 5 6" ,
         paste("--xrm Amat_b.PA",domA,sep = ""),
         "--vc-method AI",
        "--ai-maxit 30",
         "--threads 32",
         "--ignore-cove",
         "--out A_hib",sep = " "))

system(paste("hiblup",
    "--single-trait",
       "--pheno",
       phename,
        "--pheno-pos 5" ,
         paste("--xrm Amat_b.PA",domA,sep = ""),
         "--vc-method AI",
        "--ai-maxit 30",
         "--threads 32",
         "--out A_hib_ST",sep = " "))

system(paste("hiblup",
    "--single-trait",
       "--pheno",
       phename,
        "--pheno-pos 6" ,
         paste("--xrm Amat_b.PA",domA,sep = ""),
         "--vc-method AI",
        "--ai-maxit 30",
         "--threads 32",
         "--out A_hib_SC",sep = " "))


#phe_dt = fread("phe_cv.csv",sep = ",")

#x_geno  = read_geno(dir = "geno.csv",phe_dt= phe_dt,r=r,k = k)

#y_train = get_phe(phe_dt=phe_dt,r=r,k=k)

#h2 <- fread("G_hib_ST.vars", sep = "\t")$h2[1]


#effect = calEBV(x = x_geno, y = y_train, h2=round(h2,3),type = "D")

#ef_dt= makedt(effect = effect, r= r, k=k,phe_dt = phe_dt)

#fwrite(ef_dt,"dominance.csv",sep = ",")


#effect = calEBV(x = x_geno, y = y_train, h2=round(h2,3),type = "E")

#ef_dt= makedt(effect = effect, r= r, k=k,phe_dt = phe_dt)

#fwrite(ef_dt,"epistatic.csv",sep = ",")


#effect = calEBV(x = x_geno, y = y_train, h2=round(h2,3),type = "A")

#ef_dt= makedt(effect = effect, r= r, k=k,phe_dt = phe_dt)

#fwrite(ef_dt,"additive.csv",sep = ",")