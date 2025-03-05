rm(list = ls())
gc()

source("utils.r")

dir = "/home/kangziyi/poster/BaseMapHaplo/"
mapFile = "mergedMap.csv"
haploFile = "mergedHaplo.csv"

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

SP$addTraitA(nQtlPerChr = sapply(qtl_pos, length),mean=c(0,0),var = c(1,1),corA = matrix(c(1,0.3,0.3,1),nrow = 2))
SP$setVarE(h2=c(0.17,0.25)) # in the range of heritability for growth, meat yield, survival, etc
SP$addSnpChip(nSnpPerChr = sapply(snp_pos, length)) # all non-QTL SNPs saved from simulation
SP$setSexes("yes_sys") # at the time of breeding, all individuals will only be one sex

pop_founder = newPop(founderPop, simParam=SP)


nDam = 50
nSire = 25
nCrosses = 50
nProgenyPerCross = 1020
nProgeny = nProgenyPerCross
pop <- selectCross(pop_founder,
                   nFemale = nDam,nMale = nSire,
                   nCrosses = nCrosses,nProgeny = nProgenyPerCross,
                   use = "rand",
                   simParam = SP)


tv <- quantile(pop@pheno[,1], probs = 0.4)

#tv = pop@pheno[,1][order(pop@pheno[,1],decreasing =FALSE)][pop@nInd*0.4]


ped = data.table()
phe = data.table()
for(g in c(1:3)){
  
  ped = rbind(ped,data.table(ID = pop@id,father = pop@father,mother = pop@mother))
  phe = rbind(phe, data.table(ID = pop@id,env = rep(1,pop@nInd),
                              gv1 = pop@gv[,1],gv2 = pop@gv[,2],
                              phe1 = pop@pheno[,1],phe2 = pop@pheno[,2],
                              FamilyID = paste(pop@father,"_",pop@mother,sep = "")))
  
  pop = pop[(pop@pheno[,1]>tv)]
  
  pop = selectInd(pop,nInd = 1000,use = "rand")
  
  pop <- selectCross(pop,
                     nFemale = nDam,nMale = nSire,
                     nCrosses = nCrosses,nProgeny = nProgenyPerCross,
                     use = "rand",trait = 1,
                     simParam = SP)
  
}

ped = rbind(ped,data.table(ID = pop@id,father = pop@father,mother = pop@mother))
phe = rbind(phe, data.table(ID = pop@id,env = rep(1,pop@nInd),
                            gv1 = pop@gv[,1],gv2 = pop@gv[,2],
                            phe1 = pop@pheno[,1],phe2 = pop@pheno[,2],
                            FamilyID = paste(pop@father,"_",pop@mother,sep = "")))

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



system('hiblup --make-xrm --threads 32 --pedigree ped.csv --add --out Amat_b')

system("plink --file hib --geno --mind --no-sex --no-pheno --no-fid --no-parents --nonfounders --make-bed --chr-set 44 --out hib")

system("hiblup --make-xrm --threads 32 --bfile hib --add --out Gmat_b")

#system("hiblup --trans-xrm --xrm Gmat_b.GA --out Gmat")

#system("hiblup --trans-xrm --xrm Amat_b.PA --out Amat")


phe = fread("phe.csv",sep = ",")

test_ind = c()
for(f in unique(phe$FamilyID)){
  fid = phe[FamilyID==f,ID]
  test_ind = c(test_ind,sample(fid,1000/length(unique(phe$FamilyID))))
}

#Amat = fread("Amat.txt",sep = "\t")
#A_id  = fread("Amat.id.txt",header = FALSE)
#Amat = Amat[, -ncol(Amat), with = FALSE]

#Gmat = fread("Gmat.txt",sep = "\t")
#G_id  = fread("Gmat.id.txt",sep = "\t")
#Gmat = Gmat[, -ncol(Gmat), with = FALSE]



#Amat_test= as.matrix(Amat)
#rownames(Amat_test) = A_id$V1
#colnames(Amat_test) = A_id$V1

#Gmat_test= as.matrix(Gmat)
#rownames(Gmat_test) = G_id$V1
#colnames(Gmat_test) = G_id$V1

#Amat = cbind(colnames(Amat_test),Amat_test)
#Amat = as.data.table(Amat)
#Gmat = cbind(colnames(Gmat_test),Gmat_test)
#Gmat = as.data.table(Gmat)


#Amat_test = Amat_test[test_ind,test_ind]
#Amat_test = cbind(colnames(Amat_test),Amat_test)
#Amat_test = as.data.table(Amat_test)
#Gmat_test = Gmat_test[test_ind,test_ind]
#Gmat_test = cbind(colnames(Gmat_test),Gmat_test)
#Gmat_test = as.data.table(Gmat_test)

phe_test = phe[phe$ID %in% test_ind,] 
phe_trainval = phe[!phe$ID %in% test_ind,] 

val_ind = c()
for(f in unique(phe_trainval$FamilyID)){
  fid = phe_trainval[FamilyID==f,ID]
  val_ind = c(val_ind,sample(fid,1000/length(unique(phe_trainval$FamilyID))))
}

phe_val = phe_trainval[phe_trainval$ID %in% val_ind,] 




phe[ID %in% c(phe_val$ID,phe_test$ID),phe1 := NA]
phe[ID %in% c(phe_val$ID,phe_test$ID),phe2 := NA]
fwrite(phe,file = "phe_hi_continue.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
tv_2_8 <- quantile(phe[phe1==1,phe2], probs = 0.8)
tv_2_5 <- quantile(phe[phe1==1,phe2], probs = 0.5)
if(!tv_2_8>0) warnings("tv_2 is nagetive value")
phe[phe1==1 & phe2>tv_2_8,phe2:=tv_2_8]
phe[phe1==0,phe2 := NA]

#fwrite(Amat, "Amat.csv", sep = ",", col.names = FALSE)
#fwrite(Amat_test,"Amat_test.csv",sep = ",", col.names = FALSE)

#fwrite(Gmat, "Gmat.csv", sep = ",", col.names = FALSE)
#fwrite(Gmat_test,"Gmat_test.csv",sep = ",", col.names = FALSE)

fwrite(phe_test,file = "phe_test.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
fwrite(phe_val,file = "phe_val.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
fwrite(phe,file = "phe_hi.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

if(!tv_2_5>0) warnings("tv_2 is nagetive value")
phe[phe1==1 & phe2<=tv_2_5,phe2:=0]
phe[phe1==1 & phe2>tv_2_5,phe2:=1]
phe[phe1==0,phe2 := 0]
if(nrow(phe[phe1%in%c(0,1),])!=1500) warnings("not the same traning set")
fwrite(phe,file = "phe_hi_th.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
#===============================================================================================================

phename = "phe_hi.csv"
phename = "phe_hi_th.csv"
phename = "phe_hi_continue.csv"
phename = "phe_hi_soft.csv"
domG = ",Gmat_b.GD"
domA = ",Amat_b.PD"
domA = ""
domG = ""
system(paste("hiblup",
    "--multi-trait",
       "--pheno",
       phename,
        "--pheno-pos 5 6" ,
         paste("--xrm Gmat_b.GA",domG,sep=""),
         "--vc-method AI",
        "--ai-maxit 30",
         "--threads 32",
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

#system(paste("python3 performDL.py -mr Amat.csv",
  #          "-pr phe.csv",
 #           "-me Amat_test.csv",
  #          "-pe phe_test.csv",
 #           "-o /home/kangziyi/temp",
 #           "-p 5 6",
 #           "-t 5",
#            sep = " " ))

#system(paste("python3 performDL.py -mr Gmat.csv",
           # "-pr phe.csv",
           # "-me Gmat_test.csv",
           # "-pe phe_test.csv",
           # "-o /home/kangziyi/temp",
           # "-p 5 6",
           # "-t 5",
           # sep = " " ))

#===============================================================================================================
source("utils.r")
r = 10
k = 5
phe_test = fread("phe_test.csv",sep = ",")

#======================================================
#phe_dt = fread("phe_cv.csv",sep = ",")
#nm = clnm(r = 10, k = 5)
#phe_test = phe_dt[phe_dt[[nm]]!="train",]
#======================================================

#ef_dt_d = fread("dominance.csv",sep = ",")
#ef_dt_a = fread("additive.csv",sep = ",")
#ef_dt_e = fread("epistatic.csv",sep = ",")

#d = ef_dt_d$effect[match(phe_test$ID,ef_dt_d$ID)]
#a = ef_dt_a$effect[match(phe_test$ID,ef_dt_a$ID)]
#e = ef_dt_e$effect[match(phe_test$ID,ef_dt_e$ID)]



pblup_mt = fread("A_hib.phe1.rand",sep = "\t")
pblup = fread("A_hib_ST.rand",sep = "\t")
#pblup_mt = fread("A_hib.phe2.rand",sep = "\t")
#pblup= fread("A_hib_SC.rand",sep = "\t")
#ADL1 = fread("predictions_trait_1.csv",sep = ",")
#ADL2 = fread("predictions_trait_2.csv",sep = ",")

gblup_mt= fread("G_hib.phe1.rand",sep = "\t")
gblup = fread("G_hib_ST.rand",sep = "\t")
A_dt = gblup

#Ghib2 = fread("G_hib.phe2.rand",sep = "\t")
#Ghib4 = fread("G_hib_SC.rand",sep = "\t")
#GDL1 = fread("predictions_trait_1.csv",sep = ",")
#GDL2 = fread("predictions_trait_2.csv",sep = ",")

output_pac  = data.table()

output_wac = data.table()


output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = pblup,model_name = "pblup"))


output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = pblup_mt,model_name = "pblup_mt"))


output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = gblup,model_name = "gblup"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = gblup,model_name = "gblup"))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = gblup_mt,model_name = "gblup_mt"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = gblup_mt,model_name = "gblup_mt"))

output_wac 
output_pac


cnn = fread("predictions_t_cnn.csv",sep = ",")

fc= fread("predictions_t_fc.csv",sep = ",")

mul= fread("predictions_t_mul.csv",sep = ",")

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = cnn,model_name = "cnn"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = cnn,model_name = "cnn"))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = fc,model_name = "fc"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = fc,model_name = "fc"))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = mul,model_name = "mul"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = mul,model_name = "mul"))


output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = cnn,model_name = "cnn",a_dt = A_dt))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = cnn,model_name = "cnn",a_dt = A_dt))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = fc,model_name = "fc",a_dt = A_dt))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = fc,model_name = "fc",a_dt = A_dt))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = mul,model_name = "mul",a_dt = A_dt))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = mul,model_name = "mul",a_dt = A_dt))


cnn_mt= fread("predictions_c_t_cnn.csv",sep = ",")

fc_mt= fread("predictions_c_t_fc.csv",sep = ",")

mul_mt= fread("predictions_c_t_mul.csv",sep = ",")


output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = cnn_mt,model_name = "cnn_mt"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = cnn_mt,model_name = "cnn_mt"))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = fc_mt,model_name = "fc_mt"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = fc_mt,model_name = "fc_mt"))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = mul_mt,model_name = "mul_mt"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = mul_mt,model_name = "mul_mt"))


output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = cnn_mt,model_name = "cnn_mt",a_dt = A_dt))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = cnn_mt,model_name = "cnn_mt",a_dt = A_dt))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = fc_mt,model_name = "fc_mt",a_dt = A_dt))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = fc_mt,model_name = "fc_mt",a_dt = A_dt))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = mul_mt,model_name = "mul_mt",a_dt = A_dt))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = mul_mt,model_name = "mul_mt",a_dt = A_dt))




cnn_soft= fread("predictions_c_t_cnn_soft.csv",sep = ",")

fc_soft= fread("predictions_c_t_fc_soft.csv",sep = ",")

mul_soft= fread("predictions_c_t_mul_soft.csv",sep = ",")


output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = cnn_soft,model_name = "cnn_soft"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = cnn_soft,model_name = "cnn_soft"))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = fc_soft,model_name = "fc_soft"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = fc_soft,model_name = "fc_soft"))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = mul_soft,model_name = "mul_soft"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = mul_soft,model_name = "mul_soft"))


output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = cnn_soft,model_name = "cnn_soft",a_dt = A_dt))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = cnn_soft,model_name = "cnn_soft",a_dt = A_dt))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = fc_soft,model_name = "fc_soft",a_dt = A_dt))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = fc_soft,model_name = "fc_soft",a_dt = A_dt))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = mul_soft,model_name = "mul_soft",a_dt = A_dt))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = mul_soft,model_name = "mul_soft",a_dt = A_dt))


setorder(output_wac,wac)
setorder(output_pac,ac)
output_wac 
output_pac




cnn_linear= fread("predictions_c_t_cnn_linear.csv",sep = ",")

fc_linear= fread("predictions_c_t_fc_linear.csv",sep = ",")

mul_linear= fread("predictions_c_t_mul_linear.csv",sep = ",")


output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = cnn_linear,model_name = "cnn_linear"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = cnn_linear,model_name = "cnn_linear"))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = fc_linear,model_name = "fc_linear"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = fc_linear,model_name = "fc_linear"))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = mul_linear,model_name = "mul_linear"))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = mul_linear,model_name = "mul_linear"))


output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = cnn_linear,model_name = "cnn_linear",a_dt = A_dt))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = cnn_linear,model_name = "cnn_linear",a_dt = A_dt))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = fc_linear,model_name = "fc_linear",a_dt = A_dt))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = fc_linear,model_name = "fc_linear",a_dt = A_dt))

output_wac  = rbind(output_wac,calWac(phe_test = phe_test,predictions_dt = mul_linear,model_name = "mul_linear",a_dt = A_dt))

output_pac = rbind(output_pac,calac(phe_test = phe_test,predictions_dt = mul_linear,model_name = "mul_linear",a_dt = A_dt))


setorder(output_wac,wac)
setorder(output_pac,ac)
output_wac 
output_pac



# Initialize result tables
output_wac <- data.table()
output_pac <- data.table()

# Load necessary data
datasets <- list(
  pblup_mt = "A_hib.phe1.rand",
  pblup = "A_hib_ST.rand",
  gblup_mt = "G_hib.phe1.rand",
  gblup = "G_hib_ST.rand",
  cnn = "predictions_t_cnn.csv",
  fc = "predictions_t_fc.csv",
  mul = "predictions_t_mul.csv",
  cnn_mt = "predictions_c_t_cnn.csv",
  fc_mt = "predictions_c_t_fc.csv",
  mul_mt = "predictions_c_t_mul.csv",
  cnn_soft = "predictions_c_t_cnn_soft.csv",
  fc_soft = "predictions_c_t_fc_soft.csv",
  mul_soft = "predictions_c_t_mul_soft.csv"
)

# Read datasets into a named list
data <- lapply(datasets, function(file) {
  sep <- ifelse(grepl("\\.csv$", file), ",", "\t")
  fread(file, sep = sep)
})

# Define models and prediction datasets
models <- c("pblup","pblup_mt","gblup", "gblup_mt", "cnn", "fc", "mul", "cnn_mt", "fc_mt", "mul_mt","cnn_soft", "fc_soft", "mul_soft")

a_dt_models <- c("cnn", "fc", "mul", "cnn_mt", "fc_mt", "mul_mt","cnn_soft", "fc_soft", "mul_soft")

# Append results for each model
for (model in models) {
  # Add results without `a_dt`
  if(! model %in%c("pblup","pblup_mt")) output_wac <- rbind(output_wac, calWac(phe_test = phe_test, predictions_dt = data[[model]], model_name = model))

  output_pac <- rbind(output_pac, calac(phe_test = phe_test, predictions_dt = data[[model]], model_name = model))

  # Add results with `a_dt` if applicable
  if (model %in% a_dt_models) {
    output_wac <- rbind(output_wac, calWac(phe_test = phe_test, predictions_dt = data[[model]], model_name = model, a_dt = data[["gblup"]]))
    output_pac <- rbind(output_pac, calac(phe_test = phe_test, predictions_dt = data[[model]], model_name = model, a_dt = data[["gblup"]]))
  }
}

# Display results
setorder(output_wac,wac)
setorder(output_pac,ac)
output_wac
output_pac


