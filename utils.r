library(AlphaSimR)
library(data.table)

makeped = function(z){
    z[!z%in%c(0,1,2)] = -9
    z[z==2]= 22
    z[z==1]= 12
    z[z==0]= 11
    needname = rownames(z)
    z = as.data.table(cbind(needname,z))
    return(z)
  }

   #make .map file required by plink 
  makemap = function(...){
  mapck = getSnpMap()
  mapck$id = mapck$chr
  mapck$chr = rownames(mapck)
  mapck$site = mapck$pos
  mapck$pos = rep(0,nrow(mapck))
  return(mapck)
  }

clnm = function(r,k){
  return(colnm = paste("r",r,"k",k,sep = "_"))
}

calmaf = function(pop){
  x_geno = pullSnpGeno(pop)
  maf = apply(x_geno,2,function(x){
    maf = sum(x)/(2*length(x))
    if(maf>0.5) maf = 1 - maf
    return(maf)
  })
  return(mean(maf))
}

makedt = function(phe_dt,r,k,effect){
  cn = clnm(r,k)
  trainID = phe_dt$ID[phe_dt[[cn]]=="train"]
  testID = phe_dt$ID[phe_dt[[cn]]!="train"]

  ID = c(trainID, testID)
  ef_dt= data.table(ID = ID, effect = effect)
  
  
  # Remove the first column (ID column) and return only the SNP information
  return(ef_dt)
}
read_geno <- function(dir_geno, phe_dt,r,k) {
  cn = clnm(r,k)
  trainID = phe_dt$ID[phe_dt[[cn]]=="train"]
  testID = phe_dt$ID[phe_dt[[cn]]!="train"]
  # Load the genotype data
  geno <- fread(dir_geno, sep = ",", header = TRUE)
  
  # Reorder rows based on trainID and testID
  geno_ordered <- geno[match(c(trainID, testID), geno[[1]]), ]
  
  # Remove the first column (ID column) and return only the SNP information
  return(geno_ordered[, -1, with = FALSE])
}


get_phe = function(phe_dt,trainID,r,k){
  cn = clnm(r,k)

  return(phe_dt$phe1[phe_dt[[cn]]=="train"])
}

makeGDE = function(x,type = "D"){

  x = as.matrix(x)
  
  if(type=="E"){

      M = x - 1
      
      E <- 0.5 * ((M %*% t(M)) * (M %*% t(M))) - 0.5 * ((M * M) %*% t(M * M))
      
      E <- E / (sum(diag(E)) / nrow(E))
      
      A = diag(1,nrow(E))
      
      E = E * 0.99 + A * 0.01
      
      inverse <- solve(E)

  }else if(type == "D"){

      P = apply(x,2,function(col){
      
        pi = sum(col)/(2*length(col))
      
        if(pi>0.5) pi = 1-pi
      
        return(pi)
      
      })

      W = apply(x,2,function(col){
          
          pi = sum(col)/(2*length(col))
          
          if(pi>0.5){
          
            pi = 1-pi
          
            AA_bool = which(col == 0)
          
            aa_bool = which(col == 2)
          
            col[AA_bool] = 2
          
            col[aa_bool] = 0
          
          }
          
          aa = -2 * (pi^2)
          
          Aa = 2 * pi * (1-pi)
          
          AA = -2 * ((1-pi)^2)
          
          aa_bool = which(col == 0)
          
          Aa_bool = which(col == 1)
          
          AA_bool = which(col == 2)
          
          col[AA_bool] = AA
          
          col[aa_bool] = aa
          
          col[Aa_bool] = Aa
          
          return(col)
      
      })

      D <- (W %*% t(W)) / sum((2 * P * (1 - P))^2)
      
      A = diag(1,nrow(D))
      
      D = D * 0.99 + A * 0.01
      
      inverse <- solve(D)

  }else{
      P = apply(x,2,function(col){
        
        pi = sum(col)/(2*length(col))
        
        if(pi>0.5) pi = 1-pi
        
        return(pi)
     
      })
      Z = apply(x,2,function(col){
        pi = sum(col)/(2*length(col))

        if(pi>0.5){
        
          pi = 1-pi
        
          AA_bool = which(col == 0)
        
          aa_bool = which(col == 2)
        
          col[AA_bool] = 2
        
          col[aa_bool] = 0
        
        }
        
        col = col - 2*pi
        
        return(col)
      })
    
      G = (Z %*% t(Z)) / sum((2 * P * (1 - P)))
      
      A = diag(1,nrow(G))

      G = G * 0.99 + A * 0.01
      
      inverse <- solve(G)

  }
  return(inverse)
}

calEBV = function(x,y,h2,type = "D"){
  train_len = length(y)

  whole_len = nrow(x)
  
  Z = matrix(0, nrow = train_len, ncol = whole_len)
  
  # Set the diagonal elements to 1
  diag_indices = 1:train_len
  
  Z[cbind(diag_indices, diag_indices)] = 1
  
  mat = makeGDE(x = x,type = type)
  
  if(type == "D"){
  
    lamb = (1 - h2*0.1)/(h2*0.1)
  
  }else if(type=="E"){
  
    lamb = (1 - h2*0.1)/(h2*0.1)
  
  }else{
  
    lamb = (1 - h2)/h2
  
  }
  #glamb = (1 - h2)/h2
  #dlamb = (1 - h2*0.1)/(h2*0.1)
  #elamb = (1 - h2*0.1)/(h2*0.1)

  y = y - mean(y)
  
  Z_t =  t(Z)
  
  effect  = solve(Z_t %*% Z + mat * lamb) %*% Z_t %*% y

  return(effect)
}



calac = function(phe_test,predictions_dt,model_name,a_dt=NULL){

  prefixes <- c("cnn", "fc", "mul")
  suffixes <- c("_mt", "_soft", "_linear")

  combine_name <- paste(rep(prefixes, each = length(suffixes)), suffixes, sep = "")


  if(model_name %in% combine_name){
  predictions = "Prediction_t"
}else if(model_name %in% prefixes){
  predictions = "Prediction"
}else if(model_name %in% c("gblup_mt","gblup")){
  predictions = "Gmat_b.GA"
}else{
  predictions = "Amat_b.PA"
}

if(!is.null(a_dt)){
  ebv1 = predictions_dt[[predictions]][match(phe_test$ID,predictions_dt$ID)] + a_dt$Gmat_b.GA[match(phe_test$ID,a_dt$ID)]
  model_name = paste(model_name,"_a",sep = "")
}else{
  ebv1 = predictions_dt[[predictions]][match(phe_test$ID,predictions_dt$ID)]
}

  GV1 = phe_test$gv1


  acc  = cor(GV1,ebv1)

  output = data.table(mname = model_name, ac = acc,accov = cov(GV1,ebv1),va = var(ebv1),vg = var(GV1))


  return(output)

}



calWac = function(phe_test,predictions_dt,a_dt = NULL,model_name){

prefixes <- c("cnn", "fc", "mul")
suffixes <- c("_mt", "_soft", "_linear")

combine_name <- paste(rep(prefixes, each = length(suffixes)), suffixes, sep = "")


if(model_name %in% combine_name){
  predictions = "Prediction_t"
}else if(model_name %in% prefixes){
  predictions = "Prediction"
}else if(model_name %in% c("gblup_mt","gblup")){
  predictions = "Gmat_b.GA"
}else{
  predictions = "Amat_b.PA"
}

if(!is.null(a_dt)){
  ebv1 = predictions_dt[[predictions]][match(phe_test$ID,predictions_dt$ID)] + a_dt$Gmat_b.GA[match(phe_test$ID,a_dt$ID)]
}else{
  ebv1 = predictions_dt[[predictions]][match(phe_test$ID,predictions_dt$ID)]
}

GV1 = phe_test$gv1

fac = data.table()
fgv = c()
febv = c()
for( f in unique(phe_test$FamilyID)){
  phe_test_F = phe_test[FamilyID==f,]
  phe1_f = phe_test_F$phe1
  GV1_f = phe_test_F$gv1

  if(!is.null(a_dt)){

    ebv1_f = predictions_dt[[predictions]][match(phe_test_F$ID,predictions_dt$ID)] +a_dt$Gmat_b.GA[match(phe_test_F$ID,a_dt$ID)]

  }else{
    ebv1_f = predictions_dt[[predictions]][match(phe_test_F$ID,predictions_dt$ID)]
  }
  if(! model_name %in%c("pblup","pblup_mt")){
      acg = cor(GV1_f,ebv1_f)
      acp = cor(phe1_f,ebv1_f)
  }else{
      acg = 0
      acp = 0
  }

  meanp = mean(phe1_f)
  meang = mean(GV1_f)
  meanebv = mean(ebv1_f)


  ebv1[match(phe_test_F$ID,phe_test$ID)] = ebv1[match(phe_test_F$ID,phe_test$ID)] - meanebv

  GV1[match(phe_test_F$ID,phe_test$ID)] = GV1[match(phe_test_F$ID,phe_test$ID)] - meang

  
  fac = rbind(fac, data.table(FamilyID = f, acg = acg, acp = acp,phe = meanp,gv = meang, ebv = meanebv ))
  fgv = c(fgv,rep(meang,nrow(phe_test_F)))
  febv = c(febv,rep(meanebv,nrow(phe_test_F)))
}

if(!is.null(a_dt)) model_name = paste(model_name,"_a",sep = "")

  if(! model_name %in%c("pblup","pblup_mt")){
     output = data.table(mname = model_name,wac =mean(fac$acg) ,minac = min(fac$acg),maxac = max(fac$acg),sdac = sd(fac$acg),
 famacc = cor(fac$gv,fac$ebv), wfac = cor(ebv1,GV1),fcov = cov(fgv,febv),wcov = cov(ebv1,GV1),fwcov = cov(fgv,ebv1),wfcov = cov(GV1,febv),
 fva = var(febv),fvg = var(fgv),wva = var(ebv1),wvg = var(GV1),fwcova = cov(febv,ebv1),fwcovg = cov(fgv,GV1))
  }else{
     output = data.table(mname = model_name,wac =0 ,minac = 0,maxac = 0,sdac = 0,
 famacc = cor(fac$gv,fac$ebv), wfac = 0,fcov = cov(fgv,febv),wcov = 0,fwcov = 0,wfcov = cov(GV1,febv),
 fva = var(febv),fvg = var(fgv),wva = 0,wvg = var(GV1),fwcova = 0,fwcovg = cov(fgv,GV1))
  }

 return(output)

}


calcvmean = function(dt){
  result <- dt[, lapply(.SD, mean), by = mname]
  return(result)
}


extractMeanSd = function(dt){
  result_mean <- dt[, lapply(.SD, mean), by = mname]

  result_sd <- dt[, lapply(.SD, sd), by = mname]


  bool = (colnames(result_sd)!="mname")

  colnames(result_sd)[bool] = paste(colnames(result_sd)[bool],"_sd",sep = "")

  result <- merge(result_mean, result_sd, by = "mname")
  
  return(result)
}