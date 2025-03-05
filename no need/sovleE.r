source("utils.r")

r  = 10

k = 5

phe_dt = fread("phe_cv.csv",sep = ",")

x_geno  = read_geno(dir = "geno.csv",phe_dt= phe_dt,r=r,k = k)

y_train = get_phe(phe_dt=phe_dt,r=r,k=k)

h2 <- fread("G_hib_ST.vars", sep = "\t")$h2[1]

effect = calEBV(x = x_geno, y = y_train, h2=round(h2,3),type = "D")

ef_dt= makedt(effect = effect, r= r, k=k,phe_dt = phe_dt)

fwrite(ef_dt,"dominance.csv",sep = ",")