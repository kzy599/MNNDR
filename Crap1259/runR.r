rm(list = ls())
gc()

source("utils.r")
percent = ""#
#========================================================================================================

result_pac = data.table()
result_wac = data.table()

for(r in c(1:10)){

  output_pac  = data.table()

  output_wac = data.table()


    for( k in c(1:5)){

        phe = fread("phe_cv.csv",sep = ",")
        cl = clnm(r,k)
        phe_test = phe[phe[[cl]] == "test",]
        phe_val = phe[phe[[cl]] == "val",]
        fwrite(phe_test,file = "phe_test.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")
        fwrite(phe_val,file = "phe_val.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

        phe[phe2==-999,phe2 := NA]
        phe[ID %in% c(phe_val$ID,phe_test$ID),phe1 := NA]
        phe[ID %in% c(phe_val$ID,phe_test$ID),phe2 := NA]
        fwrite(phe,file = "phe_hi.csv",sep = ",",col.names = TRUE,row.names = FALSE,quote = FALSE,na = "NA")

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
                "--pheno-pos 3 4" ,
                "--dcovar 5 5",
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
                "--pheno-pos 3" ,
                "--dcovar 4,5",
                paste("--xrm Gmat_b.GA",domG,sep=""),
                "--vc-method AI",
                "--ai-maxit 30",
                "--threads 32",
                "--out G_hib_ST",sep = " "))

        #   system(paste("hiblup",
        #     "--multi-trait",
        #       "--pheno",
        #       phename,
        #         "--pheno-pos 3 4" ,
        #         paste("--xrm Amat_b.PA",domA,sep = ""),
        #         "--vc-method AI",
        #         "--ai-maxit 30",
        #         "--threads 32",
        #         "--ignore-cove",
        #         "--out A_hib",sep = " "))

        # system(paste("hiblup",
        #     "--single-trait",
        #       "--pheno",
        #       phename,
        #         "--pheno-pos 3" ,
        #         paste("--xrm Amat_b.PA",domA,sep = ""),
        #         "--vc-method AI",
        #         "--ai-maxit 30",
        #         "--threads 32",
        #         "--out A_hib_ST",sep = " "))


        system('bash -c "source activate tf-gpu && python main.py --Th"')

        # system('bash -c "source activate tf-gpu && python main.py --Th --doubleT"')
        
        system('bash -c "source activate tf-gpu && python main.py --Th --doubleT --linear"')


        # Load necessary data
        datasets <- list(
          gblup_mt = "G_hib.phe1.rand",
          gblup = "G_hib_ST.rand",
          cnn = "predictions_t_cnn.csv",
          fc = "predictions_t_fc.csv",
          mul = "predictions_t_mul.csv",
          cnn_linear = "predictions_c_t_cnn_linear.csv",
          fc_linear = "predictions_c_t_fc_linear.csv",
          mul_linear = "predictions_c_t_mul_linear.csv"
        )

        # Read datasets into a named list
        data <- lapply(datasets, function(file) {
          sep <- ifelse(grepl("\\.csv$", file), ",", "\t")
          fread(file, sep = sep)
        })

        # Define models and prediction datasets
        #models <- c("pblup","pblup_mt","gblup", "gblup_mt", "cnn", "fc", "mul", "cnn_mt", "fc_mt", "mul_mt","cnn_soft", "fc_soft", "mul_soft","cnn_linear","fc_linear","mul_linear")
        models <- c("gblup", "gblup_mt", "cnn", "fc", "mul", "cnn_linear", "fc_linear", "mul_linear")
        #a_dt_models <- c("cnn", "fc", "mul", "cnn_mt", "fc_mt", "mul_mt","cnn_soft", "fc_soft", "mul_soft","cnn_linear","fc_linear","mul_linear")
        # a_dt_models <- c("cnn", "fc", "mul", "cnn_mt", "fc_mt", "mul_mt","cnn_linear","fc_linear","mul_linear")
        # Append results for each model
        for (model in models) {
          # Add results without `a_dt`
          #if(! model %in%c("pblup","pblup_mt")) output_wac <- rbind(output_wac, calWac(phe_test = phe_test, predictions_dt = data[[model]], model_name = model))
          # output_wac <- rbind(output_wac, calWac(phe_test = phe_test, predictions_dt = data[[model]], model_name = model))
          output_pac <- rbind(output_pac, calac(phe_test = phe_test, predictions_dt = data[[model]], model_name = model))

          # # Add results with `a_dt` if applicable
          # if (model %in% a_dt_models) {
          #   output_wac <- rbind(output_wac, calWac(phe_test = phe_test, predictions_dt = data[[model]], model_name = model, a_dt = data[["gblup"]]))
          #   output_pac <- rbind(output_pac, calac(phe_test = phe_test, predictions_dt = data[[model]], model_name = model, a_dt = data[["gblup"]]))
          # }
        }

        # Display results
        # setorder(output_wac,wac)
        setorder(output_pac,ac2)
        # output_wac
        output_pac
    }

    result_pac = rbind(result_pac, calcvmean(dt = output_pac))

    # result_wac = rbind(result_wac, calcvmean(dt = output_wac))
}


pac = extractMeanSd(dt = result_pac)

# wac = extractMeanSd(dt = result_wac)

fwrite(pac,paste("pac_mean_sd",percent,".csv",sep = ""),sep = ",")

# fwrite(wac,"wac_mean_sd.csv",sep = ",")


fwrite(result_pac,paste("pac",percent,".csv",sep = ""),sep = ",")

# fwrite(result_wac,"wac.csv",sep = ",")

warnings()


