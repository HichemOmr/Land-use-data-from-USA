
# clear all existing objects  

rm(list = ls())

# Function to Install and Load R Packages

# Specify the list of required packages to be installed and load    

p = c("pdist", "fields", "TOC", "R.matlab", "nnet") 
# , "ggplot2", "maps", "reshape", "gridExtra", "Rcpp", "sp", "raster", "caret");

install_package<-function(pack)
{if(!(pack %in% row.names(installed.packages())))
{
  update.packages(ask=F)
  install.packages(pack,dependencies=F) # T
}
 require(pack,character.only=TRUE)
}


# Call the Function install_package

for(pack in p) {install_package(pack)}


# Call all needed functions used to perform land use analysis, data normalisation, data split, calibration, validation, mapping ... 

source("needed_funtions1.R")

# read_data is a function to read a dataset in MATLAB format (mat file)

read_data <- function(matfile, xvr, yvr, IndLU_t1, IndLU_t2, CodeNU, CodeU){

	   inputFile1 = paste("data/", matfile, sep="") 
         data1 = loadMat(inputFile1)
	   data1 = data1[complete.cases(data1), ]

########################################
# initialisation 
# xvr = 6:11; yvr = 14; ratio = 0.7

# m = length(xvr)

change1    =  data1[ ((data1[,IndLU_t1]==CodeNU) & (data1[,IndLU_t2]==CodeU)) ,] # 4 and 5
no_change1 =  data1[ (data1[,IndLU_t1]==CodeNU & data1[,IndLU_t2]==CodeNU) ,]
# exc_layer1 =  data1[ ( (data1[,IndLU_t1]==CodeU)  & (data1[,IndLU_t2]==CodeNU) ) | ( (data1[,IndLU_t1]==CodeU)  & (data1[,IndLU_t2]==CodeU) ),] 

# rrr = balanced_dataset(change1, no_change1, xvr, yvr) ####
# ch = rrr$change 						  ###
# no_ch = rrr$no_change			                    ### 

r = list(ch = change1, no_ch = no_change1)

return(r)

}

run_ltm_and_ltmCluster <- function(change1, no_change1, xvr, m, yvr, ratio, K) { # IndLU_t1, IndLU_t2, CodeNU, CodeU
# r$change, r$n_ch

dataset1 = rbind( cbind(change1, rep(1, nrow(change1)) ), cbind(no_change1, rep(0, nrow(no_change1) ))) ### 

dataset1N <- fn.normalize(dataset1, xvr, yvr) # dataset1 


	# data split 
	# LT = fn.splitSRS(dataset1N, yvr, ratio) # 
      LT = splitdt(dataset1N, ratio)
	L = LT$L
	T = LT$T 
      # print("size of T set: ")
      # print(nrow(T))

      # print("number of changed cells in T set: ")
      # print(nrow(T[T[,yvr]==1,]))

# colnames(L) <- c("id", "xcoord","ycoord", "LUt1", "LUt2", "x1", "x2", "x3", "x4", "x5", "x6", "var1", "var2", "Output")
# colnames(T) <- c("id", "xcoord","ycoord", "LUt1", "LUt2", "x1", "x2", "x3", "x4", "x5", "x6", "var1", "var2", "Output")

########################################
# LTM, 		without clustering prior learning 
########################################
# install.packages('TOC') # library(TOC)

##############################################################################################################################
###############################################################
# run the LTM model, without clustering prior learning 
###############################################################
##############################################################################################################################

#############
threshold = round(nrow(dataset1N[dataset1N[,yvr]==1,]) / nrow(dataset1N), 2)  # T 
print("th = ") 
print(threshold) 
#############

newT = fnLUCAnalysis_maxP(L, dataset1N, xvr, yvr, m, threshold) # dataset1N, T, rbind(L,T)
## newT = fnLUCAnalysis(L, dataset1N, xvr, yvr, m) # dataset1N, T, rbind(L,T)

###

cm = table(newT[,yvr], newT[,yvr+2]) 
print("cross-tab from the LTM model: ")
print(cm)

errors = apply(newT, 1, metrics, yvr, yvr+2)

# Total Operating Characteristic - TOC and ROC curves
# TOC performance curves for the a) LTM and b) LTM-Cluster

tocd <- TOC(newT[,yvr], newT[,yvr+2], mask=NULL, nthres = 100)
# rocd <- ROC(newT[,yvr], newT[,yvr+2], mask=NULL, nthres = 100)

pcm = pcm_evaluation(cm)
# mean(pcm)


##############################################################################################################################
# run the LTM-Cluster model, with clustering prior learning 
##############################################################################################################################

# S = dataset1N  # dataset1N # rbind(L, T) 		# the entire dataset 
# K=3
# test the model on testing set T

resk = ltm_cluster(L, dataset1N, xvr, yvr, K) # S set is the entire dataset, without the cluster number at the end of column 

newT_cl = rbind(resk[[1]]$r_k, resk[[2]]$r_k, resk[[3]]$r_k)

cm_ltm_cluster = table(newT_cl[,yvr], newT_cl[,yvr+3]) # ypred

print("cross-tab from the LTM-Cluster model: ")
print(cm_ltm_cluster)

errors_ltm_cluster = apply(newT_cl, 1, metrics, yvr, yvr+3) 

# TOC performance curves for the a) LTM and b) LTM-Cluster

tocd_ltm_cl <- TOC(newT_cl[,yvr], newT_cl[,yvr+3], mask=NULL, nthres = 100)
# rocd_ltm_cl <- ROC(newT_cl[,yvr], newT_cl[,yvr+3], mask=NULL, nthres = 100)

pcm_ltm_cluster  = pcm_evaluation(cm_ltm_cluster )
# mean(pcm_ltm_cluster)

res= list(pcm_ltm = pcm, cm_ltm = cm, errors_ltm = errors, toc_ltm =tocd,  
change_ltm = change1, no_change_ltm = no_change1, newT = newT, 
pcm_ltm_cluster = pcm_ltm_cluster, cm_ltm_cluster = cm_ltm_cluster, errors_ltm_cluster = errors_ltm_cluster, toc_ltm_cluster =tocd_ltm_cl,  
newT_cl= newT_cl)

return(res)
}


fct_result <-  function (res_case, IT) { # IT is the number of replications 

auc_ltm = c()
auc_ltm_cluster = c()
pcm_ltm = c()
pcm_ltm_cluster = c()
ER_ltm = c()
ER_ltm_cluster = c()

for (i in 1:IT) {

auc_ltm = c( auc_ltm, res_case[,i]$toc_ltm@AUC )
auc_ltm_cluster = c( auc_ltm_cluster, res_case[,i]$toc_ltm_cluster@AUC )

pcm_ltm = c( pcm_ltm, mean(res_case[,i]$pcm_ltm))
pcm_ltm_cluster = c( pcm_ltm_cluster, mean(res_case[,i]$pcm_ltm_cluster))

ER_ltm = c(ER_ltm, round((res_case[,i]$cm_ltm[1,2] + res_case[,i]$cm_ltm[2,1]) * 100 / sum(res_case[,i]$cm_ltm),3) )
ER_ltm_cluster = c(ER_ltm_cluster, round((res_case[,i]$cm_ltm_cluster[1,2] + res_case[,i]$cm_ltm_cluster[2,1]) * 100 / sum(res_case[,i]$cm_ltm_cluster),3) )

}

print(auc_ltm)
print(auc_ltm_cluster)

print(ER_ltm)
print(ER_ltm_cluster)

print(pcm_ltm)
print(pcm_ltm_cluster)

auc_ltmMSD = list(mean = mean(auc_ltm) , sd = sd(auc_ltm))
pcm_ltmMSD = list(mean = mean(pcm_ltm), sd = sd(pcm_ltm))
ER_ltmMSD = list(mean = mean(ER_ltm), sd = sd(ER_ltm))

auc_ltm_clusterMSD = list(mean = mean(auc_ltm_cluster), sd = sd(auc_ltm_cluster))
pcm_ltm_clusterMSD = list(mean = mean(pcm_ltm_cluster), sd = sd(pcm_ltm_cluster))
ER_ltm_clusterMSD = list(mean = mean(ER_ltm_cluster), sd = sd(ER_ltm_cluster))


result_mean = rbind( c(AUC = auc_ltmMSD$mean, PCM = pcm_ltmMSD$mean, ER = ER_ltmMSD$mean), 
c(auc_ltm_clusterMSD$mean, pcm_ltm_clusterMSD$mean, ER_ltm_clusterMSD$mean) ) 

result_sd =   rbind( c(AUC = auc_ltmMSD$sd, PCM = pcm_ltmMSD$sd, ER = ER_ltmMSD$sd), 
c(auc_ltm_clusterMSD$sd, pcm_ltm_clusterMSD$sd, ER_ltm_clusterMSD$sd) ) 

row.names(result_mean) <- c("res_ltm_mean", "res_ltm_cluster_mean" )
row.names(result_sd) <- c("res_ltm_sd", "res_ltm_cluster_sd" )

result = list(mean = result_mean, sd = result_sd)

return(result)
}


