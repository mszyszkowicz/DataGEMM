########################################################################
##
## Nonlinear concentraton-response function: model fitting and plotting procedure
## 
## Date: Sept 19, 2016
## Version 2.16
## 
## Purpose: This program will fit a series of nonlinear concentration-response functions using
##		Cox proportional hazards modeling to identify optimal nonlinear relationship
##
## It produces 3 outputs: (1) a graph showing nonlinear relationship based on best fitting model and an essemble model,
## 				  (2) a summary table listing coefficient, standard error, loglik, and function form of each cox model that was examined, and
##				  (3) a descriptive table showing distribution of air pollution exposure variable in the origianl and trimmed datasets 
##
##
## For more details about this modeling approach, please refer to the accompanying paper Masoud, Szyszkowicz,...,Burnett et al 2016 Air Quality, Atmosphere and Health
##
## Should there be any question with this program, please contact:
## 	Drs Hong Chen (Hong.Chen@oahpp.ca) and Rick Burnett (Rick.Burnett@hc-sc.gc.ca)
## 
##
## Please read the following notes before running the program
##
## 1. All categorical variables should be defined as factor variables before calling bestcox().
##	To do this, you may use, for example, inputdata$education <- as.factor(inputdata$education).
##
## 2. When calling bestcox(), you need to define the parameters based on the column names of your input dataset.
##
## 3. Your input dataset should be a dataframe.
##
## 4. By default, the plot will show only nonlinear curve based on best fitting model. This can be done either by leaving out
##	best_model parameter in the bestcox() or explicitly define best_model=1. 
##
##	Set "best_model"=2, the program will display the program will display concentration-response curve from ensemble models based on all models examined.
##
##    Set "best_model"=3, the program will display two concentration-response curves by overlaying best fitting model and ensemble models based on best 3 models examined.
##
##	Set "best_model"=4, the program will display two concentration-response curves by overlaying best fitting model and ensemble models based on ALL models examined..
##
## 5. By default, the bestcox() will translate air pollution exposure variable such that it has a minimum of unity.
##    
##	If users do not wish to translate data, you need to specify translate=F or FALSE
##
## 6. By default, output_dir is set as R's default working directory
##
##	If users wish to specify a different output directory, you need to ensure to have read/write permission 
##
## 7. By default, the option "lowest" is set to be "minimum". This means that the concentration-response curve covers the full range of exposure variable. 
##
##    Set "lowest"=1, the program will display concentration-response curve from 1 to max concentration of exposure variable.
##
## 8. By default, the option "x_min" is set to be 1. This means that the concentration-response plot has x axis with min=1. 
##
##    8.1. Set "x_min" as any positive numeric, x axis will start from this user-defined value.
##
##    8.2. Set "x_min" as any NON-numeric character including blank (eg, x_min=, x_min=NO, x_min=Null), x axis will start from the minimum of exposure variable.
##
## 9. By default, the option "export" is set to be FALSE. This means that this program will NOT produce any matrix for HRs derived over the range of AP exposure. 
##
##    9.1. Set "export" as TRUE or T, this program will produce two csv files: 
##			(1) a table "export_hr.csv" presenting HRs derived from ~1000 curves at a series of locations: 0% to 100% by 1% over the range of AP exposure, and
##			(2) a table "ap_exp.csv" presenting concentrations at 0% to 100% by 1% over the range of AP exposure.
##
##
## To illustrate how to use this routine, 8 examples are given below:
##
## Example 1: (fit a time-fixed cox model, without trimming and translating air pollution variable, without strata variable, c-r curve from min to max, and overlay best fitting model with an ensemble model based on 3 best models)
## out <- bestcox(data = inputdata, translate=F, timename = "time", casename = "status", lowperc = 0, upperc = 100, 
##        expo_name = "no2", cova_name = c("age","sex"), expo_unit = "ppb", best_model=3, output_dir="F:\\your output directory\\")
## out
##
## Example 2: (fit a time-fixed cox model, with trimmed and translated air pollution variable, without strata variable, c-r curve from min to max, and overlay best fitting model with an ensemble model based on 3 best models)
## out <- bestcox(data = inputdata, translate=T, timename = "time", casename = "death", lowperc = 1, upperc = 99, 
##      expo_name = c("pm25"), cova_name = c("age", "sex"), strata_name= NA, expo_unit = c("ug/m3"), best_model=3, output_dir="F:\\your output directory\\")
## out
##
## Example 3: (fit a time-fixed cox model, with trimmed and translated air pollution variable, with strata variable, and overlay best fitting model with an ensemble model based on ALL models examined)
## out <- bestcox(data = inputdata, timename = "time", casename = "status", lowperc = 1, upperc = 99, 
##        expo_name = c("no2"), cova_name = c("age"), strata_name= c("sex", "inst"), expo_unit = c("ppb"), best_model=4, output_dir="F:\\your output directory\\")
## out
##
## Example 4: (fit a time-varying cox model, with trimmed and translated air pollution variable, without strata variable, c-r curve from min to max, and only display best-fitting model)
## out <- bestcox(data = inputdata, start = "tstart", end="tstop", casename="status", lowperc = 1, upperc = 99,
##	  expo_name="no2", cova_name = c("sex", "age", "income", "bmi"), expo_unit = c("ppb"), best_model=1, output_dir="F:\\your output directory\\")
## out
##
## Example 5: (fit a time-varying cox model, with trimmed and translated air pollution variable, without strata variable, c-r curve from min to max, and only display an ensemble model based on ALL models examined)
##
## out <- bestcox(data = inputdata, start = c("T1"), end=c("T2"), casename = "Status", lowperc = 1, upperc = 99, 
##      expo_name = c("pm25"), cova_name = c("bmi"), strata_name= NA, expo_unit = c("ug/m3"), best_model=2, output_dir="F:\\your output directory\\")
## out
##
## Example 6: (fit a time-varying cox model, with trimmed and translated air pollution variable, without strata variable, c-r curve from 1 to max, and only display an ensemble model based on ALL models examined)
##
## out <- bestcox(data = inputdata, start = c("T1"), end=c("T2"), casename = "Status", lowperc = 1, upperc = 99, lowest = 1,
##      expo_name = c("pm25"), cova_name = c("bmi"), strata_name= NA, expo_unit = c("ug/m3"), best_model=2, output_dir="F:\\your output directory\\")
## out
##
## Example 7: (fit a time-varying cox model, with trimmed and translated air pollution variable, without strata variable, c-r curve from 1 to max, x axis starts at 10 unit, and only display an ensemble model based on ALL models examined)
##
## out <- bestcox(data = inputdata, start = c("T1"), end=c("T2"), casename = "Status", lowperc = 1, upperc = 99, lowest = 1,
##      expo_name = c("pm25"), cova_name = c("bmi"), strata_name= NA, expo_unit = c("ug/m3"), x_min=10, best_model=2, output_dir="F:\\your output directory\\")
## out
##
##
## Example 8: (fit a time-varying cox model, with trimmed and translated air pollution variable, without strata variable, c-r curve from 1 to max, export two csv files, x axis starts at 10 unit, and only display an ensemble model based on ALL models examined)
##
## out <- bestcox(data = inputdata, start = c("T1"), end=c("T2"), casename = "Status", lowperc = 1, upperc = 99, lowest = 1,
##      expo_name = c("pm25"), cova_name = c("bmi"), strata_name= NA, expo_unit = c("ug/m3"), x_min=10, best_model=2, export=T, output_dir="F:\\your output directory\\")
## out
########################################################################



#########################################################
## function to search for the optimal nonlinear model
#########################################################

library(survival)
library(MASS)

bestcox <- function (data,timename,translate=TRUE,start=NA,end=NA,casename,lowperc,upperc,lowest="MINIMUM",expo_name,cova_name,strata_name=NA,expo_unit,x_min=1,best_model=1,output_dir=NA, export=F){

  data$ap <- data[,expo_name]
  #data$ap[data$ap<1] <- 1						# convert any conc < 1 to 1

  # check lowest option

  if (!toupper(lowest)=="MINIMUM" & !lowest==1){
    stop("The option lowest should be either MINIMUM or 1.")    
  }

  # Parse out user's choice

  if (!is.numeric(best_model)){
    stop("Please enter a numeric value for parameter best_model.")    
  }

  if (best_model==1){
    overlay <- NA								# best fitting model
    only_ensemble_model <- 0
  } else if (best_model==2) {
    overlay <- "all"							# only ensemble model based on all models examined 
    only_ensemble_model <- 1
  } else if (best_model==3) {
    overlay <- "best3"							# overlay best fitting model with ensemble model based on best 3 models
    only_ensemble_model <- 0
  } else {
    overlay <- "all"							# overlay best fitting model with ensemble model based on all models examined 
    only_ensemble_model <- 0
  }

  # Characterize the distribution of air pollution exposure variable in the original dataset

  size_original <- length(data$ap)
  summary_original <- summary(data$ap)

  # Trim data, data translation, and tau

  data_trim <- subset(data, ap <= quantile(data$ap,upperc/100,na.rm = T) & ap >= quantile(data$ap, lowperc/100, na.rm = T))
  rm(data)									# remove original dataset to save memory space

  if (translate){
    ap.min <- min(data_trim$ap, na.rm = T)			# translate AP
    tran <- 0
  } else {
    ap.min <- 0								# no translate 
    tran <- 0
  }

  data_trim$ap_trans <- data_trim$ap - ap.min + tran 		# translate z to have the min of 0
  data_trim$ap_trans_lm <- data_trim$ap - ap.min		# translate z to have the min of 0
  
  if (translate){
  	data_trim$ap_trans_log <- log(1+data_trim$ap_trans)
  }else{
  	data_trim$ap_trans_log <- log(1+data_trim$ap_trans)
  }
 
  set_tau <- NA
  if (length(set_tau)==0 || any(is.na(set_tau)) || set_tau=='') {
	# set_tau <- 0.1
    	tau.2   <- 0.1*(max(data_trim$ap_trans, na.rm = T) - min(data_trim$ap_trans, na.rm = T))
	# set_tau <- 0.2
	tau.2[2]<- 0.2*(max(data_trim$ap_trans, na.rm = T) - min(data_trim$ap_trans, na.rm = T))

	num_tau <- 2
  }else{
  	tau <- set_tau*(max(data_trim$ap_trans, na.rm = T) - min(data_trim$ap_trans, na.rm = T))
	num_tau <- 1
  }
  	
  step_history <- NULL
  
  # Create a function to extract LL from Cox PH model

  coxmodel <- function(funcform, loca_perc, method = "breslow"){
    
    # Create f(z)*Wt

    if (loca_perc < 0){
      loca_para <- min(data_trim$ap_trans, na.rm = T) + (min(data_trim$ap_trans, na.rm = T) - quantile(data_trim$ap_trans, abs(loca_perc)/100, na.rm = T))
      loca_para_lm <- min(data_trim$ap_trans_lm, na.rm = T) + (min(data_trim$ap_trans_lm, na.rm = T) - quantile(data_trim$ap_trans_lm, abs(loca_perc)/100, na.rm = T))		# translate z to have the min of 0
    }else{
      loca_para <- quantile(data_trim$ap_trans, loca_perc/100, na.rm = T)
      loca_para_lm <- quantile(data_trim$ap_trans_lm, loca_perc/100, na.rm = T)		# translate z to have the min of 0
    }

    logit_w <- 1/(1+exp(-(data_trim$ap_trans - loca_para)/tau))
    logit_w_lm <- 1/(1+exp(-(data_trim$ap_trans_lm - loca_para_lm)/tau))			# translate z to have the min of 0	
    capture_mu <- NA												# capture mu

    if (funcform == "linear"){
      data_trim$expo <- logit_w_lm * data_trim$ap_trans_lm						# translate z to have the min of 0	
	capture_mu <- loca_para_lm										# capture mu
    }
    if (funcform == "log"){
      data_trim$expo <- logit_w * data_trim$ap_trans_log
	capture_mu <- loca_para											# capture mu
    }
    if (funcform == "pure.linear"){
      data_trim$expo <- data_trim$ap_trans_lm								# translate z to have the min of 0
	capture_mu <- NA												# capture mu
    }
    if (funcform == "pure.log"){
      data_trim$expo <- data_trim$ap_trans_log
	capture_mu <- NA												# capture mu
    }

    rm(logit_w)													# remove logit_w to save memory space
    rm(logit_w_lm)												# remove logit_w_lm to save memory space
    
    # Create coxph formula

    ## without strata variable
    if (length(strata_name)==0 || any(is.na(strata_name)) || strata_name=='') {
	   if (length(start)==0 || any(is.na(start)) || start=='') {
	   	coxformula <- paste("Surv(",timename,",",casename,")~expo+",paste(cova_name,collapse="+"),sep="")
	   } else {
	   	coxformula <- paste("Surv(",start,",",end,",",casename,")~expo+",paste(cova_name,collapse="+"),sep="")
	   }
    ## with strata variable
    } else {
	   if (length(start)==0 || any(is.na(start)) || start=='') {
	   	coxformula <- paste("Surv(",timename,",",casename,")~expo+",paste(cova_name,collapse="+"),"+","strata(",paste(strata_name, collapse=","),")",sep="")
	   } else {
	   	coxformula <- paste("Surv(",start,",",end,",",casename,")~expo+",paste(cova_name,collapse="+"),"+","strata(",paste(strata_name, collapse=","),")",sep="")
	   }
    }

    # Call coxph

    coxfit <- coxph(as.formula(coxformula), data = data_trim, method = method)
    
    est <- summary(coxfit)$coefficients[c("expo"),c("coef","se(coef)")]
    result <- data.frame(coef = est[1], se.coef=est[2], LL = coxfit$loglik[2], mu=capture_mu)	# capture mu

    rm(coxfit)														# remove coxfit to save memory space

    return(result)

  } # end of Cox PH model
  
  # Transformation Form Selection (8 or 16 Models, depending if users specify tau)

  if (length(set_tau)==0 || any(is.na(set_tau)) || set_tau=='') {

	tau <- tau.2[1]	# tau=0.1*range
  	step_a <- data.frame(funcform = c(rep("linear",4),rep("log",4)), loca_perc = rep(c(0,25,50,75),2), 
            coef = NA, se.coef = NA, LL = NA, mu=NA)
  	for (i in 1:nrow(step_a)){
    		step_a[i,c("coef","se.coef","LL", "mu")] <- coxmodel(step_a$funcform[i], step_a$loca_perc[i])
	}
	step_a$tau <- tau
	LL.step_a  <- step_a[step_a$LL == max(step_a$LL,na.rm = T),]$LL

	tau <- tau.2[2]	# tau=0.2*range
  	step_b <- data.frame(funcform = c(rep("linear",4),rep("log",4)), loca_perc = rep(c(0,25,50,75),2), 
            coef = NA, se.coef = NA, LL = NA, mu=NA)
  	for (i in 1:nrow(step_b)){
    		step_b[i,c("coef","se.coef","LL", "mu")] <- coxmodel(step_b$funcform[i], step_b$loca_perc[i])
	}
	step_b$tau <- tau
	LL.step_b  <- step_b[step_b$LL == max(step_b$LL,na.rm = T),]$LL

	if (LL.step_a > LL.step_b){
		step_0 <- step_a

		rejected.models <- step_b					# retain rejected tau and related model outputs
		rejected.tau <- step_b[1,]$tau				# retain rejected tau and related model outputs
		
	}else{
		step_0 <- step_b

		rejected.models <- step_a					# retain rejected tau and related model outputs
		rejected.tau <- step_a[1,]$tau				# retain rejected tau and related model outputs
	}

  }else{

  	step_0 <- data.frame(funcform = c(rep("linear",4),rep("log",4)), loca_perc = rep(c(0,25,50,75),2), 
                       coef = NA, se.coef = NA, LL = NA, mu=NA)
  	for (i in 1:nrow(step_0)){
    		step_0[i,c("coef","se.coef","LL", "mu")] <- coxmodel(step_0$funcform[i], step_0$loca_perc[i])
  	}
	step_0$tau <- tau	
  }

  step_history <- step_0
  step_history <- subset(step_history, select=-c(tau)) 							# drop tau
  
  step_0 <- step_0[step_0$LL == max(step_0$LL,na.rm = T),]
  funcform <- as.character(step_0$funcform)

  tau <- step_0$tau													# define best tau
  set_tau <- tau/(max(data_trim$ap_trans, na.rm = T) - min(data_trim$ap_trans, na.rm = T))		# best set_tau (0.1 or 0.2)
  step_0 <- subset(step_0, select=-c(tau)) 									# drop tau  

  set_tau_rejected <- rejected.tau/(max(data_trim$ap_trans, na.rm = T) - min(data_trim$ap_trans, na.rm = T))	# rejected set_tau

  # Transformation Location Parameter Selection

  #para0 <- step_0$loca_perc
  #
  #if(para0 == 25){
  #  step_temp <- cbind(funcform, loca_perc = 0 ,as.vector(coxmodel(funcform, 0)))
  #  step_0 <- rbind(step_0, step_temp)
  #  step_history <- rbind(step_history, step_temp)
  #}else if(para0 == 50){
  #  step_temp <- cbind(funcform, loca_perc = 75 ,as.vector(coxmodel(funcform, 75)))
  #  step_0 <- rbind(step_0, step_temp)
  #  step_history <- rbind(step_history, step_temp)
  #}
  
  #step_0 <- step_0[step_0$LL == max(step_0$LL,na.rm = T),]

  para0 <- as.numeric(step_0$loca_perc)
  
  step_temp <- rbind(cbind(funcform, loca_perc = (para0 + 5) ,as.vector(coxmodel(funcform, para0 + 5))),
                     cbind(funcform, loca_perc = (para0 - 5) ,as.vector(coxmodel(funcform, para0 - 5))))
  step_0 <- rbind(step_0, step_temp)
  step_history <- rbind(step_history, step_temp)
  
  step_1 <- step_0[step_0$LL == max(step_0$LL,na.rm = T),]
  para1 <- as.numeric(step_1$loca_perc)
  
  # STOP if reaching 15% below mu=1 or LL is not smaller	

  ##while(!(para1 %in% c(-5,100,para0))){  
  while(!(para1 %in% c(-15,100,para0))){  
    para_temp <- para1 + (para1 - para0)
    para0 <- para1
    para1 <- para_temp
    step_temp <- cbind(funcform, loca_perc = para1 ,as.vector(coxmodel(funcform, para1)))
    step_0 <- rbind(step_1, step_temp)
    step_history <- rbind(step_history, step_temp)
    step_1 <- step_0[step_0$LL == max(step_0$LL,na.rm = T),]
    para1 <- as.numeric(step_1$loca_perc)
  }
  
  # Reformat step_history and add back rejected model outputs 

  rownames(step_history) <- NULL
  step_history$tau <- set_tau			# define best tau (0.1 or 0.2)
  
  step_history_1to8 <- step_history[1:8,]
  step_history_9tolast <- step_history[9:length(step_history[,1]),]

  rejected.models$tau <- set_tau_rejected 	# define set_tau_rejected (0.2 or 0.1)

  step_history_1to8 <- rbind(step_history_1to8, rejected.models)
  step_history <- rbind(step_history_1to8, step_history_9tolast)		# re populate "step_history"

  # Add ensemble weights for 3 models around the best fit, ie., based on best mu with + and - 5th percentile 

  step_history_sort3 <- step_history
  rownames(step_history_sort3) <- NULL		# re start row number

  bestLL <- max(step_history_sort3$LL)
  step_history_sort3$best3 <- ifelse(step_history_sort3$LL==bestLL, 1, 0) 
  step_history_sort3$iteration <- rownames(step_history_sort3)
  
  best.models <- step_history_sort3[step_history_sort3$best3==1,]
  best.models.sort <-  best.models[order(best.models$iteration),]
  best.final.LL.iteration <- best.models.sort[length(best.models.sort[,1]),]$iteration
  final.mu <- best.models.sort[length(best.models.sort[,1]),]$loca_perc
  best.model.form <- best.models.sort[length(best.models.sort[,1]),]$funcform

  if (sum(step_history_sort3$best3)==2) {
	if (length(step_history_sort3[,1]) > best.final.LL.iteration){
		step_history_sort3[length(step_history_sort3[,1]),]$best3 <- 3
	} else if (length(step_history_sort3[,1]) == best.final.LL.iteration) {
		sec.final.mu <- best.models.sort[(length(best.models.sort[,1])-1),]$loca_perc
		if (final.mu > sec.final.mu) {
			step_history_sort3[step_history_sort3$loca_perc==(final.mu-10) & step_history_sort3$funcform==best.model.form,]$best3 <- 3
		} else {
			step_history_sort3[step_history_sort3$loca_perc==(final.mu+10) & step_history_sort3$funcform==best.model.form,]$best3 <- 3		
		}
	}
  } else if (sum(step_history_sort3$best3)==1) {
	##if (final.mu==-5) {
	##	step_history_sort3[step_history_sort3$loca_perc==0 & step_history_sort3$funcform==best.model.form,]$best3 <- 2
	##	step_history_sort3[step_history_sort3$loca_perc==5 & step_history_sort3$funcform==best.model.form,]$best3 <- 3
	if (final.mu==-15) {
		step_history_sort3[step_history_sort3$loca_perc==-10 & step_history_sort3$funcform==best.model.form,]$best3 <- 2
		step_history_sort3[step_history_sort3$loca_perc==-5 & step_history_sort3$funcform==best.model.form,]$best3 <- 3
	} else if (final.mu==100) {
		step_history_sort3[step_history_sort3$loca_perc==95 & step_history_sort3$funcform==best.model.form,]$best3 <- 2
		step_history_sort3[step_history_sort3$loca_perc==90 & step_history_sort3$funcform==best.model.form,]$best3 <- 3
	} else {
		step_history_sort3[step_history_sort3$loca_perc==(final.mu-5) & step_history_sort3$funcform==best.model.form,]$best3 <- 2
		step_history_sort3[step_history_sort3$loca_perc==(final.mu+5) & step_history_sort3$funcform==best.model.form,]$best3 <- 3
	}
  }

  # in rare occasion, models with set_tau_rejected may be assigned best3=2 or 3, thus need to be assigned to 0
  step_history_sort3[step_history_sort3$tau==set_tau_rejected,]$best3 <-0  

  nn <- length(step_history_sort3[step_history_sort3$best3>=1,]$LL)	# number of best fit models, normally this should be 3
  step_history_sort3$LL.diff <- NA
  step_history_sort3[step_history_sort3$best3>=1,]$LL.diff <- exp(step_history_sort3[step_history_sort3$best3>=1,]$LL - max(step_history_sort3[step_history_sort3$best3>=1,]$LL)) 
  step_history_sort3$wt.final3 <- NA
  step_history_sort3[step_history_sort3$best3>=1,]$wt.final3 <- step_history_sort3[step_history_sort3$best3>=1,]$LL.diff / sum(step_history_sort3[step_history_sort3$best3>=1,]$LL.diff, na.rm = T)
  step_history_sort3<-subset(step_history_sort3, select =-c(LL.diff, best3, iteration))

  if (nn >= 1) {
    step_history_sort <- step_history_sort3
  } else {
    step_history_sort <- step_history
    step_history_sort$wt <- NA
    step_history_sort$wt.final3 <- NA
  }

  # Add ensemble weights for all the models

  step_history_sort$wt <- exp(step_history_sort$LL-max(step_history_sort$LL))/sum(exp(step_history_sort$LL-max(step_history_sort$LL)))
  rownames(step_history_sort) <- NULL

  # Plot 

  if (nn >= 1) {
    finalmodels.best.nn <- subset(step_history_sort, !is.na(wt.final3)) 			## by default, limit to the best nn=3 models
    finalmodels.best.nn.final <- finalmodels.best.nn[order(finalmodels.best.nn$LL),]

    ##if overlay==ALL then ensemble model would include all models examined
    if (!is.na(overlay) & toupper(as.character(overlay))=='ALL'){
      nn <- length(step_history_sort[,1])
	step_history_sort.nn <- step_history_sort
	step_history_sort.nn$wt.final3 <- step_history_sort$wt
 	finalmodels.best.nn.final <- step_history_sort.nn[order(step_history_sort.nn$LL),]
    }

    ## result <- try(plot.bestmodel(ap_data=data_trim$ap, finalmodels=finalmodels.best.nn.final, expo_name=expo_name, unit=expo_unit, nn=nn, overlay=overlay, set_tau=set_tau, set_tau_reject=set_tau_rejected, tran=tran, translate=translate))
    ## result <- try(plot.bestmodel(ap_data=data_trim$ap, finalmodels=finalmodels.best.nn.final, expo_name=expo_name, unit=expo_unit, nn=nn, overlay=overlay, ensemblemodel=only_ensemble_model, set_tau=set_tau, set_tau_reject=set_tau_rejected, tran=tran, translate=translate))
    result <- try(plot.bestmodel(ap_data=data_trim$ap, finalmodels=finalmodels.best.nn.final, expo_name=expo_name, 
		unit=expo_unit, nn=nn, overlay=overlay, ensemblemodel=only_ensemble_model, set_tau=set_tau, set_tau_reject=set_tau_rejected, 
		tran=tran, translate=translate, lowest2=lowest, x_min=x_min, export=export, output_dir=output_dir))

    if (class(result)=="try-error") {
	print("Unable to find joint model! Plot is suppressed and only summary table is produced.")
	next
    }
 
  }else{
	print("Unable to find optimal model! Plot is suppressed and only summary table is produced.")
  }
 
  # Pure linear z model: exp(beta*z) (note that loca_perc = 75 is only a place holder)

  step_pure_linear <- cbind(funcform="pure.linear", loca_perc = 75 ,as.vector(coxmodel(funcform="pure.linear", 75)))
  step_pure_linear$loca_perc <- NA
  step_pure_linear$tau <- NA			# NA for tau			
  step_pure_linear$wt.final3 <- NA
  step_pure_linear$wt <- NA
  rownames(step_pure_linear) <- NULL

  # Pure log(z) model: exp(beta*log(z)) (note that loca_perc = 75 is only a place holder)

  step_pure_log <- cbind(funcform="pure.log", loca_perc = 75 ,as.vector(coxmodel(funcform="pure.log", 75)))
  step_pure_log$loca_perc <- NA
  step_pure_log$tau <- NA			# NA for tau
  step_pure_log$wt.final3 <- NA
  step_pure_log$wt <- NA
  rownames(step_pure_log) <- NULL

  # Output summary table from each step including coef, std, and wt  

  step_history_sort <- rbind(step_history_sort,step_pure_linear) 
  step_history_sort <- rbind(step_history_sort,step_pure_log) 
  ##colnames(step_history_sort)[colnames(step_history_sort)=="loca_perc"] <- "location"
  colnames(step_history_sort)[2] <- "location"
  colnames(step_history_sort)[4] <- "se"
  colnames(step_history_sort)[8] <- "finalwt"

  # re-order summary table

  step_history_sort_final <- step_history_sort[,c("funcform","location","mu","tau","coef","se","LL","wt","finalwt")]
  step_history_sort_final$LL <- format(round(step_history_sort_final$LL, 5), nsmall = 6)	## show 6 decimals of LL
  #return(step_history_sort_final)

  #if (num_tau == 1) {
  #	step_history_sort_final$tau <- set_tau
  #}else{
  	#step_history_sort_final$tau <- NA
	#step_history_sort_final[1:8,]$tau <- 0.1
	#step_history_sort_final[9:16,]$tau <- 0.2
	#step_history_sort_final[17:length(step_history_sort_final[,1])-2,]$tau <- set_tau
  #	step_history_sort_final$tau <- set_tau
  #}

#  write.table(step_history_sort_final, file="search.results.csv", sep = ",", col.names = NA)
  write.table(step_history_sort_final, file=paste(output_dir, "search.results.csv", sep=""), sep = ",", col.names = NA)

  # add 2 summary tables showing the distribution of air pollution exposure variable in the original and trimmed datasets, respectively

  #size_original <- length(data$ap)
  size_trimmed <- length(data_trim$ap)
  #summary_original <- summary(data$ap)
  summary_trimmed <- summary(data_trim$ap)
  overall_summary <- list("count_of_obs_in_original_dataset"=size_original,
				  "distr_of_exp_in_original_dataset"=summary_original, 
				  "count_of_obs_in_trimmed_dataset"=size_trimmed,
				  "distr_of_exp_in_trimmed_dataset"=summary_trimmed, 
				  "summary_model_fitting"=step_history_sort_final)  

  # output all 3 summary tables: (1) summary of model fitting and (2) distr of exp var in both original and trimmed datasets

  return(overall_summary)

}



#########################################################
## function to plot nonlinear CR relationship
#########################################################

plot.bestmodel <- function(ap_data, finalmodels, expo_name, unit, nn, overlay, ensemblemodel, set_tau, set_tau_reject, tran, translate, lowest2, x_min, export, output_dir){

  if (toupper(lowest2)=="MINIMUM") {
  	x <- seq(tran, max(ap_data)-min(ap_data)+tran, 0.1) 			# based on range of trimmed and translated ap_data
  	nx <- length(x)

  	if (finalmodels[length(finalmodels[,1]),]$funcform == "linear") { 
  		x <- x-tran									# translate z to have the min of 0
  	}

  	if (translate) { 
		x <- x 
		x.bk <- x									# save a backup copy
  	} else {
		x <- x									# if translate=FALSE, use full range of original AP data, but scale to 0:(max-min+0)
		x.bk <- seq(min(ap_data), max(ap_data), 0.1)			# if translate=FALSE, use full range of original AP data
  		nx <- length(x)
 	 }
  } else {
  	x <- seq(1, max(ap_data), 0.1) 						# from 1 to max of ap_data
	nx <- length(x)
	x.bk <- x
  }	

  # prepare x1 for use in simulation, varying depending on perc, set_tau, and funcform

  sim.x1 <- function(x_sim, perc, set_tau_sim, funcform_sim){

  	if (perc < 0){
    		mu <- min(x_sim, na.rm = T) + (min(x_sim, na.rm = T)-quantile(x_sim, abs(perc)/100, na.rm = T))
  	}else{
    		mu <- quantile(x_sim, perc/100, na.rm = T)
  	}

  	tau_sim <- set_tau_sim*(max(x_sim)-min(x_sim))

  	logit <- exp((x_sim-mu)/tau_sim)/(1+exp((x_sim-mu)/tau_sim))

   	if (funcform_sim== "linear"){
		x_sim <- x_sim-min(x_sim)		# set min(x_sim) to 0 for linear model
    		x1<-x_sim*logit
  	}
  	if (funcform_sim == "log"){
		if (translate){
    			x1<-log(1+x_sim)*logit
  		}else{
    			x1<-log(1+x_sim)*logit
		}
  	}
  
 	# x1 <- x1-min(x1)		# set min(x1) as reference, thus beta*min(x1)=0

   	return(x1)
   }

  # Consider only optimal model - simulate 1000 realizations based on se.coef alone

  nsim<-1000
  ran<-matrix(0, nsim, 1)
  rr<-matrix(0, nsim, nx)
  medRR<-matrix(0, nx, 1)
  upcl<-matrix(0, nx, 1)
  lowcl<-matrix(0, nx, 1)

  loca_perc_sim <- finalmodels[finalmodels$LL==max(finalmodels$LL, na.rm=T),]$loca_perc
  funcform_sim <- finalmodels[length(finalmodels[,1]),]$funcform
  x1 <- sim.x1(x_sim=x, perc=loca_perc_sim, set_tau_sim=set_tau, funcform_sim=funcform_sim)

  for (i in 1:nsim) {
    ran[i,]<-rnorm(1, finalmodels[length(finalmodels[,1]),]$coef, finalmodels[length(finalmodels[,1]),]$se.coef)
    for (j in 1:length(x)) {
    rr[i,j]<-exp(ran[i,1]*x1[j])
    }
  }

  for (j in 1:length(x)) {
    medRR[j] <- mean(rr[,j])
    lowcl[j] <- quantile(rr[,j], 0.025)
    upcl[j] <- quantile(rr[,j], 0.975)
  }

  # Add 3 best models or ALL models - simulate 1000 realizations based on se.coef AND weights derived from LL

  if (nn >= 2) {

    nsim<-1000
    ran.3<-matrix(0, nsim, 1)
    rr.3<-matrix(0, nsim, nx)
    medRR.3<-matrix(0, nx, 1)
    upcl.3<-matrix(0, nx, 1)
    lowcl.3<-matrix(0, nx, 1)

    pp <- 1			# position variable in the 1000 sim
    nn <- nn			# consider top 3 models
    nsim.sum <- 0		# count of nsim to ensure the last run lead to rownum of exactly 1000

    # k from 0 to nn-1: varying depending on models included and pooled

    for (k in 0:(nn-1)) {

      nsim.wt <- nsim * round(finalmodels[length(finalmodels[,1])-k,]$wt.final3, digits = 3)

  	loca_perc_sim <- finalmodels[length(finalmodels[,1])-k,]$loca_perc
  	funcform_sim <- finalmodels[length(finalmodels[,1])-k,]$funcform
      if (finalmodels[length(finalmodels[,1])-k,]$tau==set_tau_reject){ 
		x1 <- sim.x1(x_sim=x, perc=loca_perc_sim, set_tau_sim=set_tau_reject, funcform_sim=funcform_sim)
      	}else{
  		x1 <- sim.x1(x_sim=x, perc=loca_perc_sim, set_tau_sim=set_tau, funcform_sim=funcform_sim)
		}

      if (k==nn-1) {nsim.wt <- nsim - nsim.sum}
      for (i in pp:(pp+nsim.wt-1)) {
		# for models with weights~0, i may exceed 1000 thus throw an error on out of bound
		if (i<=1000){		
      		ran.3[i,]<-rnorm(1, finalmodels[length(finalmodels[,1])-k,]$coef, finalmodels[length(finalmodels[,1])-k,]$se.coef)
      		for (j in 1:length(x)){
      			rr.3[i,j]<-exp(ran.3[i,1]*x1[j])
      		}
		}
      }
      pp <- pp+nsim.wt
      nsim.sum <- nsim.sum+nsim.wt

    }

    for (j in 1:length(x)) {
      medRR.3[j] <- mean(rr.3[,j])
      lowcl.3[j] <- quantile(rr.3[,j], 0.025)
      upcl.3[j] <- quantile(rr.3[,j], 0.975)
    }
  }

  # transform x back to original scale of AP data for plotting

  if (translate) { 
    if (finalmodels[length(finalmodels[,1]),]$funcform == "linear") { 
      x <- x+min(ap_data)							# translate z to have the min of 0
    } else {
	x <- x+min(ap_data)-tran 						# translate z to have the min of 1
    }
  }else{
	x <- x.bk									# scale back x to have the min of 0
  }

  # export two csv files: one rr matrix and the other ap_exp at a series locations at 0% to 100% by 1% over ap exp

  if (export) { 

	## export first csv file

	ncol.rr.3 <- length(rr.3[1,])
	rank.rr.3 <- rep(1:ncol.rr.3,1)

	sel.rank.rr.3 <- NA
	for (p in 0:100){
		sel.rank.rr.3 <- rbind(sel.rank.rr.3, quantile(rank.rr.3, p/100 ,na.rm = T))
	}
	sel.rank.rr.3 <- na.omit(sel.rank.rr.3)
	sel.rank.rr.3 <- round(sel.rank.rr.3)
	sel.rank.rr.3 <- sel.rank.rr.3[1:length(sel.rank.rr.3)]

	export.rr <- matrix(0, length(rr.3[,1]), length(sel.rank.rr.3))

	for (pp in 1:length(sel.rank.rr.3)){
		export.rr[,pp] <- rr.3[,sel.rank.rr.3[pp]] 
	}

  	write.table(export.rr, file=paste(output_dir, "export.hr.csv", sep=""), sep = ",", col.names = NA)

	## export second csv file
	
	x.ap <- data.frame(x)
	x.ap$id <- 1:nrow(x.ap)

	rank.dataframe <- data.frame(id=sel.rank.rr.3)

	export.x <- x.ap[which( x.ap$id %in% rank.dataframe$id ),]
	export.x <- subset(export.x, select=-c(id))
	names(export.x)[names(export.x) == 'x'] <- 'air.pollutant'
 	
  	write.table(export.x, file=paste(output_dir, "ap.exp.csv", sep=""), sep = ",", col.names = NA)

  }

  # whether or not overlay ensemble model 

  if (length(overlay)==0 || any(is.na(overlay)) || overlay==''){

    # plot optimal model only
 
    par(las = 1 , cex = 1)

    # use predicted values from optimal model for the plot

    medRR.3 <-  medRR
    lowcl.3 <- lowcl
    upcl.3 <- upcl

    if (x_min>=0){ 	
    	plot(x, upcl.3, lwd=4, type="l",  col="#DEEBF7", frame.plot=T, xlim=c(x_min, max(x)), ylim=c(min(lowcl.3, na.rm=TRUE)-0.25, max(upcl.3, na.rm=TRUE)+0.25), ylab="Hazard Ratio", xlab=paste(toupper(expo_name), " (", unit, ")"))
    } else {
    	plot(x, upcl.3, lwd=4, type="l",  col="#DEEBF7", frame.plot=T, ylim=c(min(lowcl.3, na.rm=TRUE)-0.25, max(upcl.3, na.rm=TRUE)+0.25), ylab="Hazard Ratio", xlab=paste(toupper(expo_name), " (", unit, ")"))
    }

    polygon(x=c(x,  rev(x)), y=c(lowcl.3, rev(upcl.3)), col="#DEEBF7", border=NA, lty=2)
    lines(x, medRR.3, lwd=3, col="#08519C")
    lines(x, lowcl.3, lwd=2,  col="#DEEBF7")  

    nn <- 1000

  } else {

    # plot joint model

    par(las = 1 , cex = 1)

    if (nn < 2) {
      # in case if only optimal model existed, use predicted values from the optimal model for the plot
      medRR.3 <-  medRR
      lowcl.3 <- lowcl
      upcl.3 <- upcl
    }

    # in rare occasion, max(upcl.3)==Inf and/or min(lowcl.3)==Inf

    upcl.3.excluded.inf <- upcl.3[upcl.3<100]
    max.ylim <- max(upcl.3.excluded.inf, na.rm=TRUE)+0.25
    if (is.na(max.ylim)){max.ylim <- 5}

    lowcl.3.excluded.inf <- lowcl.3[lowcl.3<100]
    min.ylim <- min(lowcl.3.excluded.inf, na.rm=TRUE)-0.25
    if (is.na(min.ylim)){min.ylim <- 0}

    if (x_min>=0){ 
    	plot(x, upcl.3, lwd=4, type="l",  col="#DEEBF7", frame.plot=T, xlim=c(x_min, max(x)), ylim=c(min.ylim, max.ylim), ylab="Hazard Ratio", xlab=paste(toupper(expo_name), " (", unit, ")"))
    } else {
    	plot(x, upcl.3, lwd=4, type="l",  col="#DEEBF7", frame.plot=T, ylim=c(min.ylim, max.ylim), ylab="Hazard Ratio", xlab=paste(toupper(expo_name), " (", unit, ")"))
    } 

    polygon(x=c(x,  rev(x)), y=c(lowcl.3, rev(upcl.3)), col="#DEEBF7", border=NA, lty=2)
    lines(x, medRR.3, lwd=3, col="#08519C")
    lines(x, lowcl.3, lwd=2,  col="#DEEBF7")  

    # Overlay optimal model

    if (ensemblemodel==0) {
	lines(x, medRR, lwd=2, col="red")
 	lines(x, lowcl, lwd=1, lty=2, col="red")
    	lines(x, upcl, lwd=1, lty=2, col="red")
    }

  }

  # Add rugs, legend, and reference line

  ## add rugs=ticks at datapoints
  ####axis(side = 1 , line = -1.2 , at = jitter(x) , labels = F , tick = T , tcl = 0.8 , lwd.ticks = 0.1 , lwd = 0)
  ## rugs and labels at 1Q, median and 3Q
  ####axis(side = 1 , line = -1.0 , at = fivenum(x)[2:4], lwd = 0 , tick = T, tcl = 1.2 , lwd.ticks = 1 , col.ticks = "black" , labels = c("Quartile 1","Median","Quartile 3"), cex.axis = 0.7, col.axis = "black" , padj = -2.8)
  ####axis(side = 1 , line = 0.0 , at = fivenum(x)[2:4], lwd = 0 , tick = T, tcl = 0.2 , lwd.ticks = 1 , col.ticks = "black", labels = FALSE)

  ## add legend and RR=1 line
  ####if (nn > 100) {
  ####legend("topleft", inset=c(0,0), c("optimal model"), col=c("blue"), lty=1, lwd=3)
  ####} else if (nn >= 2) {
  ####legend("topleft", inset=c(0,0), c("ensemble model", "optimal model"), col=c("blue", "red"), lty=1, lwd=3)
  ####} else {
  ####legend("topleft", inset=c(0,0), c("optimal model", "optimal model"), col=c("blue", "red"), lty=1, lwd=3)
  ####}
  abline(1,0, col = "gray", lty=3, lwd=1)
  box(bty = "n")

}
