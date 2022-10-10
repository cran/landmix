
#' Mixture porportion group identifier
#'
#' Computes the mixture proportion group identifier.
#'
#' @param n sample size, must be at least 1.
#' @param m number of different mixture proportions, must be at least 2.
#' @param p number of populations, must be at least 2.
#' @param qvs a numeric matrix of size \code{p} by \code{m} containing all possible
#' mixture proportions (i.e., the probability of belonging to each population k, k=1,...,p.).
#' @param q.use a numeric matrix of size \code{p} by \code{n} containing the
#' mixture proportions for each person in the sample.
#'
#' @section Details:
#' The matrix \code{qvs} contains all mixture proportions. The function
#' \code{compute.uset} assigns individual \eqn{i} to mixture proportion \eqn{u},
#' for \eqn{u=1,\ldots,m}, if \eqn{q_i=qvs_{u}}.
#'
#' @return A numeric vector of length \code{n} containing the mixture proportion
#' group to which each person belongs.
#'
#' @export
compute.uset <- function(n,m,p,qvs,q.use){
  uset <- rep(0,n)          ## indicator of which subgroup qvs
  for(jj in 1:n){
    for(ii in 1:m){
      q.tmp <- all(q.use[1:p,jj]==qvs[,ii])
      if(q.tmp==TRUE){
        uset[jj] <- ii
      }
    }
  }
  return(uset)
}


make.data.set <- function(
  n,m,p,qvs,q,
  x,delta,ww,zz){

  #####################
  ## Compute u-group ##
  #####################
  uset <- compute.uset(n,m,p,qvs,q)

  true.groups <- rep(0,length(delta))
  

  data <- data.frame(x=x,delta=delta,t(q),ww,zz,uset,true.groups)
  colnames(data) <- c("x","delta",paste("q",1:p,sep=""),"w","z","uset","group")

  ## sort in increasing order by x
  data <- data[order(data$x),]
  return(data)
}


#' Sample sizes of mixture proportion groups.
#'
#' Computes the number of individuals in each mixture proportion group.
#'
#' @param n sample size, must be at least 1.
#' @param m number of different mixture proportions, must be at least 2.
#' @param p number of populations, must be at least 2.
#' @param qvs a numeric matrix of size \code{p} by \code{m} containing all possible
#' mixture proportions (i.e., the probability of belonging to each population k, k=1,...,p.).
#' @param q.use a numeric matrix of size \code{p} by \code{n} containing the
#' mixture proportions for each person in the sample.
#'
#' @export
compute.r <- function(n,m,p,qvs,q.use){
  r <- rep(0,m)		    ## number in each qvs subgroup
  for(jj in 1:n){
    for(ii in 1:m){
      q.tmp <- all(q.use[1:p,jj]==qvs[,ii])
      if(q.tmp==TRUE){
        r[ii] <- r[ii] + 1
      }
    }
  }
  return(r)
}


get_method_label <- function(run.NPNA,
                             run.NPNA_avg){

 
  ## Names for methods run
  est.names <- NULL
    if(run.NPNA==TRUE){
    est.names <- c(est.names,"NPNA")
  }

  if(run.NPNA_avg==TRUE){
    est.names <- c(est.names,"NPNA_avg")
  }

  return(est.names)
}


####################################
####################################
##
##
## Functions for the main method
##
##
####################################
####################################

landmix.estimator.wrapper <- function(
  n,m,p,qvs,q,
  x,delta,ww,zz,
  run.NPNA,
  run.NPNA_avg,
  tval,tval0,
  z.use,w.use){


  ##############################
  ## setup output for storage ##
  ##############################
  est.names <- get_method_label(run.NPNA,
                                run.NPNA_avg)

  method.label <- est.names


  ###############
  ## Compute r ##
  ###############
  r <- compute.r(n,m,p,qvs,q)



  ###################
  ## Make data set ##
  ###################
  data <- make.data.set(
    n,m,p,qvs,q,
    x,delta,ww,zz)



  #######################
  ## main computations ##
  #######################
  problem <- NULL


  ####################
  ## run estimators ##
  ####################
  data.out <- estimator.main(data,
                             n,p,m,r,qvs,
							 tval,tval0,
                             method.label,
                             z.use,w.use)

  return(data.out)
}



landmix.estimator <- function(n,m,p,qvs,q,
				x,delta,ww,zz,
				run.NPNA,
				run.NPNA_avg,
				tval,tval0,
				z.use,w.use){

	############################
	## run the main procedure ##
	############################
	estimators.out <- landmix.estimator.wrapper(
		n,m,p,qvs,q,
		x,delta,ww,zz,
		run.NPNA,
		run.NPNA_avg,
		tval,tval0,
		z.use,w.use
	)

	## check if we have a problem with NPNA estimator.
	Ft.estimate <- estimators.out$Ft.store
	St.estimate <- estimators.out$Sout.store
	
	return(list(Ft.estimate=Ft.estimate,
		 St.estimate=St.estimate))
}


get.new.qs <- function(data,n,p,m,r,qvs,tval0,method.label,update.qs=FALSE){

  ############################
  ## separate results by t0 ##
  ############################
  new.output <- vector("list",length(tval0))
  names(new.output) <- paste("t0",tval0,sep="")

  #####################################
  ## initialize values for r, qvs, q ##
  #####################################
  r.use <- vector("list",length(method.label))
  names(r.use) <- method.label
  qvs.use <- r.use
  q.use <- r.use

  for(kk in 1:length(method.label)){
    r.use[[kk]] <- r
    qvs.use[[kk]] <- qvs
    q.use[[kk]] <- data[,paste("q",1:p,sep="")]
  }

  ######################
  ## update r, qvs, q ##
  ######################
      for(tt0 in 1:length(tval0)){
      new.output[[tt0]] <- list(qvs=qvs.use,q=q.use,r=r.use)
    }
  return(new.output)
}

## Produces arrays with NA components.
## The arrays will be used to store the results for the estimated
## distribution function and prediction accuracy measures (e.g., AUC, BS).
get.null.theta <- function(theta.names=c("estimator","prediction"),
                           first.label.name,
                           tval,tval0,z.use,w.use,Ft.name="Ft",p){


  null.theta <- get.empty.list(theta.names)

  null.theta[["estimator"]] <- array(NA,dim=c(length(first.label.name),
                                              length(tval),
                                              length(tval0),
                                              length(z.use),
                                              length(w.use),p),
                                     dimnames=list(
                                       method=first.label.name,
                                       time=paste("tt",tval,sep=""),
                                       time0=paste("tt0",tval0,sep=""),
                                       zz=z.use,
                                       ww=w.use,
                                       Ft=paste(Ft.name,1:p,sep="")))

  null.theta[["prediction"]] <- array(NA,dim=c(length(first.label.name),
                                               length(tval),
                                               length(tval0),2
  ),
  dimnames=list(
    method=first.label.name,
    time=paste("tt",tval,sep=""),
    time0=paste("tt0",tval0,sep=""),
    measure=c("AUC","BS")
  ))

  return(null.theta)
}

## add dimension to each array in a list of arrays.
add.dimension.null <- function(null.theta,location,label.dim,label.name){
  null.theta0 <- null.theta
  for(u in 1:length(null.theta)){
    if(!is.null(null.theta0[[u]])){
      if(location=="first"){
        new.names <- list(iters=list(label.name),dimnames(null.theta0[[u]]))

        null.theta[[u]] <- array(NA,dim=c(label.dim,dim(null.theta0[[u]])),
                                 dimnames=unlist(new.names,recursive=FALSE))
      } else if(location=="last"){
        new.names <- list(dimnames(null.theta0[[u]]),val=list(label.name))
        null.theta[[u]] <- array(NA,dim=c(dim(null.theta[[u]]),label.dim),
                                 dimnames=unlist(new.names,recursive=FALSE))
      }
    }
  }
  return(null.theta)
}

## function to make all (z,w) combinations
get.zw.combinations <- function(z.use,w.use){
  data.evaluate.list <- list(z=z.use,w=w.use)
  data.evaluate.orig <- expand.grid(data.evaluate.list)
  return(data.evaluate.orig)
}



## get storage matrices for Ft, Sout
# dim.col=p
tmp.storage <- function(data.tmp,dim.col,dim.names="Ft"){
  out <- matrix(NA,nrow=nrow(data.tmp),ncol=dim.col)
  out.name <- paste(dim.names,1:dim.col,sep="")
  out <- cbind(data.tmp,out)
  colnames(out) <- c(colnames(data.tmp),out.name)
  return(out)
}

#############################################################################
## Store data set that will be used to estimate the distribution functions ##
## and compute the prediction accuracy measures.                           ##
#############################################################################
get.tmp.storage <- function(method,
                            z.use,w.use,data,p,m){

  if(method=="NPNA_avg" ){
    ## storage for calculation
    data.zw <- data[,c("z","w")]
  } else{
    ## storage for calculation
    data.zw <- get.zw.combinations(z.use,w.use)
  }

  ## Ft and Sout have n rows
  Ft.tmp.out <- tmp.storage(data.tmp=data.zw,dim.col=p,dim.names="Ft")


  Sout.tmp.out <- tmp.storage(data.tmp=data.zw,dim.col=m,dim.names="St")

  ## remove duplicates
  data.zw <- data.zw[!duplicated(data.zw),]

  list(data.zw=data.zw,Ft.tmp.out=Ft.tmp.out,Sout.tmp.out=Sout.tmp.out)
}

## get all list of arrays used in the analysis.
all.null.theta <- function(theta.names=c("estimator","prediction"),
                           first.label.name,
                           tval,tval0,z.use,w.use,Ft.name,p,
                           label.dim.simus,label.name.simus){

  null.theta <- get.null.theta(theta.names,
                               first.label.name,
                               tval,tval0,z.use,w.use,Ft.name,p)
  ## remove method
  method.names <- which(names(dimnames(null.theta$estimator))=="method")
  null.theta.nomethod <- array(NA,dim=dim(null.theta$estimator)[-method.names],
                               dimnames=dimnames(null.theta$estimator)[-method.names])


  ## add simus layer
  null.theta.simus <- add.dimension.null(null.theta,location="first",
                                         label.dim=label.dim.simus,label.name=label.name.simus)

  ## ci results layer
  null.theta.ci <-  add.dimension.null(null.theta,location="last",
                                       label.dim=3,label.name=c("varest","varlo","varhi"))

  ## est, ci results
  ##null.theta.est.ci <-  add.dimension.null(null.theta,location="last",
  ##              label.dim=4,label.name=c("est","varest","varlo","varhi"))

  ## simus and ci results layer to null.theta
  ##null.theta.simus.ci <-  add.dimension.null(null.theta.simus,location="last",
  ##              label.dim=3,label.name=c("varest","varlo","varhi"))

  ## simus, est, ci results
  null.theta.simus.est.ci <- add.dimension.null(null.theta.simus,location="last",
                                                label.dim=4,label.name=c("est","varest","varlo","varhi"))

  list(null.theta=null.theta,null.theta.simus=null.theta.simus,
       ##null.theta.simus.ci=null.theta.simus.ci, ## not used
       null.theta.ci=null.theta.ci,
       null.theta.simus.est.ci=null.theta.simus.est.ci,
       ##null.theta.est.ci =  null.theta.est.ci, ## not used
       null.theta.nomethod=null.theta.nomethod )
}



estimator.main <- function(data,
                           n,p,m,r,qvs,
						   tval,tval0,
                           method.label,
                           z.use,w.use){

  #############################
  ## dummy variables needed  ##
  #############################
  problem <- NULL	## for NPNA estimator
  bootvar <- FALSE  	## for kincohort estimator

  ##################
  ## store output ##
  ##################

  ## for F_p(t) data
  nsimu <- 1
  Ft.null.theta <- all.null.theta(theta.names=c("estimator","prediction"),
                                  first.label.name=method.label,
                                  tval,tval0,z.use,w.use,Ft.name="Ft",p,
                                  label.dim.simus=nsimu,
                                  label.name.simus=paste("iters",1:nsimu,sep=""))

  ## for S_m(t) data
  Sout.null.theta <- all.null.theta(theta.names=c("estimator","prediction"),
                                    first.label.name=method.label,
                                    tval,tval0,z.use,w.use,Ft.name="St",p=m,
                                    label.dim.simus=nsimu,
                                    label.name.simus=paste("iters",1:nsimu,sep=""))

    ## Run estimator evaluated at specified (z,w) pairs

    ## storage for output
    Ft.store <- Ft.null.theta$null.theta$estimator
    Sout.store <- Sout.null.theta$null.theta$estimator

	cv_folds <- 1



  for(cvv in 1:cv_folds){

	  ###########################
	  ## split data into folds ##
	  ###########################
	  if(cv_folds>1){
		num_train <- round(0.6*n)
		train_ind <- sample(1:n, num_train, replace = FALSE)
		test_ind <- setdiff(1:n, train_ind)

		data_train <- data[train_ind,]
		data_test <- data[test_ind,]
	  } else {
		data_train <- data
		data_test <- data
	  }

	  ## sort in increasing order by x
	  data_train <- data_train[order(data_train$x),]
	  data_test <- data_test[order(data_test$x),]

	  n_use <- nrow(data_train)
  	  p_use <- p
	  q_use <- t(data_train[,c("q1","q2")])
	  qvs_use <- t(as.matrix(unique(t(q_use))))
	  tmp_info <- get_new_qvs(qvs_use,qvs)
	  qvs_use <- tmp_info$qvs_new
	  m_use<- tmp_info$m_new
	  r_use <- compute.r(n_use,m_use,p_use,qvs_use,q_use)

	  #################
	  ## update q's? ##
	  #################
	  get.info.use <- get.new.qs(data_train,n_use,p_use,m_use,r_use,qvs_use,tval0,
									method.label,update.qs=FALSE)


	  ##########################
	  ## estimation procedure ##
	  ##########################
	  for(tt0 in 1:length(tval0)){
		##cat("tt0=",tt0)


		for(kk in 1:length(method.label)){
		  ##cat("kk=",kk)


		  ########################
		  ## Index used for tt0 ##
		  ########################
			tt0.use <- tt0
		  

		  ####################################################
		  ## subset of data where subjects survived past t0 ##
		  ####################################################
		  #index_train_use <- which(data$x>tval0[tt0.use] & data$delta>0)## 8/24/2016: wrong
		  index_train_use <- which(data_train$x>tval0[tt0.use]) ## 8/24/2016: this is the correct one.
		  data_train_subset <- data_train[index_train_use,]

		  index_test_use <- which(data_test$x>tval0[tt0.use])
		  data_test_subset <- data_test[index_test_use,]

		  ############################
		  ## information used at t0 ##
		  ############################
		  info.tmp <- get.info.use[[tt0.use]]

		  ## qvs.use depends on method of Ft estimator
		  qvs.tmp <- info.tmp$qvs[[method.label[kk]]]

		  ## r info for weights
		  r.tmp <- info.tmp$r[[method.label[kk]]]

		  #####################
		  ## set up the data ##
		  #####################
		  tmp.out <- get.tmp.storage(method=method.label[kk],
										z.use,w.use,
										data_test_subset,p,m)

		  data_test_zw <- tmp.out$data.zw

		  Ft_test_out <- tmp.out$Ft.tmp.out
		  Sout_test_out <- tmp.out$Sout.tmp.out

		  ## information to run kin-cohort estimators. From data_train
		  q.tmp <- t(info.tmp$q[[method.label[kk]]])
		  n.tmp <- nrow(data_train)
		  m.tmp <- length(table(data_train$q1))
		  x.tmp <- data_train[,"x"]
		  delta.tmp <- data_train[,"delta"]


	

		  for(tt in 1:length(tval)){
			##cat("tt=",tt)

			###########################################
			## Valid to compute F(t|t0) when t> t0. ##
			###########################################
			if(tval[tt] > tval0[tt0]){

			  ####################
			  ## NPMLE estimator ##
			  ####################


				####################
				## NPNA estimator ##
				####################

				#####################################
				## what(z,w,u) data to evaluate at ##
				#####################################
				data_test_zw_evaluate <- get.zw.evaluate(data_test_zw,m)

				######################################
				## get data for landmark estimation ##
				######################################
				data.land <- data_train_subset[,c("x","delta","z","w","uset")]

				## run NPNA estimator at tt, and report values at data.evaluate
				S.NPNA.zw.out <- S.NPNA.zw(t=tval[tt],
										   data=data.land,
										   newdata=data_test_zw_evaluate)

								## report output
				NPNA.out <-  as.data.frame(S.NPNA.zw.out$data.out)

				## Form F(t),St.estimate
				for(ll in 1:nrow(Ft_test_out)){

				  ## get appropriate rows for (z,w) estimates
				  z.tmp <- Ft_test_out[ll,"z"]
				  w.tmp <- Ft_test_out[ll,"w"]

				  ## get appropriate subset of NPNA.out
				  NPNA.out.tmp <- NPNA.out[which(NPNA.out$w==w.tmp & NPNA.out$z==z.tmp),]
				  Sout.tmp <- NPNA.out.tmp$survival.v.new
				  Sout_test_out[ll,paste("St",1:m,sep="")] <- Sout.tmp

				  u.index <- NPNA.out.tmp$u
				  u.tmp <- t(qvs.tmp[,u.index])

				  Ft_test_out[ll,paste("Ft",1:p,sep="")] <-
					solve(t(u.tmp) %*% diag(r.tmp) %*% u.tmp) %*%
					t(u.tmp) %*%  diag(r.tmp) %*% (1-Sout.tmp)


				}

			 


				######################################
				## output results at specific (z,w) ##
				######################################
				if(method.label[kk]!="NPNA_avg"){
				  Ft.store[kk,tt,tt0,,,] <- unflatten.array(Ft.store[kk,tt,tt0,,,],
															dim.order=c("zz","ww","Ft"),
															Ft_test_out[,paste("Ft",1:p,sep="")],
															flatten.name="Ft")

				  Sout.store[kk,tt,tt0,,,] <- unflatten.array(Sout.store[kk,tt,tt0,,,],
															  dim.order=c("zz","ww","Ft"),
															  Sout_test_out[,paste("St",1:m,sep="")],
															  flatten.name="Ft")
				} else{
				  ## take average over all (z,w) combinations
				  Ft.tmp.avg <- apply(Ft_test_out[,paste("Ft",1:p,sep="")],2,mean)
				  Ft.store[kk,tt,tt0,,,] <- repeat.zw(Ft.tmp.avg,z.use,w.use,p)

				  Sout.tmp.avg <- apply(Sout_test_out[,paste("St",1:m,sep="")],2,mean)
				  Sout.store[kk,tt,tt0,,,] <- repeat.zw(Sout.tmp.avg,z.use,w.use,m)
				}
			  }
		  }
		}
	  }
  }
  list(Ft.store=Ft.store,Sout.store=Sout.store)
}


## function to repeat Ft entries
repeat.zw <- function(x,z.use,w.use,p){
  ## entries go in as p, w, z
  ## entries go out as z by w by p
  x.tmp <- array(x,dim=c(p,length(w.use),length(z.use)))
  x.tmp <- aperm(x.tmp,perm=c(3,2,1))
  return(x.tmp)
}

## Change a matrix to an array.
unflatten.array <- function(x,dim.order,data.change,
                            flatten.name="time"){

  ## get correct permutation order of array
  xorig <- x
  if(!is.null(dim.order)){
    x <- aperm(x,dim.order)
  }

  ## get index of flatten variable
  dim.flatten <- which(names(dimnames(x))==flatten.name)

  ## find differences
  dim.use <- 1:length(dimnames(x))
  dim.use <- setdiff(dim.use,dim.flatten)

  out  <- as.matrix(data.change[, unlist(dimnames(x)[dim.flatten])])
  dim(out) <- dim(x)[c(dim.use,dim.flatten)]
  dimnames(out) <- dimnames(x)[c(dim.use,dim.flatten)]

  ## put order as original
  array.order <-  match(names(dimnames(xorig)),names(dimnames(out)))
  out <- aperm(out,array.order)

  return(out)
}




## function to get z,w,u combinations to evaluate S_j(t|z,w) = S(t|z,w,u_j)
get.zw.evaluate <- function(data.evaluate.orig,m){
  data.evaluate <- data.evaluate.orig[rep(seq_len(nrow(data.evaluate.orig)),each=m),]
  data.evaluate <- cbind(data.evaluate,rep(1:m,nrow(data.evaluate)/m))

  colnames(data.evaluate) <- c("z","w","u")
  rownames(data.evaluate) <- 1:nrow(data.evaluate)

  return(data.evaluate)
}


###################################################
##SURVIVAL USING Z and W COVARIATE INFORMATION######
###################################################
#data structure: first column should be observed event time (or censoring time), second column should be indicator of whether an event was observed (Delta), third column should be discrete covariate Z, fourth column should be continuous covariate W, fifth column should indicate the U group that each person is in (for example if there are 3 "u" groups, u_1, u_2, u_3, then the fifth column would either have "1","2", or "3"); returns the data matrix with an extra column, the extra column is the survival probability at that Z and W
#newdata is an optional n by 3 matrix where the first column is the discrete covariate Z, the second column is the continuous covariate W and the third column is the U group. Predicted survival probabilities are estimated for these data;returns the data matrix with an extra column, the extra column is the survival probability at that Z and W

S.NPNA.zw = function(t, data, newdata = NULL, weight = NULL){
  problem <- NULL
  Xi = data[,1]
  Di = data[,2]
  Zi = data[,3]
  Wi = data[,4]
  Ui = data[,5]

  if(sum(Xi > t) == 0) {message(paste("No observations past time t=",t))}
  if(is.null(weight)) {weight = rep(1,length(Xi))}


  zi.cat = unique(Zi)
  ui.cat = unique(Ui)

  ## 3/16/2017: Below will only run if newdata is NULL
  if(is.null(newdata)){
    survival.v <- rep(NA, length = dim(data)[1]) ## Changed 1/30/2017
    ##survival.v = vector(length = dim(data)[1])
    for(j in 1:length(ui.cat)) {
      ##message(j)
      for(k in 1:length(zi.cat)) {
        ##message(k)
        ui.value = ui.cat[j]
        zi.value = zi.cat[k]

        if(sum(Zi==zi.value & Ui ==ui.value) < 10) {
          message(paste("Warning: Very few individuals with covariate value = ",
                      zi.value, ",
                      and in U group = ",ui.value))}

        if(sum(Zi==zi.value & Ui == ui.value & Xi > t) < 10) {
          message(paste("Warning: Very few individuals
                      observed to survive past t=",t," with covariate value = ",
                      zi.value, ", and in U group = ",ui.value))
        }

        if(length( Wi[Zi == zi.value & Ui == ui.value] )<2){  ## Changed 1/30/2017
          problem <- TRUE
          message(paste("Warning-Error: Less than 2 individuals with W covariate for z-covariate value = ",
                      zi.value, ", and in U group = ",ui.value,
                      ". Setting estimate=0.  W is",
                      Wi[Zi == zi.value & Ui == ui.value ], sep=""))
          P.return <- 0
        } else {
          P.return <- pred.smooth.surv(w.vector = Wi[Zi == zi.value & Ui == ui.value], t=t,
                                       data.use = data, covariate.value = 	zi.value,
                                       group.value = ui.value, weight = weight)
        }
        survival.v[Zi == zi.value & Ui == ui.value] <- P.return
        }
      }
    data = cbind(data, survival.v)
}

  if(!is.null(newdata)) {
    Zi.new = newdata[,1]
    Wi.new = newdata[,2]
    Ui.new = newdata[,3]

    survival.v.new = rep(NA,length = dim(newdata)[1]) ## Changed 1/30/2017
    ##survival.v.new = vector(length = dim(newdata)[1])
    for(j in 1:length(ui.cat)) {
      for(k in 1:length(zi.cat)) {
        ui.value = ui.cat[j]
        zi.value = zi.cat[k]

        ##3/16/2017: Added warning message
        #if(sum(Zi==zi.value & Ui ==ui.value) < 10) {
        #  message(paste("Warning: Very few individuals with covariate value = ",
        #  			zi.value, ",
        #    and in U group = ",ui.value))}

        ##3/16/2017: Added warning message
        #if(sum(Zi==zi.value & Ui == ui.value & Xi > t) < 10) {
        #  message(paste("Warning: Very few individuals
        #                  observed to survive past t=",t," with covariate value = ",
        #                  zi.value, ", and in U group = ",ui.value))
        #}


        if(length( Wi[Zi == zi.value & Ui == ui.value] )<2){ ## Changed 1/30/2017
          problem <- TRUE
          message(paste("Warning-Error: Less than 2 individuals with W covariate
                      for z-covariate value = ",
                      zi.value, ", and in U group = ",ui.value,
                      ". Setting estimate=0. W is",
                      Wi[Zi == zi.value & Ui == ui.value],
                      sep=""))
          P.return <- 0
        } else {
          P.return = pred.smooth.surv(w.vector = Wi.new[Zi.new == zi.value &
                                                          Ui.new == ui.value],
                                      t=t,data.use = data, covariate.value =  zi.value, group.value = ui.value,
                                      weight = weight)
        }
        survival.v.new[Zi.new == zi.value & Ui.new == ui.value] = P.return
      }
    }
    newdata.matrix = cbind(newdata, survival.v.new)
  }
  if(is.null(newdata)) {return(list(data.out=data,problem=problem))}
  if(!is.null(newdata)) {return(list(data.out=newdata.matrix,problem=problem))}
  }


## function to get w-distribution
w.sample <- function(n){
  cag.values <- seq(30,60,by=1)
  out <- sample(x=cag.values,n,replace=TRUE,prob=rep(1,length(cag.values))/length(cag.values))
  return(out)
}

#' @import stats
gendata.zw <- function(n,simu.setting){
  zz <- rbinom(n,1,0.5)

  if(simu.setting=="1A" | simu.setting=="1B"){
    ww <- runif(n,0,1)
  } else if(simu.setting=="2A" | simu.setting=="2B"){
    ww <- w.sample(n)
    ##ww <- runif(n,30,60)
  }


  ## setting without covariates: set z=w=0
  ##zz <- rep(0,n)
  ##ww <- rep(0,n)

  list(zz=zz,ww=ww)
}


getsd <- function(p){
  sd.use <- seq(1,2,length.out=p)
  return(sd.use)
}

mult.w <- function(covariate.dependent){
  if(covariate.dependent==TRUE){
    out <- 1
  } else {
    out <- 0
  }
  return(out)
}

mult.z <- function(covariate.dependent){
  if(covariate.dependent==TRUE){
    out <- 0.5
  } else {
    out <- 0
  }
  return(out)
}

time.generate <- function(u,t.low,t.high,w,z,mu0,mu1,mu2,sigma,r,a,p,w.constant){
  out <- mu0 - w.constant*(w-mu1) - mu2*z -
       sigma * log ((u *
       	     normalizing.constant(t.low,t.high,w,z,mu0,mu1,mu2,sigma,r,a,w.constant)-
       	       	       intersection.constant(t.low,w,z,mu0,mu1,mu2,sigma,r,a,p,w.constant) )^(1/r)-1)
  return(out)
}

log.limit <- function(u,t.low,t.high,w,z,mu0,mu1,mu2,sigma,r,a,p,w.constant){
  out <- (u * normalizing.constant(t.low,t.high,w,z,mu0,mu1,mu2,sigma,r,
      	    a,w.constant)-
                       intersection.constant(t.low,w,z,mu0,mu1,mu2,sigma,r,a,
		       p,w.constant) )^(1/r)-1
  return(out)
}


get.empty.list <- function(names){
  out <-  vector("list",length(names))
  names(out) <- names
  return(out)
}

F.exp <- function(t,w,z,mu0,mu1,mu2,sigma,r,w.constant){
  out <- ( 1+exp( - ((t-mu0)+w.constant * (w-mu1)+mu2*z)/sigma)  )^r
  return(out)
}

cumsum2 <- function(mydat)     #cumsum by row, col remains the same
{
  if(is.null(dim(mydat))) return(cumsum(mydat))
  else{
    out <- matrix(cumsum(mydat), nrow=nrow(mydat))
    out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
    return(out)
  }
}

intersection.constant <- function(t.low,w,z,mu0,mu1,mu2,sigma,r,a,p,w.constant){
  if(p==2){
    out <- a*t.low - F.exp(t.low,w,z,mu0,mu1,mu2,sigma,r,w.constant)
  } else {
    out <- 0
  }
  return(out)
}


F.constants <- function(p,covariate.dependent=TRUE){
  if(p==1){
    mu1 <- 40
    mu0 <- 43
    if(covariate.dependent==TRUE){
      ## covariates impact onset ages
      mu2 <- -0.5
      w.constant <- 1
    } else {
      mu2 <- 0
      w.constant <- 0
    }
    sigma <- 7
    a <- 0
    t.low <- 0
    t.high <- 100
    r <- -0.9
  } else {
    mu1 <- 42
    mu0 <- 48
    if(covariate.dependent==TRUE){
      ## covariates impact onset ages
      w.constant <- 1
      mu2 <- -0.5
    } else {
      mu2 <- 0
      w.constant <- 0
    }
    sigma <- 10.5
    a <- 0.0007
    t.low <- 13
    t.high <- 100
    r <- -2
  }
  list(mu0=mu0,mu1=mu1,mu2=mu2,sigma=sigma,a=a,t.low=t.low,t.high=t.high,r=r,
	w.constant=w.constant)
}

## function to set w
get.w.set <-function(w,mu1,covariate.dependent){
  if(covariate.dependent==TRUE){
    ## survival times depend on covariates
    w.set <- w
  } else {
    ## survival times DO NOT depend on covariates
    w.set <- mu1
  }
 return(w.set)
}

normalizing.constant <- function(t.low,t.high,w,z,mu0,mu1,mu2,sigma,r,a,w.constant){
  out <- a*(t.low) + F.exp(t.high,w,z,mu0,mu1,mu2,sigma,r,w.constant)-F.exp(t.low,w,z,mu0,mu1,mu2,sigma,r,w.constant)
  return(out)
}

trueinvFt <- function(p,w,z,simu.setting,covariate.dependent){
  if(simu.setting =="1A" | simu.setting=="1B"){
    out <- rep(0,p)
    constant.z <- mult.z(covariate.dependent)
    constant.w <- mult.w(covariate.dependent)
    sd.use <- getsd(p)
    for(jj in 1:p){
      ## Model: log(T)=W+0.5Z+N, so T=exp(W+0.5Z+N)
      out[jj] <- exp(constant.w * w + constant.z * z + rnorm(1,0,sd.use[jj]))
    }
  } else if(simu.setting=="2A" | simu.setting=="2B"){
    out <- rep(0,p) ## p=2 (fixed)

    ## get constants needed
    constants.F1 <- F.constants(1,covariate.dependent)
    mu0.F1 <- constants.F1$mu0
    mu1.F1 <- constants.F1$mu1
    mu2.F1 <- constants.F1$mu2
    w.constant.F1 <- constants.F1$w.constant
    sigma.F1 <- constants.F1$sigma
    a.F1 <- constants.F1$a
    t.low.F1 <- constants.F1$t.low
    t.high.F1 <- constants.F1$t.high
    r.F1 <- constants.F1$r

    ## set w-value
    w.set.F1 <- get.w.set(w,mu1.F1,covariate.dependent)

    use.value <- FALSE
    while(use.value==FALSE){
      u <- runif(1)
      if( log.limit(u,t.low.F1,t.high.F1,w=w.set.F1,z,mu0.F1,
		mu1.F1,mu2.F1,sigma.F1,r.F1,a.F1,p=1,w.constant.F1) > 0 ){
        use.value <- TRUE
      }
    }
    out[1] <- time.generate(u,t.low.F1,t.high.F1,w=w.set.F1,z,mu0.F1,
  	    mu1.F1,mu2.F1,sigma.F1,r.F1,a.F1,p=1,w.constant.F1)

    ## get constants needed
    constants.F2 <- F.constants(2,covariate.dependent)
    mu0.F2 <- constants.F2$mu0
    mu1.F2 <- constants.F2$mu1
    mu2.F2 <- constants.F2$mu2
    w.constant.F2 <- constants.F2$w.constant
    sigma.F2 <- constants.F2$sigma
    a.F2 <- constants.F2$a
    t.low.F2 <- constants.F2$t.low
    t.high.F2 <- constants.F2$t.high
    r.F2 <- constants.F2$r

    ## set w-value
    w.set.F2 <- get.w.set(w,mu1.F2,covariate.dependent)

    use.value <- FALSE
    while(use.value==FALSE){
      u <- runif(1)
      if(u < a.F2 * t.low.F2 / normalizing.constant(t.low.F2,
    	       t.high.F2,w=w.set.F2,z,mu0.F2,mu1.F2,mu2.F2,sigma.F2,r.F2,
	       a.F2,w.constant.F2) |
	       log.limit(u,t.low.F2,t.high.F2,w=w.set.F2,z,mu0.F2,
                mu1.F2,mu2.F2,sigma.F2,r.F2,a.F2,p=2,w.constant.F2) > 0 ) {
        use.value <- TRUE
      }
    }

    if(u <  a.F2 * t.low.F2 / normalizing.constant(t.low.F2,
               t.high.F2,w=w.set.F2,z,mu0.F2,mu1.F2,mu2.F2,sigma.F2,r.F2,
	       a.F2,w.constant.F2)){
      out[2] <- u * normalizing.constant(t.low.F2,
               t.high.F2,w=w.set.F2,z,mu0.F2,mu1.F2,mu2.F2,sigma.F2,
	       r.F2,a.F2,w.constant.F2)/ a.F2
    } else {
      out[2] <- time.generate(u,t.low.F2,t.high.F2,w=w.set.F2,z,mu0.F2,
  	    mu1.F2,mu2.F2,sigma.F2,r.F2,a.F2,p=2,w.constant.F2)
    }
  }
  return(out)
}


genc <- function(s,censoring.rate,simu.setting){
  if(censoring.rate==0){
    out <- s + 1
  } else if(censoring.rate==20){
    if(simu.setting=="1A"){
      out <- runif(1,0,25)
    } else if(simu.setting=="1B"){
      out <- runif(1,0,22)
    } else if(simu.setting=="2A"){
      #out <- runif(1,50,90)
      ##out <- runif(1,0,80)
      out <- runif(1,0,210)
    } else if(simu.setting=="2B"){
      out <- runif(1,0,220)
    }
  } else if(censoring.rate==40){
    if(simu.setting=="1A"){
      out <- runif(1,0,8)
    } else if(simu.setting=="1B"){
      out <- runif(1,0,4)
    } else if(simu.setting=="2A"){
      #out <- runif(1,50,90)
      ##out <- runif(1,0,80)
      out <- runif(1,10,100)
    } else if(simu.setting=="2B"){
      out <- runif(1,15,100)
    }
  }
  return(out)
}





## Computes the modified nonparametric Nelson-Aalen estimator.
#' @import stats
pred.smooth.surv <- function(w.vector=NULL, t, data.use, covariate.value, group.value, weight)
{ Xi.use = data.use[,1]
Di.use = data.use[,2]
Zi.use = data.use[,3]
Wi.use = data.use[,4]
Ui.use = data.use[,5]
if(is.null(w.vector)) {w.vector = Wi.use}
h = bw.nrd(Wi.use[Zi.use == covariate.value & Ui.use == group.value])
bandwidth = h*length(Xi.use[Zi.use == covariate.value & Ui.use == group.value])^(-.10)
K = Kern.FUN(Wi.use[Zi.use == covariate.value & Ui.use == group.value],w.vector,bandwidth)
Xi.temp = data.use[Zi.use == covariate.value & Ui.use == group.value,1]
Di.temp = data.use[Zi.use == covariate.value & Ui.use == group.value,2]
Wi.temp = data.use[Zi.use == covariate.value & Ui.use == group.value,4]
tmpind = (Xi.temp<=t)&(Di.temp==1)
tj = Xi.temp[tmpind];
kerni.1 = t(weight[Zi.use == covariate.value & Ui.use == group.value]*t(K))
pihamyt0.tj.ss = sum.I(tj, "<=", Xi.temp, Vi=t(kerni.1)) ## n.tj x n.ss matrix ##
dLamhat.tj.ss = t((kerni.1[,tmpind]))/pihamyt0.tj.ss;
#dLamhat.tj.ss[is.na(dLamhat.tj.ss)] = 0
ret = apply(dLamhat.tj.ss,2,sum)
S.return  =exp(-ret)
return(S.return)
}




GenerateData <- function(n,p,m,qvs,censoring.rate,simu.setting,covariate.dependent){


  q <- array(0,dim=c(p,n))  ## mixture proportions
  x <- rep(0,n)		    ## observed event time: min(T,C)
  delta <- rep(0,n)	    ## censoring indicator
  uset <- rep(0,n)	    ## indicator of which subgroup qvs
  r <- rep(0,m)		    ## number in each qvs subgroup
  true.groups <- rep(0,n)   ## to which group each person belongs

  ###################
  ## Simulate data ##
  ###################

  ## set evenly spaced subgroup sizes
  r0 <- seq(0,1,length=m+1)

	for(i in 1:n){
      a <- runif(1)
      for(j in 1:m){
        if( (a >= r0[j]) && (a < r0[j+1]) ){
          q[,i] <- qvs[,j]
	  r[j]  <- r[j] + 1
	  uset[i] <- j
        }
      }
    }

    ## set covariates
    zw.data <- gendata.zw(n,simu.setting)
    zz <- zw.data$zz
    ww <- zw.data$ww



    for(i in 1:n){
      q1<- q[,i]
	
	  negative_value <- TRUE
	  while(negative_value==TRUE){
		t1 <- trueinvFt(p,ww[i],zz[i],simu.setting,covariate.dependent)
		if(all(t1>0)){
			negative_value <- FALSE
		}
	  }	
	  
      a <- runif(1)
      tmp <- q1[1]
      j <- 1
      while(j <= p){
        if(a <= tmp){
          s <- t1[j]
	  true.groups[i] <- j
	  j <- p+1
        } else {
          tmp <- tmp + q1[j+1]
	  j <- j+1
        }
      }

      cens <- genc(s,censoring.rate,simu.setting)
      x[i] <- min(s,cens)
      delta[i] <- which.min(c(cens,s)) - 1
    }

  list(x=x,delta=delta,q=q,ww=ww,zz=zz,true.groups=true.groups)
}


qvs.values <- function(p,m){
  qvs <- matrix(0,nrow=p,ncol=m)
  colnames(qvs) <- paste("m",1:m,sep="")
  rownames(qvs) <- paste("p",1:p,sep="")

  tmp <- 0.2
  qvs[1,1] <- 1
  qvs[1,2] <- (1+tmp)/2
  qvs[1,3] <- tmp
  qvs[1,m] <- 1+tmp*tmp/4-tmp/2-3/4

  qvs[1,] <- sort(qvs[1,])
  qvs[2,] <- 1-qvs[1,]

  return(qvs)
}

sum.I <- function(yy,FUN,Yi,Vi=NULL)   ## sum I(yy FUN Yi)Vi
{
  if(FUN=="<"|FUN=="<=") { yy <- -yy; Yi <- -Yi}
  if(substring(FUN,2,2)=="=") yy <- yy + 1e-8 else yy <- yy - 1e-8
  pos <- rank(c(yy,Yi))[1:length(yy)] - rank(yy)
  if(is.null(Vi)){return(pos)}else{
    Vi <- cumsum2(as.matrix(Vi)[order(Yi),,drop=F])
    out <- matrix(0, nrow=length(yy), ncol=dim(as.matrix(Vi))[2])
    out[pos!=0,] <- Vi[pos,]
    if(is.null(dim(Vi))) out <- c(out)
    return(out) ## n.y x p
  }
}




########################
##HELPER FUNCTION######
########################
Kern.FUN <- function(zz,zi,bw) ## returns an (n x nz) matrix
{
  out = (VTM(zz,length(zi))- zi)/bw
  norm.k = dnorm(out)/bw
  norm.k
}


########################
##HELPER FUNCTION######
########################
VTM<-function(vc, dm){
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}


## function to have qvs.boot and qvs have the same column names
## and in the correct order.
match.names <- function(qvs.boot,qvs){
  index.bucket <- NULL
  for(kk in 1:ncol(qvs.boot)){
    for(ii in 1:ncol(qvs)){
      if(all(qvs.boot[,kk]==qvs[,ii])){
        index.bucket <- c(index.bucket,ii)
      }
    }
  }
  return(index.bucket)
}



## another function to re-order qvs
get_new_qvs <- function(qvs_new,qvs){
  ## re-name qvs_new with matching names as in qvs

  qvs_new_names <- match.names(qvs_new,qvs)
  colnames(qvs_new) <- paste("m",qvs_new_names,sep="")
  qvs_new <- qvs_new[,order(qvs_new_names)]

  ## get new m
  m_new <- ncol(qvs_new)
  return(list(qvs_new=qvs_new,m_new=m_new))
}



