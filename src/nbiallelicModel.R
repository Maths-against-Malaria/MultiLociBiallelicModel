# Title        : Model
# Objective    : Implementation of the model (EM-algorithm)
# Created by   : Christian Tsoungui Obama, Kristan. A. Schneider
# Created on   : 03.04.21
# Last modified: 06.04.22

liklh <- function (Bx, nl, ppest, laest, nx){
  la <- laest
  pp <- unlist(ppest)
  nn <- nl
  likl <- 1
  Nx <- nx
  for(u in 1:nn){
    out <- 0
    for(k in 1:Bx[[u]][[3]]){                       # For each y in Ax
      p <- sum(pp[as.numeric(Bx[[u]][[1]][[k]])])   # Summing up the frequencies of each haplotypes making y
      vz <- Bx[[u]][[2]][[k]]                       # Inclusion-exclusion term (-1)^(Nx - Ny)
      lap <- la*p
      exlap <- vz*exp(lap)
      out <- out + exlap - vz                      # exlap-vz  (I can't find out the rationale for the substraction)
    }
    out <- (out/(exp(la) - 1))^Nx[u]                # Computing the likelihood
    likl <- likl*out
  }
  likl
}

varsets <- function(l, n) {
  # l = 2^2-1 (case of biallelic loci) # n = number of heterozygote loci
  B <- array(0,c(l^n,n))
  B[1:l,1] <- 0:(l-1)
  lkmo <- l
  if(n>1){
    for(k in 2:n){
      lk <- lkmo*l
      pick1 <- (lkmo+1):lk
      B[pick1,] <- B[rep(1:lkmo,l-1),]
      B[pick1,k] <- rep(1:(l-1),each=lkmo)
      lkmo <- lk
    }
  }
  B
}

nbialModel <- function(Nx, X){
  eps <- 10^-8  # Error
  N <- sum(Nx)  # Sample size
  nn <- nrow(X) # Number of different observations present in the dataset
  n <- ncol(X)  # Number of loci
  Ax <- list()
  
  for(u in 1:nn){
    xx <- array(X[u,],c(1,n)) # xx = observation
    sel <- (1:n)[xx==2]       # Identifying the loci where the 2 alleles are observed
    l <- length(sel)          # Counting the number of loci where the 2 alleles are observed
    
    if(l==0){                 # If the infection is a haplotype (only one allele per locus)
      yy <- xx
    }else{ 
      yy <- xx[rep(1,3^l),]
      yy[,sel] <- varsets(3,l) # Set of all possible observations which combinations can form xx $\mathscr{A}_{y}$
    }
    bin <- 2^((n-1):0)
    iilist <- list()
    siglist <- list()
    for(i in 1:3^l){
      y1 <- array(yy[i,],c(1,n)) # Observation {\pmb y} in the set $\mathscr{A}_{y}$
      sel <- (1:n)[y1==2]
      l1 <- length(sel)
      if(l1==0){
        ii <- y1
      }else{
        ii <- y1[rep(1,2^l1),]
        ii[,sel] <- varsets(2,l1)
      }
      iilist[[i]] <- as.character(ii%*%bin+1)
      siglist[[i]] <- (-1)^(l-l1)
    }
    Ax[[u]] <- list(iilist,siglist,3^l)
  }
  
  # list of all occuring halotypes 
  hapl1 <- c()
  for(u in 1:nn){
    hapl1 <- c(unlist(Ax[[u]][[1]]),hapl1)
  }
  hapl1 <- unique(hapl1)
  H <- length(hapl1)
  
  pp <- array(rep(1/H,H),c(H,1))   #  Initial frequency distribution (of the observed haplotypes) for the EM algorithm
  rownames(pp) <- hapl1
  
  la <- 2                          # Initial value of lambda for the EM algo.
  
  num0 <- pp*0
  cond1 <- 1                       # Initializing the condition to stop EM algorithm 
  Bcoeff <- num0
  num <- num0
  rownames(num) <- hapl1
  rownames(Bcoeff) <- hapl1
  t <- 0
  
  while(cond1>eps && t<500){
    t <- t + 1
    Ccoeff <- 0
    Bcoeff <- num0                # reset B coefficients to 0 in next iteration
    num <- num0                   # reset numerator to 0 in next iteration
    for(u in 1:nn){               # For all possible observation
      denom <- 0
      num <- num0
      CC <- 0
      for(k in 1:Ax[[u]][[3]]){   # For all h in Ay
        p <- sum(pp[Ax[[u]][[1]][[k]],]) # Be careful with this sum!!!!!
        vz <- Ax[[u]][[2]][[k]]
        lap <- la*p
        exlap <- vz*exp(lap)
        denom <- denom + exlap-vz 
        num[Ax[[u]][[1]][[k]],] <- num[Ax[[u]][[1]][[k]],]+ exlap
        CC <- CC + exlap*p
      }
      num <- num*pp
      denom <- Nx[u]/denom
      denom <- la*denom
      Ccoeff <- Ccoeff + CC*denom
      
      Bcoeff <- Bcoeff + num*denom
      
    }
    
    Ccoeff <- Ccoeff/N

      # Replacing NaN's in Ak by 0
      cnt <- 0
      for(i in seq_along(Bcoeff)){
        if (is.nan(Bcoeff[i])){
          cnt <- cnt + 1
        }
      }
      if(cnt > 0){
        break
      }else{
        ppn <- Bcoeff/(sum(Bcoeff))
      }
      ##************* 1-Dimensional Newton-Raphson to estimate the lambda parameter
      cond2 <- 1
      xt <- Ccoeff
      tau <- 0

      while(cond2 > eps &&  tau<300){
        tau <- tau+1
        ex <- exp(-xt)
        xtn <- xt + (1-ex)*(xt + Ccoeff*ex - Ccoeff)/(ex*xt+ex-1)

        if(is.nan(xtn) || (tau == 299) || xtn < 0){ #Replacing NA and negative lambda value by random values
          xtn <- runif(1, 0.1, 2.5)
        }
        cond2 <- abs(xtn-xt)
        xt <- xtn
      }
      # Filling the frequency values = 0 in the frequency vector
      cond1 <- abs(xt-la)+sqrt(sum((pp-ppn)^2))
      la <- xt
      pp <- ppn
    }

    ## Ordering the frequencies
    pp <- pp[order(as.numeric(rownames(pp))), ]

    ## Setting the frequencies of the unobserved haplotypes to 0.0
    nhapl <- 2^n

    if(length(pp)<nhapl){
       out <- t(pp)
      name <- colnames(out)
      cnt <- 0
     for (i in 1:nhapl) {
       if (is.element(as.character(i), name)){
          cnt <- cnt + 1
        }else{
          pp <- append(pp, list(x = 0.0), i-1)
        }
      }
    }

  return(c(la, pp))
}

hapl <- function(n){
  H <- array(0,c(2^n,n))
  H[1:2,1] <- c(0,1)
  for(k in 2:n){
    H[(2^(k-1)+1):2^k,1:(k-1)] <- H[1:2^(k-1),1:(k-1)]
    H[(2^(k-1)+1):2^k,k] <- 1
  }
  H <- H[,n:1]
  H
}

adhocModel <- function(X){
  nloci <- ncol(X)
  # estimate haplotype frequencies
  nhpl <- 2^nloci
  p <- rep(0, nhpl)
  
  # extract unambiguous observations
  X1 <- X[rowSums(X==2)<2,]

  if(!all(is.na(X1))){  # if there are unambiguous infections
    n1 <- nrow(X1)
    if(is.null(n1)){ # if there is only one unambiguous infection
      n1 <- 1
    }
    X <- matrix(X1, nrow = n1)
    # find indexes of multiple infections
    idx1 <- which(rowSums(X==2)==1)

    if(length(idx1)>0){
      # single infections
      s <- X[-idx1,]
    }else {
      s <- X
    }
    
    # find all the haplotypes in X
    for(i in idx1){
      y <- X[i,]
      idx2 <- which(y==2)
      h <- array(rep(y,2), c(nloci, 2))
      h[idx2,] <- c(0,1)
      # add haplotypes in s
      s <- rbind(s,t(h))
    }
    
    # binary representation
    bin <- 2^((nloci-1):0)
    pp <- rowSums(s*bin)+1
    pp <- table(pp)/length(pp)
    
    # indexes
    idx <- as.numeric(names(pp))
    p[idx] <- pp
  }
  p
}
