#newick2phylog function read tree in eNewick format. 
#The input data files will be x.tre.
newick2phylog <- function (x.tre, add.tools =FALSE, call = match.call()){ 
  lvs_num<-function(x.tree){
    x.tree<-unlist(strsplit(x.tre, ""))
    lvs<-0 
    for(i in 1:length(x.tree)){
      if(x.tree[i]==":"){
        lvs=lvs+1
           }
        }
      return((lvs+1)/2)
      }
  complete <- function(x.tre){
    if(length(x.tre) > 1) {
      w<- ""
      for(i in 1:length(x.tre)){ 
        w <- paste(w, x.tre[i], sep = "")}
        x.tre <- w
        }
    ndroite <- nchar(gsub("[^)]", "", x.tre))
    ngauche <- nchar(gsub("[^(]", "", x.tre))
    if(ndroite != ngauche) {stop(paste(ngauche, "( versus", ndroite, ")"))}
    if (regexpr(";", x.tre) == -1){stop("';' not found")}
      i    <- 0
      kint <- 0
      kext <- 0
      arret <- FALSE
      lvs<-lvs_num(x.tree)
      if(regexpr("\\[", x.tre) != -1) {
        x.tre <- gsub("\\[[^\\[]*\\]", "", x.tre, ext = FALSE)
        }
      x.tre <- gsub(" ", "", x.tre, ext = FALSE)
      while(!arret){
        i <- i + 1
        if(substr(x.tre, i, i) == ";") 
          {arret <- TRUE}
        if(substr(x.tre, i, i + 1) == "(,") {
          kext <- kext + 1
          add <- paste("Ext", kext, sep = "")
          x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre,i + 1), sep = "")
          i <- i + 1
          }
        else if (substr(x.tre, i, i + 1) == ",,") {
          kext <- kext + 1
          add <- paste("Ext", kext, sep = "")
          x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, i + 1), sep = "")
          i <- i + 1
          }
        else if (substr(x.tre, i, i + 1) == ",)") {
          kext <- kext + 1
          add <- paste("Ext", kext, sep = "")
          x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre,i + 1), sep = "")
          i <- i + 1
          }
        else if (substr(x.tre, i, i + 1) == "(:") {
          kext <- kext + 1
          add <- paste("Ext", kext, sep = "")
          x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre, i + 1), sep = "")
          i <- i + 1
          }
        else if (substr(x.tre, i, i + 1) == ",:") {
          kext <- kext + 1
          add <- paste("Ext", kext, sep = "")
          x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre,i + 1), sep = "")
          i <- i + 1
          }
        else if (substr(x.tre, i, i + 1) == "),") {
          kint <- kint + 1
          add <- kint+lvs
          x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre,i + 1), sep = "")
          i <- i + 1
          }
        else if (substr(x.tre, i, i + 1) == "))") {
          kint <- kint + 1
          add <- kint+lvs
          x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre,i + 1), sep = "")
          i <- i + 1
          }
        else if (substr(x.tre, i, i + 1) == "):") {
          kint <- kint + 1
          add <- kint+lvs
          x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre,i + 1), sep = "")
          i <- i + 1
          }
        else if (substr(x.tre, i, i + 1) == ");") {
          add <- "Root"
          x.tre <- paste(substring(x.tre, 1, i), add, substring(x.tre,i + 1), sep = "")
          i <- i + 1
          }
        }
                        
    lab.points <- strsplit(x.tre, "[(),;]")[[1]] 
    lab.points <- lab.points[lab.points != ""]
    no.long <- (regexpr(":", lab.points) == -1)
    if(all(no.long)) { 
      lab.points <- paste(lab.points, ":", c(rep("1", length(no.long) - 1), "0.0"), sep = "")
      }
    else if (no.long[length(no.long)]) {
      lab.points[length(lab.points)] <- paste(lab.points[length(lab.points)], ":0.0", sep = "")
      }
    else if (any(no.long)) {
      stop("Non convenient data leaves or nodes with and without length")
      }
            
    w <- strsplit(x.tre, "[(),;]")[[1]]       
    w <- w[w != ""]           
    leurre <- make.names(w, unique = TRUE)           
    leurre <- gsub("[.]", "_", leurre, ext = FALSE)                   
    for(i in 1:length(w)) {
      old <- paste(w[i])
      x.tre <- sub(old, leurre[i], x.tre, ext = FALSE)
      } 
           
    w <- strsplit(lab.points, ":")         
    label <- function(x){
    lab <- x[1]
    lab <- gsub("[.]", "_", lab, ext = FALSE)
    return(lab)
    }
        
    longueur <- function(x) {
      long <- x[2]
      return(long)
      }
        
   labels <- unlist(lapply(w, label))         
     longueurs <- unlist(lapply(w, longueur))                   
     labels <- make.names(labels, TRUE)          
     labels <- gsub("[.]", "_", labels, ext = FALSE)           
     w <- labels                    
     for(i in 1:length(w)) {
       new <- w[i]
       x.tre <- sub(leurre[i], new, x.tre, ext = FALSE)
       }
                       
     cat <- rep("", length(w))                   
     for(i in 1:length(w)) {
       new <- w[i]
       if(regexpr(paste(")", new, sep = ""), x.tre, ext = FALSE) != -1) 
         cat[i] <- "int"
       else if (regexpr(paste(",", new, sep = ""), x.tre, ext = FALSE) != -1) 
         cat[i] <- "ext"
       else if (regexpr(paste("(", new, sep = ""), x.tre, ext = FALSE) != -1) 
         cat[i] <- "ext"
       else cat[i] <- "unknown"
         }
          
      return(list(tre = x.tre, noms = labels, poi = as.numeric(longueurs),cat = cat))
                    
      }#end of complete function
        
    res <- complete(x.tre)
    poi <- res$poi 
    nam <- res$noms     
    names(poi) <- nam
    cat <- res$cat
    res <- list(tre = res$tre)       
    res$leaves <- poi[cat == "ext"]         
    names(res$leaves) <- nam[cat == "ext"]
    res$nodes <- poi[cat == "int"]
    names(res$nodes) <- nam[cat == "int"]
    listclass <- list()
    dnext <- c(names(res$leaves), names(res$nodes))  
    res$nde<-c(names(res$nodes))
    res$nme<-c(names(res$leaves))
    listpath <- as.list(dnext)
    names(listpath) <- dnext
    x.tre <- res$tre
    print(x.tre)
    while (regexpr("[(]", x.tre) != -1) {
      a <- regexpr("([^()]*)", x.tre, ext = FALSE)
      n1 <- a[1] + 1
      n2 <- n1 - 3 + attr(a, "match.length")
      chasans <- substring(x.tre, n1, n2)
      chaavec <- paste("(", chasans, ")", sep = "")
      nam <- unlist(strsplit(chasans, ","))
      w1 <- strsplit(x.tre, chaavec, ext = FALSE)[[1]][2]
      parent <- unlist(strsplit(w1, "[,);]", ext = FALSE))[1]
      listclass[[parent]] <- nam
      x.tre <- gsub(chaavec, "", x.tre, ext = FALSE)
      w2 <- which(unlist(lapply(listpath, function(x) any(x[1] == nam))))
      for(i in w2) {
        listpath[[i]] <- c(parent, listpath[[i]])
        }
      }#this loop print the path for each node . 
    res$parts <- listclass
    res$paths <- listpath
    dnext <- c(res$leaves, res$nodes)
    names(dnext) <- c(names(res$leaves), names(res$nodes))    
    res$droot <- unlist(lapply(res$paths, function(x) sum(dnext[x])))
    res$call <- call
    class(res) <- "phylog"
    if (!add.tools) 
        return(res)
    return(newick2phylog.addtools(res))
    }# end of function newick2phylog.

comb<-function(n,r){
  return(factorial(n)/(factorial(n-r)*factorial(r)))
  }    
    
pair_fcn<-function(tmp){#return pair for tmp sequences.
  numl=comb(length(tmp),2)   
  count=0
  posit<-array(0,c(numl))
  for(i in 1:length(tmp)){
    for(j in 1:length(tmp)){
      if(i<j){
        count=count+1
        posit[count]=paste(c(tmp[i]),c(tmp[j]),sep=",")
        }
      }
    }  
  return(posit)
  }# end of pair_fcn.

pair_array<-function(tmp){#generate pair in array format. 
  pair <-pair_fcn(tmp)
  p_arr<-matrix(c(0),nrow=length(pair),ncol=2)
  for(i in 1:length(pair)){
    p_arr[i,1]=unlist(strsplit(pair[i],","))[1]
    p_arr[i,2]=unlist(strsplit(pair[i],","))[2]
    }
  return(p_arr)
  }# end of function pair_array.
    
ord_fcn<-function(){ #this function gets the acenstor-descendants relationship.
  bmtp<-matrix(rev(res$nde),ncol=1) 
  rvpt <-rev((res$parts))
  rept<-array(0,c(length(rvpt),2))
  for(i in 1:length(rvpt)){
    rept[i,]=unlist(rvpt[i])
    }           
  cmb<-cbind(bmtp,rept)
  brnlen<-res$droot[(length(tipnames)+1):length(res$droot)]
  root<-matrix(cmb[1,],nrow=1)
  cmb<-cmb[-1,]
  brnlen<-brnlen[1:(length(brnlen)-1)]
  new_ord<-order(brnlen,decreasing=TRUE)
  cmb<-cmb[new_ord,]
  cmb<-rbind(root,cmb)
  return(cmb)
  }# end of function ord_fcn.

getntn<-function(res){# this function gets rid of unnecessarily "_" symbol. 
  size<-length(res$parts)
  relarr<-array(0,c(size,3))
  rvpt <-(res$parts)
  rept<-array(0,c(length(rvpt),2))
  for(i in 1:length(rvpt)){
    rept[i,]=unlist(rvpt[i])
    }
  for(i in 1:size){
    relarr[i,1]<-names(res$parts)[i]
    }
  relarr[,2:3]<-rept
  temp<-matrix(0,row<-size)
    
  for(j in 2:3){
    for (i in 1: size){ 
      stmp<-unlist(strsplit(relarr[,j][i], "_" ))
      temp[i]<-stmp[1]
      }
    relarr[,j]<-temp
    }
    
  ndlen<- res$droot[!(res$droot==1)]
  nam<-names(ndlen)
  ck1<-array(0,c(length(nam)))
  count<-0
  for(ele in c(nam)){
    count<-count+1
    len <- length( unlist(strsplit(ele ,"_" )))
    if( len==2 ){ck1[count]<-1}
    }
  ndlen<-ndlen[!ck1]
  new_ord<-order(ndlen)
  relarr<-relarr[new_ord,]
    
  return(relarr)
  }# end of function getntn.
    
getbrnlen<-function(res){#this function gets branch length.
  ndlen<- res$droot[!(res$droot==1)]  
  nam<-names(ndlen)
  ck1<-array(0,c(length(nam)))
  count<-0
  for(ele in c(nam)){
    count<-count+1
    len <- length( unlist(strsplit(ele ,"_" )))
    if(len==2){ck1[count]<-1}
    }
  ndlen<-ndlen[!ck1]
  ndlen<-sort(ndlen)
  ck2<-array(0,c(length(ndlen)))
  for(i in 1:(length(ndlen)-1)){
    if(abs(ndlen[i]-ndlen[i+1])<10^(-5)){ck2[i]=1}
    }  
  ndlen<-ndlen[!ck2]  
  brnlen<-array(0,c(length(ndlen)))
  tmplen<-ndlen    
  for(i in 1:(length(brnlen)-1)){
    brnlen[i]<-tmplen[i+1]-tmplen[i]
    }
  brnlen[length(brnlen)] <-1-tmplen[(length(tmplen))]
  return(brnlen)
  }# end of function getbrnlen.

### The cov_mtx function will return the covaraince matrix. 
cov_mtx<-function(bt){    
  ins_fcn<-function(ist,sqc){#finds position to insert between two parents, 
                             #for hybrdization only.
    ist<-as.numeric(unlist(strsplit(ist,"X"))[2])
    arr<-array(0,c(length(otmp)))
    for(i in 1:length(arr)){
      arr[i]<-as.numeric(unlist(strsplit(sqc[i],"X"))[2])  
      }
    insp<-which(arr==(ist-1))+1
    return(insp)
    }
  var_fcn<-function(){#return the variance.
  for(i in 1:length(otmp)){#use to fill other diagonal. 
    newi<-which(rownames(mtx)%in%otmp[i])              
    oldi<-which(rownames(omtx)%in%otmp[i])
    mtx[newi,newi]<-omtx[oldi,oldi]
    }#fill in old value from omtx exclude the new hyd.
        
  prn1<-tmp[which(tmp%in%ins)-1]#grab elements.
  prn2<-tmp[which(tmp%in%ins)+1]
  prn1<-which(rownames(omtx) %in% prn1)#grab position according to prn1.
  prn2<-which(rownames(omtx) %in% prn2)
  
  ##############################################################
  vhii<-bt^2*(omtx[prn1,prn1]+omtx[prn2,prn2]+2*omtx[prn1,prn2]) 
  ##############################################################
  
  #fill in value into matrix with formula var(R)=bt^2(var(X)+var(Y)+2cov(X,Y)). 
  hii<-which(!(tmp %in% otmp))#use to insert variance for hyd.
  mtx[hii,hii]<-vhii#fill in the diagonal hyd. 
  return(mtx)
  }#formula for insertion hyd variance.
        
        
  fillspcmtx<-function(){#fill matrix due to sepciation.       
    elm<-function(){ #use to cut one row of the pair array which the speciation happens.
    ck<-c(tmp[nsi],tmp[nsj])
    for(i in 1:dim(pn_arr)[1]){
      if(sum(pn_arr[i,]==ck)==2){break}
      }
    return(i)}
        
  pn_arr<-pair_array(tmp)
  po_arr<-pair_array(otmp)
        
  #search new speciate position.
  nsi<-which(!(tmp %in% otmp))[1]
  nsj<-which(!(tmp %in% otmp))[2]
  osii<-which(!(otmp %in% tmp))
  mtx[nsi,nsj]<- omtx[osii,osii]
  #Fill in value: the covariance for 2 speciated species equal the variance of the parent.       
  pn_arr<-pn_arr[-elm(),]#delete the data that is already used.
       
  #The following fills covaraince components by the previous matrix.
  while(length(pn_arr[,1])>0){
    newi<-which(rownames(mtx) %in% pn_arr[1,1])
    newj<-which(rownames(mtx) %in% pn_arr[1,2])           
    if(tmp[nsi] %in% pn_arr[1,]){
      otg<-which(!(pn_arr[1,] %in%  tmp[nsi]))
      oldi<- which( rownames(omtx) %in% otmp[osii])
      oldj<-which(rownames(omtx) %in% pn_arr[1,otg])
      }
            
    if(tmp[nsj] %in% pn_arr[1,] ){
      otg<-which(!(pn_arr[1,] %in%  tmp[nsj]))
      oldi<- which( rownames(omtx) %in% otmp[osii])
      oldj<-which(rownames(omtx) %in% pn_arr[1,otg])
      }
                
    if(!(tmp[nsi] %in% pn_arr[1,]) && !(tmp[nsj] %in% pn_arr[1,])){
      #detect common between omtx and mtx.   
      oldi<-which(rownames(omtx) %in% pn_arr[1,1])
      oldj<-which(rownames(omtx) %in% pn_arr[1,2])
      }
    mtx[newi,newj]<-omtx[oldi,oldj]
    pn_arr<-pn_arr[-1,]#delete row. 
    if(length(pn_arr)==2){pn_arr<-matrix(pn_arr,nrow=1)}
    }#end of while loop.
            
    mtx<-mtx+t(mtx)        
    mtx[nsi,nsi]<-omtx[osii,osii]+ branchlength[length(tmp)-1]
    mtx[nsj,nsj]<-omtx[osii,osii]+ branchlength[length(tmp)-1]
    dianew<-which(tmp %in% otmp )
    diaold<-which(otmp %in% tmp )
    for(i in 1:length(dianew)){
      mtx[dianew[i],dianew[i]]<-omtx[diaold[i],diaold[i]]+branchlength[length(tmp)-1]
      }
    return(mtx)
    }#end of fillspcmtx.
        
    fillhydmtx<-function(){#fill in value into matrix due to hybridzation.   
      pn_arr<-pair_array(tmp)    
      while(length(pn_arr[,1])>0){
        newi<-which(rownames(mtx) %in% pn_arr[1,1])
        newj<-which(rownames(mtx) %in% pn_arr[1,2])
        if(ins %in% pn_arr[1,]){#ins is the hybridized node. 
          otg<-pn_arr[1,which(!(pn_arr[1,] %in% ins ))]
          otgj<-which(rownames(omtx) %in% otg)
          #the other species, could be the hybrdized node's parent or others.
          #find the parent of ins.
          prn1<-tmp[which(tmp%in%ins)-1]#grab element.
          prn2<-tmp[which(tmp%in%ins)+1]        
          prn1<-which(rownames(omtx) %in% prn1)#grab position.
          prn2<-which(rownames(omtx) %in% prn2)
                
          mtx[newi,newj]<-bt*(omtx[prn1,otgj] +omtx[prn2,otgj]) 
          }else{#this is not hyd node, just read from previous mtx.
          #only need to read from rownames().
          oldi<-which(rownames(omtx) %in% pn_arr[1,1])
          oldj<-which(rownames(omtx) %in% pn_arr[1,2])
          mtx[newi,newj]<-omtx[oldi,oldj]
          }#end of else loop .
        pn_arr<-pn_arr[-1,] # delete data array after using it.
        if(length(pn_arr)==2){pn_arr<-matrix(pn_arr,nrow=1)}
          }#end of while loop.
          return(mtx)
      }#end of fillhydmtx.
            
  #THE MAIN PROGRAM for covariance matrix.
  ckins<-FALSE # use to check the hybrdized event.
  rept<-data[,2:3]# the descedant nodes.
  bmtp<-matrix((data)[,1],ncol=1) #the acenstor node.
    
  loop<-2
  tmp=array(0,c(loop))
  if(loop==2){tmp=rept[1,]
    otmp<-tmp
    mtx<-diag(branchlength[1],c(length(tmp)))
    rownames(mtx)<-c(tmp)
    colnames(mtx)<-c(tmp)
    omtx<-mtx
    }#end of loop==2 
  while(loop<length(bmtp)){#loaded the acenstor-descendant data. 
    loop<-loop+1#use loop to use the data
    tmp=array(0,c(length(otmp)+1))#the new seq.
    mtx<-matrix(0,nrow=length(tmp),ncol=length(tmp))
    q=loop-1#index for matching the right element: will use below. 
    op<-which(otmp==bmtp[q])#index for insertion position.
    if(length(op)!=0){#op!=0 means that we've  detected speciation.
      tmp[op:(op+1)]=rept[q,] #insertion the new speciation species.      
        if(op==1){tmp[(op+2):length(tmp)]=otmp[(op+1):(length(tmp)-1)]}
        if((op+1)==length(tmp)){tmp[1:(op-1)]=otmp[1:(op-1)] }
        if(op!=1 && (op+1)!=length(tmp)){
          tmp[(op+2):length(tmp)]=otmp[(op+1):(length(tmp)-1)] 
          tmp[1:(op-1)]=otmp[1:(op-1)]}                              
          rownames(mtx)<-c(tmp)
          colnames(mtx)<-c(tmp)
          mtx<- fillspcmtx()
          otmp<-tmp
          omtx<-mtx                
          }else{#  op = 0 means that we have detected the hybridize event.
            ins<-(bmtp[q])#grab the insertion element, ins will be used in the fillhydmtx function.
            insp<-ins_fcn(ins,otmp)#catch the position for insertion.
            tmp[insp]<-ins #insert the hyd element.
            tmp[(insp+1):length(tmp)]=otmp[insp:(length(tmp)-1)]
            tmp[1:(insp-1)]=otmp[1:(insp-1)]
            rownames(mtx)<-c(tmp)
            colnames(mtx)<-c(tmp)
            diamtx<-var_fcn()
            mtx<- fillhydmtx()
            mtx<-mtx+t(mtx)+diamtx            
            otmp<-tmp
            omtx<-mtx
            ckins<-TRUE #since we did an insertion, the next step is to replace 3 elements. 
            }#end of the length(op)!=0 if-else. 
                
       if(ckins){#replace 3 elements in tmp sequence.  
         tmp<-array(0,c(length(tmp)))
         tmp[which(otmp==ins)]<- rept[loop-1,1] # replaced with hyd element.
         tmp[which(otmp == bmtp[loop])] = rept[loop,which(rept[loop,]!=ins)]
         tmp[which(otmp == bmtp[loop+1])] = rept[loop+1,which(rept[loop+1,]!=ins)]            
         #replace 3 new nodes.
         tx1<-which(otmp==ins)
         tx2<-which(otmp == bmtp[loop])
         tx3<-which(otmp == bmtp[loop+1])
         for(i in 1:length(tmp)){
           if(i != tx1 && i!=tx2 && i!=tx3){
             tmp[i]=otmp[i]
             }
           }
             
       otmp<-tmp      
       rownames(mtx)<-c(tmp)
       colnames(mtx)<-c(tmp)      
       mtx<-mtx+diag(branchlength[length(tmp)-1],c(length(tmp)) )             
       omtx<-mtx
       ckins<-FALSE
       loop<-loop+2          
       }#end of replace 3 elements
     }#end of while loop

   if(sum(tipnames%in%tmp)!=nleaves){#catches the last speciation event.
     tmp<-tipnames
     mtx<-matrix(0,nrow=length(tmp),ncol=length(tmp))
     rownames(mtx)<-c(tmp)
     colnames(mtx)<-c(tmp)
     mtx<-fillspcmtx()
     }#end of if (sum(tipnames%in%tmp)!=nleaves).
  return(mtx)
  }#end of cov_mtx 


NegLogLike<-function(bt){
  W<-cov_mtx(bt) 
  solveW<-solve(W)  
  mu<-c((t(one)%*%solveW%*%Y)/(t(one)%*%solveW%*%one))
  Y<-Y-mu
  sigma_sq<-c(abs((t(Y)%*%solveW%*%Y)/n))
  neglogML<- n/2*log(2*pi) +  n/2*log(sigma_sq)+1/(2*sigma_sq)*t(Y)%*%solveW%*%Y+1/2*log(abs(det(W)))
  return(neglogML)
  }#end of NegLogLike
    

    
MLE<-function(bt){#MLE estimators.
  W<-cov_mtx(bt)
  solveW<-solve(W) 
  mu<-c((t(one)%*%solveW%*%Y)/(t(one)%*%solveW%*%one))
  Y<-Y-mu
  sigma_sq <- c(abs((t(Y)%*%solveW%*%Y)/n))
  return(c(bt,mu,sigma_sq))
  }#end of MLE
 
chl_decom<-function(M){#M=LL^t, returns L^{-1}. 
  decom<-t(chol(M))
  return(solve(decom))
  }

BMhydcvt<-function(bt){
  Y<-(Y-mu_hat)/sqrt(sigmasq_hat)
  Y<-chl_decom(cov_mtx(bt))%*%Y
  return(Y)
  }

####################################################################
###################### MAIN PROGRAM ################################
####################################################################

#HYD 3 taxa
#x.tre=c("((1:0.6,(2:0.6)5:0)4:0.4,(5:0,3:0.6)6:0.4 )7:0;")
#Y<-runif(3)

#5 taxa 2 HYD in EX2  
#x.tre<-c(" ((1,((2,(3)7)6)10)9,((10,(7,4)8)11,5)12)13; ")
#Y<-rnorm(5)

#x.tre<-c(" ((1:0.35,((2:0.2,(3:0.2)7:0)6:0.15)10:0)9:0.65,((10:0,(7:0,4:0.2)8:0.15)11:0.2,5:0.55)12:0.45)13:0; ")
#Y<-rnorm(5)

#5 taxa  2 HYD events nongall EX 3 Example in Tex
#x.tre=c("((1:0.3,(2:0.3)10:0)9:0.7,((10:0,(3:0.1,(4:0.1)7:0)6:0.2)11:0.3,(7:0,5:0.1)8:0.5)12:0.4)13:0; ")
#Y<-rnorm(5)

#5 taxa 1 HYD event EX 1
#x.tre<-c("(((1:0.5,(2:0.5)7:0)6:0.1,(7:0,3:0.5)8:0.1)10:0.4,(4:0.3,5:0.3)9:0.7)11:0;")
#Y<-matrix(c(1.2,3.5,1.6,2.7,4.3),ncol=1)

#6 taxa 2 HYD 1 spc:  Algorithm example
#x.tre<-c( "((1:0.65,((2:0.3,3:0.3)7:0.35)12:0)11:0.35,((12:0,(4:0.5,(5:0.5)9:0)8:0.15)13:0.25,(9:0,6:0.5)10:0.4)14:0.1)15:0;")
#Y<-rnorm(6)

#8taxa 2HYD event  Paper data EX 4
#x.tre=c("((1:0.85,((2:0.35,(3:0.15,(4:0.15)10:0)9:0.2)12:0.25,(((10:0,5:0.15)11:0.25,6:0.4)13:0.2)15:0)14:0.25)17:0.15,((15:0,7:0.6)16:0.15,8:0.75)18:0.25)19:0;")
#Y<-rnorm(8)

# 9 taxa real data
x.tre=c("(((1:0.35,2:0.35)10:0.5,(3:0.7,(4:0.45,(5:0.1,(6:0.1)13:0)12:0.35)11:0.25)16:0.15)18:0.15,((13:0,7:0.1)14:0.55,(8:0.2,9:0.2)15:0.45)17:0.35)19:0;")
Y<-matrix(c( 16.00,12.00,16.00,23.00,18.00,18.00,25.00,20.00,14.00  ),ncol=1)

#22 taxa 
#x.tre<-c("((((((1:0.62,(((2:0.185,(5:0.185)28:0)27:0.185,(16:0.37)34:0)33:0.15,((3:0.16,4:0.16)26:0.29,((28:0,6:0.185)29:0.03,(7:0.095,(12:0.095)24:0)23:0.12)30:0.235)37:0.07)38:0.1)39:0.05,8:0.67)41:0.07,((9:0.64,10:0.64)40:0.055,11:0.695)42:0.045)43:0.105,(((24:0,13:0.095)25:0.195,14:0.29)31:0.505,(15:0.755,((34:0,17:0.37)35:0.06,(18:0.335,19:0.335)32:0.095)36:0.325)44:0.04)45:0.05)46:0.08,(20:0.89,21:0.89)47:0.035)48:0.075,22:1)49:0;")
#Y<-rnorm(22)

main<-function(x.tre,Y){
res<-newick2phylog(x.tre)
data<-getntn(res)
branchlength<-getbrnlen(res)###check more examples
tipnames<-sort(names(res$droot[which(res$droot==1)]))
nleaves<-length(tipnames)
n<-nleaves
one <- matrix(1,nrow=n)
SS<-optimize(NegLogLike,c(-100,100))
mle<-MLE(SS$minimum)
mu_hat<-mle[2]
sigmasq_hat<-mle[3]
print(mle)
}
