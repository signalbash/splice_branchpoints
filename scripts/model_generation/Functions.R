get_ppt=function(i,df,branchpoint_df_x){
  dist=branchpoint_df_x$dist.2[i]
  seq=substr(df[i,2], 251, (250+dist))
  seq_vect=unlist(strsplit(seq, ""))
  pyramidines=which(seq_vect =="T" | seq_vect=="C")
  if(length(pyramidines) >1){
    pyra_dist=vector()
    for(p in 2:length(pyramidines)){
      pyra_dist[p-1]=pyramidines[p]-pyramidines[p-1]
    }
    longest_rn=0
    best_rn=0
    best_start=0
    percent_pyra=0
    start=1
    startv=vector()
    rnv=vector()
    polyperv=vector()
    
    for(p in 1:length(pyra_dist)){
      if(pyra_dist[p]==1 & longest_rn==0){
        longest_rn=1
        start=p
        seq_p=seq_vect[(pyramidines[start]):(pyramidines[start+longest_rn])]
        percent_pyra=length(which(seq_p=="C"|seq_p=="T"))/length(seq_p)
        
        if(longest_rn > best_rn){
          best_rn=longest_rn
          best_start=start
          best_percent=percent_pyra
        }
      }else if(pyra_dist[p]==1 | pyra_dist[p]==2){
        longest_rn=longest_rn+1
        seq_p=seq_vect[(pyramidines[start]):(pyramidines[start+longest_rn])]
        percent_pyra=length(which(seq_p=="C"|seq_p=="T"))/length(seq_p)
        
        if(longest_rn > best_rn){
          best_rn=longest_rn
          best_start=start
          best_percent=percent_pyra
        }
      }else{
        startv=append(startv, best_start)
        rnv=append(rnv, longest_rn)
        polyperv=append(polyperv, percent_pyra)
        longest_rn=0
        start=p+1
        best_start=start
      }
    }
    startv=append(startv, start)
    rnv=append(rnv, longest_rn)
    polyperv=append(polyperv, percent_pyra)
    
    
    pyra_df=data.frame(set=c(1:length(rnv)), rnv,startv,polyperv)
    pyra_df=arrange(pyra_df, plyr::desc(rnv))
    
    best_start=pyra_df$startv[which.max(pyra_df$rnv)]
    best_run=pyra_df$rnv[which.max(pyra_df$rnv)]
    best_polyper=pyra_df$polyperv[which.max(pyra_df$rnv)]
    
    pyra_df=pyra_df[pyra_df$rnv >=10,]
    if(dim(pyra_df)[1]!=0){
      run_py <- vector()
      ps_py <- vector()
      polyper_py <- vector()
      for(py in 1:length(pyra_df[,1])){
        seq_p=seq_vect[(pyramidines[pyra_df$startv[py]]):(pyramidines[pyra_df$startv[py]+pyra_df$rnv[py]])]
        startp=1
        
        run=length(seq_p)
        
        polyper_py[py] <-pyra_df$polyperv[py]
        while(polyper_py[py] < 0.8 & run >=10){
          pur=which(seq_p=="A"| seq_p=="G")
          purn=(length(seq_p)+1-pur)*-1
          pur2=c(pur,purn)
          minp=which.min(abs(pur2))
          
          point=pur2[minp]
          if(point < 0){
            rm=length(seq_p)-point*-1
            new_seq=seq_p[1:rm]
            run=length(new_seq)
          }else{
            new_seq=seq_p[(point+1):length(seq_p)]
            run=length(new_seq)
            startp=startp+point
          }
          polyper_py[py]=length(which(new_seq=="C"|new_seq=="T"))/length(new_seq)
          seq_p=new_seq
        }
        
        ps=((pyramidines[pyra_df$startv[py]])+startp-1)
        #tract=seq_vect[ps:(ps+run-1)]
        run_py[py]=run
        ps_py[py]=ps
      }
      
      best_start=ps_py[which.max(run_py)]
      best_run=run_py[which.max(run_py)]
      best_polyper=polyper_py[which.max(run_py)]
      
    }
    line=c(best_start,best_run,best_polyper)
  }else if(length(pyramidines)==1){
    line=c(pyramidines, 1,1)
  }else{
    line=c(0,0,0)
  }
  return(line)
}
library(stringr)