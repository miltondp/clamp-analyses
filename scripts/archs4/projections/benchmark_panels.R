# Draw projection benchmark panels for Supplementary Figure 6.
projection_benchmark_panels <- function(root) {
  method_colours <- c('CLAMPfull'='#3e348b','CLAMPbase'='#0072B2','PLIER'='#EE7733',
                      'ARCHS4 projection'='#923155','RNA-Seq'='#9c9f36',
                      'ARCHS4 projection null'='#999999')
  plot_comparisons <- function(d, pairs, key, paired=TRUE) {
    ans <- rbindlist(lapply(pairs,function(pair) {
      left <- d[method==pair[1],c(key,'value'),with=FALSE]
      right <- d[method==pair[2],c(key,'value'),with=FALSE]
      setnames(left,'value','left_value');setnames(right,'value','right_value')
      aligned <- merge(left,right,by=key)
      stopifnot(nrow(aligned)==nrow(left),nrow(aligned)==nrow(right))
      data.table(left=pair[1],right=pair[2],n=nrow(aligned),p=wilcox.test(aligned$left_value,aligned$right_value,paired=paired,exact=FALSE)$p.value)
    }))
    ans[,q:=p.adjust(p,'BH')]
    ans[,label:=vapply(q,function(value) {
      if(value>=.001)return(sprintf('italic(q) == "%.3f"',value))
      parts <- strsplit(formatC(value,format='e',digits=2),'e')[[1]]
      sprintf('italic(q) == %s %%*%% 10^displaystyle(%d)',trimws(parts[1]),as.integer(parts[2]))
    },character(1))]
    ans
  }
  distribution <- function(d,title,metric,pairs,key,out_file=NULL,levels_order=NULL,colours=method_colours,paired=TRUE) {
    stopifnot(!anyDuplicated(d[,c('method',key),with=FALSE]),all(is.finite(d$value)))
    means <- d[,.(mean_value=mean(value),max_value=max(value)),by=method][order(-mean_value)]
    levels <- if(is.null(levels_order)) means$method else levels_order
    comparisons <- plot_comparisons(d,pairs,key,paired)
    fwrite(comparisons,file.path(root,if(!is.null(out_file))out_file else if('rep'%in%key)'gtex_paired_comparisons.tsv' else 'pseudobulk_paired_comparisons.tsv'),sep='\t')
    d[,method:=factor(method,levels=levels)];means[,method:=factor(method,levels=levels)]
    ymax <- max(d$value);ymin <- min(d$value)-diff(range(d$value))*.06
    yrange <- ymax-ymin;if(yrange==0)yrange<-1
    means[,label_y:=max_value+yrange*.04]
    bracket_start <- .16
    bracket_step <- .19
    comparisons[,`:=`(x=match(left,levels),xend=match(right,levels))]


    comparisons[,`:=`(lo=pmin(x,xend),hi=pmax(x,xend),mid=(x+xend)/2)]
    comparisons[,half:=ifelse(grepl('10\\^',label),1.0,.65)]
    comparisons[,`:=`(ext_lo=pmin(lo,mid-half),ext_hi=pmax(hi,mid+half))]
    comparisons[,len:=hi-lo]
    setorder(comparisons,len,lo)
    comparisons[,level:=0L]
    for(i in seq_len(nrow(comparisons))) {
      lvl <- 1L
      repeat {
        taken <- comparisons[seq_len(i-1L)][level==lvl & ext_lo<comparisons$ext_hi[i] & ext_hi>comparisons$ext_lo[i]]
        if(!nrow(taken))break
        lvl <- lvl+1L
      }
      set(comparisons,i,'level',lvl)
    }
    comparisons[,y:=ymax+yrange*(bracket_start+bracket_step*(level-1L))]
    comparisons[,label_y:=y+yrange*.04]
    ggplot(d,aes(method,value,fill=method))+
      geom_boxplot(width=.56,outlier.shape=NA,linewidth=.3,alpha=.8)+
      geom_point(position=position_jitter(width=.2,height=0,seed=123),shape=21,size=.6,fill='white',stroke=.2)+
      geom_point(data=means,aes(method,mean_value),shape=23,size=1.6,fill='white',colour='black',inherit.aes=FALSE)+
      geom_text(data=means,aes(method,label_y,label=sprintf('%.3f',mean_value)),vjust=0,size=5/.pt,inherit.aes=FALSE)+
      geom_segment(data=comparisons,aes(x=x,xend=xend,y=y,yend=y),inherit.aes=FALSE,linewidth=.2)+
      geom_segment(data=comparisons,aes(x=x,xend=x,y=y-yrange*.025,yend=y),inherit.aes=FALSE,linewidth=.2)+
      geom_segment(data=comparisons,aes(x=xend,xend=xend,y=y-yrange*.025,yend=y),inherit.aes=FALSE,linewidth=.2)+
      geom_text(data=comparisons,aes(x=(x+xend)/2,y=label_y,label=label),inherit.aes=FALSE,size=5/.pt,vjust=0,parse=TRUE)+
      scale_fill_manual(values=colours,guide='none')+
      labs(x=NULL,y=metric)+nm_theme()+
      theme(axis.text.x=element_text(angle=45,hjust=1,size=5),plot.title=element_text(size=5.5,hjust=.5),plot.margin=margin(2,2,2,2,'mm'))+
      coord_cartesian(ylim=c(ymin,max(comparisons$label_y)+yrange*.07),clip='off')
  }
  gt <- fread(file.path(root,'gtex_full_ari.tsv'))
  gt578 <- fread(file.path(root,'gtex_sensitivity_repetitions.tsv'))[k==578 & subset==rep,.(method='ARCHS4 projection',rep,ari)]
  gt <- rbind(gt[method!='ARCHS4 projection'],gt578);setnames(gt,'ari','value')
  pb <- fread(file.path(root,'pseudobulk_adjusted.tsv'));setnames(pb,'excess_r','value')
  gt_plot <- distribution(gt,'GTEx, 578 projected GMs','ARI',
      lapply(c('CLAMPfull','CLAMPbase','RNA-Seq','PLIER','ARCHS4 projection null'),function(m)c('ARCHS4 projection',m)), 'rep', paired=FALSE)
  pb_plot <- distribution(pb,'Pseudobulk, 47 targets','Adjusted max Pearson r',
      lapply(c('CLAMPfull','CLAMPbase','PLIER','ARCHS4 projection null'),function(m)c('ARCHS4 projection',m)),c('dataset','cell_type'))
  gt_plot <- gt_plot+scale_y_continuous(breaks=c(.6,.7,.8))
  list(gtex=gt_plot,pseudobulk=pb_plot)
}
