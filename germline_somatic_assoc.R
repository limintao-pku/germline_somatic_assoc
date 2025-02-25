#Load some packages
#If any of the following packages are not available, please install them before proceeding
#The installation of deconstructSigs and BSgenome.Hsapiens.UCSC.hg19 may require BiocManager::install
library(data.table)
library(stringr)
requireNamespace("readxl")
requireNamespace("perm")
requireNamespace("deconstructSigs")
requireNamespace("BSgenome.Hsapiens.UCSC.hg19")
requireNamespace("survival")

#Some user-defined functions
#identical_use_cl and match.stage will be used to judge whether two records correspond to the same TCGA patient based on clinical data
#get_cn_index will be used to obtain certain CNV indexes
#perm_test_cont will be used to perform the permutation test on continuous indexes
#get_sbs_pattern, get_sbs_avg, get_sig_pattern, and perm_test_sig will be used in the mutational signature analysis
match.stage<-function(stage1,stage2,exact=F){
  if(is.na(stage1)|is.na(stage2)){return(T)}
  if(!exact){
    stage1<-str_remove_all(stage1,"(Stage )|A|B|C")
    stage2<-str_remove_all(stage2,"(Stage )|A|B|C")
  }
  stage1==stage2
}
identical_use_cl<-function(sex1,sex2,cancer_type1,cancer_type2,
                           stage1,stage2,match_stage=T,match_exact=F,
                           days_to_birth1,days_to_birth2,days_to_birth_tol=1,
                           age_at_diagnosis1,age_at_diagnosis2,age_at_diagnosis_tol=1,
                           year_of_diagnosis1,year_of_diagnosis2,year_of_diagnosis_tol=0){
  if(sex1!=sex2){return(F)}
  if(cancer_type1!=cancer_type2){return(F)}
  if(match_stage){
    if(!match.stage(stage1,stage2,match_exact)){
      return(F)
    }
  }
  if(abs(days_to_birth1-days_to_birth2)>days_to_birth_tol){return(F)}
  if(abs(age_at_diagnosis1-age_at_diagnosis2)>age_at_diagnosis_tol){return(F)}
  if(abs(year_of_diagnosis1-year_of_diagnosis2)>year_of_diagnosis_tol){return(F)}
  return(T)
}
get_cn_index<-function(f_id,path,erbb2=F){
  gcnc_avg<-function(x){
    #Average gene copy number change
    x<-na.omit(x)
    mean(abs(x-mean(x)))
  }
  
  f_path<-paste(path,f_id,sep="/")
  n<-length(f_path)
  out<-rep(NaN,n)
  for(i in 1:n){
    f_name_i<-dir(f_path[i])
    f_name_i<-f_name_i[str_detect(f_name_i,"copy_number")]
    stopifnot(length(f_name_i)==1)
    
    path_i<-paste(f_path[i],f_name_i,sep="/")
    file_i<-fread(path_i)
    
    if(erbb2){
      out[i]<-file_i$copy_number[file_i$gene_name=="ERBB2"]
    }else{
      cn_autosome<-file_i$copy_number[!file_i$chromosome%in%c("chrX","chrY")]
      out[i]<-gcnc_avg(cn_autosome)
    }
    
    if(i==n){
      cat("iter",i,"/",n,"\r\n")
    }else{
      cat("iter",i,"/",n,"\r")
    }
  }
  out
}
perm_test_cont<-function(withgm_id,nogm_id,data,num_sim=5e5,stat="mean",na.rm=F,...){
  if(na.rm){
    id<-data$sampleID[!is.na(data$value)]
    withgm_id<-withgm_id[withgm_id%in%id]
    nogm_id<-nogm_id[nogm_id%in%id]
    stopifnot(length(withgm_id)>0&length(nogm_id)>0)
  }
  x<-rep.int(c(1,0),times=c(length(withgm_id),length(nogm_id)))
  id<-c(withgm_id,nogm_id)
  y<-data$value[match(id,data$sampleID)]
  stopifnot(!anyNA(y))
  
  y0<-y
  if(stat!="mean"){y<-rank(y)}
  p_perm<-perm::permTS(y[x==1],y[x==0],method="exact.mc",control=perm::permControl(nmc=num_sim-1,setSEED=F,tsmethod="abs"),...)$p.value
  y<-y0
  
  v1<-y[x==1]
  v2<-y[x==0]
  
  list(
    p.value=p_perm,
    dat_sum=list(
      group1=c(as.list(summary(v1)),list(sd=sd(v1))),
      group2=c(as.list(summary(v2)),list(sd=sd(v2)))
    )
  )
}
get_sbs_pattern<-function(id,mut_dat,count_lower=1){
  stopifnot(!anyDuplicated(id))
  n<-length(id)
  out_list<-list()
  for(i in 1:n){
    if(i==n){
      cat("iter",i,"/",n,"\r\n")
    }else{
      cat("iter",i,"/",n,"\r")
    }
    
    dat_i<-mut_dat[mut_dat$sampleID==id[i],]
    if(nrow(dat_i)<count_lower){next}
    
    sigs.input<-suppressWarnings(
      mut.to.sigs.input(
        mut.ref=dat_i, 
        sample.id="sample_one", 
        chr="chr", 
        pos="pos", 
        ref="ref", 
        alt="alt"
      )
    )
    
    list_i<-list(as.numeric(sigs.input))
    names(list_i)<-id[i]
    out_list<-c(out_list,list_i)
  }
  out_list
}
get_sbs_avg<-function(sbs_list,lower=1,upper=Inf,cut=0,k=1){
  x<-rep(0,96)
  for(i in 1:length(sbs_list)){
    x_i<-sbs_list[[i]]
    s_x_i<-sum(x_i)
    if(s_x_i>=lower&s_x_i<=upper){
      if(s_x_i>=cut){
        x<-x+x_i/s_x_i*k
      }else{
        x<-x+x_i/s_x_i
      }
    }
  }
  x/sum(x)
}
get_sig_pattern<-function(sigs.input,signatures.ref=signatures.cosmic,scale=NULL,...){
  if(!is.null(scale)){
    sigs.input<-sigs.input*scale
    sigs.input<-sigs.input/sum(sigs.input)
  }
  dim(sigs.input)<-c(1,96)
  colnames(sigs.input)<-colnames(signatures.cosmic)
  rownames(sigs.input)<-"1"
  sigs.input<-as.data.frame(sigs.input)
  sig_out<-whichSignatures(tumor.ref=sigs.input, 
                           signatures.ref=signatures.ref,
                           signature.cutoff=0,contexts.needed=T,
                           sample.id=1,...)
  weights<-as.matrix(sig_out$weights)
  
  weights/sum(weights)
}
perm_test_sig<-function(sbs_list1,sbs_list2,n_cores=30,n_each=335,
                        signatures.ref=signatures.cosmic,scale=NULL,...){
  n1<-length(sbs_list1)
  n2<-length(sbs_list2)
  n_all<-n1+n2
  loc<-1:n_all
  sbs_list<-c(sbs_list1,sbs_list2)
  
  idx1<-get_sig_pattern(get_sbs_avg(sbs_list1,k=1),signatures.ref=signatures.ref,scale=scale,...)
  idx2<-get_sig_pattern(get_sbs_avg(sbs_list2,k=1),signatures.ref=signatures.ref,scale=scale,...)
  diff_obs<-as.numeric(idx2-idx1)
  names(diff_obs)<-rownames(signatures.ref)
  
  my_task<-function(n_perm){
    t0<-Sys.time();pid<-0
    cp<-seq(0,1,by=0.2);progress_previous<-0
    diff_perm<-matrix(NaN,nrow=nrow(signatures.ref),ncol=n_perm)
    for(i in 1:n_perm){
      progress_now<-cp[rev(which(i/n_perm>=cp))[1]]
      if(progress_now>progress_previous){
        system(sprintf("echo \"%s\"",paste("Worker",pid,"|",paste0(progress_now*100,"%"),"|",capture.output(Sys.time()-t0))))
      }
      progress_previous<-progress_now
      
      loc2<-sample.int(n_all,n2)
      idx1<-get_sig_pattern(get_sbs_avg(sbs_list[-loc2],k=1),signatures.ref=signatures.ref,scale=scale,...)
      idx2<-get_sig_pattern(get_sbs_avg(sbs_list[loc2],k=1),signatures.ref=signatures.ref,scale=scale,...)
      diff_perm[,i]<-idx2-idx1
    }
    diff_perm
  }
  
  out_prl<-parallel::mclapply(rep(n_each,n_cores),my_task,mc.cores=n_cores)
  
  diff_perm<-matrix(NaN,nrow=nrow(signatures.ref),ncol=n_cores*n_each)
  for(i in 1:n_cores){
    diff_perm[,((i-1)*n_each+1):(i*n_each)]<-out_prl[[i]]
  }
  rownames(diff_perm)<-rownames(signatures.ref)
  
  return(list(diff_obs=diff_obs,diff_perm=diff_perm))
}

#Prepare the data of P/LP GMs according to https://doi.org/10.1016/j.cell.2018.03.039
##file_path1 here should point to the file of Table S2 by Huang et al.
file_path1<-"./1-s2.0-S0092867418303635-mmc2.xlsx";stopifnot(file.exists(file_path1))
if(T){
  ###Considering the sample size, both Pathogenic_variants and truncating Prioritized_VUSs will be treated as P/LP
  cell_plp_gm<-readxl::read_excel(file_path1,sheet=1)
  cell_pvus_gm<-readxl::read_excel(file_path1,sheet=2)
  cell_plp_gm<-as.data.frame(cell_plp_gm);cell_pvus_gm<-as.data.frame(cell_pvus_gm)
  dim(cell_plp_gm);dim(cell_pvus_gm)
  identical(colnames(cell_plp_gm),colnames(cell_pvus_gm))#true
  
  table(cell_pvus_gm$Variant_Classification)
  loc<-cell_pvus_gm$Variant_Classification%in%c("frameshift_variant","stop_gained")
  cell_gm<-rbind(cell_plp_gm,cell_pvus_gm[loc,])
  table(cell_gm$cancer)
  cell_gm<-cell_gm[cell_gm$cancer%in%c("COAD","READ"),]
  dim(cell_gm)#51, 173
  
  ###Gender
  cell_gm$gender<-tolower(cell_gm$gender)
  
  ###Stage
  table(cell_gm$ajcc_pathologic_tumor_stage);anyNA(cell_gm$ajcc_pathologic_tumor_stage)
  loc<-which(!str_detect(cell_gm$ajcc_pathologic_tumor_stage,"Stage"))
  cell_gm$ajcc_pathologic_tumor_stage[loc]<-NA_character_
  table(cell_gm$ajcc_pathologic_tumor_stage);sum(is.na(cell_gm$ajcc_pathologic_tumor_stage))
  
  ###Other variables required for matching
  suppressWarnings({
    cell_gm$birth_days_to<-as.numeric(cell_gm$birth_days_to)
    cell_gm$age_at_initial_pathologic_diagnosis<-as.numeric(cell_gm$age_at_initial_pathologic_diagnosis)
    cell_gm$initial_pathologic_dx_year<-as.numeric(cell_gm$initial_pathologic_dx_year)
  })
  loc<-is.na(cell_gm$birth_days_to)|is.na(cell_gm$age_at_initial_pathologic_diagnosis)|is.na(cell_gm$initial_pathologic_dx_year)
  cell_gm<-cell_gm[!loc,]
  dim(cell_gm)#51, 173
}

#TCGA clinical and somatic data folder
file_folder1<-"path to TCGA clinical and somatic data folder"
#Please prepare the following data
# file_folder1
#      |___COAD&READ
#            |___clinical
#                  |___clinical.tsv ##https://portal.gdc.cancer.gov/ (obtained together with ASCAT2 gene-level copy number data)
#            |___CNV
#                  |___xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx ##https://portal.gdc.cancer.gov/ (ASCAT2 gene-level copy number data folders)
#                  |___... ##a total of 680 file folders
#            |___ERBB2
#                  |___data_cna.txt ##https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018
#            |___sig&survival
#                  |___Colorectum_clean_somatic_mutations_for_signature_analysis.txt ##https://ftp.sanger.ac.uk/pub/cancer/AlexandrovEtAl/somatic_mutation_data/Colorectum/
#                  |___data_clinical_patient.txt ##https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018
#                  |___data_clinical_sample.txt ##https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018
#                  |___data_mutations.txt ##https://www.cbioportal.org/study/summary?id=coadread_tcga_pan_can_atlas_2018
#            |___gdc_sample_sheet.year-month-day.tsv ##https://portal.gdc.cancer.gov/
#
#      |___ESCA
#            |___clinical
#                  |___clinical.tsv ##https://portal.gdc.cancer.gov/ (obtained together with ASCAT2 gene-level copy number data)
#            |___CNV
#                  |___xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx ##https://portal.gdc.cancer.gov/ (ASCAT2 gene-level copy number data folders)
#                  |___... ##a total of 183 file folders
#            |___ERBB2
#                  |___data_cna.txt ##https://www.cbioportal.org/study/summary?id=esca_tcga_pan_can_atlas_2018
#            |___gdc_sample_sheet.year-month-day.tsv ##https://portal.gdc.cancer.gov/
#
#      |___STAD
#            |___clinical
#                  |___clinical.tsv ##https://portal.gdc.cancer.gov/ (obtained together with ASCAT2 gene-level copy number data)
#            |___CNV
#                  |___xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx ##https://portal.gdc.cancer.gov/ (ASCAT2 gene-level copy number data folders)
#                  |___... ##a total of 435 file folders
#            |___ERBB2
#                  |___data_cna.txt ##https://www.cbioportal.org/study/summary?id=stad_tcga_pan_can_atlas_2018
#            |___gdc_sample_sheet.year-month-day.tsv ##https://portal.gdc.cancer.gov/
#
#      |___perm_out.Rdata ##optional; downloaded from the GitHub repository; you can generate the perm_out data manually with this script
#
#If you have troubles in obtaining these data, you can contact the corresponding author or the first author

#Prepare the tumor gene-level copy number data and clinical data of the TCGA CRC patients
#ASCAT2 is selected as it results in minimal sample loss
file_path2<-paste(file_folder1,"COAD&READ","clinical/clinical.tsv",sep="/");stopifnot(file.exists(file_path2))
file_path3<-paste(file_folder1,"COAD&READ","CNV",sep="/");length(dir(file_path3))==681
files1<-dir(paste(file_folder1,"COAD&READ",sep="/"))
files1<-files1[str_detect(files1,"gdc_sample_sheet")];length(files1)==1
file_path4<-paste(file_folder1,"COAD&READ",files1,sep="/")
if(T){
  tcga_cl<-fread(file_path2)
  tcga_cl<-as.data.frame(tcga_cl)
  dim(tcga_cl)
  
  tcga_cl<-tcga_cl[!duplicated(tcga_cl$case_id),]
  dim(tcga_cl)#621, 197
  
  loc<-tcga_cl$days_to_birth=="'--"|tcga_cl$age_at_diagnosis=="'--"|tcga_cl$year_of_diagnosis=="'--"
  tcga_cl<-tcga_cl[!loc,]
  dim(tcga_cl)#517, 197
  
  table(tcga_cl$gender)
  table(tcga_cl$project_id)
  
  ##To enhance the contrast between the TCGA and our cohorts, TCGA patients with a reported Asian race are excluded
  table(tcga_cl$race);anyNA(tcga_cl$race)
  tcga_cl<-tcga_cl[tcga_cl$race!="asian",]
  dim(tcga_cl)#506, 197
  
  table(tcga_cl$ajcc_pathologic_stage);anyNA(tcga_cl$ajcc_pathologic_stage)
  loc<-tcga_cl$ajcc_pathologic_stage=="'--"|tcga_cl$ajcc_pathologic_stage=="Unknown"
  tcga_cl$ajcc_pathologic_stage[loc]<-NA_character_
  
  tcga_cl$days_to_birth<-as.numeric(tcga_cl$days_to_birth)
  tcga_cl$age_at_diagnosis<-as.numeric(tcga_cl$age_at_diagnosis)
  tcga_cl$year_of_diagnosis<-as.numeric(tcga_cl$year_of_diagnosis)
  tcga_cl$project_id<-str_remove(tcga_cl$project_id,"TCGA-")
}

#Match the samples from tcga_cl (clinical data) to cell_gm (GM data)
cell_gm$sampleID<-NA_character_
if(T){
  for(i in 1:nrow(cell_gm)){
    res_i<-logical(nrow(tcga_cl))
    for(j in 1:nrow(tcga_cl)){
      res_i[j]<-identical_use_cl(
        tcga_cl$gender[j],cell_gm$gender[i],
        tcga_cl$project_id[j],cell_gm$cancer[i],
        tcga_cl$ajcc_pathologic_stage[j],cell_gm$ajcc_pathologic_tumor_stage[i],T,F,#Removing 'stage' does not lead to more matched patients 
        tcga_cl$days_to_birth[j],cell_gm$birth_days_to[i],days_to_birth_tol=0,
        tcga_cl$age_at_diagnosis[j]/365.25,cell_gm$age_at_initial_pathologic_diagnosis[i],age_at_diagnosis_tol=1,
        tcga_cl$year_of_diagnosis[j],cell_gm$initial_pathologic_dx_year[i],year_of_diagnosis_tol=0
      )
    }
    if(any(res_i)){
      loc<-which(res_i)
      if(length(loc)>1){message("multi-hit: ",i)}
      cell_gm$sampleID[i]<-tcga_cl$case_submitter_id[loc[1]]
    }
  }
  cell_gm$sampleID
  cell_gm<-cell_gm[!is.na(cell_gm$sampleID),]
  nrow(cell_gm)#37
}

#P/LP GM distribution
if(T){
  gene_list<-list(
    MMR=c('EXO1','HMGB1','LIG1','MLH1','MLH3','MSH2','MSH3','MSH4','MSH5','MSH6','PCNA','PMS1','PMS2',
          'POLD1','POLD2','POLD3','POLD4','RFC1','RFC2','RFC3','RFC4','RFC5','RPA1','RPA2','RPA3'),
    HRR=c('BLM','BRCA1','BRCA2','DMC1','EME1','EME2','GEN1','HFM1','MRE11A','MUS81','NBN','PPP4C',
          'PPP4R1','PPP4R2','PPP4R4','RAD50','RAD51','RAD51B','RAD51C','RAD51D','RAD52','RAD54B',
          'RAD54L','RAD54L2','RDM1','RECQL','RECQL4','RECQL5','RMI1','RMI2','RPA1','RPA2','RPA3',
          'SEM1','SLX1A','SLX4','PPP4R3A','PPP4R3B','SPO11','TOP3A','TOP3B','WRN','XRCC2','XRCC3'),
    FA=c('CENPS','BLM','BRCA1','BRCA2','BRIP1','FAAP100','FAAP24','FAN1','FANCA','FANCC','FANCD2',
         'FANCE','FANCF','FANCG','FANCI','FANCL','FANCM','HES1','PALB2','RAD51','RAD51C','RMI1',
         'RMI2','CENPX','TELO2','TOP3A','TOP3B','UBE2T','USP1','WDR48'),
    CPF=c('AEN','ATM','ATR','ATRIP','CHEK1','CHEK2','HUS1','HUS1B','PER1','PER2','PER3','RAD1','RAD17',
          'RAD9A','RAD9B','RFC2','RFC3','RFC4','RFC5','TIMELESS','TIPIN','TP53'),
    BER=c('APEX1','APLF','APTX','CCNO','FEN1','HMGB1','LIG1','LIG3','MBD4','MPG','MUTYH','NEIL1',
          'NEIL2','NEIL3','NTHL1','OGG1','PARP1','PARP2','PARP3','PARP4','PCNA','PNKP','POLB','POLD1',
          'POLD2','POLD3','POLD4','POLE','POLE2','POLE3','POLE4','POLL','SMUG1','TDG','TDP1','UNG','XRCC1'),
    RecQ=c('RECQL','WRN','BLM','RECQL4','RECQL5')
  )
  
  ##Obtain cell_gm$path
  cell_gm$path<-NA_character_
  for(i in 1:length(cell_gm$HUGO_Symbol)){
    g_i<-cell_gm$HUGO_Symbol[i]
    p_i<-rep(NA,5)
    for(j in 1:5){
      p_i[j]<-g_i%in%gene_list[[j]]
    }
    p_i<-paste(names(gene_list)[which(p_i)],collapse=", ")
    if(p_i==""){
      p_i<-"Other"
    }
    cell_gm$path[i]<-p_i
  }
  table(cell_gm$path)
  cell_gm$path<-factor(cell_gm$path,levels=c("MMR","HRR","HRR, FA","FA","CPF","Other"))
  
  ##Obtain cell_gm$Variant_Classification2
  table(cell_gm$Variant_Classification)
  cell_gm$Variant_Classification2<-c("frameshift","missense","splice","nonsense")[
    match(
      cell_gm$Variant_Classification,
      c("frameshift_variant","missense_variant","splice_donor_variant","stop_gained")
    )
  ]
  table(cell_gm$Variant_Classification2);anyNA(cell_gm$Variant_Classification2)#false
  cell_gm$Variant_Classification2<-factor(cell_gm$Variant_Classification2,
                                          levels=c("frameshift","nonsense","splice","missense"))
  
  ##Obtain cell_gm$gene_name
  x<-names(sort(table(cell_gm$HUGO_Symbol),decreasing=T))
  cell_gm$gene_name<-factor(cell_gm$HUGO_Symbol,levels=x)
  
  ##GM distribution
  dist_mat<-as.matrix(table(cell_gm$gene_name,cell_gm$Variant_Classification2))
  t(dist_mat)
  table(cell_gm$path)
  
  ##Codes for plotting are not included in this script
}

#Use gdc_sample_sheet to locate the CNV file for each patient
if(T){
  ss<-fread(file_path4)
  dim(ss)#681, 8
  
  ##Remove the comma
  ss$`Case ID`[1:5]
  ss$`Case ID`<-vapply(strsplit(ss$`Case ID`,", ",fixed=T),function(x){x[1]},"")
  all(tcga_cl$case_submitter_id%in%ss$`Case ID`)#true
  
  ##Obtain file id
  tcga_cl$f_id<-ss$`File ID`[match(tcga_cl$case_submitter_id,ss$`Case ID`)]
  anyNA(tcga_cl$f_id)#false
  anyDuplicated(tcga_cl$f_id)>0#false
}

#Calculate the average gene copy number change (sex chromosomes are excluded) and inferred ERBB2 CN
tcga_cl$gcnc_avg<-get_cn_index(tcga_cl$f_id,file_path3,erbb2=F)
tcga_cl$erbb2_cn<-get_cn_index(tcga_cl$f_id,file_path3,erbb2=T)

#Compare the average gene copy number changes between different GM groups
#Comparisons of the TMB between groups will be presented later 
if(T){
  dat0<-data.frame(value=tcga_cl$gcnc_avg,sampleID=tcga_cl$case_submitter_id)
  all(cell_gm$sampleID%in%tcga_cl$case_submitter_id)#true
  
  id_list<-list(
    MMR=unique(cell_gm$sampleID[cell_gm$HUGO_Symbol%in%gene_list$MMR]),
    HRR=unique(cell_gm$sampleID[cell_gm$HUGO_Symbol%in%gene_list$HRR]),
    FA=unique(cell_gm$sampleID[cell_gm$HUGO_Symbol%in%gene_list$FA]),
    CPF=unique(cell_gm$sampleID[cell_gm$HUGO_Symbol%in%gene_list$CPF]),
    BER=unique(cell_gm$sampleID[cell_gm$HUGO_Symbol%in%gene_list$BER]),
    RecQ=unique(cell_gm$sampleID[cell_gm$HUGO_Symbol%in%gene_list$RecQ]),
    GM=unique(cell_gm$sampleID),
    noGM=tcga_cl$case_submitter_id[!tcga_cl$case_submitter_id%in%cell_gm$sampleID]
  )
  
  ##With GM; the permutation test generally finishes within 45s
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20250206)
  perm_test_cont(withgm_id=id_list$GM,
                 nogm_id=id_list$noGM,
                 data=dat0,
                 num_sim=5e5,stat="mean",alternative="two.sided")
  ##0.84699
  
  ##MMR
  set.seed(20250206)
  perm_test_cont(withgm_id=id_list$MMR,
                 nogm_id=id_list$noGM,
                 data=dat0,
                 num_sim=5e5,stat="mean",alternative="less")
  ##0.005036
  
  ##HRR
  set.seed(20250206)
  perm_test_cont(withgm_id=id_list$HRR,
                 nogm_id=id_list$noGM,
                 data=dat0,
                 num_sim=5e5,stat="mean",alternative="greater")
  ##0.018298
  
  ##RecQ
  set.seed(20250206)
  perm_test_cont(withgm_id=id_list$RecQ,
                 nogm_id=id_list$noGM,
                 data=dat0,
                 num_sim=5e5,stat="mean",alternative="greater")
  ##0.003696
}

#Mutational signature analysis, correlation analysis, and survival analysis
file_folder2<-paste(file_folder1,"COAD&READ","sig&survival",sep="/")
file_path5<-paste(file_folder2,"Colorectum_clean_somatic_mutations_for_signature_analysis.txt",sep="/");stopifnot(file.exists(file_path5))
file_path6<-paste(file_folder2,"data_mutations.txt",sep="/");stopifnot(file.exists(file_path6))
file_path7<-paste(file_folder2,"data_clinical_patient.txt",sep="/");stopifnot(file.exists(file_path7))
file_path8<-paste(file_folder2,"data_clinical_sample.txt",sep="/");stopifnot(file.exists(file_path8))
##Compare signature structures between different GM groups
if(T){
  ###dat_AlexandrovEtAl and sbs_AlexandrovEtAl
  dat_AlexandrovEtAl<-fread(file_path5);head(dat_AlexandrovEtAl)
  table(dat_AlexandrovEtAl$V2)
  dat_AlexandrovEtAl<-dat_AlexandrovEtAl[dat_AlexandrovEtAl$V2=="subs",]
  summary(dat_AlexandrovEtAl$V5-dat_AlexandrovEtAl$V4)
  
  dat_AlexandrovEtAl2<-data.frame(
    sampleID=dat_AlexandrovEtAl$V1,
    ref=dat_AlexandrovEtAl$V6,
    alt=dat_AlexandrovEtAl$V7,
    pos=dat_AlexandrovEtAl$V5,
    chr=paste0("chr",dat_AlexandrovEtAl$V3),
    sample_one=1
  )
  dim(dat_AlexandrovEtAl2)#204785, 6
  
  ###Get the SBS pattern for each patient; may take several minutes
  library(deconstructSigs)
  sbs_AlexandrovEtAl<-get_sbs_pattern(unique(dat_AlexandrovEtAl2$sampleID),dat_AlexandrovEtAl2)
  
  ###tcga_sm and sbs_tcga
  tcga_sm<-fread(file_path6)
  tcga_sm$sampleID<-substr(tcga_sm$Tumor_Sample_Barcode,1,12)
  
  tcga_sm2<-data.frame(
    sampleID=tcga_sm$sampleID,
    ref=tcga_sm$Reference_Allele,
    alt=tcga_sm$Tumor_Seq_Allele2,
    pos=tcga_sm$End_Position,
    chr=paste0("chr",tcga_sm$Chromosome),
    type=tcga_sm$Variant_Type,
    sample_one=1
  )
  tcga_sm2<-tcga_sm2[tcga_sm2$type=="SNP",]
  tcga_sm2$var_name<-paste(tcga_sm2$chr,tcga_sm2$pos,tcga_sm2$ref,tcga_sm2$alt,tcga_sm2$sampleID,sep=",")
  tcga_sm2<-tcga_sm2[!duplicated(tcga_sm2$var_name),]
  dim(tcga_sm2)#293562, 8
  
  ###Get the SBS pattern for each patient; may take several minutes
  sbs_tcga<-get_sbs_pattern(unique(tcga_sm2$sampleID),tcga_sm2)
  
  ###Different research may adopt different data filtering criteria for somatic mutations
  ###As many signatures are defined according to the studies by Alexandrov et al., we use sbs_AlexandrovEtAl as the standard
  adj<-get_sbs_avg(sbs_AlexandrovEtAl)/get_sbs_avg(sbs_tcga)
  dat_sig<-rbind(
    get_sig_pattern(get_sbs_avg(sbs_tcga[names(sbs_tcga)%in%id_list$noGM]),
                    tri.counts.method="exome2genome",
                    signatures.ref=signatures.cosmic,
                    scale=adj),
    get_sig_pattern(get_sbs_avg(sbs_tcga[names(sbs_tcga)%in%id_list$GM]),
                    tri.counts.method="exome2genome",
                    signatures.ref=signatures.cosmic,
                    scale=adj),
    get_sig_pattern(get_sbs_avg(sbs_tcga[names(sbs_tcga)%in%id_list$MMR]),
                    tri.counts.method="exome2genome",
                    signatures.ref=signatures.cosmic,
                    scale=adj),
    get_sig_pattern(get_sbs_avg(sbs_tcga[names(sbs_tcga)%in%id_list$HRR]),
                    tri.counts.method="exome2genome",
                    signatures.ref=signatures.cosmic,
                    scale=adj),
    get_sig_pattern(get_sbs_avg(sbs_tcga[names(sbs_tcga)%in%id_list$RecQ]),
                    tri.counts.method="exome2genome",
                    signatures.ref=signatures.cosmic,
                    scale=adj)
  )
  rownames(dat_sig)<-c("Without GM","With GM","MMR","HRR","RecQ")
  dat_sig
  
  ###Permutation test
  ###Computationally intensive; we performed the tests on a Linux platform with 30 logical cores; it takes about three hours
  ###Reduce n_each can accelerate the tests
  ###The results can be directly downloaded from the GitHub repository
  if(T){
    if(F){
      sbs_list1<-sbs_tcga[names(sbs_tcga)%in%id_list$noGM]
      RNGkind("L'Ecuyer-CMRG")
      set.seed(20250207)
      perm_out<-list(
        GM=perm_test_sig(sbs_list1,sbs_tcga[names(sbs_tcga)%in%id_list$GM],
                         n_cores=30,n_each=1000,
                         signatures.ref=signatures.cosmic,scale=adj,
                         tri.counts.method="exome2genome"),
        MMR=perm_test_sig(sbs_list1,sbs_tcga[names(sbs_tcga)%in%id_list$MMR],
                          n_cores=30,n_each=1000,
                          signatures.ref=signatures.cosmic,scale=adj,
                          tri.counts.method="exome2genome"),
        HRR=perm_test_sig(sbs_list1,sbs_tcga[names(sbs_tcga)%in%id_list$HRR],
                          n_cores=30,n_each=1000,
                          signatures.ref=signatures.cosmic,scale=adj,
                          tri.counts.method="exome2genome"),
        RecQ=perm_test_sig(sbs_list1,sbs_tcga[names(sbs_tcga)%in%id_list$RecQ],
                           n_cores=30,n_each=1000,
                           signatures.ref=signatures.cosmic,scale=adj,
                           tri.counts.method="exome2genome")
      )
    }else{
      load(paste(file_folder1,"perm_out.Rdata",sep="/"))
    }
  }
  
  get_p_one_sided<-function(perm_out,loc,loc0=1:1e4){
    out<-matrix(NaN,nrow=2,ncol=length(perm_out))
    colnames(out)<-names(perm_out)
    rownames(out)<-c("greater","less")
    for(i in 1:length(perm_out)){
      perm_v<-colSums(perm_out[[i]]$diff_perm[loc,loc0,drop=F])
      obs_v<-sum(perm_out[[i]]$diff_obs[loc])
      out[,i]<-c(mean(perm_v>=obs_v),mean(perm_v<=obs_v))
    }
    out
  }
  ###Signature 3
  get_p_one_sided(perm_out,3,1:3e4)
  ###Signature 6
  get_p_one_sided(perm_out,6,1:3e4)
  ###Signature 12&26
  get_p_one_sided(perm_out,c(12,26),1:3e4)
}
##investigate the correlation between gcnc_avg and erbb2_cn
if(T){
  cor.test(tcga_cl$gcnc_avg,tcga_cl$erbb2_cn)
  cor.test(tcga_cl$gcnc_avg,log2(tcga_cl$erbb2_cn))
  ###plot(tcga_cl$gcnc_avg,tcga_cl$erbb2_cn,ylim=c(0,40))
}
##Investigate whether average gene copy number change affects prognosis
if(T){
  ###Extract necessary data from tcga_surv for the survival analysis
  tcga_surv<-suppressWarnings(readLines(file_path7))
  tcga_surv<-read.table(text=tcga_surv[-(1:4)],sep="\t",header=T)
  tcga_cl$OS_MONTHS<-tcga_surv$OS_MONTHS[match(tcga_cl$case_submitter_id,tcga_surv$PATIENT_ID)]
  tcga_cl$OS_STATUS<-tcga_surv$OS_STATUS[match(tcga_cl$case_submitter_id,tcga_surv$PATIENT_ID)]
  
  table(tcga_cl$OS_STATUS);anyNA(tcga_cl$OS_STATUS)
  tcga_cl$OS_STATUS<-as.numeric(substr(tcga_cl$OS_STATUS,1,1))
  
  tcga_cl$age_diagnosis<-tcga_cl$age_at_diagnosis/365.25
  tcga_cl$cancer_type<-tcga_cl$project_id
  ###We will focus on gcnc_avg_bin
  tcga_cl$gcnc_avg_bin<-tcga_cl$gcnc_avg>median(tcga_cl$gcnc_avg)
  
  ###We use msi_score as a covariate and extract it from tcga_msi
  tcga_msi<-suppressWarnings(readLines(file_path8))
  tcga_msi<-read.table(text=tcga_msi[-(1:4)],sep="\t",header=T)
  summary(tcga_msi$MSI_SENSOR_SCORE)
  tcga_cl$msi_score<-tcga_msi$MSI_SENSOR_SCORE[match(tcga_cl$case_submitter_id,tcga_msi$PATIENT_ID)]
  
  library(survival)
  mm<-model.matrix(~OS_MONTHS+OS_STATUS+gcnc_avg_bin+gcnc_avg+msi_score+cancer_type+gender+age_diagnosis,
                   data=tcga_cl)
  head(mm);dim(mm)#478, 9
  mm<-as.data.frame(mm[,-1])
  
  fit0<-coxph(
    Surv(OS_MONTHS,OS_STATUS)~gcnc_avg+msi_score+cancer_typeREAD+gendermale+age_diagnosis,
    data=mm
  );summary(fit0);cox.zph(fit0)
  
  fit<-coxph(
    Surv(OS_MONTHS,OS_STATUS)~gcnc_avg_binTRUE+msi_score+cancer_typeREAD+gendermale+age_diagnosis,
    data=mm
  )
  summary(fit)
  
  fit_curve<-survfit(
    fit,conf.int=0.95,conf.type="log-log",
    newdata=data.frame(
      gcnc_avg_binTRUE=c(0,1),msi_score=rep(mean(mm$msi_score),2),
      cancer_typeREAD=rep(mean(mm$cancer_typeREAD),2),gendermale=rep(mean(mm$gendermale),2),
      age_diagnosis=rep(mean(mm$age_diagnosis),2)
    )
  )
  summary(fit_curve)
  #plot(fit_curve)
}

#Compare the TMB between different GM groups
if(T){
  tcga_cl$tmb<-tcga_msi$TMB_NONSYNONYMOUS[match(tcga_cl$case_submitter_id,tcga_msi$PATIENT_ID)]
  summary(tcga_cl$tmb)
  dat0<-data.frame(value=tcga_cl$tmb,sampleID=tcga_cl$case_submitter_id)
  
  ##With GM
  RNGkind("L'Ecuyer-CMRG")
  set.seed(20250210)
  perm_test_cont(withgm_id=id_list$GM,
                 nogm_id=id_list$noGM,
                 data=dat0,
                 num_sim=5e5,stat="loc",alternative="greater",na.rm=T)
  ##0.015766
  
  ##MMR
  set.seed(20250210)
  perm_test_cont(withgm_id=id_list$MMR,
                 nogm_id=id_list$noGM,
                 data=dat0,
                 num_sim=5e5,stat="loc",alternative="greater",na.rm=T)
  ##0.00177
  
  ##HRR
  set.seed(20250210)
  perm_test_cont(withgm_id=id_list$HRR,
                 nogm_id=id_list$noGM,
                 data=dat0,
                 num_sim=5e5,stat="loc",alternative="two.sided",na.rm=T)
  ##0.925526
  
  ##RecQ
  set.seed(20250210)
  perm_test_cont(withgm_id=id_list$RecQ,
                 nogm_id=id_list$noGM,
                 data=dat0,
                 num_sim=5e5,stat="loc",alternative="two.sided",na.rm=T)
  ##0.998654
}

#Examine the HHR-ERBB2 association in the ESCA&STAD cohort
if(T){
  ##prepare the cell_gm data set
  if(T){
    ###Considering the sample size, both Pathogenic_variants and truncating Prioritized_VUSs will be treated as P/LP
    cell_plp_gm<-readxl::read_excel(file_path1,sheet=1)
    cell_pvus_gm<-readxl::read_excel(file_path1,sheet=2)
    cell_plp_gm<-as.data.frame(cell_plp_gm);cell_pvus_gm<-as.data.frame(cell_pvus_gm)
    dim(cell_plp_gm);dim(cell_pvus_gm)
    identical(colnames(cell_plp_gm),colnames(cell_pvus_gm))#true
    
    table(cell_pvus_gm$Variant_Classification)
    loc<-cell_pvus_gm$Variant_Classification%in%c("frameshift_variant","stop_gained")
    cell_gm<-rbind(cell_plp_gm,cell_pvus_gm[loc,])
    table(cell_gm$cancer)
    cell_gm<-cell_gm[cell_gm$cancer%in%c("ESCA","STAD"),]
    dim(cell_gm)#83, 173
    
    ###Gender
    cell_gm$gender<-tolower(cell_gm$gender)
    
    ###Stage
    table(cell_gm$ajcc_pathologic_tumor_stage);anyNA(cell_gm$ajcc_pathologic_tumor_stage)
    loc<-which(!str_detect(cell_gm$ajcc_pathologic_tumor_stage,"Stage"))
    cell_gm$ajcc_pathologic_tumor_stage[loc]<-NA_character_
    table(cell_gm$ajcc_pathologic_tumor_stage);sum(is.na(cell_gm$ajcc_pathologic_tumor_stage))
    
    ###Other variables required for matching
    suppressWarnings({
      cell_gm$birth_days_to<-as.numeric(cell_gm$birth_days_to)
      cell_gm$age_at_initial_pathologic_diagnosis<-as.numeric(cell_gm$age_at_initial_pathologic_diagnosis)
      cell_gm$initial_pathologic_dx_year<-as.numeric(cell_gm$initial_pathologic_dx_year)
    })
    loc<-is.na(cell_gm$birth_days_to)|is.na(cell_gm$age_at_initial_pathologic_diagnosis)|is.na(cell_gm$initial_pathologic_dx_year)
    cell_gm<-cell_gm[!loc,]
    dim(cell_gm)#79, 173
  }
  
  ##Prepare tcga_cl
  file_path9<-paste(file_folder1,"ESCA","clinical/clinical.tsv",sep="/");stopifnot(file.exists(file_path9))
  file_path10<-paste(file_folder1,"STAD","clinical/clinical.tsv",sep="/");stopifnot(file.exists(file_path10))
  if(T){
    tcga_cl<-rbind(fread(file_path9),fread(file_path10))
    tcga_cl<-as.data.frame(tcga_cl)
    dim(tcga_cl)
    
    tcga_cl<-tcga_cl[!duplicated(tcga_cl$case_id),]
    dim(tcga_cl)#619, 197
    
    loc<-tcga_cl$days_to_birth=="'--"|tcga_cl$age_at_diagnosis=="'--"|tcga_cl$year_of_diagnosis=="'--"
    tcga_cl<-tcga_cl[!loc,]
    dim(tcga_cl)#499, 197
    
    table(tcga_cl$gender)
    table(tcga_cl$project_id)
    
    ##Reported Asian patients are excluded
    table(tcga_cl$race);anyNA(tcga_cl$race)
    tcga_cl<-tcga_cl[tcga_cl$race!="asian",]
    dim(tcga_cl)#386, 197
    
    table(tcga_cl$ajcc_pathologic_stage);anyNA(tcga_cl$ajcc_pathologic_stage)
    loc<-tcga_cl$ajcc_pathologic_stage=="'--"|tcga_cl$ajcc_pathologic_stage=="Unknown"
    tcga_cl$ajcc_pathologic_stage[loc]<-NA_character_
    
    tcga_cl$days_to_birth<-as.numeric(tcga_cl$days_to_birth)
    tcga_cl$age_at_diagnosis<-as.numeric(tcga_cl$age_at_diagnosis)
    tcga_cl$year_of_diagnosis<-as.numeric(tcga_cl$year_of_diagnosis)
    tcga_cl$project_id<-str_remove(tcga_cl$project_id,"TCGA-")
  }
  
  ##Match the samples from tcga_cl with cell_gm
  cell_gm$sampleID<-NA_character_
  if(T){
    for(i in 1:nrow(cell_gm)){
      res_i<-logical(nrow(tcga_cl))
      for(j in 1:nrow(tcga_cl)){
        res_i[j]<-identical_use_cl(
          tcga_cl$gender[j],cell_gm$gender[i],
          tcga_cl$project_id[j],cell_gm$cancer[i],
          tcga_cl$ajcc_pathologic_stage[j],cell_gm$ajcc_pathologic_tumor_stage[i],T,F,
          tcga_cl$days_to_birth[j],cell_gm$birth_days_to[i],days_to_birth_tol=0,
          tcga_cl$age_at_diagnosis[j]/365.25,cell_gm$age_at_initial_pathologic_diagnosis[i],age_at_diagnosis_tol=1,
          tcga_cl$year_of_diagnosis[j],cell_gm$initial_pathologic_dx_year[i],year_of_diagnosis_tol=0
        )
      }
      if(any(res_i)){
        loc<-which(res_i)
        if(length(loc)>1){message("multi-hit: ",i)}
        cell_gm$sampleID[i]<-tcga_cl$case_submitter_id[loc[1]]
      }
    }
    cell_gm$sampleID
    cell_gm<-cell_gm[!is.na(cell_gm$sampleID),]
    nrow(cell_gm)#51
  }
  
  ##prepare the ERBB2 CNV data
  file_path11<-paste(file_folder1,"ESCA","ERBB2/data_cna.txt",sep="/");stopifnot(file.exists(file_path11))
  file_path12<-paste(file_folder1,"STAD","ERBB2/data_cna.txt",sep="/");stopifnot(file.exists(file_path12))
  if(T){
    esca_cnv<-fread(file_path11)
    stad_cnv<-fread(file_path12)
    colnames(esca_cnv)[1:5]
    colnames(stad_cnv)[1:5]
    tcga_cnv<-rbind(
      t(esca_cnv[esca_cnv$Hugo_Symbol=="ERBB2",-c(1:2)]),
      t(stad_cnv[stad_cnv$Hugo_Symbol=="ERBB2",-c(1:3)])
    )
    tcga_cnv<-as.data.frame(tcga_cnv)
    rownames(tcga_cnv)<-substr(rownames(tcga_cnv),1,12)
    dim(tcga_cnv)#620, 1
    anyDuplicated(rownames(tcga_cnv))#0
    
    tcga_cl$erbb2_cnv<-tcga_cnv$V1[match(tcga_cl$case_submitter_id,rownames(tcga_cnv))]
    anyNA(tcga_cl$erbb2_cnv)#false
    tcga_cl$erbb2_amp<-tcga_cl$erbb2_cnv==2
    mean(tcga_cl$erbb2_amp)#0.1450777
  }
  
  ##compare the prevalence of ERBB2 amp between different GM groups
  if(!require("BNDtest",quietly=T)){
    if(!requireNamespace("remotes",quietly=T)){
      install.packages("remotes")
    }
    ###We need this package to perform the Barnard's exact test
    remotes::install_github("limintao-pku/BNDtest")
  }
  if(T){
    gene_list<-list(
      HRR=c('BLM','BRCA1','BRCA2','DMC1','EME1','EME2','GEN1','HFM1','MRE11A','MUS81','NBN','PPP4C',
            'PPP4R1','PPP4R2','PPP4R4','RAD50','RAD51','RAD51B','RAD51C','RAD51D','RAD52','RAD54B',
            'RAD54L','RAD54L2','RDM1','RECQL','RECQL4','RECQL5','RMI1','RMI2','RPA1','RPA2','RPA3',
            'SEM1','SLX1A','SLX4','PPP4R3A','PPP4R3B','SPO11','TOP3A','TOP3B','WRN','XRCC2','XRCC3'),
      FA=c('CENPS','BLM','BRCA1','BRCA2','BRIP1','FAAP100','FAAP24','FAN1','FANCA','FANCC','FANCD2',
           'FANCE','FANCF','FANCG','FANCI','FANCL','FANCM','HES1','PALB2','RAD51','RAD51C','RMI1',
           'RMI2','CENPX','TELO2','TOP3A','TOP3B','UBE2T','USP1','WDR48'),
      RecQ=c('RECQL','WRN','BLM','RECQL4','RECQL5'),
      BRCA=c("BRCA1","BRCA2")
    )
    id_list<-list(
      HRR_FA=unique(cell_gm$sampleID[cell_gm$HUGO_Symbol%in%c(gene_list$HRR,gene_list$FA)]),
      HRR=unique(cell_gm$sampleID[cell_gm$HUGO_Symbol%in%gene_list$HRR]),
      FA=unique(cell_gm$sampleID[cell_gm$HUGO_Symbol%in%gene_list$FA]),
      RecQ=unique(cell_gm$sampleID[cell_gm$HUGO_Symbol%in%gene_list$RecQ]),
      BRCA=unique(cell_gm$sampleID[cell_gm$HUGO_Symbol%in%gene_list$BRCA]),
      noGM=tcga_cl$case_submitter_id[!tcga_cl$case_submitter_id%in%cell_gm$sampleID]
    )
    
    tcga_cl$erbb2_amp<-factor(tcga_cl$erbb2_amp,levels=c(T,F))
    
    ###ESCA and STAD are combined together
    if(T){
      #HRR/FA
      mat<-rbind(
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$HRR_FA]),
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$noGM])
      );mat
      p<-mat[,1]/rowSums(mat)
      list(rr=p[1]/p[2],p.value=Barnard_test(mat,alternative="greater"))
      
      #HRR
      mat<-rbind(
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$HRR]),
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$noGM])
      );mat
      p<-mat[,1]/rowSums(mat)
      list(rr=p[1]/p[2],p.value=Barnard_test(mat,alternative="greater"))
      
      #FA
      mat<-rbind(
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$FA]),
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$noGM])
      );mat
      p<-mat[,1]/rowSums(mat)
      list(rr=p[1]/p[2],p.value=Barnard_test(mat,alternative="greater"))
      
      #BRCA1/2
      mat<-rbind(
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$BRCA]),
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$noGM])
      );mat
      p<-mat[,1]/rowSums(mat)
      list(rr=p[1]/p[2],p.value=Barnard_test(mat,alternative="greater"))
      
      #RecQ
      mat<-rbind(
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$RecQ]),
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$noGM])
      );mat
      p<-mat[,1]/rowSums(mat)
      list(rr=p[1]/p[2],p.value=Barnard_test(mat,alternative="greater"))
    }
    
    ###STAD only
    if(T){
      loc<-tcga_cl$project_id=="STAD"
      
      #HRR/FA
      mat<-rbind(
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$HRR_FA&loc]),
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$noGM&loc])
      );mat
      p<-mat[,1]/rowSums(mat)
      list(rr=p[1]/p[2],p.value=Barnard_test(mat,alternative="greater"))
      
      #HRR
      mat<-rbind(
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$HRR&loc]),
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$noGM&loc])
      );mat
      p<-mat[,1]/rowSums(mat)
      list(rr=p[1]/p[2],p.value=Barnard_test(mat,alternative="greater"))
      
      #FA
      mat<-rbind(
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$FA&loc]),
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$noGM&loc])
      );mat
      p<-mat[,1]/rowSums(mat)
      list(rr=p[1]/p[2],p.value=Barnard_test(mat,alternative="greater"))
      
      #BRCA1/2
      mat<-rbind(
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$BRCA&loc]),
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$noGM&loc])
      );mat
      p<-mat[,1]/rowSums(mat)
      list(rr=p[1]/p[2],p.value=Barnard_test(mat,alternative="greater"))
      
      #RecQ
      mat<-rbind(
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$RecQ&loc]),
        table(tcga_cl$erbb2_amp[tcga_cl$case_submitter_id%in%id_list$noGM&loc])
      );mat
      p<-mat[,1]/rowSums(mat)
      list(rr=p[1]/p[2],p.value=Barnard_test(mat,alternative="greater"))
    }
    
    ###Due to limited sample size, ESCA was not investigated separately
    table(tcga_cl$project_id)#93, 293
  }
}
