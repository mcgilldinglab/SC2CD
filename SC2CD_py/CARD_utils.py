import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()

def R_load(local,py):
    ro.globalenv['path_local']=local
    ro.globalenv['path_py']=py
    r('''setwd(path_local)''')
    r_run=r('''
            knitr::opts_chunk$set(echo = TRUE,warning=FALSE, message=FALSE,dpi=600)
            library(dplyr)
            library(wrMisc)
            library(Rcpp)
            library(fields)
            library(gtools)
            library(scatterpie)
            library(ggplot2)
            library(ggforce)
            library(RcppArmadillo)
            library(SummarizedExperiment)
            library(SingleCellExperiment)
            library(ape)
            library(reticulate)
            library(NMF)
    
            use_python(path_py)
            np<-import("numpy")
            pd<-import('pandas')
            sc<-import('scanpy')
            ad<-import('anndata')
            h5py<-import('h5py')
            random<-import('random')
            torch<-import('torch')
            sk<-import('sklearn')
            ot<-import('ot')
            
            source('R/CARD.imputation.R')
            source('R/CARD.prop.R')
            source('R/CARD.refFree.R')
            source('R/CARD.SCMapping.R')
            source('R/RcppExports.R')
            source('R/utilities.R')
            source('R/visualization.R')
            source('R/measure.R')
            source('MuSiC-master/R/analysis.R')
            source('MuSiC-master/R/music2.R')
            source('MuSiC-master/R/construct.R')
            source('MuSiC-master/R/plot.R')
            source('MuSiC-master/R/simulation.R')
            source('MuSiC-master/R/utils.R')
            source('MuSiC-master/R/zzz.R')
            source('R/CARD.prop.R')
            
            sourceCpp('src/CARDfree.cpp')
            sourceCpp('src/CARDref.cpp')
            ''')

def Gen_Imput(path_data,adata_var_index):
    ro.globalenv['path']=path_data
    ro.globalenv['adata_var_index']=adata_var_index
    r_run=r('''
            setwd(path)
            
            adata<-sc$read_visium(path='',count_file ='filtered_feature_bc_matrix.h5')
            obs<-adata$obs
            sp_loc_ini<-data.frame(x=obs$array_row,y=obs$array_col,row.names = rownames(obs))
            
            sp_count<-np$transpose(adata$X)
            colnames(sp_count)<-rownames(adata$obs)
            rownames(sp_count)<-adata_var_index
            
            scrna<-sc$read("./scRNA.h5ad")
            scrna$var_names_make_unique()
            obs_sc<-scrna$obs
            sc_count<-data.frame(t(scrna$X))
            rownames(sc_count)<-rownames(scrna$var)
            sc_count<-as(as.matrix(sc_count), "sparseMatrix")
            
            sc_meta<-data.frame(cellID=colnames(sc_count),cellType=obs_sc$cell_type,sampleInfo='sample1')
            rownames(sc_meta)<-colnames(sc_count)
             
            setwd(path_local)
            ''')



def Gen_Imput_LIBD(path_data,adata_var_index):
    ro.globalenv['path']=path_data
    ro.globalenv['adata_var_index']=adata_var_index
    r_run=r('''
            setwd(path)
            
            adata<-sc$read_visium(path='',count_file ='filtered_feature_bc_matrix.h5')
            obs<-adata$obs
            sp_loc_ini<-data.frame(x=obs$array_row,y=obs$array_col,row.names = rownames(obs))

            sp_count<-np$transpose(adata$X)
            colnames(sp_count)<-rownames(obs)
            setwd(path_local)
            gene_name<-read.csv('Data/1.DLPFC/151673/gene_name.csv',header = T,row.names = 1)
            sc_aa<-read.csv('Data/1.DLPFC/151673/joint_sc.csv',header=T)
            setwd(path)
            rownames(sp_count)<-adata_var_index
            scrna<-sc$read("./scRNA.h5ad")
            scrna$var_names_make_unique()
            obs_sc<-scrna$obs
            sc_count<-data.frame(t(scrna$X))
            rownames(sc_count)<-rownames(scrna$var)
            sc_count<-sc_count[!is.na(sc_aa$gene_ids),]
            rownames(sc_count)<-sc_aa$gene_ids[!is.na(sc_aa$gene_ids)]
            obs1<-scrna$obs
            ngene<-length(obs1$cell_type)
            colnames(sc_count)<-seq(1,ngene)
            
            sc_meta<-data.frame(cellID=seq(1,ngene),cellType=obs1$cell_type,sampleInfo='sample1')
            sc_meta<-as.matrix(sc_meta)
            rownames(sc_meta)<-seq(1,ngene)
            sc_meta<-data.frame(sc_meta)

            setwd(path_local)
            ''')


def Set_Obj(cluster_pred):
    ro.globalenv['domain']=cluster_pred
    r_run=r('''
        obs$pred<-as.factor(domain)
        sum_info<-obs%>%group_by(pred)%>%summarise(mean_x=mean(array_row),mean_y=mean(array_col),
                                        sd_x=sd(array_row),sd_y=sd(array_col))
        sp_loc_mean<-left_join(obs,sum_info,by='pred')%>%select(mean_x,mean_y)
        rownames(sp_loc_mean)<-rownames(sp_loc_ini)
        sp_loc<-cbind(sp_loc_ini,sp_loc_mean)
        CARD_obj<- createCARDObject_imp(sc_count = sc_count, sc_meta = sc_meta,
                                                      spatial_count = sp_count, spatial_location = sp_loc,
                                	                  ct.varname = "cellType", 
                                                      ct.select = unique(sc_meta$cellType),
                                                      sample.varname = "sampleInfo",
                                                      minCountGene = 100, minCountSpot = 5) 

        ''')



def Get_Dec(lam_c=0,sigma1=0.1,sigma2=0.1):
    ro.globalenv['lam_CARD']=lam_c
    ro.globalenv['sigma1']=sigma1
    ro.globalenv['sigma2']=sigma2
    r_run=r('''CARD_obj_new<-CARD_deconvolution_imp(
        CARD_object = CARD_obj,
        sigma1=sigma1,sigma2=sigma2,lambda=lam_CARD)''')
    cell_dec_py=pd.DataFrame(r('''CARD_obj_new@Proportion_CARD'''))
    cell_dec_py.index=r('''rownames(CARD_obj_new@Proportion_CARD)''')
    cell_dec_py.columns=r('''colnames(CARD_obj_new@Proportion_CARD)''')
    return cell_dec_py
    

def Set_Obj_CARDfree(path_data, cluster_pred):
    ro.globalenv['domain']=cluster_pred
    ro.globalenv['path']=path_data
    r_run=r('''
        setwd(path)
        markerList <- readRDS("markerList_example.rds")
        qc <- clean_markerList_for_CARDfree(
                      markerList = markerList,
                      spatial_count = sp_count,
                      min_markers_per_type = 5
                    )
        markerList <- qc$markerList        
        setwd(path_local)
            
        obs$pred<-as.factor(domain)
        sum_info<-obs%>%group_by(pred)%>%summarise(mean_x=mean(array_row),mean_y=mean(array_col),
                                        sd_x=sd(array_row),sd_y=sd(array_col))
        sp_loc_mean<-left_join(obs,sum_info,by='pred')%>%select(mean_x,mean_y)
        rownames(sp_loc_mean)<-rownames(sp_loc_ini)
        sp_loc<-cbind(sp_loc_ini,sp_loc_mean)
        CARDfree_obj<- createCARDfreeObject_imp(markerList = markerList,
                                                    spatial_count = sp_count,
                                                    spatial_location = sp_loc,
                                                    minCountGene = 100,
                                                    minCountSpot = 5) 
                                                    
        ''')



def Get_Dec_CARDfree(lam_c=0,sigma1=0.1,sigma2=0.1):
    ro.globalenv['lam_CARD']=lam_c
    ro.globalenv['sigma1']=sigma1
    ro.globalenv['sigma2']=sigma2
    r_run=r('''CARDfree_obj_new<-CARD_refFree_imp(
        CARDfree_object = CARDfree_obj,
        sigma1=sigma1,sigma2=sigma2,lambda=lam_CARD)''')
    cell_dec_py=pd.DataFrame(r('''CARDfree_obj_new@Proportion_CARD'''))
    cell_dec_py.index=r('''rownames(CARDfree_obj_new@Proportion_CARD)''')
    cell_dec_py.columns=r('''colnames(CARDfree_obj_new@Proportion_CARD)''')
    return cell_dec_py





