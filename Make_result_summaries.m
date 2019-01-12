%% Find ICN pairs based on nchoosek and add the anatomical label
% ICNs, anatomical regions, pvals and beta coefficients for all states 

% load your GICA_dfnc.mat
ComponentList = dfncInfo.comps; % order of anatomical regions 1-100

% Enter the anatomical regions/ artifact for the 100 independent components
% I'm sure there's a way to code this
L = { ...
'Art'; ...                           
'Art'; ...                           
'Putamen R+L'; ...                           
'Art'; ...                           
'Cerebellum R+L'; ...          
'Art'; ...         
'Postcentral Gyrus R+L'; ...        
'Art'; ... 
'Postcentral Gyrus L'; ...                    
'Postcentral Gyrus R'; ...      
'Superior Temporal Gyrus R+L'; ...                          
'Postcentral Gyrus L'; ...                           
'Paracentral Lobule medial'; ...   
'Planum polare R+L'; ...            
'Cerebellum R+L'; ...                       
'Cerebellum R+L'; ...                
'Lingual Gyrus R+L'; ...             
'Putamen R+L'; ...                  
'Precentral Gyrus R+L'; ...                 
'Temporal Pole R+L'; ...                 
'Art'; ...                     
'Precenrtal Gyrus R+L'; ...            
'Fusiform Gyrus R+L'; ...                          
'Superior Parietal Lobule R+L'; ...                           
'Art'; ...    
'Cerebellum R+L'; ...                           
'Temporal Pole R+L'; ...   
'Art'; ...                          
'Art'; ...                           
'Art'; ...  
'Art'; ...                          
'Art'; ...                        
'Art'; ...          
'Art'; ... 
'Fusiform Gyrus R+L'; ...    
'rACC medial'; ...                          
'preSMA medial'; ...                           
'Art'; ...                           
'Posterior Middle Temporal Gyrus R+L'; ...              
'Precuneus medial'; ...    
'Supramarginal Gyrus L'; ...      
'Lateral Inferior Occipital Gyrus R+L'; ...      
'Inferior Parietal lobule R+L'; ...                
'Art'; ...                    
'SMA medial'; ...                
'Art'; ...                  
'Art'; ...                           
'Art'; ... 
'Art'; ...                
'Inferior Frontal Gyrus L'; ...          
'Cuneus R+L'; ...                          
'Art'; ...                          
'Posterior Middle Temporal Gyrus R+L'; ...                        
'Middle Frontal Gyrus R+L'; ...  
'Art'; ...    
'Art'; ...                          
'Superior Parietal lobule R+L'; ...                      
'Inferior Temporal Gyrus R+L'; ...              
'Inferior Parietal lobule L'; ...                          
'dACC medial'; ...                          
'Superior Frontal Gyrus medial'; ...                          
'rACC medial'; ...                          
'Frontal Pole L'; ...            
'Inferior Occipital Gyrus R+L'; ...  
'Inferior Frontal Gyrus R'; ...     
'Superior Temporal Gyrus R+L'; ...
'Art'; ...                           
'Art'; ... 
'Art'; ...   
'Insular Cortex R+L'; ...                        
'Precuneus medial'; ...                         
'Art'; ...  
'Inferior Occipital Gyrus R+L'; ...                          
'Hippocampus R+L'; ...     
'Inferior Parietal lobule R+L'; ...   
'Art'; ...                         
'Cerebellum R+L'; ...                          
'Ventral Striatum R+L'; ...                 
'Art'; ...                           
'Art'; ...                          
'Art'; ...                 
'Art'; ... 
'Art'; ...                        
'Art'; ...                         
'Superior Frontal Gyrus R+L'; ...                       
'Art'; ...                           
'Inferior Occipital Gyrus R+L'; ...     
'Art'; ...                          
'IPL & MFG R+L'; ...    
'Lingual Gyrus R+L'; ...                          
'Art'; ...            
'Art'; ...                         
'Anterior Insula R+L'; ...                           
'Angular Gyrus R+L'; ...                         
'Anterior Insular R+L'; ...               
'Angular Gyrus R'; ...                          
'Posterior Hippocampus R+L'; ...                           
'Middle Frontal Gyrus R+L'; ...                           
'Art'; ...             
'Art'; ...            
};

ActualRegions = cell(size(ComponentList));
for ii = 1:size(ActualRegions,2)
    ActualRegions{ii} = L{ComponentList(ii)};
end
% ActualRegions is now in order of the FNC matrix

Cx = nchoosek(1:size(ComponentList,2),2); % Indexing system that represents the FNC matrix

%% 'clean' states, i.e., same subjects that have been used in mStepwise (no
% NaNs, no missing MatEd
state1_clean = state1_all(goodID_state1_total,:);
state2_clean = state2_all(goodID_state2_total,:);
state3_clean = state3_all(goodID_state3_total,:);
state4_clean = state4_all(goodID_state4_total,:);
state5_clean = state5_all(goodID_state5_total,:);

%% State 1 - model Dx, sex, age, Dx * sex
% Dx effects
sig_state1_Dx = find(p_state1_FD_mat(:,1) < icatb_fdr(p_state1_FD_mat,0.05));
ICNs_state1_Dx = cell(size(sig_state1_Dx,1),8); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta, for categorical vars means for groups

for ii = 1:size(sig_state1_Dx,1)
    ICNs_state1_Dx{ii,1} = Cx(sig_state1_Dx(ii),1);
    ICNs_state1_Dx{ii,2} = Cx(sig_state1_Dx(ii),2);
    ICNs_state1_Dx{ii,3} = ActualRegions{Cx(sig_state1_Dx(ii),1)};
    ICNs_state1_Dx{ii,4} = ActualRegions{Cx(sig_state1_Dx(ii),2)};
    ICNs_state1_Dx{ii,5} = p_state1_FD_mat(sig_state1_Dx(ii),1);
    ICNs_state1_Dx{ii,6} = t_state1_FD_mat(sig_state1_Dx(ii),1);
    ICNs_state1_Dx{ii,7} = nanmean(state1_clean(Dx(goodID_state1_total) == 0,sig_state1_Dx(ii)));
    ICNs_state1_Dx{ii,8} = nanmean(state1_clean(Dx(goodID_state1_total) == 1,sig_state1_Dx(ii)));
end

clear ii 

Summary_state1_Dx = cell2table(ICNs_state1_Dx, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval','mean_TD','mean_PS'});
writetable(Summary_state1_Dx,'Summary_state1_Dx_withFD.xlsx')

% Sex effects
sig_state1_sex = find(p_state1_FD_mat(:,2) < icatb_fdr(p_state1_FD_mat,0.05));
ICNs_state1_sex = cell(size(sig_state1_sex,1),8); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta, for categorical vars means for groups

for ii = 1:size(sig_state1_sex,1)
    ICNs_state1_sex{ii,1} = Cx(sig_state1_sex(ii),1);
    ICNs_state1_sex{ii,2} = Cx(sig_state1_sex(ii),2);
    ICNs_state1_sex{ii,3} = ActualRegions{Cx(sig_state1_sex(ii),1)};
    ICNs_state1_sex{ii,4} = ActualRegions{Cx(sig_state1_sex(ii),2)};
    ICNs_state1_sex{ii,5} = p_state1_FD_mat(sig_state1_sex(ii),2);
    ICNs_state1_sex{ii,6} = t_state1_FD_mat(sig_state1_sex(ii),2);
    ICNs_state1_sex{ii,7} = nanmean(state1_clean(sexOrder(goodID_state1_total) == 1,sig_state1_sex(ii)));
    ICNs_state1_sex{ii,8} = nanmean(state1_clean(sexOrder(goodID_state1_total) == 2,sig_state1_sex(ii)));
end

clear ii 

Summary_state1_sex = cell2table(ICNs_state1_sex, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval','mean_female','mean_male'});
writetable(Summary_state1_sex,'Summary_state1_sex_withFD.xlsx')

% Age effects
sig_state1_age = find(p_state1_FD_mat(:,3) < icatb_fdr(p_state1_FD_mat,0.05));
ICNs_state1_age = cell(size(sig_state1_age,1),6); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state1_age,1)
    ICNs_state1_age{ii,1} = Cx(sig_state1_age(ii),1);
    ICNs_state1_age{ii,2} = Cx(sig_state1_age(ii),2);
    ICNs_state1_age{ii,3} = ActualRegions{Cx(sig_state1_age(ii),1)};
    ICNs_state1_age{ii,4} = ActualRegions{Cx(sig_state1_age(ii),2)};
    ICNs_state1_age{ii,5} = p_state1_FD_mat(sig_state1_age(ii),3);
    ICNs_state1_age{ii,6} = t_state1_FD_mat(sig_state1_age(ii),3);
end

clear ii 

Summary_state1_age = cell2table(ICNs_state1_age, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval'});
writetable(Summary_state1_age,'Summary_state1_age_withFD.xlsx')

% Dx * sex effects
sig_state1_DxSex = find(p_state1_FD_mat(:,5) < icatb_fdr(p_state1_FD_mat,0.05));
ICNs_state1_DxSex = cell(size(sig_state1_DxSex,1),10); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta, mean_TD_fem, mean_TDS_mal, mean_PS_fem, mean_PS_mal,

for ii = 1:size(sig_state1_DxSex,1)
    ICNs_state1_DxSex{ii,1} = Cx(sig_state1_DxSex(ii),1);
    ICNs_state1_DxSex{ii,2} = Cx(sig_state1_DxSex(ii),2);
    ICNs_state1_DxSex{ii,3} = ActualRegions{Cx(sig_state1_DxSex(ii),1)};
    ICNs_state1_DxSex{ii,4} = ActualRegions{Cx(sig_state1_DxSex(ii),2)};
    ICNs_state1_DxSex{ii,5} = p_state1_FD_mat(sig_state1_DxSex(ii),5);
    ICNs_state1_DxSex{ii,6} = t_state1_FD_mat(sig_state1_DxSex(ii),5);
    ICNs_state1_DxSex{ii,7} = nanmean(state1_clean(sexOrder(goodID_state1_total) == 1 & Dx(goodID_state1_total) == 0,sig_state1_DxSex(ii)));
    ICNs_state1_DxSex{ii,8} = nanmean(state1_clean(sexOrder(goodID_state1_total) == 2 & Dx(goodID_state1_total) == 0,sig_state1_DxSex(ii)));
    ICNs_state1_DxSex{ii,9} = nanmean(state1_clean(sexOrder(goodID_state1_total) == 1 & Dx(goodID_state1_total) == 1,sig_state1_DxSex(ii)));
    ICNs_state1_DxSex{ii,10} = nanmean(state1_clean(sexOrder(goodID_state1_total) == 2 & Dx(goodID_state1_total) == 1,sig_state1_DxSex(ii)));
end

clear ii 

Summary_state1_DxSex = cell2table(ICNs_state1_DxSex, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval', 'mean_TD_fem','mean_TD_mal','mean_PS_fem','mean_PS_mal'});
writetable(Summary_state1_DxSex,'Summary_state1_DxSex_withFD.xlsx')


%% State 2 - model sex, age, matEd
% sex effects
sig_state2_sex = find(p_state2_FD_mat(:,1) < icatb_fdr(p_state2_FD_mat,0.05));
ICNs_state2_Sex = cell(size(sig_state2_sex,1),8); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta, for categorical vars means for groups

for ii = 1:size(sig_state2_sex,1)
    ICNs_state2_Sex{ii,1} = Cx(sig_state2_sex(ii),1);
    ICNs_state2_Sex{ii,2} = Cx(sig_state2_sex(ii),2);
    ICNs_state2_Sex{ii,3} = ActualRegions{Cx(sig_state2_sex(ii),1)};
    ICNs_state2_Sex{ii,4} = ActualRegions{Cx(sig_state2_sex(ii),2)};
    ICNs_state2_Sex{ii,5} = p_state2_FD_mat(sig_state2_sex(ii),1);
    ICNs_state2_Sex{ii,6} = t_state2_FD_mat(sig_state2_sex(ii),1);
    ICNs_state2_Sex{ii,7} = nanmean(state2_clean(sexOrder(goodID_state2_total) == 1,sig_state2_sex(ii)));
    ICNs_state2_Sex{ii,8} = nanmean(state2_clean(sexOrder(goodID_state2_total) == 2,sig_state2_sex(ii)));
end

clear ii 

Summary_state2_sex = cell2table(ICNs_state2_Sex, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval','mean_female','mean_male'});
writetable(Summary_state2_sex,'Summary_state2_sex_withFD.xlsx')

% Age effects
sig_state2_age = find(p_state2_FD_mat(:,2) < icatb_fdr(p_state2_FD_mat,0.05));
ICNs_state2_age = cell(size(sig_state2_age,1),6); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state2_age,1)
    ICNs_state2_age{ii,1} = Cx(sig_state2_age(ii),1);
    ICNs_state2_age{ii,2} = Cx(sig_state2_age(ii),2);
    ICNs_state2_age{ii,3} = ActualRegions{Cx(sig_state2_age(ii),1)};
    ICNs_state2_age{ii,4} = ActualRegions{Cx(sig_state2_age(ii),2)};
    ICNs_state2_age{ii,5} = p_state2_FD_mat(sig_state2_age(ii),2);
    ICNs_state2_age{ii,6} = t_state2_FD_mat(sig_state2_age(ii),2);
end

clear ii 

Summary_state2_age = cell2table(ICNs_state2_age, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval'});
writetable(Summary_state2_age,'Summary_state2_age_withFD.xlsx')

% Maternal education effects
sig_state2_MatEd = find(p_state2_FD_mat(:,3) < icatb_fdr(p_state2_FD_mat,0.05));
ICNs_state2_MatEd = cell(size(sig_state2_MatEd,1),6); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state2_MatEd,1)
    ICNs_state2_MatEd{ii,1} = Cx(sig_state2_MatEd(ii),1);
    ICNs_state2_MatEd{ii,2} = Cx(sig_state2_MatEd(ii),2);
    ICNs_state2_MatEd{ii,3} = ActualRegions{Cx(sig_state2_MatEd(ii),1)};
    ICNs_state2_MatEd{ii,4} = ActualRegions{Cx(sig_state2_MatEd(ii),2)};
    ICNs_state2_MatEd{ii,5} = p_state2_FD_mat(sig_state2_MatEd(ii),3);
    ICNs_state2_MatEd{ii,6} = t_state2_FD_mat(sig_state2_MatEd(ii),3);
end

clear ii 

Summary_state2_MatEd = cell2table(ICNs_state2_MatEd, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval'});
writetable(Summary_state2_MatEd,'Summary_state2_MatEd_withFD.xlsx')

%% State 3 - Dx, sex, age, MatEd
% Dx effects
sig_state3_Dx = find(p_state3_FD_mat(:,1) < icatb_fdr(p_state3_FD_mat,0.05));
ICNs_state3_Dx = cell(size(sig_state3_Dx,1),8); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta, for categorical vars means for groups

for ii = 1:size(sig_state3_Dx,1)
    ICNs_state3_Dx{ii,1} = Cx(sig_state3_Dx(ii),1);
    ICNs_state3_Dx{ii,2} = Cx(sig_state3_Dx(ii),2);
    ICNs_state3_Dx{ii,3} = ActualRegions{Cx(sig_state3_Dx(ii),1)};
    ICNs_state3_Dx{ii,4} = ActualRegions{Cx(sig_state3_Dx(ii),2)};
    ICNs_state3_Dx{ii,5} = p_state3_FD_mat(sig_state3_Dx(ii),1);
    ICNs_state3_Dx{ii,6} = t_state3_FD_mat(sig_state3_Dx(ii),1);
    ICNs_state3_Dx{ii,7} = nanmean(state3_clean(Dx(goodID_state3_total) == 0,sig_state3_Dx(ii)));
    ICNs_state3_Dx{ii,8} = nanmean(state3_clean(Dx(goodID_state3_total) == 1,sig_state3_Dx(ii)));
end

clear ii 

Summary_state3_Dx = cell2table(ICNs_state3_Dx, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval','mean_TD','mean_PS'});
writetable(Summary_state3_Dx,'Summary_state3_Dx_withFD.xlsx')

% sex effects
sig_state3_sex = find(p_state3_FD_mat(:,2) < icatb_fdr(p_state3_FD_mat,0.05));
ICNs_state3_sex = cell(size(sig_state3_sex,1),8); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state3_sex,1)
    ICNs_state3_sex{ii,1} = Cx(sig_state3_sex(ii),1);
    ICNs_state3_sex{ii,2} = Cx(sig_state3_sex(ii),2);
    ICNs_state3_sex{ii,3} = ActualRegions{Cx(sig_state3_sex(ii),1)};
    ICNs_state3_sex{ii,4} = ActualRegions{Cx(sig_state3_sex(ii),2)};
    ICNs_state3_sex{ii,5} = p_state3_FD_mat(sig_state3_sex(ii),2);
    ICNs_state3_sex{ii,6} = t_state3_FD_mat(sig_state3_sex(ii),2);
    ICNs_state3_sex{ii,7} = nanmean(state3_clean(sexOrder(goodID_state3_total) == 1,sig_state3_sex(ii)));
    ICNs_state3_sex{ii,8} = nanmean(state3_clean(sexOrder(goodID_state3_total) == 2,sig_state3_sex(ii)));
end

clear ii 

Summary_state3_sex = cell2table(ICNs_state3_sex, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval','mean_female','mean_male'});
writetable(Summary_state3_sex,'Summary_state3_sex_withFD.xlsx')

% age effects
sig_state3_age = find(p_state3_FD_mat(:,3) < icatb_fdr(p_state3_FD_mat,0.05));
ICNs_state3_age = cell(size(sig_state3_age,1),6); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state3_age,1)
    ICNs_state3_age{ii,1} = Cx(sig_state3_age(ii),1);
    ICNs_state3_age{ii,2} = Cx(sig_state3_age(ii),2);
    ICNs_state3_age{ii,3} = ActualRegions{Cx(sig_state3_age(ii),1)};
    ICNs_state3_age{ii,4} = ActualRegions{Cx(sig_state3_age(ii),2)};
    ICNs_state3_age{ii,5} = p_state3_FD_mat(sig_state3_age(ii),3);
    ICNs_state3_age{ii,6} = t_state3_FD_mat(sig_state3_age(ii),3);
end

clear ii 

Summary_state3_age = cell2table(ICNs_state3_age, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval'});
writetable(Summary_state3_age,'Summary_state3_age_withFD.xlsx')

% MatEd effects
sig_state3_MatEd = find(p_state3_FD_mat(:,4) < icatb_fdr(p_state3_FD_mat,0.05));
ICNs_state3_MatEd = cell(size(sig_state3_MatEd,1),6); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state3_MatEd,1)
    ICNs_state3_MatEd{ii,1} = Cx(sig_state3_MatEd(ii),1);
    ICNs_state3_MatEd{ii,2} = Cx(sig_state3_MatEd(ii),2);
    ICNs_state3_MatEd{ii,3} = ActualRegions{Cx(sig_state3_MatEd(ii),1)};
    ICNs_state3_MatEd{ii,4} = ActualRegions{Cx(sig_state3_MatEd(ii),2)};
    ICNs_state3_MatEd{ii,5} = p_state3_FD_mat(sig_state3_MatEd(ii),4);
    ICNs_state3_MatEd{ii,6} = t_state3_FD_mat(sig_state3_MatEd(ii),4);
end

clear ii 

Summary_state3_MatEd = cell2table(ICNs_state3_MatEd, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval'});
writetable(Summary_state3_MatEd,'Summary_state3_MatEd_withFD.xlsx')

%% State 4 - Dx, sex, age, MatEd, Dx*Age, Sex * MatEd
% Dx effects
sig_state4_Dx = find(p_state4_FD_mat(:,1) < icatb_fdr(p_state4_FD_mat,0.05));
ICNs_state4_Dx = cell(size(sig_state4_Dx,1),8); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta, for categorical vars means for groups

for ii = 1:size(sig_state4_Dx,1)
    ICNs_state4_Dx{ii,1} = Cx(sig_state4_Dx(ii),1);
    ICNs_state4_Dx{ii,2} = Cx(sig_state4_Dx(ii),2);
    ICNs_state4_Dx{ii,3} = ActualRegions{Cx(sig_state4_Dx(ii),1)};
    ICNs_state4_Dx{ii,4} = ActualRegions{Cx(sig_state4_Dx(ii),2)};
    ICNs_state4_Dx{ii,5} = p_state4_FD_mat(sig_state4_Dx(ii),1);
    ICNs_state4_Dx{ii,6} = t_state4_FD_mat(sig_state4_Dx(ii),1);
    ICNs_state4_Dx{ii,7} = nanmean(state4_clean(Dx(goodID_state4_total) == 0,sig_state4_Dx(ii)));
    ICNs_state4_Dx{ii,8} = nanmean(state4_clean(Dx(goodID_state4_total) == 1,sig_state4_Dx(ii)));
end

clear ii 

Summary_state4_Dx = cell2table(ICNs_state4_Dx, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval','mean_TD','mean_PS'});
writetable(Summary_state4_Dx,'Summary_state4_Dx_withFD.xlsx')

% sex effects
sig_state4_sex = find(p_state4_FD_mat(:,2) < icatb_fdr(p_state4_FD_mat,0.05));
ICNs_state4_sex = cell(size(sig_state4_sex,1),8); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state4_sex,1)
    ICNs_state4_sex{ii,1} = Cx(sig_state4_sex(ii),1);
    ICNs_state4_sex{ii,2} = Cx(sig_state4_sex(ii),2);
    ICNs_state4_sex{ii,3} = ActualRegions{Cx(sig_state4_sex(ii),1)};
    ICNs_state4_sex{ii,4} = ActualRegions{Cx(sig_state4_sex(ii),2)};
    ICNs_state4_sex{ii,5} = p_state4_FD_mat(sig_state4_sex(ii),2);
    ICNs_state4_sex{ii,6} = t_state4_FD_mat(sig_state4_sex(ii),2);
    ICNs_state4_sex{ii,7} = nanmean(state4_clean(sexOrder(goodID_state4_total) == 1,sig_state4_sex(ii)));
    ICNs_state4_sex{ii,8} = nanmean(state4_clean(sexOrder(goodID_state4_total) == 2,sig_state4_sex(ii)));
end

clear ii 

Summary_state4_sex = cell2table(ICNs_state4_sex, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval','mean_female','mean_male'});
writetable(Summary_state4_sex,'Summary_state4_sex_withFD.xlsx')

% age effects
sig_state4_age = find(p_state4_FD_mat(:,3) < icatb_fdr(p_state4_FD_mat,0.05));
ICNs_state4_age = cell(size(sig_state4_age,1),6); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state4_age,1)
    ICNs_state4_age{ii,1} = Cx(sig_state4_age(ii),1);
    ICNs_state4_age{ii,2} = Cx(sig_state4_age(ii),2);
    ICNs_state4_age{ii,3} = ActualRegions{Cx(sig_state4_age(ii),1)};
    ICNs_state4_age{ii,4} = ActualRegions{Cx(sig_state4_age(ii),2)};
    ICNs_state4_age{ii,5} = p_state4_FD_mat(sig_state4_age(ii),3);
    ICNs_state4_age{ii,6} = t_state4_FD_mat(sig_state4_age(ii),3);
end

clear ii 

Summary_state4_age = cell2table(ICNs_state4_age, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval'});
writetable(Summary_state4_age,'Summary_state4_age_withFD.xlsx')

% MatEd effects
sig_state4_MatEd = find(p_state4_FD_mat(:,4) < icatb_fdr(p_state4_FD_mat,0.05));
ICNs_state4_MatEd = cell(size(sig_state4_MatEd,1),6); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state4_MatEd,1)
    ICNs_state4_MatEd{ii,1} = Cx(sig_state4_MatEd(ii),1);
    ICNs_state4_MatEd{ii,2} = Cx(sig_state4_MatEd(ii),2);
    ICNs_state4_MatEd{ii,3} = ActualRegions{Cx(sig_state4_MatEd(ii),1)};
    ICNs_state4_MatEd{ii,4} = ActualRegions{Cx(sig_state4_MatEd(ii),2)};
    ICNs_state4_MatEd{ii,5} = p_state4_FD_mat(sig_state4_MatEd(ii),4);
    ICNs_state4_MatEd{ii,6} = t_state4_FD_mat(sig_state4_MatEd(ii),4);
end

clear ii 

Summary_state4_MatEd = cell2table(ICNs_state4_MatEd, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval'});
writetable(Summary_state4_MatEd,'Summary_state4_MatEd_withFD.xlsx')

% Dx * Age effects
sig_state4_DxAge = find(p_state4_FD_mat(:,5) < icatb_fdr(p_state4_FD_mat,0.05));
ICNs_state4_DxAge = cell(size(sig_state4_DxAge,1),6); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state4_DxAge,1)
    ICNs_state4_DxAge{ii,1} = Cx(sig_state4_DxAge(ii),1);
    ICNs_state4_DxAge{ii,2} = Cx(sig_state4_DxAge(ii),2);
    ICNs_state4_DxAge{ii,3} = ActualRegions{Cx(sig_state4_DxAge(ii),1)};
    ICNs_state4_DxAge{ii,4} = ActualRegions{Cx(sig_state4_DxAge(ii),2)};
    ICNs_state4_DxAge{ii,5} = p_state4_FD_mat(sig_state4_DxAge(ii),5);
    ICNs_state4_DxAge{ii,6} = t_state4_FD_mat(sig_state4_DxAge(ii),5);
end

clear ii 

Summary_state4_DxAge = cell2table(ICNs_state4_DxAge, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval'});
writetable(Summary_state4_DxAge,'Summary_state4_DxAge_withFD.xlsx')

% to further explore interaction pull out correlations for TD and PS
Corr_state4_DxAge = zeros(size(state4_allClean,1), size(sig_state4_DxAge,1));
for ii = 1:size(sig_state4_DxAge,1)
    Corr_state4_DxAge(:,ii) = state4_allClean(:,sig_state4_DxAge(ii));
end

Corr_state4_DxAge_tbl = table(Corr_state4_DxAge(:,1),Corr_state4_DxAge(:,2),Corr_state4_DxAge(:,3),ageCon_cleaned,Dx_cleaned, 'VariableNames',{'Pair1','Pair2','Pair3','age','Dx'});
writetable(Corr_state4_DxAge_tbl,'Correlation_state4_DxAge_sigPairs_NEW.xlsx')

% Sex * MatEd
sig_state4_SexMatEd = find(p_state4_FD_mat(:,6) < icatb_fdr(p_state4_FD_mat,0.05));
ICNs_state4_SexMatEd = cell(size(sig_state4_SexMatEd,1),6); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state4_SexMatEd,1)
    ICNs_state4_SexMatEd{ii,1} = Cx(sig_state4_SexMatEd(ii),1);
    ICNs_state4_SexMatEd{ii,2} = Cx(sig_state4_SexMatEd(ii),2);
    ICNs_state4_SexMatEd{ii,3} = ActualRegions{Cx(sig_state4_SexMatEd(ii),1)};
    ICNs_state4_SexMatEd{ii,4} = ActualRegions{Cx(sig_state4_SexMatEd(ii),2)};
    ICNs_state4_SexMatEd{ii,5} = p_state4_FD_mat(sig_state4_SexMatEd(ii),7);
    ICNs_state4_SexMatEd{ii,6} = t_state4_FD_mat(sig_state4_SexMatEd(ii),7);
end

clear ii 

Summary_state4_SexMatEd = cell2table(ICNs_state4_SexMatEd, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval'});
writetable(Summary_state4_SexMatEd,'Summary_state4_SexMatEd_withFD.xlsx')

%% State 5 - Dx, sex, age, MatEd, Dx * MatEd
% Dx effects
sig_state5_Dx = find(p_state5_FD_mat(:,1) < icatb_fdr(p_state5_FD_mat,0.05));
ICNs_state5_Dx = cell(size(sig_state5_Dx,1),8); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta, for categorical vars means for groups

for ii = 1:size(sig_state5_Dx,1)
    ICNs_state5_Dx{ii,1} = Cx(sig_state5_Dx(ii),1);
    ICNs_state5_Dx{ii,2} = Cx(sig_state5_Dx(ii),2);
    ICNs_state5_Dx{ii,3} = ActualRegions{Cx(sig_state5_Dx(ii),1)};
    ICNs_state5_Dx{ii,4} = ActualRegions{Cx(sig_state5_Dx(ii),2)};
    ICNs_state5_Dx{ii,5} = p_state5_FD_mat(sig_state5_Dx(ii),1);
    ICNs_state5_Dx{ii,6} = t_state5_FD_mat(sig_state5_Dx(ii),1);
    ICNs_state5_Dx{ii,7} = nanmean(state5_clean(Dx(goodID_state5_total) == 0,sig_state5_Dx(ii)));
    ICNs_state5_Dx{ii,8} = nanmean(state5_clean(Dx(goodID_state5_total) == 1,sig_state5_Dx(ii)));
end

clear ii 

Summary_state5_Dx = cell2table(ICNs_state5_Dx, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval','mean_TD','mean_PS'});
writetable(Summary_state5_Dx,'Summary_state5_Dx_withFD.xlsx')

% sex effects
sig_state5_sex = find(p_state5_FD_mat(:,2) < icatb_fdr(p_state5_FD_mat,0.05));
ICNs_state5_sex = cell(size(sig_state5_sex,1),8); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state5_sex,1)
    ICNs_state5_sex{ii,1} = Cx(sig_state5_sex(ii),1);
    ICNs_state5_sex{ii,2} = Cx(sig_state5_sex(ii),2);
    ICNs_state5_sex{ii,3} = ActualRegions{Cx(sig_state5_sex(ii),1)};
    ICNs_state5_sex{ii,4} = ActualRegions{Cx(sig_state5_sex(ii),2)};
    ICNs_state5_sex{ii,5} = p_state5_FD_mat(sig_state5_sex(ii),2);
    ICNs_state5_sex{ii,6} = t_state5_FD_mat(sig_state5_sex(ii),2);
    ICNs_state5_sex{ii,7} = nanmean(state5_clean(sexOrder(goodID_state5_total) == 1,sig_state5_sex(ii)));
    ICNs_state5_sex{ii,8} = nanmean(state5_clean(sexOrder(goodID_state5_total) == 2,sig_state5_sex(ii)));
end

clear ii 

Summary_state5_sex = cell2table(ICNs_state5_sex, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval','mean_female','mean_male'});
writetable(Summary_state5_sex,'Summary_state5_sex_withFD.xlsx')

% age effects
sig_state5_age = find(p_state5_FD_mat(:,3) < icatb_fdr(p_state5_FD_mat,0.05));
ICNs_state5_age = cell(size(sig_state5_age,1),6); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state5_age,1)
    ICNs_state5_age{ii,1} = Cx(sig_state5_age(ii),1);
    ICNs_state5_age{ii,2} = Cx(sig_state5_age(ii),2);
    ICNs_state5_age{ii,3} = ActualRegions{Cx(sig_state5_age(ii),1)};
    ICNs_state5_age{ii,4} = ActualRegions{Cx(sig_state5_age(ii),2)};
    ICNs_state5_age{ii,5} = p_state5_FD_mat(sig_state5_age(ii),3);
    ICNs_state5_age{ii,6} = t_state5_FD_mat(sig_state5_age(ii),3);
end

clear ii 

Summary_state5_age = cell2table(ICNs_state5_age, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval'});
writetable(Summary_state5_age,'Summary_state5_age_withFD.xlsx')

% MatEd effects
sig_state5_MatEd = find(p_state5_FD_mat(:,4) < icatb_fdr(p_state5_FD_mat,0.05));
ICNs_state5_MatEd = cell(size(sig_state5_MatEd,1),6); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state5_MatEd,1)
    ICNs_state5_MatEd{ii,1} = Cx(sig_state5_MatEd(ii),1);
    ICNs_state5_MatEd{ii,2} = Cx(sig_state5_MatEd(ii),2);
    ICNs_state5_MatEd{ii,3} = ActualRegions{Cx(sig_state5_MatEd(ii),1)};
    ICNs_state5_MatEd{ii,4} = ActualRegions{Cx(sig_state5_MatEd(ii),2)};
    ICNs_state5_MatEd{ii,5} = p_state5_FD_mat(sig_state5_MatEd(ii),4);
    ICNs_state5_MatEd{ii,6} = t_state5_FD_mat(sig_state5_MatEd(ii),4);
end

clear ii 

Summary_state5_MatEd = cell2table(ICNs_state5_MatEd, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval'});
writetable(Summary_state5_MatEd,'Summary_state5_MatEd_withFD.xlsx')

% Dx * MatEd effects
sig_state5_DxMatEd = find(p_state5_FD_mat(:,5) < icatb_fdr(p_state5_FD_mat,0.05));
ICNs_state5_DxMatEd = cell(size(sig_state5_DxMatEd,1),6); % ICN1, ICN2, anat. region1, anat, region2, p-val, beta,

for ii = 1:size(sig_state5_DxMatEd,1)
    ICNs_state5_DxMatEd{ii,1} = Cx(sig_state5_DxMatEd(ii),1);
    ICNs_state5_DxMatEd{ii,2} = Cx(sig_state5_DxMatEd(ii),2);
    ICNs_state5_DxMatEd{ii,3} = ActualRegions{Cx(sig_state5_DxMatEd(ii),1)};
    ICNs_state5_DxMatEd{ii,4} = ActualRegions{Cx(sig_state5_DxMatEd(ii),2)};
    ICNs_state5_DxMatEd{ii,5} = p_state5_FD_mat(sig_state5_DxMatEd(ii),5);
    ICNs_state5_DxMatEd{ii,6} = t_state5_FD_mat(sig_state5_DxMatEd(ii),5);
end

clear ii 

Summary_state5_DxMatEd = cell2table(ICNs_state5_DxMatEd, 'VariableNames',{'ICN1','ICN2','anat1','anat2','pval','tval'});
writetable(Summary_state5_DxMatEd,'Summary_state5_DxMatEd_withFD.xlsx')

% to further explore interaction pull out correlations for TD and PS
Corr_state5_DxMatEd = zeros(size(state5_allClean,1), size(sig_state5_DxMatEd,1));
for ii = 1:size(sig_state5_DxMatEd,1)
    Corr_state5_DxMatEd(:,ii) = state5_allClean(:,sig_state5_DxMatEd(ii));
end

Corr_state5_DxMatEd_tbl = table(Corr_state5_DxMatEd(:,1),MatEdu_cleaned,Dx_cleaned, 'VariableNames',{'Pair1','MatEd','Dx'});
writetable(Corr_state5_DxMatEd_tbl,'Correlation_state5_DxMatEd_sigPairs_NEW.xlsx')

