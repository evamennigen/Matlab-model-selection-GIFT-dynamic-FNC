%% This code performas backward model selection using mStepwise.m and then mancovan.m (all from mancovan toolbox in GIFT) on dynamic functional connectivity states
% by E. Mennigen, 02/2018, UCLA
% Note: states should have no NaNs, just use subjects who actually occupied that state
% Models are fitted to each state separately; 
% Syntax for mStepwise.m: First group variable - categorical, then covariates - continuous variables

% Load states you created with GIFT - dynamic FNC toolbox
% Here, I concatenate states in a slightly different order
state1_all = [aTD_state1;aPS_state1];
state2_all = [aTD_state2;aPS_state2];
state3_all = [aTD_state3;aPS_state3];
state4_all = [aTD_state4;aPS_state4];
state5_all = [aTD_state5;aPS_state5];

% Replace 0 with NaNs
state1_all(state1_all == 0) = NaN;

% Here, I added mean frame-wise displacement (mean FD) to the model to further correct for motion
load('/project/data/meanFD_perSub.mat')
meanFD_GICAorder = meanFD_perSub';

% Load grouping variable
meanFD_ordered = [meanFD_GICAorder(TD);meanFD_GICAorder(PS)];

Dx = [zeros(452,1); ones(129,1)]; % states are concatenated in this way first 452 rows -> TD, then 129 rows -> PS
TD=(find(group==0)); 
PS=(find(group==1));

ageTD = age_cnb(TD); % bring variable in order of dFNC states
agePS = age_cnb(PS);
ageCon = [ageTD;agePS];

sexTD = complete_data.sex(TD);
sexPS = complete_data.sex(PS);
sexOrder = [sexTD;sexPS];


MatEdu_TD = complete_data.Maternal_education(TD);
MatEdu_PS = complete_data.Maternal_education(PS);
MatEduOrder = [MatEdu_TD;MatEdu_PS];

%% Indexing so there are no NaNs in variables 
goodID_MatEdu = isnan(MatEduOrder);
goodID_MatEdu = find(goodID_MatEdu == 0);

covariates = [ageCon, MatEduOrder];

GroupVar = [Dx, sexOrder];

goodID_state1 = find(isnan(state1_all(:,1)) == 0);
goodID_state1_total = intersect(goodID_MatEdu,goodID_state1);

%% mSTEPWISE - backward model selection. The final model is printed to the screen and saved in stats_AllVars_state
% you have to adapt these lines after running model selection + univariate
% tests for each state (or copy this chunk of code and adapt it accordingly
% - but you will have to make some changes on the go)

[T_AllVars_state1, p_AllVars_state1, stats_AllVars_state1] = mStepwise(state1_all(goodID_state1_total,:), GroupVar(goodID_state1_total,:), covariates(goodID_state1_total,:), 0.05, {'group-group' 'group-covariate' 'verbose' 'SVD'})
% which terms have been dropped is printed to the screen and the final model
% See mStepwise.m for further info on covariates, group variables and their interaction 

Terms_for_mancovan_state1 = stats_AllVars_state1.Terms; 
DesignMat_for_mancovan_state1 = stats_AllVars_state1.X; 

%% Univariate tests
% find the significant terms in the model 

[sig_terms_state1, I, J]= mUnique(Terms_for_mancovan_state1);
sig_terms_state1_FD = sig_terms_state1;
sig_terms_state1_FD{1,7} = 5; % this is adding an extra term for mean FD (which was bypassed from the model selection); cave! the index for sig_terms_state1_FD depends on your number of significant variables

Terms_for_mancovan_state1_FD = Terms_for_mancovan_state1;
Terms_for_mancovan_state1_FD{1,7} = 5; % this is adding an extra term for mean FD (which was bypassed from the model selection); cave! the index for sig_terms_state1_FD depends on your number of significant variables

DesignMat_for_mancovan_state5_FD = DesignMat_for_mancovan_state1;
DesignMat_for_mancovan_state5_FD(:,7) = meanFD_ordered(goodID_state1_total); % adding mean FD to the model

t_state1_FD = cell(1,length(sig_terms_state1_FD)-1); % just creating empty structures to store t- and p-values and stats from univariate tests
p_state1_FD = cell(1,length(sig_terms_state1_FD)-1);
stats_state1_FD = cell(1,length(sig_terms_state1_FD)-1);


for ii = 2:length(sig_terms_state1_FD)  % First entry is a constant that you don't want to test
     fprintf('Working on term %d of %d\n', ii-1, length(sig_terms_state1_FD)-1)
     [t_tmp, p_tmp, stats_tmp] = mT(state1_all(goodID_state1_total,:), DesignMat_for_mancovan_state1_FD, Terms_for_mancovan_state1_FD, sig_terms_state1_FD{ii}, {'verbose'}); 
     t_state1_FD{ii-1} = t_tmp;
     p_state1_FD{ii-1} = p_tmp;
     stats_state1_FD{ii-1} = stats_tmp;
     clear t_tmp p_tmp stats_tmp Rsq
end

t_state1_FD_mat = reshape(cell2mat(t_state1_FD),[1711 6]); % Just reshaping the t- and p-values to get dimensionality of connectivity matrix
p_state1_FD_mat = reshape(cell2mat(p_state1_FD),[1711 6]);

fdr_state1 = icatb_fdr(p_state1_FD_mat(:),0.05); % Calculating the FDR threshold for this state

% The variable stats.B stores the regression coefficients; first row is constant term
state1_Dx_coeff = stats_state1_FD{1,1}.B(2,sig_state1_Dx)'; 

%% Summary metrics from dynamic FNC analysis like Fraction of Time (FT) and Mean Dwell Time (MDT)
% mSTEPWISE 
covariates = [ageCon, MatEduOrder, meanFD_ordered];
covariates_clean = covariates(goodID_MatEdu,:);

GroupVar = [Dx, sexOrder];
GroupVar_clean = GroupVar(goodID_MatEdu,:);

% Fraction of time
[T_AllVars_GoodSubs_FT, p_AllVars_GoodSubs_FT, stats_AllVars_GoodSubs_FT] = mStepwise(FT_all(goodID_MatEdu), GroupVar_clean, covariates_clean, 0.05, {'group-group' 'group-covariate' 'verbose'})
% which terms have been dropped is printed to the screen and the final model

Terms_for_mancovan_FT = stats_AllVars_GoodSubs_FT.Terms; % output from mstepwise
DesignMat_for_mancovan_FT = stats_AllVars_GoodSubs_FT.X; % output from mstepwise

%% Univariate tests

% FT
[sig_terms_FT, I, J]= mUnique(Terms_for_mancovan_FT);
sig_terms_FT{1,5} = 5; 
Terms_for_mancovan_FT{1,5} = 5;
DesignMat_for_mancovan_FT(:,5) = meanFD_ordered(goodID_MatEdu,:);

t_FT = cell(1,length(sig_terms_FT)-1);
p_FT = cell(1,length(sig_terms_FT)-1);
stats_FT = cell(1,length(sig_terms_FT)-1);

for ii = 2:length(sig_terms_FT)  % no test for constant
     fprintf('Working on term %d of %d\n', ii-1, length(sig_terms_FT)-1)
     [t_tmp, p_tmp, stats_tmp] = mT(FT_all(goodID_MatEdu,:), DesignMat_for_mancovan_FT, Terms_for_mancovan_FT, sig_terms_FT{ii}, {'verbose'}); 
     t_FT{ii-1} = t_tmp;
     p_FT{ii-1} = p_tmp;
     stats_FT{ii-1} = stats_tmp;
     clear t_tmp p_tmp stats_tmp Rsq
end

t_stateFT_mat = reshape(cell2mat(t_FT),[5 4]); % dimensionality is number of states x number of significant variables
p_stateFT_mat = reshape(cell2mat(p_FT),[5 4]);

FDR_FT = icatb_fdr(p_stateFT_mat(:),0.05);

%% MDT
[T_AllVars_GoodSubs_MDT, p_AllVars_GoodSubs_MDT, stats_AllVars_GoodSubs_MDT] = mStepwise(MDT_all(goodID_MatEdu), GroupVar_clean, covariates_clean, 0.05, {'group-group' 'group-covariate' 'verbose'})
% which terms have been dropped is printed to the screen and the final model

Terms_for_mancovan_MDT = stats_AllVars_GoodSubs_MDT.Terms; % output from mstepwise
DesignMat_for_mancovan_MDT = stats_AllVars_GoodSubs_MDT.X; % output from mstepwise

[sig_terms_MDT, I, J]= mUnique(Terms_for_mancovan_MDT);
sig_terms_MDT{1,5} = 5; 
Terms_for_mancovan_MDT{1,5} = 5;
DesignMat_for_mancovan_MDT(:,5) = meanFD_ordered(goodID_MatEdu,:);

t_MDT = cell(1,length(sig_terms_MDT)-1);
p_MDT = cell(1,length(sig_terms_MDT)-1);
stats_MDT = cell(1,length(sig_terms_MDT)-1);

for ii = 2:length(sig_terms_MDT)  % no test for constant
     fprintf('Working on term %d of %d\n', ii-1, length(sig_terms_MDT)-1)
     [t_tmp, p_tmp, stats_tmp] = mT(MDT_all(goodID_MatEdu,:), DesignMat_for_mancovan_MDT, Terms_for_mancovan_MDT, sig_terms_MDT{ii}, {'verbose'}); 
     t_MDT{ii-1} = t_tmp;
     p_MDT{ii-1} = p_tmp;
     stats_MDT{ii-1} = stats_tmp;
     clear t_tmp p_tmp stats_tmp Rsq
end

t_stateMDT_mat = reshape(cell2mat(t_MDT),[5 4]); % dimensionality is number of states x number of significant variables
p_stateMDT_mat = reshape(cell2mat(p_MDT),[5 4]);

FDR_MDT = icatb_fdr(p_stateMDT_mat(:),0.05);


