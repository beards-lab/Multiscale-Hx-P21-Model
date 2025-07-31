%{ 
This script runs the optimization. 
%}

clear;

addpath data\
addpath model\
addpath parameters\

figure (1000)
clf 

n_cores = 2; 
delete(gcp('nocreate'))
% p = parpool('local',n_cores); 

%% Initialize 
% Select animal 
animal_id = 11;

% Parameters for optimization 
INDMAP = [5, 6, 11, 12, 13, 14, 15, 16, 17, 18] % top 10

% Animal information 
all_nx_animal_ids = [11 12 51 52 54 55 56];
all_hx_animal_ids = [1, 4, 5, 7, 8, 9, 10, 57, 58, 59, 61, 62];  

% Load data structure 
table = readtable('P21_data_input.xlsx','PreserveVariableNames',true); % read excel data 
data = make_datastructure_P21(animal_id,table); % create the data structure 
data.MgATP_cytoplasm = 8; 
data.MgADP_cytoplasm = 0.05; 
data.Pi_cyto         = 1.3; 

% Get filenames, get paraemters, and set hx flag 
if ismember(animal_id,all_nx_animal_ids) == 1 % check if Nx
    filename = sprintf('opt_pars_Nx%d.mat',animal_id); % to save the optimization results 
    exp_filename = sprintf('Nx%d_data.mat',animal_id); % to load the experimental data 
    [pars,UB,LB,data] = parameters_Nx(data);
elseif ismember(animal_id,all_hx_animal_ids) == 1 % check if Hx
    filename = sprintf('opt_pars_Hx%d.mat',animal_id); 
    exp_filename = sprintf('Hx%d_data.mat',animal_id);
    data.hx_flag = 1; 
    [pars,UB,LB,data] = parameters_Hx(data);
end

% Load experimental data
load (exp_filename,'V_LV_mean','P_LV_mean','V_RV_mean','P_RV_mean')
data.V_LV_avg = V_LV_mean; 
data.V_RV_avg = V_RV_mean; 
data.P_LV_avg = P_LV_mean; 
data.P_RV_avg = P_RV_mean;

%% Set Global 
ALLPARS  = pars;
ODE_TOL  = 1e-6; 
DIFF_INC = sqrt(ODE_TOL);

gpars.INDMAP   = INDMAP;
gpars.ALLPARS  = ALLPARS;
gpars.ODE_TOL  = ODE_TOL;
gpars.DIFF_INC = DIFF_INC;

data.gpars = gpars;

%% Optimization 

optx   = pars(INDMAP); 
opthi  = UB(INDMAP);
optlow = LB(INDMAP);

fun = @(x) model_wrap(x,data); 

% options = optimoptions('lsqnonlin','Display','iter','UseParallel',true);
options = optimoptions('lsqnonlin','Display','iter','UseParallel',true,'FiniteDifferenceStepSize',1e-2,'MaxFunctionEvaluations',2000);

% local optimization 
xopt = lsqnonlin(fun,optx,optlow,opthi,options);

% % multistart 
% optx_ms = ones(length(optx),20); 
% optx_ms(:,1) = optx;
% size_optx_ms = size(optx_ms); 
% i_length = size_optx_ms(2); 
% j = 1; 
% 
% % initialize
% xopt_mat = ones(10,20);
% J_mat = ones(1,20);
% 
% for i = 2:i_length
%     optx_ms(:,i) = unifrnd(optlow,opthi); 
%     [xopt,J] = lsqnonlin(fun,optx_ms(:,i),optlow,opthi,options);
%     xopt_mat(:,j) = xopt;
% %     rout_mat(:,j) = rout;
%     J_vec(j) = J;
%     j = j+1;
% end
% 
% loc = find(J_vec == min(J_vec)); 
% xopt = xopt_mat(:,loc);

%% Only want to change parameters in INDMAP

pars(INDMAP) = xopt; 

% run model with optimized parameters 
[outputs,rout,J] = model_sol(pars,data);

optpars = exp(pars);
disp('optimized parameters')
disp([INDMAP' optpars(INDMAP)])
disp([INDMAP' pars(INDMAP)])

% save (filename)

elapsed_time = toc;
elapsed_time = elapsed_time/60

% iternation number you're at, function value you're at, f(x) = value of
% cost function J

%% Plot

addpath data\
addpath model\
addpath data_struct\
addpath parameters\

params =  {'$C_{SA}$' '$C_{SV}$' '$C_{PA}$'  '$C_{PV}$' ...       
    '$R_{SA}$' '$R_{PA}$' ...
    '$R_m$' '$R_a$' '$R_t$' '$R_p$'   ... 
    '$Amref_{LV}$' '$Amref_{SEP}$' '$Amref_{RV}$' ...
    '$Vw_{LV}$' '$Vw_{SEP}$' '$Vw_{RV}$' ...
    '$k_{TS}$','$k_{TR}$'};     

titlename = erase(filename,'pars');
titlename = erase(titlename,'.mat');
titlename = erase(titlename,'_');

figure(100)
clf 
figure(101)
clf 
[~,index] = sort(J_mat);
for i = 1:size(xopt_mat,2)
%     pars_opt = pars;    
    current_index = index(i);
    xopt = xopt_mat(:,current_index);
    pars(INDMAP) = xopt;
    [outputs,rout,J] = model_sol(pars,data);
    if isempty(outputs) == 1
        V_LV = 0; 
        V_RV = 0; 
        P_LV = 0; 
        P_RV = 0; 
    else 
        V_LV = outputs.volumes.V_LV;
        V_RV = outputs.volumes.V_RV;
        P_LV = outputs.pressures.P_LV;
        P_RV = outputs.pressures.P_RV;
    end
    V_LV_avg = data.V_LV_avg/1000;
    V_RV_avg = data.V_RV_avg/1000;
    P_LV_avg = data.P_LV_avg;
    P_RV_avg = data.P_RV_avg;
    
    hfig100 = figure(100); 
    sgtitle(titlename)
    subplot(5,4,i)
    plot(V_LV,P_LV,'r')
    hold on 
    plot(V_LV_avg,P_LV_avg,'k')
    hold on 
    plot(V_RV,P_RV,'b')
    hold on 
    plot(V_RV_avg,P_RV_avg,'k')
    if isempty(outputs) == 1
        title('N/A')
    else 
        title('Cost: ',num2str(round(J,3)))
    end    
    set(gca,'FontSize',10)

end

hfig101 = figure(101); 
hold on
title(titlename)
boxplot((xopt_mat(:,INDMAP)'))
xticklabels(params(INDMAP))
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',15)

if ~exist('Figures', 'dir')
    mkdir('Figures')
end
if ~exist('Figures/Appendix', 'dir')
    mkdir('Figures/Appendix')
end

% Change animal number 
print(hfig100,'-dpng',strcat('Figures/','/Appendix','/Fa_Nx5.png'))
print(hfig101,'-dpng',strcat('Figures/','/Appendix','/Fb_Nx5.png'))

print(hfig100,'-depsc2',strcat('Figures/','/Appendix','/Fa_Nx5.eps'))
print(hfig101,'-depsc2',strcat('Figures/','/Appendix','/Fb_Nx5.eps'))