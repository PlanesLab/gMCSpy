clear all
close all
clc
addpath(genpath('/scratch/a905383/Software/cplex1210/cplex/matlab/x86-64_linux'))
addpath(genpath('/scratch/a905383/gmcspy/matlab_codes/cobratoolbox'))
addpath(genpath('/scratch/a905383/gmcspy/matlab_codes/GMCS.mat/src'))
addpath(genpath('/scratch/a905383/gmcspy/matlab_codes/GMCS.mat/src_old'))
initCobraToolbox(0)



model_array = {'e_coli_core', 'yeast8','iML1515', 'iJN1463', 'Human-GEM-1.16.0_CultureMedia', 'Human-GEM-1.16.0_raw', 'Recon3D'};
separate_isoform_array = {'', '', '', '', '', '', ''};

%model_array = {'yeast8'};
%separate_isoform_array = {''};

%model_array = model_array([1]);
%separate_isoform_array = separate_isoform_array([1]);
forceLength_array = {'False', 'True'};
forceLength_array = forceLength_array([2]);
nIter_array = 0:9;
solver_array = {'bioinformatics'};

TT = combvec(1:length(model_array), 1:length(nIter_array), 1:length(solver_array), 1:length(forceLength_array), nan, nan);
TT = cell2table(num2cell(TT)');
TT.Properties.VariableNames = {'model', 'iter', 'solver', 'forceLength', 'time', 'size'};
TT.separate_isoform = reshape(separate_isoform_array(TT.model),[],1);
TT.model = reshape(model_array(TT.model),[],1);
TT.iter = reshape(nIter_array(TT.iter),[],1);
TT.solver = reshape(solver_array(TT.solver),[],1);
TT.forceLength = reshape(forceLength_array(TT.forceLength),[],1);
TT.filename_res = strcat('GMCS_mcs_matlab_', TT.model, '_', TT.solver, '_', TT.forceLength, '_', cellfun(@num2str, num2cell(TT.iter)), '.csv');
TT.filename_times = strcat('GMCS_times_matlab_', TT.model, '_', TT.solver, '_', TT.forceLength, '_', cellfun(@num2str, num2cell(TT.iter)), '.csv');
TT
TT = sortrows(TT);
TT

numWorkers = 16;
printLevel = 1;
max_gmcs = 10000;
max_len_gmcs = 3;
timeLimit = 1000;
checkGMCS = 0;


for i = 1:size(TT,1)
    % for i = 79:80
    
    TT(i,:)
    
    model = ['models' filesep TT.model{i} '.mat'];
    model = readCbModel(model);

    G_mat_name = [TT.model{i} '_' TT.solver{i}];
    
    % define different ways to store the matrix
    if strcmp(TT.solver{i}, 'bioinformatics')
        changeCobraSolver('ibm_cplex', 'all');
        tictoc = tic;
        [gmcs, gmcs_time] = calculateGeneMCS_old(G_mat_name, model, max_gmcs, max_len_gmcs,...
            'printLevel', printLevel,...
            'timeLimit', timeLimit, ...
            'forceLength', strcmp(TT.forceLength{i}, 'True'),...
            'numWorkers', numWorkers);
        delete G_* tmp.mat
        TT.time(i) = toc(tictoc);
        TT.size(i) = numel(gmcs);
    elseif strcmp(TT.solver{i}, 'ibm_cplex') || strcmp(TT.solver{i}, 'gurobi')
        addpath(genpath('/scratch/a905383/gmcspy/matlab_codes/GMCS.mat/src'))
        tictoc = tic;
        [gmcs, gmcs_time] = calculateGeneMCS(G_mat_name, model, max_gmcs, max_len_gmcs,...
            'printLevel', printLevel, ...
            'timeLimit', timeLimit, ...
            'forceLength', strcmp(TT.forceLength{i}, 'True'), ...
            'solver', TT.solver{i}, ...
            'numWorkers', numWorkers);
        TT.time(i) = toc(tictoc);
        TT.size(i) = numel(gmcs);
        delete G_* tmp.mat
    end
    
    TT2 = cell(numel(gmcs), max(cellfun(@numel, gmcs)));
    for jj = 1:length(gmcs)
        TT2(jj, 1:numel(gmcs{jj})) =  gmcs{jj};
    end
    TT2 = cell2table(TT2);
    writetable(TT2, ['results' filesep TT.filename_res{i}],...
        'Delimiter', ',', 'WriteRowNames', 1, 'WriteVariableNames', 1);
    
    writetable(cell2table(gmcs_time), ['results' filesep TT.filename_times{i}],...
        'Delimiter', ',', 'WriteRowNames', 1, 'WriteVariableNames', 1);
    
end

aux = datestr(datetime('now'), 'yy-mm-dd-HH-MM-SS');
writetable(TT, ['results' filesep 'results_matlab_' aux '.csv'], 'Delimiter', ',', 'WriteRowNames', 1, 'WriteVariableNames', 1);