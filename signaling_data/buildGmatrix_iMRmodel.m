function [G, G_ind, G_rxns, related, n_genes_KO, G_time, iMRmodel_info] = buildGmatrix_iMRmodel(model_name, model_struct, regulatoryNetwork, nLayers, pos, numWorkers)
% Build the G matrix in eGPR networks required for the calculation of
% genetic Minimal Cut Sets (gMCSs) on integrated metabolic and regulatory
% models.
%
% USAGE:
%
%    [G, G_ind, related, n_genes_KO, G_time, eGPRinfo] = buildGmatrix_eGPR(model_name, model_struct, regulatoryNetwork, nLayers)
%
% INPUTS:
%    model_name:            Name of the metabolic-regulatory model under study.
%    model_struct:          Metabolic model structure (COBRA Toolbox format).
%    regulatoryNetwork:     Regulatory information (Source gene, target gene, interaction (+/-)).
%    nLayers:               Number of regulatory layers that are wanted to be added to the integrated model.
%
% OUTPUTS:
%    G:                 G matrix.
%    G_ind:             Gene knockouts associated with each row in the G matrix.
%    related:           Relationships among rows of the G matrix.
%    n_genes_KO:        Number of genes whose knockout is added by each row of the G matrix.
%    G_time:            Calculation times for each of the steps.
%    iMRmodel_info:     Structure with the information of the integration of the metabolic and regulatory model: 
%                       
%                       * .layersPer_eGPR - number of regulatory layers
%                       added in each metabolic GPR rule.
%                       * .reason - The reason why that number of layers
%                       has been added:
%                           - 1: It is not possible to add all the wanted 
%                           regulatory layers because it is needed more 
%                           regulatory information.
%                           - 2: The addition of the last regulatory layer
%                           creates a cycle, so the regulatory layer is not added. 
%                           - 3: All the wanted regulatory layers are added
%                           to the model.
%                       * .genes_eGPR - all the genes (metabolic and
%                       regulatory) of each GPR rule.
% 
%
% EXAMPLE:
%
%    [G, G_ind, related, n_genes_KO, G_time, eGPRinfo] = buildGmatrix_eGPR(Human1_TRRUST1, HumanGEM, TRRUST_interactions, 1)
%
% .. Authors:
%       - Naroa Barrena, 09/12/2022, University of Navarra, TECNUN School of Engineering
%       - Luis V. Valcarcel, 09/12/2022, University of Navarra, TECNUN School of Engineering.
%       - Francisco J. Planes, 09/12/2022, University of Navarra, TECNUN School of Engineering.

tic
% global CBTDIR
% tmpFolderName = [CBTDIR filesep '.tmp'];
% if ~exist(tmpFolderName,'dir')  % Create directories if needed
%     mkdir(tmpFolderName)
% end
% if ~exist([tmpFolderName filesep 'iMR_models'],'dir')
%     mkdir([tmpFolderName filesep 'iMR_models'])
% end
% 
% search_filename = [tmpFolderName filesep 'iMR_models' filesep 'iMR_' model_name '_' num2str(nLayers) '_layers.mat'];
% if exist(search_filename, 'file')
%     load(search_filename);
% else
if exist('pos','var') == 0
    pos = [];
end
search_filename = ['iMR_' model_name '_' num2str(nLayers) '_layers.mat'];
n_mcs = 100000000000;

% Unique GPR rules are calculated. The gMCSs are computed once per each
% unique GPR rule.
uniqueGPR = unique(model_struct.grRules, 'stable');
metGenes = cellfun(@unique, regexp(uniqueGPR, '([^\(\)\+\!\|\&\=\s\or\(and)]+)', 'match'), 'Uniform', 0);

iMRmodel_info = struct('layersPer_eGPR', zeros(1, length(uniqueGPR)), 'reason', zeros(1, length(uniqueGPR)), 'genes_eGPR', []);

for i = 1:length(uniqueGPR)
    disp(i)
    if uniqueGPR(i)~=""
        eGPR_network = {};
        BooleanRules = {};
        BooleanRules(1,1) = strcat('target = ' ,{' '}, uniqueGPR(i));
        if isempty(find(pos == i))
            if  ~isempty(find(ismember(table2cell(regulatoryNetwork(:,2)), metGenes{i,1}),1)) %if any ot the genes in the GPR rule is in the signaling network as targets
                for layer = 1:nLayers
                    pre_eGPR_network = eGPR_network;
                    if layer == 1
                        genes_sel = metGenes{i,1};
                        eGPR_new = table2cell(regulatoryNetwork(ismember(table2cell(regulatoryNetwork(:,2)), genes_sel),:));
                    else
                        genes_sel = eGPR_network(:,1);
                        eGPR_new = table2cell(regulatoryNetwork(ismember(table2cell(regulatoryNetwork(:,2)), genes_sel),:));
                    end
                    eGPR_network = [eGPR_network; eGPR_new];
                    eGPR_network_new = cell2table(eGPR_network);
                    eGPR_network = table2cell(unique(eGPR_network_new, 'rows'));
                    if size(pre_eGPR_network,1) == size(eGPR_network,1)
                        iMRmodel_info.layersPer_eGPR(i) = layer - 1;
                        if size(eGPR_network,1) == 1
                            iMRmodel_info.genes_eGPR{i} = unique([unique([eGPR_network(:,1) eGPR_network(:,2)]'); metGenes{i,1}']);
                        else
                            iMRmodel_info.genes_eGPR{i} =  unique([unique([eGPR_network(:,1) eGPR_network(:,2)]); metGenes{i,1}']);
                        end
                        iMRmodel_info.reason(i) = 1;
                        eGPR_network = pre_eGPR_network;
                        break
                    else
                        iMRmodel_info.layersPer_eGPR(i) = layer;
                    end

                    [cycle, eGPR_network, iMRmodel_info] = checkCycles(eGPR_network, pre_eGPR_network, iMRmodel_info, metGenes{i,1}, i, layer);
                    if cycle == 1
                        break
                    end
                end
            else
                iMRmodel_info.reason(i) = 3;
            end
        end
        if size(eGPR_network,1)==1
            iMRmodel_info.genes_eGPR{i} = unique([unique([eGPR_network(:,1) eGPR_network(:,2)]'); metGenes{i,1}']);
        elseif size(eGPR_network,1)==0 && ~isempty(metGenes{i,1})
            iMRmodel_info.genes_eGPR{i} = metGenes{i,1}';
        else
            iMRmodel_info.genes_eGPR{i} =  unique([unique([eGPR_network(:,1) eGPR_network(:,2)]); metGenes{i,1}']);
        end
        for j = 1:size(eGPR_network,1)
            if cell2mat(eGPR_network(j,3)) == 1
                BooleanRules(j+1,1) = strcat(eGPR_network(j,2), ' = ', {' '}, eGPR_network(j,1));
            elseif cell2mat(eGPR_network(j,3)) == -1
                BooleanRules(j+1,1) = strcat(eGPR_network(j,2), ' = !', eGPR_network(j,1));
            end
        end
        [RxnFormulas, auxKO_Nodes, inputs] = generateFormulas(BooleanRules); 
        target = 'target';
        [gMCSs] = generateModel(RxnFormulas, auxKO_Nodes, target, n_mcs, inputs, 1, numWorkers);
        allgMCSs{i} = gMCSs;
    else
        allgMCSs{i} = {};
    end 

end
time = toc;
G_time{1, 1} = '------ TIMING ------';
G_time{1, 2} = '--- MCSs in eGPR networks ---';
G_time{2, 1} = 'Total';
G_time{2, 2} = time;
save(search_filename, 'allgMCSs', 'iMRmodel_info', 'G_time');
% end
tic;

[~, ~, posGPR] = unique(model_struct.grRules, 'stable');
% Once the gMCSs per each eGPR are computed, the G matrix is built.
n_rxns_or_and = length(allgMCSs);
lengths = cellfun(@length, allgMCSs, 'Uniform', 0);

G_ind = vertcat(allgMCSs{:});
G_ind = cellfun(@transpose, G_ind, 'UniformOutput', false);
G_ind = cellfun(@sort, G_ind, 'UniformOutput', false);
G_mat = cellfun(@cell2mat, G_ind, 'UniformOutput', false);

[duplicated, b, a] = unique(G_mat, 'stable');
G = spalloc(length(duplicated), length(posGPR), 5*length(duplicated));
k = 1;
for i = 1:length(lengths)
    disp(['1_' num2str(i) ' of ' num2str(n_rxns_or_and)]) 
    if lengths{i}>0
        reaction = posGPR == i;
        k = k:lengths{i}+k-1;
        G(a(k),reaction) = 1;
        k = k(end)+1;
    end
end

G_ind = G_ind(b,:);

n_genes_KO = cellfun(@length, G_ind);
[n_genes_KO, ind] = sort(n_genes_KO, 'ascend');
G_ind = G_ind(ind);
G = G(ind, :);
n_G_ind = length(G_ind);

related = zeros(0,2);
% generate matrix that relates genes and G_ind
Gind2genes_genes = unique([G_ind{:}]);
Gind2genes_mat = spalloc(length(G_ind),length(Gind2genes_genes), sum(cellfun(@length,G_ind)));
for i = 1:length(G_ind)
    Gind2genes_mat(i,ismember(Gind2genes_genes, G_ind{i})) = 1;
end
% use matrix to search G_inds that contains lower order G_inds
for i = 1:n_G_ind
    act_G_ind = G_ind{i};
    n_act_G_ind = length(act_G_ind);
    pos = find(n_genes_KO > n_act_G_ind);
    % if a G_ind contains a smaller order one, all columns for these genes
    % should be one
    pos = pos(mean(Gind2genes_mat(pos,ismember(Gind2genes_genes,act_G_ind)),2)==1);
    
    for j = 1:length(pos)
        % Increase the G matrix
        G(pos(j), :) = G(pos(j), :) + G(i, :);
        
        % Add the relationships between G_inds
        related(end+1,:) = [pos(j), i];
    end
end
G = double(G>0);

if size(related)>0
    un_related = unique(related(:, 1));
    n_un_related = length(un_related);
    for i = 1:n_un_related
        act_KO = un_related(i);
        ind = related(:, 1) == act_KO;
        act_related = related(ind, 2);
        all_genes = [G_ind{act_related}];
        un_all_genes = unique(all_genes);
        n_un_all_genes = length(un_all_genes);
        n_genes_KO(act_KO) = n_genes_KO(act_KO)-n_un_all_genes;
    end
else
    related = NaN;
end

time = toc;

G_time{3, 1} = '------ TIMING ------';
G_time{3, 2} = '--- G matrix ---';
G_time{4, 1} = 'Total G matrix';
G_time{4, 2} = time;

final_filename = [pwd filesep 'G_matrices/G_' model_name '_' num2str(nLayers) '_layers.mat'];
G_rxns = model_struct.rxns;
save(final_filename, 'G', 'G_ind', 'G_rxns', 'related', 'n_genes_KO', 'G_time', 'iMRmodel_info');

end

function [cycle, network_2, iMRmodel_info] = checkCycles(network_2, network_pre, iMRmodel_info, metGenes, pos, layer)
reactions = {};
abr = {};
names = {};

for j = 1:size(network_2,1)
    reactions(j,1) = strcat(network_2(j,1), ' -> ', {' '}, network_2(j,2));
    abr{j,1} = ['R' num2str(j)];
    names{j,1} = ['Reaction_' num2str(j)];
end

%Check of the presence of cycles
submodel = createModel(abr, names, reactions, 'printLevel', 0);
S = submodel.S;
[nrow, ncol] = size(S);

A = sparse(nrow+1, ncol);
A(1:nrow,:) = S;
A(nrow+1,:) = ones(1, ncol);

rhs = zeros(nrow+1,1);
lhs = zeros(nrow+1,1);
rhs(1:nrow,1) = 0;
lhs(1:nrow,1)= 0;
rhs(nrow+1,1) = inf;
lhs(nrow+1,1)= 1;

ub = zeros(ncol, 1);
lb = zeros(ncol, 1);

ub(1:ncol,1) = inf;
lb(1:ncol,1) = 0;
ctype = '';
ctype(1:ncol) = 'C';
obj = ones(1, ncol);
sense = 'minimize';
cplex = Cplex('Prueba');
Model = struct();
[Model.A, Model.rhs, Model.lhs, Model.ub, Model.lb, Model.obj, Model.ctype, Model.sense] = deal(A, rhs, lhs, ub, lb, obj, ctype, sense);
cplex.Model = Model;
cplex.Param.threads.Cur = 12;
cplex.DisplayFunc = [];

cplex.solve()
if cplex.Solution.status~=103
    iMRmodel_info.layersPer_eGPR(pos)=layer-1;
    if size(network_pre,1)==1
        iMRmodel_info.genes_eGPR{pos} = unique([unique([network_pre(:,1) network_pre(:,2)]'); metGenes']);
    elseif size(network_pre,1)==0
        iMRmodel_info.genes_eGPR{pos} = {};
    else
        iMRmodel_info.genes_eGPR{pos} =  unique([unique([network_pre(:,1) network_pre(:,2)]); metGenes']);
    end
    iMRmodel_info.reason(pos) = 2;
    network_2 = network_pre;
    cycle = 1;
else
    iMRmodel_info.layersPer_eGPR(pos) = layer;
    cycle = 0;
end

end

function [RxnFormulas, auxKO_Nodes, inputs] = generateFormulas(BooleanRules)
% The Boolean rules of the network are parsed to RxnFormula format.
% Auxiliary nodes and auxKO nodes (KO nodes) are added to the rules.

BooleanRulesParsed = firstPreParse(BooleanRules); %Boolean rules are preparsed.

% Input nodes are found.
% leftNodes = cellfun(@unique, regexp(extractAfter(BooleanRulesParsed,' = '), '([^\(\)\+\!\|\&\=\s]+)', 'match'), 'Uniform', 0);
leftNodes = cellfun(@unique, regexp(extractAfter(BooleanRulesParsed,' = '), '([^\(\)\!\|\&\=\s]+)', 'match'), 'Uniform', 0);
leftNodes = unique([leftNodes{~cellfun(@isempty,leftNodes)}])';

% rigthNodes = cellfun(@unique, regexp(extractBefore(BooleanRulesParsed,' = '), '([^\(\)\+\!\|\&\=\s]+)', 'match'), 'Uniform', 0);
rigthNodes = cellfun(@unique, regexp(extractBefore(BooleanRulesParsed,' = '),  '([^\(\)\!\|\&\=\s]+)', 'match'), 'Uniform', 0);
rigthNodes = unique([rigthNodes{~cellfun(@isempty,rigthNodes)}])';

inputs = setdiff(leftNodes, rigthNodes);
outputs = setdiff(rigthNodes, leftNodes);

% Aux KO nodes, which are the ones for doing the knockout of the node,
% are added to each Boolean rule.
[BooleanRulesParsed, auxKO_Nodes] = AddingAuxKO_Nodes(BooleanRulesParsed, inputs, outputs);

BooleanRulesParsed2 = BooleanRulesParsed(~strcmp(regexprep(BooleanRulesParsed, ' =.*$', ''), outputs));

auxKO_Node_off = sort(strcat(regexprep(BooleanRulesParsed2, ' =.*$', ''), '_auxKO_off'));
auxKO_Nodes = [strcat(auxKO_Nodes, '_on'); auxKO_Node_off];

leftNodes = regexprep(BooleanRulesParsed, ' =.*$', ''); %genes before =

BoolRulesRight = extractAfter(BooleanRulesParsed,'= ');

fp = FormulaParserBoolean();
model = cell(length(BoolRulesRight),1);

zz = num2str(length(BoolRulesRight));
for i = 1:length(BoolRulesRight)
%     if printLevel == 1
        disp(['Parsing Boolean ON rules ' num2str(i) ' / ' zz])
%     end
    head = fp.parseFormulaBoolean_on(BoolRulesRight{i});
    model{i} = modelParser(head,['GPR_ON_', leftNodes{i,1}, '_N_1'], [leftNodes{i,1} '_on']);
    model{i} = ReduceModel(model{i});
end

% The RxnFormulas of the model are saved and reduced in order to eliminate
% the redundant nodes.
RxnFormulas = cell(size(model));
for i = 1:length(model)
    RxnFormulas{i} = reshape(model{i},1, []);
end
RxnFormulas = reshape([RxnFormulas{:}],[],1);
RxnFormulas = unique(RxnFormulas);

%Negated rules are calculated
fp = FormulaParserBoolean();
model = cell(length(BoolRulesRight),1);

zz = num2str(length(BoolRulesRight));
for i = 1:length(BoolRulesRight)
%     if printLevel == 1
        disp(['Parsing Boolean OFF rules ' num2str(i) ' / ' zz])
%     end
    head = fp.parseFormulaBoolean_off(BoolRulesRight{i});
    model{i} = modelParser(head,['GPR_OFF_', leftNodes{i,1}, '_N_1'], [leftNodes{i,1} '_off']);
    model{i} = ReduceModel(model{i});
end

% The RxnFormulas of the model are saved and reduced in order to eliminate
% the redundant nodes.
RxnFormulas_off = cell(size(model));
for i = 1:length(model)
    RxnFormulas_off{i} = reshape(model{i},1, []);
end
RxnFormulas_off = reshape([RxnFormulas_off{:}],[],1);
RxnFormulas_off = unique(RxnFormulas_off);

RxnFormulas = [RxnFormulas; RxnFormulas_off];

end

function [preParsedBoolRules] = firstPreParse(BoolRules)
% This function preparses Boolean rules. Then, it moves all the inhibitory
% relationships into the rigth part of the rule and joins all the OR
% rules for the same node in the same line by | operator. 

auxPreParsedBoolRules = regexprep(BoolRules, '[\]\}]',')'); % replace other brackets by parenthesis
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '[\[\{]','('); % replace other brackets by parenthesis
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '*=','='); % replace *= by =
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '"',''); % eliminate "
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '/','_'); % replace / by _
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '-','_'); % replace - by _
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '([\)]\s?|\s)\s*(?i)(and)\s*?(\s?[\(]|\s)\s*', '$1&$3'); % replace all ands
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '([\)]\s?|\s)\s*(?i)(+)\s*?(\s?[\(]|\s)\s*', '$1&$3'); % replace all ands
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '([\)]\s?|\s)\s*(?i)(or)\s*?(\s?[\(]|\s)\s*', '$1|$3'); % replace all ors
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '([\)]\s?|\s)\s*(?i)(not)\s*?(\s?[\(]|\s)\s*', '$1!$3'); % replace all nots
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '[\s]?&[\s]?', ' & '); % introduce spaces around ands
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '[\s]?\|[\s]?', ' | '); % introduce spaces around ors
auxPreParsedBoolRules = regexprep(auxPreParsedBoolRules, '! ', '!'); % eliminate spaces after !

notGenes = find(contains(extractBefore(auxPreParsedBoolRules,' ='),'!')); % nodes in the left side that have !.
% For those nodes, the ! is moved to the other side and parenthesis are
% added to the expression, so as ! affects to all the rule and not only to a
% unique node.
% Example:              !A = B | C ==> A = !(B | C)

if ~isempty(notGenes)
    leftGenes = extractBefore(auxPreParsedBoolRules,' =');
    leftGenes = regexprep(leftGenes, '!', '');
    BoolRuleRight = extractAfter(auxPreParsedBoolRules,'= ');
    for i = 1:length(notGenes)
        BoolRuleRight{notGenes(i),1} = insertBefore(BoolRuleRight{notGenes(i),1}, 1, '!(');
        BoolRuleRight{notGenes(i),1} = insertAfter(BoolRuleRight{notGenes(i),1}, BoolRuleRight{notGenes(i),1}(end), ')');
    end
    auxPreParsedBoolRules = strcat(leftGenes, ' = ', BoolRuleRight);
    auxPreParsedBoolRules = insertAfter(auxPreParsedBoolRules, '=', ' ');
end

% The rules are sorted to see if different OR rules are written in different
% lines.
auxPreParsedBoolRules = sort(auxPreParsedBoolRules);
leftGenes = strtrim(extractBefore(auxPreParsedBoolRules, '='));
leftGenesUnique = unique(leftGenes);

preParsedBoolRules = cell(size(leftGenesUnique));
if length(leftGenesUnique)<length(leftGenes)
    ind = 1;
    i = 1;
    while i <= size(auxPreParsedBoolRules, 1) % times a gene is in the 
        % left side of the rule. 
        % if it is more than one time, all its rules are joined in one 
        % line by | operators.
        Index = find(strcmp(strtrim(extractBefore(auxPreParsedBoolRules, '=')),strtrim(extractBefore(auxPreParsedBoolRules{i,1}, '='))));
        preParsedBoolRules(ind,1) = auxPreParsedBoolRules(i,1);
        if length(Index) > 1
            for j = 2:length(Index)
                preParsedBoolRules(ind,1) = strcat(preParsedBoolRules{ind,1}, {' '}, '|', extractAfter(auxPreParsedBoolRules{Index(j),1}, '='));
                i = i + 1;
            end
        end
        i = i + 1;
        ind = ind + 1;
    end
else
    preParsedBoolRules = auxPreParsedBoolRules;
end
end

function [BooleanRules, auxKO_Node] = AddingAuxKO_Nodes(BooleanRules, inputs, outputs)
% This function adds the corresponding aux KO node to each node. It makes
% possible the knockout of the nodes.
% In the case of the input nodes of the model, in order to facilitate the
% indexing of the auxKO_ nodes, it is added a new rule for each input.
%
BooleanRules2 = BooleanRules(~strcmp(regexprep(BooleanRules, ' =.*$', ''), outputs));
BooleanRules_out = BooleanRules(strcmp(regexprep(BooleanRules, ' =.*$', ''), outputs));
auxKO_Node = strcat(regexprep(BooleanRules2, ' =.*$', ''), '_auxKO');
BooleanRules_new = strcat(regexprep(BooleanRules2 ,'= ', '= ('), ') & ', {' '}, auxKO_Node);

if ~isempty(inputs)
    if iscell(inputs{1})
        inputs = reshape([inputs{:}],[],1);
    else
        inputs = reshape(inputs,[],1);
    end
    BooleanRulesInput = strcat(inputs, ' = ', {' '} , inputs, '_auxKO');
    auxKO_NodeInput = strcat(inputs, '_auxKO');
    BooleanRules = [BooleanRules_out; BooleanRules_new; BooleanRulesInput];
    auxKO_Node = unique(sort([auxKO_Node; auxKO_NodeInput]));
end

end

function model = modelParser(head, nodeParent, parent)
% This function builds auxiliary nodes for each layer of the Boolean rules.
% For example:
%     g1 & g2 & g3 is one layer in the model,   RXN = g1 + g2 + g3
%     (g1 | g2) & g3 is a two layer model:      RXN = (g1 | g2) + g3
%                                               (g1 | g2) = g1
%                                               (g1 | g2) = g2

reactions = struct();
reactions.rxns = {'aux'; 'aux'};
reactions.formulas = {'aux'; 'aux'};
reactions = modelParser_aux(reactions, head, nodeParent, parent);

reactions.rxns = reactions.rxns(3:end);
reactions.formulas = reactions.formulas(3:end);

model = reactions.formulas;

end

function reactions = modelParser_aux(reactions, head, nodeParent, parent)
% Name nodes according to position and layer inside the Boolean rule

node = cell(size(head.children));
if length(head.children) == 1 %If it has an only child, it is named as its parent.
    node{1} = parent;
else %If the node has more than one child, they will be named as : nodeparent_1, nodeparent_2...
    node = strcat(nodeParent,'_', strsplit(num2str(1:length(head.children))));
    if isa(head,'OrNode') %if it is an OrNode:
        for child = 1:length(head.children)
            reactions.rxns{end+1} = node{child};
            reactions.formulas{end+1} = [node{child},' -> ',parent];
        end
    else %if it is an AndNode
        reactions.rxns{end+1} = strjoin(node,' & ');
        reactions.formulas{end+1} = [strjoin(node,' + '),' -> ',parent];
    end
end

%  Expand model with every node, using OR/AND rules
for i = 1:length(head.children)
    child = head.children(i);
    % check if the child is final layer
    if isa(child,'LiteralNode')
        reactions.rxns{end+1} = [child.id, '_', node{i}];
        reactions.formulas{end+1} = [child.id, ' -> ', node{i}];
    else %If the child is not the last layer.
        reactions = modelParser_aux(reactions, child,node{i}, node{i});
    end
end

end

function [rxnFormulasReduced] = ReduceModel(rxnFormulas)
% As the algorithm generates redundant nodes, this function is used to
% remove these nodes.

LeftHand = regexprep(rxnFormulas, ' ->.*$', '');
RightHand = regexprep(rxnFormulas, '^.* -> ', '');

idx = contains(LeftHand, '+');
LeftHandCell(~idx) = cellfun(@(x)({{x}}), LeftHand(~idx));
if sum(idx)>0
    LeftHandCell(idx) = cellfun(@strsplit, LeftHand(idx), repmat({' + '}, sum(idx),1), 'UniformOutput', 0);
end

idx = contains(RightHand, '+');
RightHandCell(~idx) = cellfun(@(x)({{x}}), RightHand(~idx));
if sum(idx)>0
    RightHandCell(idx) = cellfun(@strsplit, RightHand(idx), repmat({' + '}, sum(idx),1), 'UniformOutput', 0);
end

nodes = [reshape(LeftHandCell,1,[]), reshape(RightHandCell,1,[])];
nodes = reshape(nodes,1,[]);
nodes = reshape(unique([nodes{:}]),[],1);

intermediateNodesRaw = nodes(startsWith(nodes, 'GPR_'));
intermediateNodes = false(size(intermediateNodesRaw));

for i = 1:length(intermediateNodes)
    if sum(cellfun(@length, strfind(RightHand, intermediateNodesRaw{i}))~=0)==1
        if sum(cellfun(@length, strfind(LeftHand, intermediateNodesRaw{i}))~=0)==1
            intermediateNodes(i) = 1;
        end
    end
end
intermediateNodes = intermediateNodesRaw(intermediateNodes);

% generate synonims
intermediateNodesTranslation = cell(numel(intermediateNodes),2);
intermediateNodesTranslation(:,1) = intermediateNodes;
for i = 1:length(intermediateNodes)
    if strcmp(LeftHand,intermediateNodes{i})
        intermediateNodesTranslation{i,2} = RightHand{strcmp(LeftHand, intermediateNodes{i})};
    else
        intermediateNodesTranslation{i,2} = LeftHand{strcmp(RightHand, intermediateNodes{i})};
    end
end

for i = 1:length(intermediateNodes)
    idx = find(cellfun(@length,strfind(LeftHand, intermediateNodesTranslation{i,1}))>0);
    for j = 1:1:length(idx)
        LeftHandCell{idx(j)} = regexprep(LeftHandCell{idx(j)}, intermediateNodesTranslation{i,1}, intermediateNodesTranslation{i,2});
    end
    
    idx = find(cellfun(@length,strfind(RightHand, intermediateNodesTranslation{i,1}))>0);
    for j=1:1:length(idx)
        RightHandCell{idx(j)} = regexprep(RightHandCell{idx(j)}, intermediateNodesTranslation{i,1}, intermediateNodesTranslation{i,2});
    end
end

LeftHandCell = cellfun(@sort, LeftHandCell, 'UniformOutput', 0);
RightHandCell = cellfun(@sort, RightHandCell, 'UniformOutput', 0);

idx = cellfun(@length, LeftHandCell)>1;
LeftHand(~idx) = [LeftHandCell{~idx}];

if sum(idx)>0
    pos = find(idx>0);
    for i = 1:sum(idx)
        LeftHand(pos(i)) = cellfun(@strjoin, LeftHandCell(pos(i)), repmat({' + '}, sum(idx(pos(i))),1), 'UniformOutput', 0);
    end
end

idx = cellfun(@length, RightHandCell)>1;
RightHand(~idx) = [RightHandCell{~idx}];
if sum(idx)>0
    RightHand(idx) = cellfun(@strjoin, RightHandCell(idx), repmat({' + '}, sum(idx),1), 'UniformOutput', 0);
end

idx = ~strcmp(LeftHand, RightHand);
LeftHand = LeftHand(idx);
RightHand = RightHand(idx);

rxnFormulasReduced = strcat(LeftHand, {'  -> '}, RightHand, {' '});
end

function [gMCSs] = generateModel(RxnFormulas, auxKO_Nodes, target, n_mcs, inputs, forceLength, numWorkers)
for i = 1:size(RxnFormulas,1)
	abr{i,1} = ['R' num2str(i)];
	names{i,1} = ['Reaction_' num2str(i)];
end
model = createModel(abr, names, RxnFormulas, 'printLevel', 0);
rxn_set = cell(length(auxKO_Nodes), 1);

%input reaction to every aux KO node
for i = 1:length(auxKO_Nodes)
    model = addReaction(model, strcat('r_', auxKO_Nodes{i,1}), ...
        'reactionName', strcat('r_', auxKO_Nodes{i,1},'_input'), 'metaboliteList', {strcat(auxKO_Nodes{i,1}, '[c]')},...
        'lowerBound', 0, 'upperBound', 10000, 'stoichCoeffList', 1, 'printLevel', 0);
    rxn_set{i,1} = strcat('r_', auxKO_Nodes{i,1});
    
end

for i = 1:length(inputs)
    model = addReaction(model, strcat('r_', inputs{i,1}, '_off_input'), ...
        'reactionName', strcat('r_', inputs{i,1}, '_off_input'), 'metaboliteList', {strcat(inputs{i,1}, '_off[c]')},...
        'lowerBound', 0, 'upperBound', 100, 'stoichCoeffList', 1, 'printLevel', 0);    
end

%output reaction for the targer node
model = addReaction(model, char(strcat('r_', target, '_on')), ...
        'reactionName', char(strcat('r_', target, '_output')), 'metaboliteList', {char(strcat(target, '_on[c]'))},...
        'lowerBound', 0, 'upperBound', 10000, 'stoichCoeffList', -1, 'printLevel', 0);

    
[~, auxKO_on_pos] = unique(extractBefore(auxKO_Nodes, '_auxKO'), 'stable' );
auxKO_off_pos = setdiff(1:numel(auxKO_Nodes), auxKO_on_pos);
    
model = changeObjective(model, model.rxns(end));
% integrality_tolerance = 1e-5;
% M = 1e3;    % Big Value
alpha = 1;  % used to relate the lower bound of v variables with z variables
c = 1e-3;   % used to activate w variable
b = 1e-3;   % used to activate KnockOut constraint
% phi = 1000; % b/c;

% Build the K Matrix
[~, n_ini_rxns] = size(model.S);
K = speye(n_ini_rxns);

% Splitting
S = model.S;
% K_ind = model.rxns;
% n_K_ind = length(K_ind);
[n_mets, n_rxns] = size(S);
nbio = model.c;
t = zeros(n_rxns, 1);
t(find(nbio)) = 1;

% Permit only KOs in  exchanges of auxKO_ nodes.
rxn_set = unique(rxn_set);
tmp_set = cellfun(@ismember, model.rxns, repmat({rxn_set}, n_ini_rxns, 1), 'UniformOutput', false);
pos_set = find(cell2mat(tmp_set));
K = K(pos_set, :);
K_ind = model.rxns(pos_set);
n_K_ind = length(K_ind);


% if isempty(KO)
% ENUMERATE MCSs
% Define variables
var.u = 1:n_mets;
var.vp = var.u(end)+1:var.u(end)+n_K_ind;
var.w = var.vp(end)+1:var.vp(end)+1;
var.zp = var.w(end)+1:var.w(end)+n_K_ind;
var.zw = var.zp(end)+1:var.zp(end)+1;
n_vars = var.zw(end);
var_group.v = [var.vp var.w];
var_group.z = [var.zp var.zw];

% Define constraints
cons.Ndual = 1:size(S, 2);
cons.forceBioCons = cons.Ndual(end)+1:cons.Ndual(end)+1;
cons.auxKO = cons.forceBioCons(end)+1:cons.forceBioCons(end)+length(auxKO_on_pos);
cons.forceLength = cons.auxKO(end)+1:cons.auxKO(end)+1;
n_cons = cons.forceLength(end);

% Cplex - A matrix
A = sparse(zeros(n_cons, n_vars));
A(cons.Ndual, var.u) = S';
A(cons.Ndual, var.vp) = K';
A(cons.Ndual, var.w) = -t;
A(cons.forceBioCons, var.w) = -b;
A(cons.auxKO, var.zp(auxKO_on_pos)) = speye(length(auxKO_on_pos));
A(cons.auxKO, var.zp(auxKO_off_pos)) = speye(length(auxKO_off_pos));
if forceLength == 1
    A(cons.forceLength, var.zp(auxKO_on_pos)) = 1;
end

% Cplex - rhs and lhs vectors
rhs = zeros(n_cons, 1);
rhs(cons.Ndual, 1) = inf;
rhs(cons.forceBioCons) = -c;
rhs(cons.auxKO) = 1;
if forceLength == 1
    rhs(cons.forceLength) = 1;
end

lhs = zeros(n_cons, 1);
lhs(cons.Ndual, 1) = 0;
lhs(cons.forceBioCons) = -inf;
lhs(cons.auxKO) = 1;
if forceLength == 1
    lhs(cons.forceLength) = 1;
end

% Cplex - ub and lb vectors
ub(var.u, 1) = inf;
ub(var.vp) = inf;
ub(var.w) = inf;
ub(var.zp) = 1;
ub(var.zw) = 1;
lb(var.u, 1) = -inf;
lb(var.vp) = 0;
lb(var.w) = 0;
lb(var.zp) = 0;
lb(var.zw) = 0;

% Cplex - obj vector
obj(var.u, 1) = 0;
obj(var.vp) = 0;
obj(var.w) = 0;
obj(var.zp(auxKO_on_pos)) = 1;
obj(var.zp(auxKO_off_pos)) = 0;
obj(var.zw) = 0;

% Cplex - ctype vector
ctype(var.u) = 'C';
ctype(var.vp) = 'C';
ctype(var.w) = 'C';
ctype(var.zp) = 'B';
ctype(var.zw) = 'B';

% Cplex - sense of the optimization
sense = 'minimize';

% Cplex - Introduce all data in a Cplex structure
cplex = Cplex('MCS');
Model = struct();
[Model.A, Model.rhs, Model.lhs, Model.ub, Model.lb, Model.obj, Model.ctype, Model.sense] = deal(A, rhs, lhs, ub, lb, obj, ctype, sense);
cplex.Model = Model;

% Cplex Indicators
% z = 1  -->  v >= alpha
for ivar = 1:length(var_group.z)
    a = zeros(n_vars, 1);
    a(var_group.v(ivar)) = 1;
    cplex.addIndicators(var_group.z(ivar), 0, a, 'G', alpha);
end

% Cplex Indicators
% z = 0  -->  v <= 0
for ivar = 1:length(var_group.z)
    a = zeros(n_vars, 1);
    a(var_group.v(ivar)) = 1;
    cplex.addIndicators(var_group.z(ivar), 1, a, 'L', 0);
end

k = 0;
sP = struct();
timelimit = 5*60;
% numWorkers = 2;
integrality_tolerance = 1e-5;
[sP.mip.tolerances.integrality, sP.mip.strategy.heuristicfreq, sP.mip.strategy.rinsheur] = deal(integrality_tolerance, 1000, 50);
[sP.emphasis.mip, sP.output.clonelog, sP.timelimit, sP.threads] = deal(4, -1, max(10, timelimit), numWorkers);
[sP.preprocessing.aggregator, sP.preprocessing.boundstrength, ...
    sP.preprocessing.coeffreduce, sP.preprocessing.dependency, ...
    sP.preprocessing.dual, sP.preprocessing.fill,...
    sP.preprocessing.linear, sP.preprocessing.numpass, ...
    sP.preprocessing.presolve, sP.preprocessing.reduce,..., ...
    sP.preprocessing.relax, sP.preprocessing.symmetry] = deal(50, 1, 2, 1, 1, 50, 1, 50, 1, 3, 1, 1);
cplex = setCplexParam(cplex, sP);

% Calculation of gMCSs:
largest_mcs = 0;
% max_len_mcs = length(auxKO_on_pos);
%     max_len_mcs = 100;
max_len_mcs = 5;

while largest_mcs <= max_len_mcs && k < n_mcs && cplex.Model.rhs(cons.forceLength) <= max_len_mcs
    cplex.Param.mip.limits.populate.Cur = 40;
    cplex.Param.threads.Cur = 12;
    cplex.Param.mip.pool.relgap.Cur = 0.1;
    cplex.populate();
    n_pool = size(cplex.Solution.pool.solution, 1);
    if n_pool ~= 0
        solution = cplex.Solution.pool.solution;
        for j = 1:n_pool
            k = k+1;
            gMCSs{k, 1} = extractBetween(K_ind((solution(j).x(var.zp(auxKO_on_pos)))>0.9), "r_", "_auxKO");
            n_cons = n_cons+1;
            sol = solution(j).x(var.zp(auxKO_on_pos))>0.9;
            cplex.Model.A(n_cons, var.zp(auxKO_on_pos)) = sparse(double(sol));
            cplex.Model.rhs(n_cons) = sum(sol)-1;
            cplex.Model.lhs(n_cons) = 0;
        end
    else
        if forceLength == 1
            cplex.Model.rhs(cons.forceLength) = cplex.Model.rhs(cons.forceLength)+1;
            cplex.Model.lhs(cons.forceLength) = cplex.Model.lhs(cons.forceLength)+1;
        else
            return;
        end
    end
    try largest_mcs = max(cellfun(@length, gMCSs)); end
end
if ~exist('gMCSs', 'var')
    gMCSs = {};
end
end


