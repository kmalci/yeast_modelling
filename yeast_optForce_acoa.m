% OptForce workflow to increase Acetyl-Coa or GGPP 
% production in yeast

% author: Koray Malci

initCobraToolbox    % to start the Cobra Toolbox

changeCobraSolver('gurobi', 'all');   %to change the solver


% Introduction the model and the constraints 
model = readCbModel('yeast-GEM.mat');  
biomass = 'r_2111';
model.c(strcmp(model.rxns, biomass)) = 1;   %biomass is the objective function



% Adding exchange reaction for the target molecule
model = addExchangeRxn(model, {'s_0373[c]'}, 0, 1000);   % adding acetyl-Coa (s_0373) exchange reaction
%model = addExchangeRxn(model, {'s_0189[c]'}, 0, 1000);   % adding GGPP (s_0189) exchange reaction



% SETTING SPECIFIC CONSTRAINTS
% prespecified amount of carbon source uptake 10 mmol/grDW*hr

model = changeRxnBounds(model, 'r_1714', -10, 'b');   %glucose exchange
%model = changeRxnBounds(model, 'r_1710', -10, 'b');  %galactose exchange



% Unconstrained uptake routes for inorganic phosphate, sulfate,
% ammonia, oxygen, 

model = changeRxnBounds(model, 'r_2005', -1000, 'l');   %phosphate
model = changeRxnBounds(model, 'r_2060', -1000, 'l');   %sulphate
model = changeRxnBounds(model, 'r_1654', -1000, 'l');   %ammonium
model = changeRxnBounds(model, 'r_1992', -1000, 'l');   %oxygen
model = changeRxnBounds(model, 'r_2049', -1000, 'l');   %sodium
model = changeRxnBounds(model, 'r_2020', -1000, 'l');   %potassium
model = changeRxnBounds(model, 'r_4593', -1000, 'l');   %Cl (cloride)
model = changeRxnBounds(model, 'r_4594', -1000, 'l');   %Cu (copper)
model = changeRxnBounds(model, 'r_4595', -1000, 'l');   %Mn (mangane)
model = changeRxnBounds(model, 'r_4596', -1000, 'l');   %Zn (zinc)
model = changeRxnBounds(model, 'r_4597', -1000, 'l');   %Mg (magnesium)
model = changeRxnBounds(model, 'r_4600', -1000, 'l');   %Ca (calcium)
model = changeRxnBounds(model, 'r_1861', -1000, 'l');   %Fe (iron)
model = changeRxnBounds(model, 'r_1832', -1000, 'l');   %H (hydrogen)


% Secretion routes  for acetate, carbon dioxide, ethanol, glycolaldehyde,
% diphosphate, water, glycerol and acetaldehyde are enabled
model = changeRxnBounds(model, 'r_1634', 1000, 'u');    %acetate
model = changeRxnBounds(model, 'r_1672', 1000, 'u');    %co2
model = changeRxnBounds(model, 'r_1761', 1000, 'u');    %ethanol
model = changeRxnBounds(model, 'r_1814', 1000, 'u');    %glycolaldehyde
model = changeRxnBounds(model, 'r_4527', 1000, 'u');    %diphopshate
model = changeRxnBounds(model, 'r_2100', 1000, 'u');    %water
model = changeRxnBounds(model, 'r_1808', 1000, 'u');    %glycerol
model = changeRxnBounds(model, 'r_1631', 1000, 'u');    %acetaldehyde


% Calculating the growth rate
growthRate = optimizeCbModel(model);
fprintf('The maximum growth rate is %1.5f', growthRate.f);


% Calculating the maximum production of the molecule of interest
model = changeObjective(model, 'EX_s_0373[c]');
max_acoa = optimizeCbModel(model);
fprintf('The maximum production rate of acetyl-coA is %1.5f', max_acoa.f);


% Constraining the WT strain to maximum biomass (using the max. growth rate)         
constrWT = struct('rxnList', {{'r_2111'}}, 'rxnValues', 0.86, 'rxnBoundType', 'b');

% Constraining the MUtant strain to maximum production (using the max. production rate)
constrMU = struct('rxnList', {{'r_2111', 'EX_s_0373[c]'}}, 'rxnValues', [0, 1.70], 'rxnBoundType', 'bb');

% Flux variability analysis (FVA) for both constrained strains 
[minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ~, ~] = FVAOptForce(model, constrWT, constrMU);

% Defining the run ID
runID = 'OptForce_acetylcoA_glucose';

% Defining the constraints (the constraints from the previous processes)
constrOpt = struct('rxnList', {{'r_1714', 'r_2111', 'EX_s_0373[c]'}}, 'values', [-10, 0, 1.70]);


% Finding mustLSet
[mustLSet, pos_mustL] = findMustL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, 'runID', runID, 'outputFolder', 'OutputsFindMustL', ...
    'outputFileName', 'MustL' , 'printExcel', 1,  'printText', 1, 'printReport', 1, 'keepInputs', 1, 'verbose', 0);


% Finding mustUSet
[mustUSet, pos_mustU] = findMustU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, 'runID', runID, 'outputFolder', 'OutputsFindMustU', ...
    'outputFileName', 'MustU' , 'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, 'verbose', 0);


% Finding second order must sets

% First, the reactions found in the previous step and the exchange
% reactions are excluded
exc=findExcRxns(model);
exc_rxns=transpose({model.rxns{exc}});
added_exc_rxns = model.rxns(cellfun(@isempty, strfind(model.rxns, 'EX_')) == 0);
excludedRxns = unique([mustUSet; mustLSet; exc_rxns; added_exc_rxns]);

[mustUU, pos_mustUU, mustUU_linear, pos_mustUU_linear] = findMustUU(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, 'excludedRxns', excludedRxns,'runID', runID, ...
               'outputFolder', 'OutputsFindMustUU', 'outputFileName', 'MustUU', 'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, 'verbose', 1);
           
[mustLL, pos_mustLL, mustLL_linear, pos_mustLL_linear] = findMustLL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, 'excludedRxns', excludedRxns,'runID', runID, ...
    'outputFolder', 'OutputsFindMustLL', 'outputFileName', 'MustLL', 'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, 'verbose', 1);

% Optional set
% [mustUL, pos_mustUL, mustUL_linear, pos_mustUL_linear] = findMustUL(model, minFluxesW, maxFluxesW, 'constrOpt', constrOpt, 'excludedRxns', excludedRxns,'runID', runID, ...
%     'outputFolder', 'OutputsFindMustUL', 'outputFileName', 'MustUL', 'printExcel', 1, 'printText', 1, 'printReport', 1, 'keepInputs', 1, 'verbose', 1);



% optForce run

mustU = unique(union(mustUSet, mustUU));
mustL = unique(union(mustLSet, mustLL));

targetRxn = 'EX_s_0373[c]';
biomassRxn = 'r_2111';

k = 1;
nSets = 1;
constrOpt = struct('rxnList', {{'r_1714','r_2111'}}, 'values', [-10, 0]);

[optForceSets, posOptForceSets, typeRegOptForceSets, flux_optForceSets] = optForce(model, targetRxn, biomassRxn, mustU, mustL, minFluxesW, maxFluxesW, minFluxesM, maxFluxesM, ...
    'k', k, 'nSets', nSets, 'constrOpt', constrOpt, 'runID', runID, 'outputFolder', 'OutputsOptForce', 'outputFileName', 'OptForce', 'printExcel', 1, 'printText', 1, ...
    'printReport', 1, 'keepInputs', 1, 'verbose', 1);


fprintf('The reaction(s) found by optForce:')
disp(optForceSets)




