% OptKnock workflow to increase GGPP production in yeast
% author: Koray Malci

initCobraToolbox    % to start the Cobra Toolbox
%%
%changeCobraSolver('gurobi', 'all');   %to change the solver

model = readCbModel('yeast-GEM.mat');  
biomass = 'r_2111';

model = addExchangeRxn(model, {'s_0189[c]'}, 0, 1000);   % adding GGPP (s_0189) exchange reaction

%SETTING SPECIFIC CONSTRAINTS
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

%%

%FBA Solutions

objectiveCoeff = 1.0;
model = changeObjective(model,'r_2111',objectiveCoeff);
opt_WT = optimizeCbModel(model);
gr_WT=opt_WT.f;
solution_max = optimizeCbModel(model,'max');
max_biomass=solution_max.f;
solution_min = optimizeCbModel(model,'min');
min_biomass=solution_min.f;
printFluxVector(model,opt_WT.x,'true','false');
fprintf('The minimum and maximum solution for growth is %.5f and %.5f respectively\n', min_biomass, max_biomass);
fprintf('The growth rate is %.5f \n', gr_WT);


%Change the objective to GGPP porudction

model = changeObjective(model, 'EX_s_0189[c]', objectiveCoeff); 
fbaWT = optimizeCbModel(model);
fbaWTMin = optimizeCbModel(model, 'min');
fbaWTMax = optimizeCbModel(model, 'max');
min_ggpp_FluxWT = fbaWTMin.f;
max_ggpp_FluxWT = fbaWTMax.f;
printFluxVector(model,fbaWT.x,'true','false');
fprintf('The minimum and maximum production of GGPP before optKnock is %.10f and %.10f respectively\n', min_ggpp_FluxWT, max_ggpp_FluxWT);

%%
% Selecting the target genes to be deleted. After single gene 
% deletion, the growth rate between WT and KO strain is set to "1"
% indication the same growt rate (without negative effect)

gr=singleRxnDeletion(model,'FBA');

selectedRxnList={model.rxns{gr==1}};
%%

% Set optKnock options
% The exchange of GGPP  will be the objective of the outer problem
options = struct('targetRxn', 'EX_s_0189[c]', 'numDel', 3);

% We will impose that biomass should be at least 50% of the biomass of wild-type
% 'sense' takes one of 3 values G;greater, L;lower, E;equal

constrOpt = struct('rxnList', {{biomass}}, 'values', 0.2*fbaWT.f, 'sense', 'G');        

% We will try to find 10 optKnock sets of a maximun length of 3

threshold=5;
previousSolutions = cell(100, 1);
contPreviousSolutions = 1;

nIter = 1;
while nIter < threshold
    fprintf('...Performing optKnock analysis...')
    if isempty(previousSolutions{1})
        optKnockSol = OptKnock(model, selectedRxnList, options, constrOpt);
        %optKnockSol = OptKnock(model, selectedRxnList, options);
    else
        optKnockSol = OptKnock(model, selectedRxnList, options, constrOpt, previousSolutions);
        %optKnockSol = OptKnock(model, selectedRxnList,options);
    end
    
    % determine GGPP production and growth rate after optimization
    ggpp_FluxM1 = optKnockSol.fluxes(strcmp(model.rxns, 'EX_s_0189[c]'));
    
    growthRateM1 = optKnockSol.fluxes(strcmp(model.rxns, biomass));
    
    setM1 = optKnockSol.rxnList;
    
    if ~isempty(setM1)
        previousSolutions{contPreviousSolutions} = setM1;
        contPreviousSolutions=contPreviousSolutions + 1;
        %printing results
        fprintf('optKnock found a optKnock set of large %d composed by ',length(setM1));
        for j = 1:length(setM1)
            if j == 1
                fprintf('%s',setM1{j});
            elseif j == length(setM1)
                fprintf(' and %s',setM1{j});
            else
                fprintf(', %s',setM1{j});
            end
        end
        fprintf('\n');
        fprintf('The production of GGPP  after optimization is %.10f \n', ggpp_FluxM1);
        fprintf('The growth rate after optimization is %.10f \n', growthRateM1);
        
        fprintf('...Performing coupling analysis...\n');
        [type, maxGrowth, maxProd, minProd] = analyzeOptKnock(model, setM1, 'EX_s_0189[c]');
        fprintf('The solution is of type: %s\n', type);
        fprintf('The maximun growth rate given the optKnock set is %.10f\n', maxGrowth);
        fprintf('The maximun and minimun production of GGPP  given the optKnock set is %.10f and %.10f, respectively \n\n', minProd, maxProd);
        if strcmp(type,'growth coupled')
            singleProductionEnvelope(model, setM1, 'EX_s_0189[c]', biomass, 'savePlot', 1, 'showPlot', 1, 'fileName', ['acoA_ex2_' num2str(nIter)], 'outputFolder', 'OptKnockResults');
        end
    else
        if nIter == threshold
            fprintf('optKnock was not able to found an optKnock set\n');
        else
            fprintf('optKnock was not able to found additional optKnock sets\n');
        end
        break;
    end
    nIter = nIter + 1;
end
%%
%publish('Yeast_optKnock_ggpp.m','pdf')     %enable it when using .m file