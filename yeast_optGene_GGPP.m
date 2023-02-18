% OptGene workflow to increase GGPP production in yeast
% author: Koray Malci

initCobraToolbox    % to start the Cobra Toolbox

changeCobraSolver('gurobi', 'all');   %to change the solver
%%
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

% objectiveCoeff = 1.0;
% model = changeObjective(model,'r_2111',objectiveCoeff);
% opt_WT = optimizeCbModel(model);
% gr_WT=opt_WT.f;
% solution_max = optimizeCbModel(model,'max');
% max_biomass=solution_max.f;
% solution_min = optimizeCbModel(model,'min');
% min_biomass=solution_min.f;
% printFluxVector(model,opt_WT.x,'true','false');
% fprintf('The minimum and maximum solution for growth is %.5f and %.5f respectively\n', min_biomass, max_biomass);
% fprintf('The growth rate is %.5f \n', gr_WT);


%Change the objective to GGPP porudction

% model = changeObjective(model, 'EX_s_0189[c]', objectiveCoeff); 
% fbaWT = optimizeCbModel(model);
% fbaWTMin = optimizeCbModel(model, 'min');
% fbaWTMax = optimizeCbModel(model, 'max');
% min_ggpp_FluxWT = fbaWTMin.f;
% max_ggpp_FluxWT = fbaWTMax.f;
% printFluxVector(model,fbaWT.x,'true','false');
% fprintf('The minimum and maximum production of GGPP is %.5f and %.5f respectively\n', min_ggpp_FluxWT, max_ggpp_FluxWT);

%%
% Selecting the target genes to be deleted. After single gene 
% deletion, the growth rate between WT and KO strain is set to "1"
% indicating the same growt rate (without negative effect)

gr=singleGeneDeletion(model,'FBA');

target_genes={model.genes{gr==1}};
%%
% Finding optGene sets

fprintf('\n Finding optGene sets\n\n')
previousSolutions = cell(10, 1);
contPreviousSolutions = 1;
nIter = 0;
threshold=3;
while nIter < threshold
    fprintf('...Performing optGene analysis...\n')
    %optGene algorithm is run with the following options:
    [~, ~, ~, optGeneSol] = optGene(model, 'EX_s_0189[c]', 'r_1714', target_genes, 'MaxKOs', 3, 'Generations', 100);
    
    SET_M1 = optGeneSol.geneList;
    
    if ~isempty(SET_M1)
        previousSolutions{contPreviousSolutions} = SET_M1;
        contPreviousSolutions = contPreviousSolutions + 1;
        %printing results
        fprintf('optGene found a knockout set of large %d composed by ', length(SET_M1));
        for j = 1:length(SET_M1)
            if j == 1
                fprintf('%s ',SET_M1{j});
            elseif j == length(SET_M1)
                fprintf('and %s',SET_M1{j});
            else
                fprintf(', %s ',SET_M1{j});
            end
        end
        fprintf('\n');
        fprintf('...Performing coupling analysis...\n');
        [type, maxGrowth, maxProd, minProd] = analyzeOptKnock(model, optGeneSol.geneList, 'EX_s_0189[c]', biomass, 1);
        fprintf('The solution is of type: %s\n',type);
        fprintf('The maximun growth rate after optimizacion is %.5f\n', maxGrowth);
        fprintf('The maximun and minimun production of GGPP after optimization is %.5f and %.5f, respectively \n\n', minProd, maxProd);
        
    else
        if nIter  ==  1
            fprintf('optGene was not able to found an optGene set\n');
        else
            fprintf('optGene was not able to found additional optGene sets\n');
        end
        break;
    end
    nIter = nIter + 1;
    
end