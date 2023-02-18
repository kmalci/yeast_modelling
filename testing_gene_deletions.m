
% Gene list to be deleted from the model one by one

% author: Koray Malci

%Gene deletions for optGene acetyl-coA using galactose

geneList= {'YER019W','YFR047C','YPR127W','YHR002W','YMR303C','YIL006W','YKL024C','YPR062W'};


% Load the model to be used

model = readCbModel('yeast-GEM.mat');  
biomass = 'r_2111';

model = addExchangeRxn(model, {'s_0373[c]'}, 0, 1000);   % adding acetyl-Coa (s_0373) exchange reaction
%model = addExchangeRxn(model, {'s_0189[c]'}, 0, 1000);   % adding GGPP (s_0189) exchange reaction

%SETTING SPECIFIC CONSTRAINTS
% prespecified amount of carbon source uptake 10 mmol/grDW*hr

model = changeRxnBounds(model, 'r_1714', 0, 'b');   %glucose exchange
model = changeRxnBounds(model, 'r_1710', -10, 'b');  %galactose exchange

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

model = changeObjective(model, biomass);
fbaWT = optimizeCbModel(model);
growthRateWT = fbaWT.f;

model = changeObjective(model, 'EX_s_0373[c]'); 
fbaWTMin = optimizeCbModel(model, 'min');
fbaWTMax = optimizeCbModel(model, 'max');
min_acoA_FluxWT = fbaWTMin.f;
max_acoA_FluxWT = fbaWTMax.f;
%printFluxVector(model,fbaWT.x,'true','false');

fprintf('The minimum and maximum production of acetyl-coA before deletion is %.10f and %.10f respectively\n', min_acoA_FluxWT, max_acoA_FluxWT);
fprintf('The growth rate before deletion is %.10f \n', growthRateWT);

%Deletion loop

for i = 1:length(geneList)
    
    genename=geneList{i};
    [model, hasEffect, constrRxnNames, ~] = deleteModelGenes(model, geneList{i});
    model = changeObjective(model, biomass);
    fbaWT = optimizeCbModel(model);
    growthRateWT = fbaWT.f;

    model = changeObjective(model, 'EX_s_0373[c]'); 
    fbaWTMin = optimizeCbModel(model, 'min');
    fbaWTMax = optimizeCbModel(model, 'max');
    min_acoA_FluxWT = fbaWTMin.f;
    max_acoA_FluxWT = fbaWTMax.f;
    %fprintf('\n\n The exchange routes \n')
    %printFluxVector(model,fbaWT.x,'true','false');
    fprintf('..................');
    fprintf('\n\n The target gene is %s\n', genename);
    if hasEffect==1
        fprintf('\n\n Deletion of %s has an effect on the model \n', genename);
    else
        fprintf('\n\n Deletion of %s has NO effect on the model\n', genename);
    end
    
    if ~isempty(constrRxnNames)
        fprintf('the constrained reactions:')
        display(constrRxnNames);
        printFluxBounds(model, constrRxnNames);
    else
        fprintf('There is NO constrained reactions')
    end
   
    
    fprintf('the min and max production of acetyl-coA after deletion is %.10f and %.10f respectively\n',min_acoA_FluxWT, max_acoA_FluxWT);
    fprintf('The growth rate after deletion of is %.10f \n\n\n', growthRateWT);
    
    
    %reload the model
    
    model = readCbModel('yeast-GEM.mat');
    biomass = 'r_2111';
    
    model = addExchangeRxn(model, {'s_0373[c]'}, 0, 1000);
    
    %SETTING SPECIFIC CONSTRAINTS
    % prespecified amount of carbon source uptake 10 mmol/grDW*hr
    
    model = changeRxnBounds(model, 'r_1714', 0, 'b');   %glucose exchange
    model = changeRxnBounds(model, 'r_1710', -10, 'b');  %galactose exchange
    
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

    model = changeObjective(model, biomass);
    fbaWT = optimizeCbModel(model);
    
end