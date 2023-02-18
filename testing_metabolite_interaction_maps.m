%initCobraToolbox
%changeCobraSolver ('gurobi', 'all');
model = readCbModel('yeast-GEM.mat');  
biomass = 'r_2111';

model = addExchangeRxn(model, {'s_0373[c]'}, 0, 1000);   % adding acetyl-Coa (s_0373) exchange reaction
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

mets={'s_0373[c]','s_0367[c]','s_0218[c]','s_0028[c]','s_0018[c]','s_0019[c]','s_0943[c]',...
    's_1376[c]','s_0745[c]','s_0190[c]','s_0831[m]','s_0189[c]'}; % metabolites of interest.

%graphObj=createMetIntrcNetwork(model,mets);

%%
%FBAsolution = optimizeCbModel(model);
%fluxes=FBAsolution.x;
%objectiveCoeff = 1.0;
%model = changeObjective(model, 'EX_s_0373[c]', objectiveCoeff); 

ibm1 = find(ismember(model.rxns, biomass));
model.lb(ibm1)=0.05;
model.c(:)=0;
fbaWT = optimizeCbModel(model);
fluxesWT = fbaWT.x;


graphObjWT=createMetIntrcNetwork(model,mets,'fluxes',fluxesWT);
set(gcf,'Visible','on'); % produce figure as pop up since live editor does

%%

[del_model, ~, deleted_reactions, ~] = deleteModelGenes(model, {'YBR208C','YOR155C','YPL092W'});

display(deleted_reactions); %the following reaction list will be knocked-out if these genes are deleted

%del_model = changeObjective(del_model, 'EX_s_0373[c]', objectiveCoeff);
ibm2 = find(ismember(del_model.rxns, biomass));
del_model.lb(ibm2)=0.05;
del_model.c(:)=0;
fbaKO = optimizeCbModel(del_model);
fluxesKO = fbaKO.x;

graphObjKO=createMetIntrcNetwork(del_model,mets,'fluxes',fluxesKO);
set(gcf,'Visible','on'); % produce figure as pop up since live editor does

%%

% interaction map for all metabolites

% graphObj=createMetIntrcNetwork(model,model.mets);

%%

%test section 

% ibm1 = find(ismember(model.rxns, biomass));
% model.lb(ibm1)=0.05;
% model.c(:)=0;
% 
% if ismember(1, model.c(:)) 
%     fprintf('\n yes \n')
% else
%     
%     fprintf('\n no \n ')
% end
% 
% objectiveCoeff = 1.0;
% 
% model = changeObjective(model, biomass, objectiveCoeff);
% 
% if ismember(1, model.c(:)) 
%     fprintf('\n now yes \n')
% else
%     
%     fprintf('\n still no \n')
% end
% 
