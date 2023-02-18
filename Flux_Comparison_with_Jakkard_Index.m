% comparison of fluxes in two models using Jakkard Index
% author: Koray Malci

initCobraToolbox
changeCobraSolver('gurobi', 'all');

%%

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

%%
% objectiveCoeff = 1.0;
% model = changeObjective(model,'EX_s_0189[c]',objectiveCoeff);
% opt_WT = optimizeCbModel(model);
% %gr_WT=opt_WT.f;
% solution_max = optimizeCbModel(model,'max');
% max_ggpp=solution_max.f;
% solution_min = optimizeCbModel(model,'min');
% min_ggpp=solution_min.f;
% printFluxVector(model,opt_WT.x,'true','false');
% fprintf('The minimum and maximum production of ggpp is %.10f and %.10f respectively\n', min_ggpp, max_ggpp);
%fprintf('The growth rate is %.5f \n', gr_WT);

%%
ibm = find(ismember(model.rxns, biomass));
model.lb(ibm)=0.1;
model.c(:)=0;
%%
[del_model, ~, deleted_reactions, ~] = deleteModelGenes(model, {'YBR208C','YOR155C','YPL092W'});

display(deleted_reactions); %the following reaction list will be knocked-out if these genes are deleted   
%%
[minWT, maxWT] = fluxVariability(model);
[minDel, maxDel] = fluxVariability(del_model);

fprintf('Max. biomass of WT: %.10f/h.\n', maxWT(ibm));
fprintf('Max. biomass of KO: %.10f/h.\n\n', maxDel(ibm));

J = fvaJaccardIndex([minWT, minDel], [maxWT, maxDel]);
fprintf('Mean Jaccard index = %.4f.\n', mean(J));
%%
E = [(maxWT - minWT)/2 (maxDel - minDel)/2];
Y = [minWT minDel] + E;
X = [(1:length(Y)) - 0.1; (1:length(Y)) + 0.1]';

[~, xj] = sort(J);

f1 = figure;
if strcmp(version('-release'), '2016b')
    errorbar(X, Y(xj, :), E(xj, :), 'linestyle', 'none', 'linewidth', 2, 'capsize', 0);
else
    %errorbar(X, Y(xj, :), E(xj, :), 'linestyle', 'none', 'linewidth', 2);
    hold on
    errorbar(X(:,1), Y(xj, 1), E(xj, 1), 'linestyle', 'none', 'linewidth', 2,'Color','b');
    errorbar(X(:,2), Y(xj, 2), E(xj, 2), 'linestyle', 'none', 'linewidth', 2,'Color','r');
end
set(gca, 'xlim', [0, length(Y) + 1])

xlabel('Reaction')
ylabel('Flux range (mmol/gDW/h)')
ylim([-1500,1500])
xlim([-20,5000])
yyaxis right
plot(J(xj),'linewidth', 2)
legend('WT', 'KO', 'Jaccard','location', 'northoutside', ...
       'orientation', 'horizontal')
ylabel('Jaccard index')