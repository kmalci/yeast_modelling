% construction metabolite interaction maps
% author: Koray Malci

initCobraToolbox
changeCobraSolver ('gurobi', 'all');
model = readCbModel('yeast-GEM.mat');  
biomass = 'r_2111';

model = addExchangeRxn(model, {'s_0373[c]'}, 0, 1000);   % adding acetyl-Coa (s_0373) exchange reaction
model = addExchangeRxn(model, {'s_0189[c]'}, 0, 1000);   % adding GGPP (s_0189) exchange reaction


%SETTING SPECIFIC CONSTRAINTS
% prespecified amount of carbon source uptake 10 mmol/grDW*hr

model = changeRxnBounds(model, 'r_1714', 0, 'b');   %glucose exchange
model = changeRxnBounds(model, 'r_1710', -20, 'b');  %galactose exchange

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
    's_1376[c]','s_0745[c]','s_0190[c]','s_0831[m]','s_0189[c]','s_1311[c]'}; % metabolites of interest.

%graphObj=createMetIntrcNetwork(model,mets);

%%
% FBAsolution = optimizeCbModel(model);
% fluxes=FBAsolution.x;
% objectiveCoeff = 1.0;
% model = changeObjective(model, 'EX_s_0189[c]', objectiveCoeff); 
% 
%removing the objective function to minimise the bias
ibm1 = find(ismember(model.rxns, biomass));
model.lb(ibm1)=0.1;
model.c(:)=0;
fbaWT = optimizeCbModel(model);
fluxesWT = fbaWT.x;
% % 
Graphtitle1 = 'Metabolite Interaction map of WT';
Graphtitle1 = [Graphtitle1 newline 'NO Objective'];
graphObjWT=createMetIntrcNetwork(model,mets,'fluxes',fluxesWT, 'Graphtitle',Graphtitle1);
set(gcf,'Visible','on'); % produce figure as pop up since live editor does

%%
% of_list={'YKR009C','YKL060C','YGR032W','YMR306W','YLR342W','YCR034W','YCR028C','YBR291C','YGR244C','YOR142W','YOR245C','YDL078C','YML042W','YMR246W','YOR317W','YMR241W','YGL125W','YPL023C','YPR021C','YDR178W','YJL045W','YKL141W','YLL041C','YKL148C','YLR164W','YOL059W','YGR015C','YMR205C','YGR240C','YML120C','YFR015C','YJL137C','YKR058W','YLR258W','YNL106C','YOR109W','YOR120W','YPL028W','YEL063C','YNL270C','YLR027C','YLL052C','YLL053C','YPR192W','YIL160C','YEL047C','YML059C','YBR145W','YOL086C'};
% ok_list={'YMR303C', 'YLR303W', 'YJL097W', 'YBR249C', 'YDR035W', 'YER065C', 'YFR019W', 'YPR021C', 'YGL256W', 'YMR083W', 'YKL146W', 'YNL101W', 'YIL006W', 'YHR043C', 'YHR044C', 'YDL052C', 'YKR089C', 'YOR081C', 'YDR400W', 'YLR209C', 'YNR013C', 'YBL015W', 'YAL044C', 'YDR019C', 'YFL018C', 'YMR189W', 'YDR148C', 'YIL125W', 'YDL078C', 'YAL054C', 'YFL054C', 'YLL043W', 'YHR068W', 'YIR031C', 'YNL117W', 'YML042W', 'YLR153C', 'YBR058C-A', 'YDR062W', 'YMR296C', 'YOR283W', 'YKL152C', 'YHR067W', 'YKL212W', 'YLR284C', 'YIL160C', 'YMR205C', 'YGR240C', 'YJL071W',  'YHR100C', 'YGR256W', 'YHR183W', 'YPR128C', 'YGL080W', 'YGR243W', 'YHR162W', 'YKL060C', 'YDL080C', 'YOL059W', 'YBR011C', 'YOL157C', 'YMR207C', 'YGR192C', 'YJL052W', 'YJR009C', 'YLR174W', 'YLR299W', 'YHL018W', 'YDR305C', 'YBR084W'};
% og_list={'YCR028C', 'YDL142C', 'YER183C', 'YGR202C', 'YLR157C', 'YDR403W', 'YPR140W', 'YJR130C', 'YJL005W', 'YDR040C', 'YDR001C', 'YJR001W', 'YMR083W', 'YDR284C', 'YNL037C', 'YNL065W', 'YOR126C', 'YOR241W', 'YGR209C', 'YBR208C', 'YOR155C', 'YPL092W', 'YML054C', 'YPL268W', 'YDR400W', 'YML022W', 'YKL215C', 'YMR226C', 'YPR069C', 'YAL060W', 'YNL202W', 'YER163C', 'YER019W', 'YPR127W', 'YHR002W', 'YMR303C', 'YIL006W', 'YPR062W', 'YPL028W', 'YPL206C', 'YDL120W', 'YIL155C', 'YKL141W', 'YMR008C', 'YKL148C', 'YDL166C', 'YBR069C', 'YLR151C', 'YJR010W', 'YNL104C', 'YGR019W', 'YER010C', 'YDR503C', 'YBR006W', 'YNR057C', 'YBR281C', 'YDL040C', 'YDL198C', 'YDR173C', 'YCL038C', 'YIL099W', 'YJL200C', 'YBR180W', 'YBR291C', 'YLL041C', 'YDR196C', 'YKL212W', 'YJL196C', 'YOL064C', 'YFR019W', 'YLR153C', 'YKL103C', 'YFL030W', 'YOR163W', 'YLR245C', 'YER086W', 'YFR044C', 'YBR011C', 'YMR207C', 'YHR123W', 'YML035C', 'YLR020C', 'YOL103W', 'YDR272W', 'YPL147W', 'YDR148C', 'YPL087W', 'YHR144C', 'YGL077C', 'YGR015C', 'YMR289W','YKL055C', 'YGR012W'};
% 
test_list={'YDR284C',...
    'YKL055C',...
    'YDL078C',...
    'YBL015W',...
    'YNL117W',...
    'YDR403W',...
    'YDR503C',...
    'YER019W',...
    'YBR180W'};


%try
 
for i = 1:length(test_list)
    currentGene=test_list{i};
    [del_model, ~, deleted_reactions, ~] = deleteModelGenes(model, currentGene);

    display(currentGene);
    display(deleted_reactions); %the following reaction list will be knocked-out if these genes are deleted

    % del_model = changeObjective(del_model, 'EX_s_0189[c]', objectiveCoeff);
    ibm2 = find(ismember(del_model.rxns, biomass));
    del_model.lb(ibm2)=0.1;
    del_model.c(:)=0;
    fbaKO = optimizeCbModel(del_model);
    fluxesKO = fbaKO.x;


    Graphtitle2 = sprintf('Metabolite Interaction Map of %s', currentGene);
    %Graphtitle2 = [Graphtitle2 newline 'deleted gene(s):YML120C'];

    graphObjKO=createMetIntrcNetwork(del_model,mets,'fluxes',fluxesKO, 'Graphtitle',Graphtitle2);
    set(gcf,'Visible','on'); % produce figure as pop up since live editor does
    


        
end

%catch err
    



%end
%%

% interaction map for all metabolites

% graphObj=createMetIntrcNetwork(model,model.mets);

%%

%test section 

upreg_sin_genes={'YLR017W',...
    'YGR088W',...
    'YGR277C'};

for i=1:length(upreg_sin_genes)
   
    %removing the objective function to minimise the bias
    
    ibm = find(ismember(model.rxns, biomass));
    model.lb(ibm)=0.1;
    model.c(:)=0;
    
    model1 = model;
    model2 = model;
    
    
    del_gene=upreg_sin_genes{i};
    [delModel, ~, deleted_reactions, ~] = deleteModelGenes(model, del_gene);
    
    fprintf('\n________________________________________________________');
    fprintf('________________________________________________________\n');
    
    display(del_gene);
    
    display(deleted_reactions); %the following reaction list will be knocked-out if these genes are deleted
    
    for j=1:length(deleted_reactions)
        del_rxn=deleted_reactions{j};

        %fprintf('\nNow Being Deleted: %s\n', del_rxn)

        model1 = changeRxnBounds(model1, del_rxn, 0.2, 'l');
        model1 = changeRxnBounds(model1, del_rxn, 1, 'u');



        model2 = changeRxnBounds(model2, del_rxn, 1, 'l');
        model2 = changeRxnBounds(model2, del_rxn, 5, 'u');           
            
    end  

    
    try
        
    
    fbaWT = optimizeCbModel(model1);
    fluxesWT = fbaWT.x;  
    
    
    fbaUP = optimizeCbModel(model2);
    fluxesUP = fbaUP.x;
    
    Graphtitle1 = sprintf('Metabolite Interaction Map of WT gene %s', del_gene);
    %Graphtitle1 = [Graphtitle1 newline 'deleted gene(s):YML120C'];

    graphObjKO=createMetIntrcNetwork(model1 ,mets,'fluxes',fluxesWT, 'Graphtitle',Graphtitle1);
    set(gcf,'Visible','on'); % produce figure as pop up since live editor does
    
    catch
        
    end
    
    
    try   
      
    Graphtitle2 = sprintf('Metabolite Interaction Map of %s', del_gene);
    
    graphObjUP=createMetIntrcNetwork(model2 ,mets,'fluxes',fluxesUP, 'Graphtitle',Graphtitle2);
    
    set(gcf,'Visible','on'); % produce figure as pop up since live editor does
    
    catch
        
    end

        
end
