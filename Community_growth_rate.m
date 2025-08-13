
clc,clear;


models_name={'organism_Q_sbml','organism_P_sbml'};
%models_name={'MM','DV'};

%% Create community model


communitymodel=readCbModel([SAVEDIR filesep append(models_name{1,1},'.xml')]);

[m,n]=size(communitymodel.S);
communitymodel.rxns(:,1)=strcat(communitymodel.rxns(:,1),'_species1');
communitymodel.mets(:,1)=strcat(communitymodel.mets(1:m),'_species1');

for i=2:length(models_name)
  
  model=readCbModel([SAVEDIR filesep append(models_name{1,i},'.xml')]);
  % 
  % biomass_auto_id=[2392;2393;2394];
  % biomass_auto=['Biomass_Chlamy_auto' ,'Biomass_Chlamy_mixo','Biomass_Chlamy_hetero'];
  model.speciesTagName=strcat('_species',num2str(i));
  communitymodel=CreateCommunityModel(model,communitymodel);
end


%%
[m,n]=size(communitymodel.S);
% convert Exchange reaction to two irrevesible reaction

expression = '^EX_.*?(?<!_\[Env\])$'; %'^EX_\w*+(?<!_\[Env\])$'; %'EX_\w*(?<!_[Env])$';
Ex_index=find(~cellfun('isempty', regexp(communitymodel.rxns,expression)));

for i=1:length(Ex_index)
    if(communitymodel.lb(Ex_index(i)) < 0 & communitymodel.ub(Ex_index(i)) > 0)
        met_index=find(communitymodel.S(:,Ex_index(i))~=0);

        if(communitymodel.S(met_index,Ex_index(i))<0)
            name_rxn=extractAfter(communitymodel.rxns(Ex_index(i),1),'EX_');

            communitymodel.rxns(Ex_index(i),1)=strcat('EXCom_export_',name_rxn);
            communitymodel.rxns(n+1,1)=strcat('EXCom_uptake_',name_rxn);
            communitymodel.S(met_index,n+1)=communitymodel.S(met_index,Ex_index(i))*-1;
            communitymodel.lb(n+1,1)=0;
            communitymodel.ub(n+1,1)=communitymodel.lb(Ex_index(i))*-1;

            communitymodel.lb(Ex_index(i),1)=0;
            %communitymodel.ub(Ex_index(i),1)=1000;
            
            n=n+1;

        elseif(communitymodel.S(met_index,Ex_index(i))>0)
            communitymodel.rxns(Ex_index(i),1)=strcat('EXCom_uptake_',name_rxn);

            communitymodel.rxns(n+1,1)=strcat('EXCom_export_',name_rxn);
            communitymodel.S(met_index,n+1)=communitymodel.S(met_index,Ex_index(i))*-1;
            communitymodel.lb(n+1,1)=0;
            communitymodel.ub(n+1,1)=communitymodel.lb(Ex_index(i),1)*-1;
            communitymodel.lb(Ex_index(i),1)=0;

            n=n+1;

        end
    else
        met_index=find(communitymodel.S(:,Ex_index(i))~=0);

        if(communitymodel.S(met_index,Ex_index(i))<0)
        communitymodel.rxns(Ex_index(i),1)=strcat('EXCom_export_',extractAfter(communitymodel.rxns(Ex_index(i),1),'EX_'));
        else
        communitymodel.rxns(Ex_index(i),1)=strcat('EXCom_uptake_',extractAfter(communitymodel.rxns(Ex_index(i),1),'EX_'));
        end
    end  
end % for

%%
[m, n] = size(communitymodel.S);
S = communitymodel.S(1:m, 1:n);

% Set very small values to zero (adjust the threshold if needed)
threshold = 1e-6;
S(abs(S) < threshold) = 0;

% Convert to cell array with formatted numbers
S_cell = cell(m+1, n+1);
S_cell(1, 2:n+1) = communitymodel.rxns(1:n,1);
S_cell(2:m+1, 1) = communitymodel.mets(1:m,1);

for i = 2:m+1
    for j = 2:n+1
        if S(i-1,j-1) == 0
            S_cell{i,j} = 0;  % Write zero
        else
            S_cell{i,j} = S(i-1,j-1);  % Keep original value
        end
    end
end
xlswrite('model\S.xlsx',S_cell,1);

header={'rxns'};
output(1,1)=header;
output(2:n+1,1)=communitymodel.rxns;
xlswrite('model\S.xlsx',output,2);

clear output;
header={'mets'};
output(1,1)=header;
output(2:m+1,1)=communitymodel.mets;
xlswrite('model\S.xlsx',output,3);

clear output;
header={'lb'};
output(1,1)=header;
output(2:n+1,1)=num2cell(communitymodel.lb);
xlswrite('model\S.xlsx',output,4);

clear output;
header={'ub'};
output(1,1)=header;
output(2:n+1,1)=num2cell(communitymodel.ub) ;
xlswrite('model\S.xlsx',output,5);



