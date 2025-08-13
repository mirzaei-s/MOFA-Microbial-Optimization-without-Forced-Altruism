clc,clear;
global CBT_LP_SOLVER
if isempty(CBT_LP_SOLVER)
    initCobraToolbox(false)
end



models_name={'organism_Q_sbml','organism_P_sbml'};
models_name={'MM','DV'};

for i=1:length(models_name)
    model=readCbModel([SAVEDIR filesep append(models_name{1,i},'.xml')]);
    EXs=find(contains(model.rxns,'EX_'));


    for j=1:length(EXs)

        [m,n]=size(model.S);
        met_index=find(model.S(:,EXs(j))~=0);
        model.S(m+1,EXs(j))=model.S(met_index,EXs(j))*-1;

        model.S(m+1,n+1)=1;
        model.rxns(n+1,1)=strcat(model.mets(met_index,1),'_com_up');
        model.mets(m+1,1)=strcat(model.mets(met_index,1),'_[Env]');
        model.lb(n+1,1)=0;
        model.ub(n+1,1)=1000;
        model.c(n+1,1)=0;
        model.b(m+1,1)=0
        model.csense(m+1,1)='E'
        model.rules{n+1,1}=''

        model.S(m+1,n+2)=-1;
        model.rxns(n+2,1)=strcat(model.mets(met_index,1),'_com_ex');
        model.lb(n+2,1)=0;
        model.ub(n+2,1)=1000;
        model.c(n+2,1)=0;
        model.rules{n+2,1}=''
        model.genes=model.rxns;
        model.geneNames=model.rxns
        model.metCharges(m+1,1)=0
        model.metNames(m+1,1)=strcat(model.mets(met_index,1),'_[Env]');
        model.proteins=model.rxns
        model.rxnNames=model.rxns
        model.rxnNotes=model.rxns


    end
    %writeCbModel(model,strcat('D:\PhD_thesis\Implementation\CommunityGrowthRate\models\inputModels','\necom_',models_name{i},'.xml'))

    [m,n]=size(model.S);
    S=cell(m+1,n+1);

    S(1, 2:n+1) = model.rxns(1:n,1);
    S(2:m+1, 1) = model.mets(1:m,1);
    S(2:m+1,2:n+1)=num2cell(model.S(1:m,1:n));
    xlswrite(strcat('C:\model\Necom\MM_DV\S_',models_name{i}),S,1);

    header={'rxns'};
    output(1,1)=header;
    output(2:n+1,1)=model.rxns;
    xlswrite(strcat('C:\model\Necom\MM_DV\S_',models_name{i}),output,2);

    clear output;
    header={'mets'};
    output(1,1)=header;
    output(2:m+1,1)=model.mets;
    xlswrite(strcat('C:\model\Necom\MM_DV\S_',models_name{i}),output,3);

    clear output;
    header={'lb'};
    output(1,1)=header;
    output(2:n+1,1)=num2cell(model.lb);
    xlswrite(strcat('C:\model\Necom\MM_DV\S_',models_name{i}),output,4);

    clear output;
    header={'ub'};
    output(1,1)=header;
    output(2:n+1,1)=num2cell(model.ub) ;
    xlswrite(strcat('C:\model\Necom\MM_DV\S_',models_name{i}),output,5);
end
%%
%models_name={'necom_organism_Q_sbml','necom_organism_P_sbml'};

model1=readCbModel([SAVEDIR filesep append(models_name{1,1},'.xml')]);
model2=readCbModel([SAVEDIR filesep append(models_name{1,2},'.xml')]);

[m1,n1]=size(model1.S);
[m2,n2]=size(model2.S);


EX_1=find(contains(model1.rxns,'EX_'));
EX_2=find(contains(model2.rxns,'EX_'));

[index,~]=find(model1.S(:,EX_1)~=0);
k=1;
for i=1:length(index)
    if(~contains(model1.mets(index(i),1),'_[Env]'))
        i1(k)=index(i);
        k=k+1;
    end
end

clear index

[index,~]=find(model2.S(:,EX_2)~=0);
k=1;
for i=1:length(index)
    if(~contains(model2.mets(index(i),1),'_[Env]'))
        i2(k)=index(i);
        k=k+1;
    end
end

mets_1=model1.mets(i1);
mets_2=model2.mets(i2);

shared_mets=intersect(mets_1,mets_2);