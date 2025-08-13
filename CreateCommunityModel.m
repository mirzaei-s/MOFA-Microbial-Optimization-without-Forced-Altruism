function [communitymodel] = CreateCommunityModel(model,communitymodel)
%%
[m1,n1]=size(communitymodel.S);
[m2,n2]=size(model.S);


EX_1=find(contains(model.rxns,'EX_'));
EX_2=find(contains(communitymodel.rxns,'EX_'));

[i1,~]=find(model.S(:,EX_1)~=0);
[i2,~]=find(communitymodel.S(:,EX_2)~=0);

mets_1=model.mets(i1);
mets_2=extractBefore(communitymodel.mets(i2),'_species');
%
% mets_1=extractBefore(model.mets(i1),'_sp2');
% mets_2=extractBefore(communitymodel.mets(i2),'_sp1');


%erase(temp_mets_2,'_[Env]');
shared_mets=intersect(mets_1,mets_2);
clear mets_1 mets_2;
%%
for i=1:length(shared_mets)
    met1_index=find(strcmp(model.mets(:,1),shared_mets(i)));
    temp=find(model.S(met1_index,:)~=0);  %find related reaction of shared metabolite
    rxn1_index=temp(find(contains(model.rxns(temp),'EX_')));
    %model.rxns(rxn1_index)

    if(model.lb(rxn1_index) < 0 & model.ub(rxn1_index) > 0) % if exchange reaction is reversible

        if(model.S(met1_index,rxn1_index)<0) % reaction is export in S matrix

            model.S(:,n2+1)=zeros(m2,1);
            model.S(met1_index,n2+1)=-1*model.S(met1_index,rxn1_index);
            model.rxns(n2+1)=strcat( 'EX_uptake_',extractAfter(model.rxns(rxn1_index,1),'EX_'));
            model.lb(n2+1)=0;
            model.ub(n2+1)=model.lb(rxn1_index)*-1;
            model.c(n2+1)=0;

            model.rxns(rxn1_index,1)=strcat( 'EX_export_',extractAfter(model.rxns(rxn1_index,1),'EX_'));
            model.lb(rxn1_index)=0;

        else
            model.S(:,n2+1)=zeros(m2,1);
            model.S(met1_index,n2+1)=-1*model.S(met1_index,rxn1_index);
            model.rxns(n2+1)=strcat( 'EX_export_',extractAfter(model.rxns(rxn1_index,1),'EX_'));
            model.lb(n2+1)=0;
            model.ub(n2+1)=model.lb(rxn1_index)*-1;
            model.c(n2+1)=0;
            model.rxns(rxn1_index,1)=strcat( 'EX_uptake_',extractAfter(model.rxns(rxn1_index,1),'EX_'));
            model.lb(rxn1_index)=0;
        end
        n2=n2+1;

    elseif(model.S(met1_index,rxn1_index)<0) %% model.lb < 0 & model.ub < 0
        model.rxns(rxn1_index,1)=strcat( 'EX_export_',extractAfter(model.rxns(rxn1_index,1),'EX_'));

    else
        model.rxns(rxn1_index,1)=strcat( 'EX_uptake_',extractAfter(model.rxns(rxn1_index,1),'EX_'));


    end %% model.lb < 0 & model.ub > 0
    clear temp;
    %--------------------------------------------------------------------------------------------------
    %                                             second model
    %---------------------------------------------------------------------------------------------------

    met2_index=find(strcmp(extractBefore(communitymodel.mets(:,1),'_species'),shared_mets(i)));
    temp=find(communitymodel.S(met2_index,:)~=0);  %find related reaction of shared metabolite
    rxn2_index=temp(find(contains(communitymodel.rxns(temp),'EX_')));
    %communitymodel.rxns(rxn2_index)

    if(~any(contains(communitymodel.rxns(rxn2_index),'_[Env]')))

        if(communitymodel.lb(rxn2_index) < 0 & communitymodel.ub(rxn2_index) > 0) % if exchange reaction is reversible
            if(communitymodel.S(met2_index,rxn2_index)<0) % reaction is export in S matrix

                %communitymodel.S(:,n1+1)=zeros(m1,1);
                communitymodel.S(met2_index,n1+1)=-1*communitymodel.S(met2_index,rxn2_index);
                communitymodel.rxns(n1+1)=strcat( 'EX_uptake_',extractAfter(communitymodel.rxns(rxn2_index,1),'EX_'));
                communitymodel.lb(n1+1,1)=0;
                communitymodel.ub(n1+1,1)=communitymodel.lb(rxn2_index)*-1;
                communitymodel.c(n1+1,1)=0;
                communitymodel.rxns(rxn2_index,1)=strcat( 'EX_export_',extractAfter(communitymodel.rxns(rxn2_index,1),'EX_'));
                communitymodel.lb(rxn2_index)=0;


            else
                communitymodel.S(:,n1+1)=zeros(m1,1);
                communitymodel.S(met2_index,n1+1)=-1*communitymodel.S(met2_index,rxn2_index);
                communitymodel.rxns(n1+1)=strcat('EX_export_',extractAfter(communitymodel.rxns(rxn2_index,1),'EX_'));
                communitymodel.lb(n1+1)=0;
                communitymodel.ub(n1+1)=communitymodel.lb(rxn2_index)*-1;
                communitymodel.c(n1+1)=0;

                communitymodel.rxns(rxn2_index,1)=strcat( 'EX_uptake_',extractAfter(communitymodel.rxns(rxn2_index,1),'EX_'));
                communitymodel.lb(rxn2_index)=0;

            end
            n1=n1+1;

        elseif(communitymodel.S(met2_index,rxn2_index)>0) %% model.lb < 0 & model.ub < 0
            communitymodel.rxns(rxn2_index,1)=strcat( 'EX_uptake_',extractAfter(communitymodel.rxns(rxn2_index,1),'EX_'));
        elseif(communitymodel.S(met2_index,rxn2_index)<0)
            communitymodel.rxns(rxn2_index,1)=strcat( 'EX_export_',extractAfter(communitymodel.rxns(rxn2_index,1),'EX_'));

        end % communitymodel.lb(rxn2_index) < 0 & communitymodel.ub(rxn2_index) > 0
    end   % any(contains(communitymodel.rxns(rxn2_index),'_[Env]')

    clear met1_inde rxn1_index met2_index rxn2_index
end %for

%% Merge models

communitymodel.S(m1+1:m2+m1,n1+1:n1+n2)=model.S(1:m2,1:n2);
communitymodel.rxns(n1+1:n1+n2,1)=strcat(model.rxns(1:n2),model.speciesTagName);
communitymodel.mets(m1+1:m1+m2,1)= strcat(model.mets(1:m2),model.speciesTagName); %model.mets(1:m2,1) ; % % %
communitymodel.lb(n1+1:n1+n2,1)=model.lb(1:n2,1);
communitymodel.ub(n1+1:n1+n2,1)=model.ub(1:n2);
communitymodel.c(n1+1:n1+n2,1)=model.c(1:n2,1);

n=n1+n2;
m=m1+m2;

for i=1:length(shared_mets)
   
    met_index=find(strcmp(extractBefore(communitymodel.mets(1:m1+m2),'_species'),shared_mets(i)));
    [~,temp]=find(communitymodel.S(met_index(:,1),:)~=0);
    rxn_index=temp(find(contains(communitymodel.rxns(temp),'EX_')));

    if(~any(contains(communitymodel.rxns(rxn_index),'_[Env]')))

        communitymodel.mets(m+1,1)=shared_mets(i);

        communitymodel.rxns(n+1,1)=strcat('EX_export_',shared_mets(i),'_[Env]');
        communitymodel.S(m+1,n+1)=-1;
        communitymodel.lb(n+1,1)=0;
        communitymodel.ub(n+1,1)=1000;
        n=n+1;

        communitymodel.rxns(n+1,1)=strcat('EX_uptake_',shared_mets(i),'_[Env]');
        communitymodel.S(m+1,n+1)=1;
        communitymodel.lb(n+1,1)=0;
        communitymodel.ub(n+1,1)=1000;
        m=m+1;
        n=n+1;
        
        communitymodel.rxns(rxn_index,1)=strcat('EXCom_',extractAfter(communitymodel.rxns(rxn_index,1),'EX_'));

        for j=1:length(rxn_index)
            jj=find(communitymodel.S(:,rxn_index(j,1))~=0);
            communitymodel.S(m,rxn_index(j,1))=communitymodel.S(jj,rxn_index(j,1))*-1;

            if(~contains(communitymodel.mets(jj,1),'_species'))
                communitymodel.mets(jj,1)=strcat(communitymodel.mets(jj,1),'_species',extractAfter(communitymodel.rxns(rxn_index(j,1)),'_species'))
            end
        end

    else

        rxn=find(contains(communitymodel.rxns(rxn_index),'species'));
        communitymodel.rxns(rxn_index,1)=strcat('EXCom_',extractAfter(communitymodel.rxns(rxn,1),'EX_'));
        ii=find(communitymodel.S(:,rxn(1,1))~=0)

        temp=find(contains(communitymodel.rxns(rxn_index),'_[Env]'));
        jj=find(communitymodel.S(:,temp(1,1))~=0);

        for j=1:length(rxn)
            communitymodel.S(jj(1,1),rxn(j,1))=communitymodel.S(ii,rxn(j,1))*-1;
        end

    end % any(contains( ','_[Env]' )

end


end