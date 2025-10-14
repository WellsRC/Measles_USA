clear;
clc;


Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')

R_0=zeros(length(County_Data.State),length(County_Data.beta_j));

for ii=1:length(County_Data.State)
    State_temp=County_Data.State{ii};
    State_temp(State_temp==' ')='_';
    M=readtable([pwd '/State_Contact_Matrix/United_States_subnational_' State_temp '_M_overall_contact_matrix_18']);
    M=table2array(M);
    
    Tot_Pop=County_Data.Total_Population(ii);
    Pop=table2array(County_Data.Population(ii,:));
    Pop=Tot_Pop.*Pop;

    % Equation 13
    % (https://pmc.ncbi.nlm.nih.gov/articles/PMC7088810/pdf/11538_2010_Article_9623.pdf)
    
    for bb=1:length(County_Data.beta_j)
        A_0=County_Data.beta_j(bb).*repmat(Pop,18,1).*M./repmat(Pop,18,1);
        A_0(repmat(Pop,18,1)==0)=0;
        R_0(ii,bb)=max(eig(A_0));
    end
end

County_Data.R_0=R_0;

save([Vaccine '_Immunity.mat'],'County_Data')
