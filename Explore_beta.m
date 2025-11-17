clear;
clc;
close all;

Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')

Known_Ind_Cases=zeros(length(County_Data.County),1);

Measles_Cases=readtable('County_Level_Measles_Cases_Adjusted.csv');

Imported_Case=zeros(length(County_Data.County),1);

% Known imported cases

beta_Cap=zeros(length(County_Data.GEOID),1);
for cc=1:length(County_Data.GEOID)
    t_f=str2double(County_Data.GEOID{cc})==Measles_Cases.GEOID & strcmp(Measles_Cases.type,'imported') & ~isnan(Measles_Cases.case_count);
    if(sum(t_f)>0)
        Imported_Case(cc)=Measles_Cases.case_count(t_f);
    end

    t_f=str2double(County_Data.GEOID{cc})==Measles_Cases.GEOID & strcmp(Measles_Cases.type,'local') & ~isnan(Measles_Cases.case_count);
    if(sum(t_f)>0)
        Known_Ind_Cases(cc)=Measles_Cases.case_count(t_f);
    end
    beta_Cap(cc)=pchip(County_Data.R_0(cc,:),County_Data.beta_j,18.1);
end


f_out=find(Known_Ind_Cases>=3);
b_temp=zeros(length(f_out),1);
R0=zeros(length(b_temp),1);


for jj=1:length(b_temp)
    xx=County_Data.Final_Size_Est(f_out(jj),:);
    bb=County_Data.beta_j;
    bb=bb(xx>0);
    xx=xx(xx>0);
    b_temp(jj)=pchip(xx,bb,Known_Ind_Cases(f_out(jj)));
    R0(jj)=interp1(County_Data.beta_j',County_Data.R_0(f_out(jj),:)',b_temp(jj))';
end
M_beta=median(beta_Cap);
N_temp=1-County_Data.Total_Immunity(f_out);
opts=optimoptions('fmincon','FunctionTolerance',10^(-9),'MaxFunctionEvaluations',10^6,'MaxIterations',10^6);
[X,fval]=lsqnonlin(@(x) Objective_Function_Transmission(x,N_temp,b_temp,County_Data,Known_Ind_Cases,f_out,true),[log10(200)	-1.07901083234047	-0.665808554320118	0.0474997453437891 log10(54.9)],[2 log10(0) log10(0.15) 0 -4],[3 log10(0.15) log10(0.5) log10(1.5) 3]);
[X,fval]=fmincon(@(x) Objective_Function_Transmission(x,N_temp,b_temp,County_Data,Known_Ind_Cases,f_out,false),[X(1:4) log10(sqrt(fval./66))],[],[],[],[],[2 log10(0) log10(0.15) 0 -4],[3 log10(0.15) log10(0.5) log10(1.5) 3],[],opts);


Xt=10.^X;
Xt(1)=-Xt(1);
figure(1)
scatter((N_temp),b_temp);
hold on
ii=linspace(0,0.5,1001);

    plot(ii,Transmission_Relation(ii,Xt(1:4)),'k');

scatter((N_temp),b_temp,20,'r','filled');
xlabel('Susceptible')
ylabel('beta')




X_Full=X;
L_Full=-fval;

for jj=1:10
    lb=[2. -2.5 -0.75 0 -1.4];
    ub=[2.7 -0.85 -0.45 log10(1.5) -0.8];
    
    XS=lb+(ub-lb).*rand(10^5,5);
    LS=zeros(10^5,1);
    parfor ii=1:10^5
        LS(ii)=-Objective_Function_Transmission(XS(ii,:),N_temp,b_temp,County_Data,Known_Ind_Cases,f_out,false);
    end

    XS=XS(~isnan(LS),:);
    LS=LS(~isnan(LS));

    XS=XS(~isinf(LS),:);
    LS=LS(~isinf(LS));

    X_Full=[X_Full;XS];
    L_Full=[L_Full;LS];
end



w=cumsum(exp(L_Full-max(L_Full)))./sum(exp(L_Full-max(L_Full)));

r_samp=rand(5.*10^3,1);

X_Samp=zeros(5.*10^3,5);
L_Samp=zeros(5.*10^3,1);
for ii=1:5.*10^3
    indx=find(r_samp(ii)<=w,1,"first");
    X_Samp(ii,:)=X_Full(indx,:);
    L_Samp(ii)=L_Full(indx);
end

X=10.^X;
X(1)=-X(1);

X_Samp=10.^X_Samp;
X_Samp(:,1)=-X_Samp(:,1);


beta_j=Transmission_Relation(N_temp,X(1:4));

FS=zeros(length(beta_j),1);
for jj=1:length(beta_j)
    FS(jj)=interp1(County_Data.beta_j',County_Data.Final_Size_Est(f_out(jj),:)',beta_j(jj))';
end


figure(2)
scatter((N_temp),b_temp);
hold on
ii=linspace(0,0.5,1001);
for jj=1:5000
    plot(ii,Transmission_Relation(ii,X_Samp(jj,1:4)),'k');
end
scatter((N_temp),b_temp,20,'r','filled');
xlabel('Susceptible')
ylabel('beta')

K=[FS Known_Ind_Cases(f_out)];

save('Explore_Beta_Relation_Parameters.mat','X_Samp','L_Samp');