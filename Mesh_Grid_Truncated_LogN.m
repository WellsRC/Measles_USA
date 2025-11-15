clear;
clc;

avg_fs=linspace(0,8,500);
std_lgn=linspace(-5,3,750);
[avg_fs,std_lgn]=meshgrid(avg_fs,std_lgn);
old_size=size(avg_fs);
mu_lgn=zeros(size(avg_fs));
avg_fs=avg_fs(:);
std_lgn=std_lgn(:);
mu_lgn=mu_lgn(:);
p0=linspace(0,1-10^(-3),10^4);
opts=optimoptions('fmincon','FunctionTolerance',10^(-9),'MaxFunctionEvaluations',10^3,'StepTolerance',10^(-9));
parfor kk=1:length(mu_lgn)
    temp_std=10.^std_lgn(kk);
    afs=10.^avg_fs(kk);

    

    J=(temp_std.*(1-(p0))-afs.*(p0).*(1-(p0).^temp_std)).^2;
    x0=log10(p0(J==min(J)));
    mu_lgn(kk)=fmincon(@(z) (temp_std.*(1-(10.^z))-afs.*(10.^z).*(1-(10.^z).^temp_std)).^2,x0,[],[],[],[],-16,0,[],opts);
end
mu_lgn=10.^mu_lgn;
avg_fs=reshape(avg_fs,old_size);
std_lgn=reshape(std_lgn,old_size);
mu_lgn=reshape(mu_lgn,old_size);

save('Turncated_log_Normal_Parameter.mat',"mu_lgn","std_lgn","avg_fs");
