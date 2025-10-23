clear;
clc;

avg_fs=linspace(0,6,500);
kv=linspace(-3.5,1,750);
[avg_fs,kv]=meshgrid(avg_fs,kv);
old_size=size(avg_fs);
pv=zeros(size(avg_fs));
avg_fs=avg_fs(:);
kv=kv(:);
pv=pv(:);
p0=linspace(0,1-10^(-3),10^4);
opts=optimoptions('fmincon','FunctionTolerance',10^(-9),'MaxFunctionEvaluations',10^3,'StepTolerance',10^(-9));
parfor kk=1:length(pv)
    k=10.^kv(kk);
    afs=10.^avg_fs(kk);
    J=(k.*(1-(p0))-afs.*(p0).*(1-(p0).^k)).^2;
    x0=log10(p0(J==min(J)));
    pv(kk)=fmincon(@(z) (k.*(1-(10.^z))-afs.*(10.^z).*(1-(10.^z).^k)).^2,x0,[],[],[],[],-16,0,[],opts);
end
pv=10.^pv;
avg_fs=reshape(avg_fs,old_size);
kv=reshape(kv,old_size);
pv=reshape(pv,old_size);

save('Turncated_Negative_Binomial_Parameter.mat',"pv","kv","avg_fs");
