N_0_4=(112+248+178+486);
H_0_4=(112+178);

H_5_19=71;
N_5_19=(71+540);

H_20p=112;
N_20p=(112+170);

NV_0_4=(20+71+7);
HV_0_4=20;

NV_5_19=(9+72);
HV_5_19=9;

NV_20p=(38+161);
HV_20p=38;

pH=surrogateopt(@(x)-(log(binopdf(HV_0_4,NV_0_4,x(1).*(1-x(4))))+log(binopdf(H_0_4,N_0_4,x(1)))+log(binopdf(HV_5_19,NV_5_19,x(2).*(1-x(4))))+log(binopdf(H_5_19,N_5_19,x(2)))+log(binopdf(HV_20p,NV_20p,x(3).*(1-x(4))))+log(binopdf(H_20p,N_20p,x(3)))),[0 0 0 0],[1 1 1 1]);

pH=fmincon(@(x)-(log(binopdf(HV_0_4,NV_0_4,x(1).*(1-x(4))))+log(binopdf(H_0_4,N_0_4,x(1)))+log(binopdf(HV_5_19,NV_5_19,x(2).*(1-x(4))))+log(binopdf(H_5_19,N_5_19,x(2)))+log(binopdf(HV_20p,NV_20p,x(3).*(1-x(4))))+log(binopdf(H_20p,N_20p,x(3)))),pH,[],[],[],[],[0 0 0 0],[1 1 1 1]);

