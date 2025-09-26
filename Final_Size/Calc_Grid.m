clear;

vac_cov=0:0.001:1;
R0=1.01:0.01:20;

[vac_cov,R0]=ndgrid(vac_cov,R0);

Z=0.*vac_cov;


parfor ii=1:length(Z(:))
    if((1-vac_cov(ii)).*R0(ii)>1)
        Z(ii) = 10.^fmincon(@(x) 10.^16.*(Calc_Final_Size_All_or_Nothing_RHS(vac_cov(ii),R0(ii),10.^x)-10.^x).^2,log10(1-vac_cov(ii))-0.001,[],[],[],[],-32,log10(1-vac_cov(ii)));
    end
end


save('Grid_Lookup_Final_Size_All_or_Nothing.mat','vac_cov','R0','Z');