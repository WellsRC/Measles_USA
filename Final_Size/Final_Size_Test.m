% clear;


load('Grid_Lookup_Final_Size_All_or_Nothing.mat')
Final_Size = griddedInterpolant(vac_cov,R0,Z);

v=linspace(0,1,10001);
R0=1.1;

plot(v,Final_Size(v',1.1.*ones(size(v'))));
