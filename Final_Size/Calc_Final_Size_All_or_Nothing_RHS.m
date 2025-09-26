function RHS = Calc_Final_Size_All_or_Nothing_RHS(vac_cov,R0,zt)

RHS=(1-vac_cov).*(1-exp(-R0.*zt));

end

