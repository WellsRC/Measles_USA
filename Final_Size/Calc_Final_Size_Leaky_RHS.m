function RHS = Calc_Final_Size_Leaky_RHS(vac_up,eps_v,R0,zt)

RHS=(1-vac_up).*(1-exp(-R0.*zt))+vac_up.*(1-exp(-R0.*zt.*(1-eps_v)));

end

