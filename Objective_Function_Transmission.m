function J=Objective_Function_Transmission(X,N_temp,Known_Ind_Cases,County_Data,f_out,lsq)


X=10.^X;
X(1)=-X(1);
err=X(5);
beta_j=Transmission_Relation(N_temp,X(1:4));

FS=zeros(length(beta_j),1);
for jj=1:length(beta_j)
    FS(jj)=interp1(County_Data.beta_j',County_Data.Final_Size_Est(f_out(jj),:)',beta_j(jj))';
end

if(lsq)
     J=Known_Ind_Cases(:)-FS(:);
else
    J=-sum(log(normpdf(Known_Ind_Cases(:),FS(:),err)));
end

end

