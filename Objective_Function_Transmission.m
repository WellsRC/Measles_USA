function J=Objective_Function_Transmission(X,N_temp,b_temp,County_Data,Known_Ind_Cases,f_out,lsq)

w=Known_Ind_Cases(f_out);
w=w./sum(w); % Allowing as much uncertainty so it can be teased out in the fitting of the data
X=10.^X;
X(1)=-X(1);
err=X(5);
beta_j=Transmission_Relation(N_temp,X(1:4));

FS=zeros(length(beta_j),1);
for jj=1:length(beta_j)
    FS(jj)=interp1(County_Data.beta_j',County_Data.Final_Size_Est(f_out(jj),:)',beta_j(jj))';
end

if(lsq)
     J=w(:).*(b_temp(:)-FS(:));
else
    J=-sum(w(:).*log(normpdf(b_temp(:),beta_j(:),err)));
end

end

