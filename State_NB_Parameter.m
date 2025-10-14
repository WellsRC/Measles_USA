function [p_inf_state,k_nb_state,p_nb_state,L_NaN] = State_NB_Parameter(state_matrix,p_inf_County,mu_County)

p_inf_state=zeros(size(state_matrix,1),1);
k_nb_state=zeros(size(state_matrix,1),1);
p_nb_state=zeros(size(state_matrix,1),1);

L_NaN=false;

p_temp=p_inf_County(:)+(1-p_inf_County(:)).*poisspdf(0,mu_County(:));
for State_Num=1:size(state_matrix,1)

    mu_County_NB=state_matrix(State_Num,:)*((1-p_inf_County(:)).*mu_County(:)); 
    var_County_NB=state_matrix(State_Num,:)*((1-p_inf_County(:)).*(mu_County(:)+p_inf_County(:).*mu_County(:).^2)); 
    p_zero_county=prod(p_temp(state_matrix(State_Num,:)==1));
    if(p_zero_county==0)
        p_zero_county=10^(-64);
    end
    
    [p_inf_state(State_Num),k_nb_state(State_Num),p_nb_state(State_Num),L_NaNt]=Estimate_ZIF_NB(mu_County_NB,var_County_NB,p_zero_county);
    L_NaN=L_NaN|L_NaNt;
end


end

