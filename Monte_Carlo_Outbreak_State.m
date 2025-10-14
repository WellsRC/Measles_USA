function [Outbreak_State]=Monte_Carlo_Outbreak_State(p_inf_state,k_nb_state,p_nb_state,NS)

Outbreak_State=zeros(length(p_inf_state),NS);
for ss=1:length(p_inf_state)
    r_z=rand(NS,1);
    os=nbinrnd(k_nb_state(ss),p_nb_state(ss),NS,1);
    Outbreak_State(ss,r_z>p_inf_state(ss))=os(r_z>p_inf_state(ss));
end

end