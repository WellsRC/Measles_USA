function [beta_county] = County_Transmission(beta_transmission,X)

temp_beta=X*beta_transmission(:);

beta_county=exp(temp_beta);

end