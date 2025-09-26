function [beta_x,beta_insurance,hyp_par,hyp_par_gam_a] = Parameters(x)

% Intercept based on the type of insurance
beta_insurance.Uninsured=x(1);
beta_insurance.Private=x(2);
beta_insurance.Public=x(3);

% Coefficients for the logisitic regression model
beta_x=x(4:end-4);
beta_x=beta_x(:);

hyp_par=10.^x(end-3:end-2);
hyp_par_gam_a=10.^x(end-1:end);
end

