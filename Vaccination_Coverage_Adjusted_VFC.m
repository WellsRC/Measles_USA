function [Estimated_Vaccination_Coverage] = Vaccination_Coverage_Adjusted_VFC(x,beta_x,beta_insurance,q_insurance,dZ_County)

% Uninsured population
z=beta_insurance.Uninsured+x*beta_x+dZ_County;
Estimated_Vaccination_Coverage.Uninsured=1./(1+exp(-z));

% Privately insured population
z=beta_insurance.Private+x*beta_x;
Estimated_Vaccination_Coverage.Private=1./(1+exp(-z));

% Public insured population
z=beta_insurance.Public+x*beta_x+dZ_County;
Estimated_Vaccination_Coverage.Public=1./(1+exp(-z));

% Overall

Estimated_Vaccination_Coverage.Overall=(q_insurance.Uninsured.*Estimated_Vaccination_Coverage.Uninsured+q_insurance.Private.*Estimated_Vaccination_Coverage.Private+q_insurance.Public.*Estimated_Vaccination_Coverage.Public)./(q_insurance.Uninsured+q_insurance.Private+q_insurance.Public);


Estimated_Vaccination_Coverage.Overall(isnan(q_insurance.Uninsured+q_insurance.Private+q_insurance.Public))=(Estimated_Vaccination_Coverage.Uninsured(isnan(q_insurance.Uninsured+q_insurance.Private+q_insurance.Public))+Estimated_Vaccination_Coverage.Private(isnan(q_insurance.Uninsured+q_insurance.Private+q_insurance.Public))+Estimated_Vaccination_Coverage.Public(isnan(q_insurance.Uninsured+q_insurance.Private+q_insurance.Public)))./3;
end

