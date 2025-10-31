function [Productivity_Cost_Under_15_Case,Productivity_Cost_15_plus_Case,Productivity_Cost_Under_15_Contact,Productivity_Cost_15_plus_Contact,Cost_per_Contact,Cost_per_Vaccine_dose_Private,Cost_per_Vaccine_dose_VFC,Cost_per_Hospitalization,Cost_per_Non_Hospitalization,Tests_per_Contact,Cost_per_Test]=Measles_Outbreak_Cost()
Cost_per_Contact=208.67;


% 44/318 proportion of contacts to receive MMR PEP after exposure
Cost_per_Vaccine_dose_Private=2.*(95.20+25.80).*(44/318); 
Cost_per_Vaccine_dose_VFC=2.*(26.33+13.25).*(44/318);  

Cost_per_Hospitalization=9522.11;
Cost_per_Non_Hospitalization=408.08;

Tests_per_Contact=684/4011;
Cost_per_Test=(191.92+29.45);

% https://www.tandfonline.com/doi/full/10.1080/13696998.2018.1542520#abstract
% https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0105153
Productivity_Cost_Under_15_Case=round(1.339403059*57324/366*7.3,2); 
Productivity_Cost_15_plus_Case=round(1.339403059*57324/366*10.1,2);


% Quarantine 21 days
% https://www.nj.gov/health/cd/documents/topics/measles/measles_exposure_guidance_public.pdf
Productivity_Cost_Under_15_Contact=round(1.339403059*57324/366*21,2)./3; 
Productivity_Cost_15_plus_Contact=round(1.339403059*57324/366*21,2);
end