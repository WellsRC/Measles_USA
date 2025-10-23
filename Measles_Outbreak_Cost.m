function [Cost_per_Contact,Cost_per_Vaccine_dose_Private,Cost_per_Vaccine_dose_VFC,Cost_per_Hospitalization,Cost_per_Non_Hospitalization]=Measles_Outbreak_Cost()
Cost_per_Contact=223;
Cost_per_Vaccine_dose_Private=2.*(95.20+25.80).*0.38; % 38% contacts received PEP
Cost_per_Vaccine_dose_VFC=2.*(26.33+13.25).*0.38;  % 38% contacts received PEP
Cost_per_Hospitalization=7438;
Cost_per_Non_Hospitalization=325;
end