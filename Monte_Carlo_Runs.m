function Monte_Carlo_Runs(National_Reduction,NS,Scenario_Plot,Age_0_to_6)

if(~Age_0_to_6)
    Age_Reduction=[true(1,1) false(1,4)];
else
    Age_Reduction=[true(1,2) false(1,3)];
end

load('Turncated_Negative_Binomial_Parameter.mat');
F_NB = scatteredInterpolant(kv(:),avg_fs(:),log(pv(:)./(1-pv(:))));
[Total_Cases_County,Unvaccinated_Cases_County_Baseline,Vaccinated_Cases_County_Baseline,Total_Contacts_Baseline,Unvaccinated_Contacts_Baseline,Imported_Case]=Monte_Carlo_Incidence(F_NB,National_Reduction,Age_Reduction,NS,Scenario_Plot,Age_0_to_6);

if(Age_0_to_6)
    save(['Monte_Carlo_Run_' Scenario_Plot '_National_Reduction=' num2str(National_Reduction.*100) '_Ages_0_to_6.mat'],'Total_Cases_County','Unvaccinated_Cases_County_Baseline','Vaccinated_Cases_County_Baseline','Total_Contacts_Baseline','Unvaccinated_Contacts_Baseline','Imported_Case');
else
    save(['Monte_Carlo_Run_' Scenario_Plot '_National_Reduction=' num2str(National_Reduction.*100) '_Ages_0_to_4.mat'],'Total_Cases_County','Unvaccinated_Cases_County_Baseline','Vaccinated_Cases_County_Baseline','Total_Contacts_Baseline','Unvaccinated_Contacts_Baseline','Imported_Case');
end

end