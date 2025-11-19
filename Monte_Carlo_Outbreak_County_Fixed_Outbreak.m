function [Total_Cases_County,Unvaccinated_Cases_County,Vaccinated_Cases_County]=Monte_Carlo_Outbreak_County_Fixed_Outbreak(Outbreak_Size,Proportion_Size_Age_Unvaccinated,Proportion_Size_Age_Vaccinated,NS)
rng(20251009)
Total_Cases_County=Outbreak_Size.*ones(size(Proportion_Size_Age_Unvaccinated,1),NS);
Unvaccinated_Cases_County=zeros(size(Proportion_Size_Age_Unvaccinated,1),size(Proportion_Size_Age_Unvaccinated,2),NS);
Vaccinated_Cases_County=zeros(size(Proportion_Size_Age_Unvaccinated,1),size(Proportion_Size_Age_Vaccinated,2),NS);
for ss=1:size(Proportion_Size_Age_Unvaccinated,1)
    p_vec=[Proportion_Size_Age_Unvaccinated(ss,:) Proportion_Size_Age_Vaccinated(ss,:) ];
    pd = makedist('Multinomial','Probabilities',p_vec);
    for nn=1:NS
        r = random(pd,1,Total_Cases_County(ss,nn));
        for jj=1:size(Proportion_Size_Age_Unvaccinated,2)
            Unvaccinated_Cases_County(ss,jj,nn)=sum(r==jj);
            Vaccinated_Cases_County(ss,jj,nn)=sum(r==(jj+size(Proportion_Size_Age_Unvaccinated,2)));
        end
    end
end
end