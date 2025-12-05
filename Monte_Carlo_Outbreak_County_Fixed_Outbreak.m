function [Total_Cases_County,Unvaccinated_Cases_County,Vaccinated_Cases_County,Uninsured_Unvaccinated_Cases_County,Uninsured_Vaccinated_Cases_County,Public_Unvaccinated_Cases_County,Public_Vaccinated_Cases_County]=Monte_Carlo_Outbreak_County_Fixed_Outbreak(Outbreak_Size,Proportion_Size_Age_Unvaccinated,Proportion_Size_Age_Vaccinated,Proportion_Age_Unvaccinated_Uninsured,Proportion_Age_Unvaccinated_Public,Proportion_Age_Unvaccinated_Private,Proportion_Age_Vaccinated_Uninsured,Proportion_Age_Vaccinated_Public,Proportion_Age_Vaccinated_Private,NS)
rng(20251009)
Total_Cases_County=Outbreak_Size.*ones(size(Proportion_Size_Age_Unvaccinated,1),NS);
Unvaccinated_Cases_County=zeros(size(Proportion_Size_Age_Unvaccinated,1),size(Proportion_Size_Age_Unvaccinated,2),NS);
Vaccinated_Cases_County=zeros(size(Proportion_Size_Age_Unvaccinated,1),size(Proportion_Size_Age_Vaccinated,2),NS);

Uninsured_Unvaccinated_Cases_County=zeros(size(Proportion_Size_Age_Unvaccinated,1),size(Proportion_Size_Age_Unvaccinated,2),NS);
Uninsured_Vaccinated_Cases_County=zeros(size(Proportion_Size_Age_Unvaccinated,1),size(Proportion_Size_Age_Vaccinated,2),NS);

Public_Unvaccinated_Cases_County=zeros(size(Proportion_Size_Age_Unvaccinated,1),size(Proportion_Size_Age_Unvaccinated,2),NS);
Public_Vaccinated_Cases_County=zeros(size(Proportion_Size_Age_Unvaccinated,1),size(Proportion_Size_Age_Vaccinated,2),NS);

for ss=1:size(Proportion_Size_Age_Unvaccinated,1)
    p_vec=[Proportion_Size_Age_Unvaccinated(ss,:) Proportion_Size_Age_Vaccinated(ss,:) ];
    pd = makedist('Multinomial','Probabilities',p_vec);
    for nn=1:NS
        if(Total_Cases_County(ss,nn)>0)
            r = random(pd,1,Total_Cases_County(ss,nn));
            for jj=1:size(Proportion_Size_Age_Unvaccinated,2)
                Unvaccinated_Cases_County(ss,jj,nn)=sum(r==jj);
                w_ins=[Proportion_Age_Unvaccinated_Uninsured(ss,jj) Proportion_Age_Unvaccinated_Public(ss,jj) Proportion_Age_Unvaccinated_Private(ss,jj)];
                w_ins=w_ins./sum(w_ins);
                for zz=1:sum(r==jj)
                    ins_typ=rand(1);
                    if(ins_typ<=w_ins(1))
                        Uninsured_Unvaccinated_Cases_County(ss,jj,nn)=Uninsured_Unvaccinated_Cases_County(ss,jj,nn)+1;
                    elseif(ins_typ<=w_ins(2))
                        Public_Unvaccinated_Cases_County(ss,jj,nn)=Public_Unvaccinated_Cases_County(ss,jj,nn)+1;
                    end
                end

                Vaccinated_Cases_County(ss,jj,nn)=sum(r==(jj+size(Proportion_Size_Age_Unvaccinated,2)));
                w_ins=[Proportion_Age_Vaccinated_Uninsured(ss,jj) Proportion_Age_Vaccinated_Public(ss,jj) Proportion_Age_Vaccinated_Private(ss,jj)];
                w_ins=w_ins./sum(w_ins);
                for zz=1:sum(r==jj)
                    ins_typ=rand(1);
                    if(ins_typ<=w_ins(1))
                        Uninsured_Vaccinated_Cases_County(ss,jj,nn)=Uninsured_Vaccinated_Cases_County(ss,jj,nn)+1;
                    elseif(ins_typ<=w_ins(2))
                        Public_Vaccinated_Cases_County(ss,jj,nn)=Public_Vaccinated_Cases_County(ss,jj,nn)+1;
                    end
                end
            end
        end
    end
end
end