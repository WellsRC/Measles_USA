clear;
Age_Reduction=[true(1,2) false(1,3)];
Age_0_to_6=true;
for ii=0:0.005:0.05
    National_Reduction_Vaccine_Uptake(ii,Age_Reduction,Age_0_to_6);
    Monte_Carlo_Runs(ii,2500,'Sample_2025',Age_0_to_6);
end

