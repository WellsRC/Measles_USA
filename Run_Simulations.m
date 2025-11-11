clear;

for ii=0:0.01:0.05
    if(ii==0)
        National_Reduction_Vaccine_Uptake(ii,0);
        Monte_Carlo_Runs(ii,2500,'Sample_2025',0);
    else
        for yy=1:5
            National_Reduction_Vaccine_Uptake(ii,yy);
            Monte_Carlo_Runs(ii,2500,'Sample_2025',yy);
        end
    end
    
end
