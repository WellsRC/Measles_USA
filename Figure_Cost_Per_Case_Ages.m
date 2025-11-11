function Figure_Cost_Per_Case_Ages()
close all;
[~,~,~,~,~,~,~,Cost_per_Hospitalization,Cost_per_Non_Hospitalization,~,~]=Measles_Outbreak_Cost();

[p_H_Unvaccinated,p_H_Vaccinated]=Hospitalization_Probability();

load('Monte_Carlo_Run_Sample_2025_National_Reduction=0_Ages_0_to_6.mat');

Unvaccinated_Cases_County=squeeze(sum(Unvaccinated_Cases_County_Baseline,1));
Vaccinated_Cases_County=squeeze(sum(Vaccinated_Cases_County_Baseline,1));

Hospitalizations_Baseline=zeros(size(Unvaccinated_Cases_County));
Cases_Baseline=Unvaccinated_Cases_County+Vaccinated_Cases_County;

for aa=1:length(p_H_Unvaccinated)
    Hospitalizations_Baseline(aa,:)=p_H_Unvaccinated(aa)*Unvaccinated_Cases_County(aa,:)+p_H_Vaccinated(aa)*Vaccinated_Cases_County(aa,:);
end

Hospitalizations_Baseline=[Hospitalizations_Baseline(1,:); Hospitalizations_Baseline(2,:); sum(Hospitalizations_Baseline(3:4,:),1); sum(Hospitalizations_Baseline(5:end,:),1)];
Cases_Baseline=[Cases_Baseline(1,:); Cases_Baseline(2,:); sum(Cases_Baseline(3:4,:),1); sum(Cases_Baseline(5:end,:),1)];


Direct_Medical_Cost=Cost_per_Hospitalization.*Hospitalizations_Baseline+Cost_per_Non_Hospitalization.*(Cases_Baseline-Hospitalizations_Baseline);
Cost_per_Case_Baseline=Direct_Medical_Cost./Cases_Baseline;

load('Monte_Carlo_Run_Sample_2025_National_Reduction=1_Ages_0_to_6.mat');

Unvaccinated_Cases_County=squeeze(sum(Unvaccinated_Cases_County_Baseline,1));
Vaccinated_Cases_County=squeeze(sum(Vaccinated_Cases_County_Baseline,1));

Hospitalizations_Reduction=zeros(size(Unvaccinated_Cases_County));
Cases_Reduction=Unvaccinated_Cases_County+Vaccinated_Cases_County;

for aa=1:length(p_H_Unvaccinated)
    Hospitalizations_Reduction(aa,:)=p_H_Unvaccinated(aa)*Unvaccinated_Cases_County(aa,:)+p_H_Vaccinated(aa)*Vaccinated_Cases_County(aa,:);
end

Hospitalizations_Reduction=[Hospitalizations_Reduction(1,:); Hospitalizations_Reduction(2,:); sum(Hospitalizations_Reduction(3:4,:),1); sum(Hospitalizations_Reduction(5:end,:),1)];
Cases_Reduction=[Cases_Reduction(1,:); Cases_Reduction(2,:); sum(Cases_Reduction(3:4,:),1); sum(Cases_Reduction(5:end,:),1)];


Direct_Medical_Cost=Cost_per_Hospitalization.*Hospitalizations_Reduction+Cost_per_Non_Hospitalization.*(Cases_Reduction-Hospitalizations_Reduction);
Cost_per_Case_Reduction_1=Direct_Medical_Cost./Cases_Reduction;

load('Monte_Carlo_Run_Sample_2025_National_Reduction=2.5_Ages_0_to_6.mat');

Unvaccinated_Cases_County=squeeze(sum(Unvaccinated_Cases_County_Baseline,1));
Vaccinated_Cases_County=squeeze(sum(Vaccinated_Cases_County_Baseline,1));

Hospitalizations_Reduction=zeros(size(Unvaccinated_Cases_County));
Cases_Reduction=Unvaccinated_Cases_County+Vaccinated_Cases_County;

for aa=1:length(p_H_Unvaccinated)
    Hospitalizations_Reduction(aa,:)=p_H_Unvaccinated(aa)*Unvaccinated_Cases_County(aa,:)+p_H_Vaccinated(aa)*Vaccinated_Cases_County(aa,:);
end

Hospitalizations_Reduction=[Hospitalizations_Reduction(1,:); Hospitalizations_Reduction(2,:); sum(Hospitalizations_Reduction(3:4,:),1); sum(Hospitalizations_Reduction(5:end,:),1)];
Cases_Reduction=[Cases_Reduction(1,:); Cases_Reduction(2,:); sum(Cases_Reduction(3:4,:),1); sum(Cases_Reduction(5:end,:),1)];


Direct_Medical_Cost=Cost_per_Hospitalization.*Hospitalizations_Reduction+Cost_per_Non_Hospitalization.*(Cases_Reduction-Hospitalizations_Reduction);
Cost_per_Case_Reduction_2_5=Direct_Medical_Cost./Cases_Reduction;

load('Monte_Carlo_Run_Sample_2025_National_Reduction=5_Ages_0_to_6.mat');

Unvaccinated_Cases_County=squeeze(sum(Unvaccinated_Cases_County_Baseline,1));
Vaccinated_Cases_County=squeeze(sum(Vaccinated_Cases_County_Baseline,1));

Hospitalizations_Reduction=zeros(size(Unvaccinated_Cases_County));
Cases_Reduction=Unvaccinated_Cases_County+Vaccinated_Cases_County;

for aa=1:length(p_H_Unvaccinated)
    Hospitalizations_Reduction(aa,:)=p_H_Unvaccinated(aa)*Unvaccinated_Cases_County(aa,:)+p_H_Vaccinated(aa)*Vaccinated_Cases_County(aa,:);
end

Hospitalizations_Reduction=[Hospitalizations_Reduction(1,:); Hospitalizations_Reduction(2,:); sum(Hospitalizations_Reduction(3:4,:),1); sum(Hospitalizations_Reduction(5:end,:),1)];
Cases_Reduction=[Cases_Reduction(1,:); Cases_Reduction(2,:); sum(Cases_Reduction(3:4,:),1); sum(Cases_Reduction(5:end,:),1)];


Direct_Medical_Cost=Cost_per_Hospitalization.*Hospitalizations_Reduction+Cost_per_Non_Hospitalization.*(Cases_Reduction-Hospitalizations_Reduction);
Cost_per_Case_Reduction_5=Direct_Medical_Cost./Cases_Reduction;

figure('units','normalized','outerposition',[0.125 0.125 0.5 0.5]);
baseline_b=zeros(4,1);
reduction_1=zeros(4,1);
reduction_2_5=zeros(4,1);
reduction_5=zeros(4,1);
for jj=1:4
    baseline_b(jj)=mean(Cost_per_Case_Baseline(jj,~isnan(Cost_per_Case_Baseline(jj,:))));
    reduction_1(jj)=mean(Cost_per_Case_Reduction_1(jj,~isnan(Cost_per_Case_Reduction_1(jj,:))));
    reduction_2_5(jj)=mean(Cost_per_Case_Reduction_2_5(jj,~isnan(Cost_per_Case_Reduction_2_5(jj,:))));
    reduction_5(jj)=mean(Cost_per_Case_Reduction_5(jj,~isnan(Cost_per_Case_Reduction_5(jj,:))));
end
bar([1:4],[ reduction_1(:)-baseline_b(:)  reduction_2_5(:)-baseline_b(:) reduction_5(:)-baseline_b(:)]);
x_L={['Ages 0' char(8211) '4']; ['Ages 5' char(8211) '9'];['Ages 10' char(8211) '19'];['Ages 20+'];};

legend({'1% reduction';'2.5% reduction';'5% reduction'},'FontSize',16);
set(gca,'LineWidth',2,'TickDir','out','FontSize',16,'XTick',[1:4],'XTickLabel',x_L);
xlabel('Age group','FontSize',20);
ylabel('Additional cost per case','FontSize',20);
box off;
ytickformat('$%,.2f')
end

