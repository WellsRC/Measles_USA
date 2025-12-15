clear;
clc;

load('MMR_Immunity.mat','County_Data');
S=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],"UseGeoCoords",true);
V_Baseline=NaN.*zeros(length(S),6);

Pop_0_6=NaN.*zeros(length(S),1);
for cc=1:length(S)
    tf=strcmp(County_Data.GEOID,S(cc).GEOID);
    if(sum(tf)>0)
        for jj=1:5
            V_Baseline(cc,jj)=table2array(County_Data.Vaccine_Uptake(tf,jj));
        end
        V_Baseline(cc,6)=(table2array(County_Data.Vaccine_Uptake(tf,1)).*table2array(County_Data.Population(tf,1))+(2/5).*table2array(County_Data.Vaccine_Uptake(tf,2)).*table2array(County_Data.Population(tf,2)))./(table2array(County_Data.Population(tf,1))+(2/5).*table2array(County_Data.Population(tf,2)));

        Pop_0_6(cc)=(table2array(County_Data.Population(tf,1))+(2/5).*table2array(County_Data.Population(tf,2)));
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%Plot Vaccine Uptake
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
 Cb=flip({'#ffffff';
     '#fff7ec';
'#fee8c8';
'#fdd49e';
'#fdbb84';
'#fc8d59';
'#ef6548';
'#d7301f';
'#b30000';
'#7f0000';});

 for jj=1:6
     Measure_Baseline=V_Baseline(:,jj);
    
    x_baseline=linspace(0.65,0.95,length(Cb));
    C_Baseline=hex2rgb(Cb);
    
    if(jj==6)
        X_Label_Baseline=['MMR coverage among children 0' char(8211) '6 years of age'];
    else
        X_Label_Baseline=['MMR coverage among children ' num2str(5.*(jj-1)) char(8211) num2str(5.*jj-1) ' years of age'];
    end
    
    prct_label=true;
    monitary_label=false;
    
    text_v=[0.65:0.1:0.95];
    
    
    inq_txt_baseline=zeros(size(text_v));
    inq_txt_baseline(1)=-1;
    inq_txt_baseline(end)=1;
    
    Gen_Figure_Single_County(Measure_Baseline,x_baseline,C_Baseline,text_v,inq_txt_baseline,X_Label_Baseline,prct_label,monitary_label,S,jj);

    print(gcf,['Figure_S1' char(64+jj) '.png'],'-dpng','-r300');

 end
