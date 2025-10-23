function [Importation_Cases_County] = Case_Importation_Sample(Type,NS)

Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')

if(strcmp('Baseline',Type))
    Importation_Cases_County=zeros(length(County_Data.County),1);
    Measles_Cases=readtable('County_Level_Measles_Cases_Adjusted.csv');    
    for cc=1:length(Importation_Cases_County)
        t_f=str2double(County_Data.GEOID{cc})==Measles_Cases.GEOID & strcmp(Measles_Cases.type,'imported') & ~isnan(Measles_Cases.case_count);
        if(sum(t_f)>0)
            Importation_Cases_County(cc)=Measles_Cases.case_count(t_f);
        end
    end
else
    Importation_Cases_County=zeros(length(County_Data.County),NS);
    State_Importaton=readtable('NNDSS_Weekly_Data_Meales_Importation_2023_to_2025.csv');
    U_State=unique(State_Importaton.ReportingArea);
    for ss=1:length(U_State)    
        f_import=strcmp(State_Importaton.ReportingArea,U_State{ss});
        w_import=State_Importaton.Weekly_Importation(f_import);
        samp_import=w_import(randi(length(w_import),52,NS));
        samp_import=sum(samp_import,1);
        if(strcmp(U_State{ss},'New York City'))
            f_county=strcmp(County_Data.State,'New York') & (strcmp(County_Data.County,'New York') | strcmp(County_Data.County,'Kings') | strcmp(County_Data.County,'Bronx')  | strcmp(County_Data.County,'Richmond')  | strcmp(County_Data.County,'Queens'));
        elseif(strcmp(U_State{ss},'New York'))
            f_county=strcmp(County_Data.State,U_State{ss}) & ~(strcmp(County_Data.County,'New York') | strcmp(County_Data.County,'Kings') | strcmp(County_Data.County,'Bronx')  | strcmp(County_Data.County,'Richmond')  | strcmp(County_Data.County,'Queens'));
        else
            f_county=strcmp(County_Data.State,U_State{ss});
        end
        import_weight=County_Data.Total_Population(f_county)./sum(County_Data.Total_Population(f_county));
        temp_county_import=zeros(sum(f_county),NS);
        w_temp=import_weight;
        if(sum(w_temp)~=1)
            w_temp=w_temp./sum(w_temp);
        end
        if(sum(w_temp)~=1)
            dx=sum(w_temp)-1;
            w_temp(max(w_temp)==w_temp)=w_temp(max(w_temp)==w_temp)-dx;
        end
        if(sum(w_temp)~=1)
            w_temp=round(w_temp,16)./sum(round(w_temp,16));
        end
        pd = makedist('Multinomial','Probabilities',w_temp);
        for nn=1:NS
            if(samp_import(nn)>0)
                r = random(pd,1,samp_import(nn));
                for cc=1:size(temp_county_import,1)
                    temp_county_import(cc,nn)=sum(r==cc);
                end
            end
        end
        Importation_Cases_County(f_county,:)=temp_county_import;
    end

end

end

