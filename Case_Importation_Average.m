function [Importation_Cases_County] = Case_Importation_Average(Type)

Vaccine='MMR';
load([Vaccine '_Immunity.mat'],'County_Data')

    Importation_Cases_County=zeros(length(County_Data.County),1);
    State_Importaton=readtable('NNDSS_Weekly_Data_Meales_Importation_2023_to_2025.csv');
    if(strcmp('Sample_2023',Type))
        State_Importaton=State_Importaton(State_Importaton.CurrentMMWRYear==2023 & State_Importaton.MMWRWEEK==52,:);
    elseif(strcmp('Sample_2024',Type))
        State_Importaton=State_Importaton(State_Importaton.CurrentMMWRYear==2024 & State_Importaton.MMWRWEEK==52,:);
    end
    U_State=unique(State_Importaton.ReportingArea);
    for ss=1:length(U_State)    
        f_import=strcmp(State_Importaton.ReportingArea,U_State{ss});
        annual_import=State_Importaton.CumulativeYTDCurrentMMWRYear(f_import);
        if(isnan(annual_import))
            annual_import=0;
        end
        if(strcmp(U_State{ss},'New York City'))
            f_county=strcmp(County_Data.State,'New York') & (strcmp(County_Data.County,'New York') | strcmp(County_Data.County,'Kings') | strcmp(County_Data.County,'Bronx')  | strcmp(County_Data.County,'Richmond')  | strcmp(County_Data.County,'Queens'));
        elseif(strcmp(U_State{ss},'New York'))
            f_county=strcmp(County_Data.State,U_State{ss}) & ~(strcmp(County_Data.County,'New York') | strcmp(County_Data.County,'Kings') | strcmp(County_Data.County,'Bronx')  | strcmp(County_Data.County,'Richmond')  | strcmp(County_Data.County,'Queens'));
        else
            f_county=strcmp(County_Data.State,U_State{ss});
        end
        import_weight=County_Data.Total_Population(f_county)./sum(County_Data.Total_Population(f_county));
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

        Importation_Cases_County(f_county)=w_temp.*annual_import;
    end

end


