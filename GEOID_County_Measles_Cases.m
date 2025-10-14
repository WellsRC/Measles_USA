clear;
clc;

County_Data=readtable([pwd '/County_Data.xlsx'],'Sheet',['Year_2023']);
Measles_Cases=readtable('County_Level_Measles_Cases_Adjusted.csv');

Measles_Cases.GEOID=cell(height(Measles_Cases),1);
for cc=1:height(Measles_Cases)
    tf=strcmp(Measles_Cases.county{cc},County_Data.County) & strcmp(Measles_Cases.state{cc},County_Data.State);
    if(sum(tf)>0)
        Measles_Cases.GEOID{cc}=County_Data.GEOID(tf);
    end
end

writetable(Measles_Cases,'County_Level_Measles_Cases_Adjusted.csv');