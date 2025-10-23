clear;
clc;

County_T=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],'UseGeoCoords',true);

Lat_C=zeros(height(County_T),1);
Lon_C=zeros(height(County_T),1);

for cc=1:length(Lat_C)
    polyin = polyshape(County_T(cc).Lon,County_T(cc).Lat);
    [Lon_C(cc),Lat_C(cc)] = centroid(polyin);
end

temp_Distance_Matrix=zeros(length(Lat_C));
for cc=1:length(Lat_C)
    temp_Distance_Matrix(cc,:)=deg2sm(distance(Lat_C,Lon_C,Lat_C(cc),Lon_C(cc)));
end

County_GEOID={County_T.GEOID};



County_Data=readtable([pwd '/County_Data.xlsx'],'Sheet',['Year_2023']);

Population_i=zeros(length(County_Data.County));
Population_j=zeros(length(County_Data.County));
Distance_Matrix_ij=zeros(length(County_Data.County));
for cc=1:size(Population_i,1) 
    Population_i(cc,:)=County_Data.Total_Population(cc).*ones(1,length(County_Data.County));
    Population_j(:,cc)=County_Data.Total_Population(cc).*ones(1,length(County_Data.County));
    t_cc=find(strcmp(County_Data.GEOID(cc),County_GEOID));
    parfor kk=1:size(Population_i,2) 
        t_kk=find(strcmp(County_Data.GEOID(kk),County_GEOID));
        Distance_Matrix_ij(cc,kk)=temp_Distance_Matrix(t_cc,t_kk);
    end
end

save('County_Matrix_Gravity_Covariates.mat',"Distance_Matrix_ij",'Population_j','Population_i','County_GEOID')
