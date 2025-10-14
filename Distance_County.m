clear;
clc;

County_T=shaperead([pwd '\Shapefile\cb_2023_us_county_20m.shp'],'UseGeoCoords',true);

Lat_C=zeros(height(County_T),1);
Lon_C=zeros(height(County_T),1);

for cc=1:length(Lat_C)
    polyin = polyshape(County_T(cc).Lon,County_T(cc).Lat);
    [Lon_C(cc),Lat_C(cc)] = centroid(polyin);
end

Distance_Matrix=zeros(length(Lat_C));
for cc=1:length(Lat_C)
        Distance_Matrix(cc,:)=deg2sm(distance(Lat_C,Lon_C,Lat_C(cc),Lon_C(cc)));
end

County_GEOID={County_T.GEOID};



County_Data=readtable([pwd '/County_Data.xlsx'],'Sheet',['Year_2023']);

Gravity_Model=zeros(length(County_Data.County));

for cc=1:size(Gravity_Model,1) 
    parfor kk=1:size(Gravity_Model,2) 
        if(cc~=kk)
            t_cc=find(strcmp(County_Data.GEOID(cc),County_GEOID));
            t_kk=find(strcmp(County_Data.GEOID(kk),County_GEOID));
            d_ij=Distance_Matrix(t_cc,t_kk);
            Gravity_Model(cc,kk)=log(d_ij/(County_Data.Total_Population(cc).*County_Data.Total_Population(kk)));
        end
    end
end

save('County_Matrix_Gravity_Covariate.mat',"Gravity_Model",'County_GEOID')
