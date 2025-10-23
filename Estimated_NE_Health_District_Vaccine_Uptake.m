function v = Estimated_NE_Health_District_Vaccine_Uptake(v_county,pop_county,GEOID_county,Health_District_GEOIDs)

v=zeros(length(Health_District_GEOIDs),1);
for ii=1:length(Health_District_GEOIDs)
    temp_c=Health_District_GEOIDs{ii}{:};
    tf=ismember(GEOID_county,temp_c);
    v(ii)=sum(v_county(tf).*pop_county(tf))./sum(pop_county(tf));
end

end