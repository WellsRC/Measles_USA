function [p_H_Unvaccinated,p_H_Vaccinated]=Severe_Disease_Probability()
% https://academic.oup.com/cid/article/80/3/663/7756619?login=false#508422746
% Supplemental Table 3 and 4
    p_H_Unvaccinated=zeros(1,18);
    p_H_Vaccinated=zeros(1,18);

    p_H_Unvaccinated(1)=(121+198)/(121+349+198+677);
    p_H_Unvaccinated(2:4)=84/(84+904);
    p_H_Unvaccinated(5:end)=116/(116+350);

    p_H_Vaccinated(1)=(1+21)/(21+95+1+6);
    p_H_Vaccinated(2:4)=10/(10+94);
    p_H_Vaccinated(5:end)=40/(40+208);
end