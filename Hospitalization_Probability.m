function [p_H_Unvaccinated,p_H_Vaccinated,duration_hospitalization]=Hospitalization_Probability()
% https://academic.oup.com/cid/article/80/3/663/7756619?login=false#508422746
% Supplemental Table 3 and 4
    p_H_Unvaccinated=zeros(1,18);
    p_H_Vaccinated=zeros(1,18);
    duration_hospitalization=zeros(1,18);
% https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0231329
    duration_hospitalization(1:2)=3.84;
    duration_hospitalization(3:4)=5.52;
    duration_hospitalization(5:8)=5.33;
    duration_hospitalization(9:end)=5.44;

    p_H_Unvaccinated(1)=(112+178)/(112+248+178+486);
    p_H_Unvaccinated(2:4)=71/(71+540);
    p_H_Unvaccinated(5:end)=112/(112+170);

    p_H_Vaccinated(1)=20/(20+71+7);
    p_H_Vaccinated(2:4)=9/(9+72);
    p_H_Vaccinated(5:end)=38/(38+161);
end