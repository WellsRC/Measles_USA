function [p_H_Unvaccinated,p_H_Vaccinated]=Hospitalization_Probability()
    p_H_Unvaccinated=zeros(1,18);
    p_H_Vaccinated=zeros(1,18);

    p_H_Unvaccinated(1)=178/(178+486);
    p_H_Unvaccinated(2:4)=71/(71+540);
    p_H_Unvaccinated(5:end)=112/(112+170);

    p_H_Vaccinated(1)=20/(20+71);
    p_H_Vaccinated(2:4)=9/(9+72);
    p_H_Vaccinated(5:end)=38/(38+161);
end