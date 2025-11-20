clear;

for ii=25:25:200
        Monte_Carlo_Runs(0,2500,['Scenario_' num2str(ii)],0);
        Monte_Carlo_Runs(0.0025,2500,['Scenario_' num2str(ii)],1);
        Monte_Carlo_Runs(0.005,2500,['Scenario_' num2str(ii)],1);
        Monte_Carlo_Runs(0.0075,2500,['Scenario_' num2str(ii)],1);
        Monte_Carlo_Runs(0.01,2500,['Scenario_' num2str(ii)],1);
end