%% Define the parameters

rand('state',sum(100*clock));
close all
clear all

trials   =50;                           % N of trials per conditions
nsubjects=1000;                         % N of virtual subjects
                           
contingenciesRe=[.25];              % Probability of +1 for the risk option, risk bad
contingenciesRi=[.50];              % Probability of +1 for the risk option, risk bad
contingenciesRo=[.75];              % Probability of +1 for the risk option, risk bad

paramsub =[0.1 0.1 0.1];                % beta / alpha(+) / alpha(-), 
paramsub2=[0.1 0.2 0.2];                % beta / alpha(+) / alpha(-), 
paramsub3=[0.1 0.3 0.3];                % beta / alpha(+) / alpha(-), 

%% Run the simulations
for n=1:nsubjects;
    
    [choicesRe(n,:),  outcomesRe(n,:), probaRe(n,:)]  = optimisticRisk(paramsub,trials,contingenciesRe);
    [choicesRe2(n,:), outcomesRe2(n,:),probaRe2(n,:)] = optimisticRisk(paramsub2,trials,contingenciesRe);
    [choicesRe3(n,:), outcomesRe3(n,:),probaRe2(n,:)] = optimisticRisk(paramsub3,trials,contingenciesRe);
    
    [choicesRi(n,:),  outcomesRi(n,:), probaRi(n,:)]  = optimisticRisk(paramsub,trials,contingenciesRi);
    [choicesRi2(n,:), outcomesRi2(n,:),probaRi2(n,:)] = optimisticRisk(paramsub2,trials,contingenciesRi);
    [choicesRi3(n,:), outcomesRi3(n,:),probaRi3(n,:)] = optimisticRisk(paramsub3,trials,contingenciesRi);
    
    [choicesRo(n,:),  outcomesRo(n,:), probaRo(n,:)]  = optimisticRisk(paramsub,trials,contingenciesRo);
    [choicesRo2(n,:), outcomesRo2(n,:),probaRo2(n,:)] = optimisticRisk(paramsub2,trials,contingenciesRo);
    [choicesRo3(n,:), outcomesRo3(n,:),probaRo3(n,:)] = optimisticRisk(paramsub3,trials,contingenciesRo);

end

%% Plot the simulations
figure;
subplot(1,3,1);
SurfaceCurvePlot(choicesRe'-1,[0.25 0.5 0.75],[0 0 0.5],1,0.2,0,1,14,'','','');%
hold on
SurfaceCurvePlot(choicesRe2'-1,[0.25 0.5 0.75],[0.5 0 0],1,0.2,0,1,14,'','','');
hold on
SurfaceCurvePlot(choicesRe3'-1,[0.25 0.5 0.75],[0.5 0.5 0.5],1,0.2,0,1,14,'Risky bad','trial','p(risky)');
subplot(1,3,2);
SurfaceCurvePlot(choicesRi'-1,[0.25 0.5 0.75],[0 0 0.5],1,0.2,0,1,14,'','','');%
hold on
SurfaceCurvePlot(choicesRi2'-1,[0.25 0.5 0.75],[0.5 0 0],1,0.2,0,1,14,'','','');
hold on
SurfaceCurvePlot(choicesRi3'-1,[0.25 0.5 0.75],[0.5 0.5 0.5],1,0.2,0,1,14,'Risky neutral','trial','p(risky)');
subplot(1,3,3)
SurfaceCurvePlot(choicesRo'-1,[0.25 0.5 0.75],[0 0 0.5],1,0.2,0,1,14,'','','');%
hold on
SurfaceCurvePlot(choicesRo2'-1,[0.25 0.5 0.75],[0.5 0 0],1,0.2,0,1,14,'','','');
hold on
SurfaceCurvePlot(choicesRo3'-1,[0.25 0.5 0.75],[0.5 0.5 0.5],1,0.2,0,1,14,'Risky good','trial','p(risky)');

%%
figure
SurfaceCurvePlot(choicesRi'-1,[0.25 0.5 0.75],[0 0 0.5],1,0.2,0,1,14,'','','');%
hold on
SurfaceCurvePlot(choicesRi2'-1,[0.25 0.5 0.75],[0.5 0 0],1,0.2,0,1,14,'','','');
hold on
SurfaceCurvePlot(choicesRi3'-1,[0.25 0.5 0.75],[0.5 0.5 0.5],1,0.2,0,1,14,'Risky neutral','trial','p(risky)');

%%
figure;plot(choicesRo3')

