% This script runs the simulations with 2 sessions
close all
clear all

addpath './simulation'
addpath './fit'
addpath './utils'

data = load('data/fit_exp/data');
fit_params = data.data('parameters_lpp');

%% -------------------- Define parameters ------------------------------ %

models = [2];
%% both 
% Define array of conditions
conds = repelem(1:8, 24);
tmax = length(conds);  % 8 * 24 trials
noptions = 2;
nagent = 10; % n agent per subject
ncond = length(unique(conds));
rewards = cell(ncond, 1, 1);
probs = cell(ncond, 1, 1);

% AB
for i = [1, 1+4]
    rewards{i} = {[0, 1], [0, 1]};
    probs{i} = {[0.5, 0.5], [0.5, 0.5]};
end
% CD and EF
for i = [2, 3, 2+4, 3+4]
    rewards{i} = {[0, 1], [0, 1]};
    probs{i} = {[0.75, 0.25], [0.25, 0.75]};
end
% GH
for i = [4, 4+4]
    rewards{i} = {[0,  1], [0, 1]};
    probs{i} = {[0.17, 0.83], [0.83, 0.17]};
end

r_risk = repelem({{}}, tmax);
p_risk = repelem({{}}, tmax);

reversals = [84:96, 84*2:192];

for t = 1:tmax
    if ismember(t, reversals) 
        r_risk{t} = fliplr(rewards{conds(t)});
        p_risk{t} = fliplr(probs{conds(t)});
    else
        r_risk{t} = rewards{conds(t)};
        p_risk{t} = probs{conds(t)};
    end
end


%% ----------------  run simulations ----------------------------------- % 
sim_params = containers.Map(...
    {'noptions', 'tmax', 'nagent', 'p', 'r', 'conds' 'name', 'show_window', 'models'},...
    {noptions, tmax, nagent, p_risk, r_risk, conds, 'risk', true, models}...
);

data = simulation(fit_params, sim_params);
save('data/data_sim/conf', 'data');


