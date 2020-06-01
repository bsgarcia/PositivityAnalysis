% This script runs the simulations with 2 sessions
close all
clear all

addpath './simulation'
addpath './fit'
addpath './utils'

data = load('data/fit_exp/data.mat');
fit_params = data.data('parameters_lpp');

%% -------------------- Define parameters ------------------------------ %
%% both 
% Define array of conditions
models = [2];
conds = repelem(1:9, 20);
tmax = length(conds);  % 8 * 24 trials
noptions = 2;
nagent = 10; % n agent per subject
ncond = length(unique(conds));
rewards = cell(ncond, 1, 1);
probs = cell(ncond, 1, 1);

j = 1;
for i = [0.75, 0.5, 0.25]
    rewards{j} = {[0, 0], [-1, 1]};
    probs{j} = {[0.5, 0.5], [1-i, i]};
    j = j + 1;
end
for i = [0.75, 0.5, 0.25]
    rewards{j} = {[-1, -1], [-2, 0]};
    probs{j} = {[0.5, 0.5], [1-i, i]};
    j = j + 1;
end
for i = [0.75, 0.5, 0.25]
    rewards{j} = {[1, 1], [0, 2]};
    probs{j} = {[0.5, 0.5], [1-i, i]};
    j = j + 1;
end

r_risk = repelem({{}}, tmax);
p_risk = repelem({{}}, tmax);

for t = 1:tmax
    r_risk{t} = rewards{conds(t)};
    p_risk{t} = probs{conds(t)};
end

%% ----------------  run simulations ----------------------------------- % 
sim_params = containers.Map(...
    {'noptions', 'tmax', 'nagent', 'p', 'r', 'conds' 'name', 'show_window', 'models'},...
    {noptions, tmax, nagent, p_risk, r_risk, conds, 'risk', true, models}...
);

data = simulation(fit_params, sim_params);
save('data/data_sim/online_exp', 'data');


