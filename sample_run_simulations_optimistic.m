% This script runs the simulations with 2 sessions
close all
clear all

addpath './simulation'
addpath './fit'
addpath './utils'

data = load('data/fit_exp/data.mat');
fit_params = data.data('parameters_lpp');
models = data.data('whichmodel');

%% -------------------- Define parameters ------------------------------ %
%% both 
% Define array of conditions
conds = repelem(1:8, 24);
tmax = length(conds);  % 8 * 24 trials
noptions = 2;
nagent = 1; % n agent per subject
ncond = length(unique(conds));
rewards = cell(ncond, 1, 1);
probs = cell(ncond, 1, 1);

% AB
for i = [1, 1+4]
    rewards{i} = {[-1, 1], [-1, 1]};
    probs{i} = {[0.75, 0.25], [0.75, 0.25]};
end
% CD and EF
for i = [2, 2+4]
    rewards{i} = {[-1, 1], [-1, 1]};
    probs{i} = {[0.75, 0.25], [0.25, 0.75]};
end
% GH
for i = [3, 3+4]
    rewards{i} = {[-1,  1], [-1, 1]};
    probs{i} = {[0.25,0.75], [0.75, 0.25]};
end
for i = [4, 4+4]
    rewards{i} = {[-1,  1], [-1, 1]};
    probs{i} = {[0.25,0.75], [0.25, 0.75]};
end

r_risk = repelem({{}}, tmax);
p_risk = repelem({{}}, tmax);

for t = 1:tmax
    r_risk{t} = rewards{conds(t)};
    p_risk{t} = probs{conds(t)};
end

%reversals = [84:96, 84*2:192];
%% ----------------  run simulations ----------------------------------- % 
data = cell(1);
w = waitbar(0, 'Computing simulations');
for i = 1:100
    waitbar(i/100, w, sprintf('Computing simulations %d %%', i));
    sim_params = containers.Map(...
        {'noptions', 'tmax', 'nagent', 'p', 'r', 'conds' 'name', 'show_window', 'models'},...
        {noptions, tmax, nagent, p_risk, r_risk, conds, 'risk', false, models}...
    );

    idx = 1:85;
    selected_fit_params = fit_params(idx, :, :);
    data{i} = simulation(selected_fit_params, sim_params);
end
close(w);
save('data/data_sim/35_opt', 'data');


