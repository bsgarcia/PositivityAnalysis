% This script runs the simulations with 2 sessions
close all
clear all

addpath './simulation'
addpath './fit'
addpath './utils'

data = load('data/fit_exp/parameters.mat');
fit_params = data.parameters;

%% -------------------- Define common parameters --------------------- %
%% both 
% Define array of conditions
conds = repelem(1:8, 24);
tmax = length(conds);  % 8 * 24 trials

%% -------------------- Define params status quo conditions -----------%
% when reversal occurs
t_reversal = containers.Map(...
    {'statusquo1', 'statusquo2'},...
    {[36, 56, 64, 78, 84, 90, 132, 154, 162, 176, 182, 188],...
    [30, 60, 90, 126, 156, 186]});

% rewards are always the same
r = containers.Map(...
    {'statusquo1', 'statusquo2'},...
    {repelem({{[-1, 1], [-1, 1]}}, tmax), ...
    repelem({{[-1, 1], [-1, 1]}}, tmax)});

% define probabilities
p = containers.Map({'statusquo1', 'statusquo2'},...
    {repelem({{[0.25, 0.75], [0.75, 0.25]}}, tmax),...
    repelem({{[0.25, 0.75], [0.75, 0.25]}}, tmax)});

% apply reversals to probabilities (rewards do not change)
k = keys(t_reversal);
for i = 1:length(k)
    temp_t_reversal = t_reversal(k{i});
    for t = 1:tmax/2
        if ismember(t, temp_t_reversal)
            ptemp = p(k{i});
            for t2 = t:tmax
                ptemp{t2} = flip(ptemp{t2});
            end
            p(k{i}) = ptemp;
        end
    end
    clear ptemp;
    for t = tmax/2:tmax
        if ismember(t, temp_t_reversal)
            ptemp = p(k{i});
            for t2 = t:tmax
                ptemp{t2} = flip(ptemp{t2});
            end
            p(k{i}) = ptemp;
        end
    end
    clear ptemp;
end

%% -------------------- Define params risk condition --------------------%
ncond = length(unique(conds));
rewards = cell(ncond, 1, 1);
probs = cell(ncond, 1, 1);

% AB
for i = [1, 1+4]
    rewards{i} = {[0, 0], [-1, 1]};
    probs{i} = {[0.75, 0.25], [0.25, 0.75]};
end
% CD
for i = [2, 2+4]
    rewards{i} = {[0, 0], [-1, 1]};
    probs{i} = {[0.75, 0.25], [0.75, 0.25]};
end
% EF
for i = [3, 3+4]
    rewards{i} = {[0, 0], [-1, 1]};
    probs{i} = {[0.75, 0.25], [0.5, 0.5]};
end
% GH
for i = [4, 4+4]
    % GH
    rewards{i} = {[0,  -1], [0, 1]};
    probs{i} = {[0.5, 0.5], [0.5, 0.5]};
end

%rewards{i} = {[0,  -1], [-1, 1]};
   % probs{i} = {[0.75, 0.25], [0.75, 0.25]};

r_risk = repelem({{}}, tmax);
p_risk = repelem({{}}, tmax);

for t = 1:tmax
    r_risk{t} = rewards{conds(t)};
    p_risk{t} = probs{conds(t)};
end

% add risk rewards and probs to map object
r('risk') = r_risk;
p('risk') = p_risk;

%% ----------------  run simulations ---------------------------- % 
% get new keys
k = keys(r);
for i = 1:length(k)
    
    if strcmp(k{i}, 'risk')
    
        sim_params = containers.Map(...
            {'tmax', 'p', 'r', 'conds' 'name'},...
            {tmax, p(k{i}), r(k{i}), conds, k{i}}...
        );

        data = simulation(fit_params, sim_params);
        save(sprintf('data/data_sim_2sessions/%s_new', k{i}), 'data');
        clear data;
    end
end

