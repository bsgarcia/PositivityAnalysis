% This script runs the simulations
close all
clear all

addpath './simulation'
addpath './fit'
addpath './utils'

data = load('data/fit_exp/parameters.mat');
fit_params = data.parameters;


%% both 
% Define array of conditions
conds = repelem(1, 48*4);
tmax = length(conds);  % 4 * 48 trials
nmodel = 6;
whichmodel = 1:nmodel;
models = {'QLearning', 'Asymmetric', 'AsymmetricPessimistic',...
        'Perseveration', 'Priors', 'Full'};

ncond = length(unique(conds));
rewards = cell(ncond, 1, 1);
probs = cell(ncond, 1, 1);

rewards{1} = {[0, 0], [-1, 1]};
probs{1} = {[0.75, 0.25], [0.5, 0.5]};
r_risk = repelem({{}}, tmax);
p_risk = repelem({{}}, tmax);

for t = 1:tmax
    r_risk{t} = rewards{conds(t)};
    p_risk{t} = probs{conds(t)};
end

% add risk rewards and probs to map object
% rewards are always the same
r = containers.Map(...
    {'risk'},...
    {r_risk});

p = containers.Map(...
    {'risk'},...
    {p_risk});

%% ----------------  run simulations ---------------------------- % 
try
    tdata = load('data/analyse/partial_risk');
    data = tdata.data;
catch
    sim_params = containers.Map(...
            {'tmax', 'p', 'r', 'conds' 'name'},...
            {tmax, p('risk'), r('risk'), conds, 'risk'}...
     );

    data = simulation(fit_params, sim_params);

    save('data/analyse/partial_risk', 'data');
end
%% --------------------- Reformat data in order to analyse ------ %
% Format probabilitites and correct choices
% into matrices
nsub = length(data);
nmodel = 6;
tempnewp = zeros(nsub, tmax, nmodel);
tempnewcorrect = zeros(nsub, tmax, nmodel);
tempqvalues = zeros(nsub, tmax, nmodel);
tempout = zeros(nsub, tmax, nmodel);
for sub = 1:nsub
    %             if ~(mod(sub, 100))
    %                 progressbar(sub/nsub, 50, sub==1,...
    %                     sprintf('Reformatting data for cond %s', condlabels{i}));
    %             end
    tempnewp(sub, :, :) = data{sub}(:, 4, :);
    tempnewcorrect(sub, :,  :) = data{sub}(:, 5, :);
    tempqvalues(sub, :, :) = data{sub}(:, 6, :);
    tempout(sub, :, :) = data{sub}(:, 2, :);
end
% Average agents with identical cognitive parameters
newp = zeros(nsub/50, tmax, nmodel);
newcorrect = zeros(nsub/50, tmax, nmodel);
newqvalues = zeros(nsub/50, tmax, nmodel);
newout = zeros(nsub/50, tmax, nmodel);
poolsub = [1, 50];
for sub = 1:nsub/50
    for model = whichmodel
        for t = 1:tmax
            newp(sub, t, model) = mean(...
                tempnewp(poolsub(1):poolsub(2), t, model));
            newcorrect(sub, t, model) = mean(...
                tempnewcorrect(poolsub(1):poolsub(2), t, model));
            newqvalues(sub, t, model) = mean(...
                tempqvalues(poolsub(1):poolsub(2), t, model));
            newout(sub, t, model) = sum(...
                tempout(poolsub(1):poolsub(2), t, model));
        end
    end
    poolsub = poolsub + 50;
end

%% ------------------------ Plots ----------------------------- %
figure('Renderer', 'painters', 'Position', [10 10 1600 1200]);
suptitle('risk');
for model = whichmodel
    subplot(3, 2, model);
    reversalplot(...
        newp(:, :, model)',... %data
        [],... %time when reversal occurs
        [0, 4*48],... %time when cond changes
        ones(3) * 0.5,... % chance lvl
        'k',... %curve color
        0.7,... %linewidth
        0.7,... % alpha
        0, 1,... % ylims
        9,... %fontsize,
        models{model},... %title,
        'trials',... %xlabel
        'p of option 2'... % ylabel
    );
    
end

figure('Renderer', 'painters', 'Position', [10 10 1600 1200]);
suptitle('risk');
for model = whichmodel
    subplot(3, 2, model);
    reversalplot(...
        newqvalues(:, :, model)',... %data
        [],... %time when reversal occurs
        [0, 4*48],... %time when cond changes
        zeros(3),... % chance lvl
        'k',... %curve color
        0.7,... %linewidth
        0.7,... % alpha
        -1, 1,... % ylims
        9,... %fontsize,
        models{model},... %title,
        'trials',... %xlabel
        'Q of option 2'... % ylabel
        );
end

figure('Renderer', 'painters', 'Position', [10 10 1600 1200]);
suptitle('risk');
for model = whichmodel
    subplot(3, 2, model);
    reversalplot(...
        newout(:, :, model)',... %data
        [],... %time when reversal occurs
        [0, 4*48],... %time when cond changes
        zeros(3),... % chance lvl
        'k',... %curve color
        0.7,... %linewidth
        0.7,... % alpha
        0, 5,... % ylims
        9,... %fontsize,
        models{model},... %title,
        'trials',... %xlabel
        'Outcome'... % ylabel
        );
end