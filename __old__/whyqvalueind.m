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
model = 1;
nagent = 50;

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
tempdeltaq = zeros(nsub, tmax, nmodel);
tempout = zeros(nsub, tmax, nmodel);
for sub = 1:nsub
    %             if ~(mod(sub, 100))
    %                 progressbar(sub/nsub, nagent, sub==1,...
    %                     sprintf('Reformatting data for cond %s', condlabels{i}));
    %             end
    tempnewp(sub, :, :) = data{sub}(:, 4, :);
    tempdeltaq(sub, :, :) = data{sub}(:, 7, :);
    tempout(sub, :, :) = data{sub}(:, 2, :);
end
% Average agents with identical cognitive parameters
newp = zeros(nsub/nagent, tmax, nmodel);
newdeltaq = zeros(nsub/nagent, tmax, nmodel);
newout = zeros(nsub/nagent, tmax, nmodel);
poolsub = [1, nagent];
for sub = 1:nsub/nagent
    mysub = randsample(poolsub(1):poolsub(2), 1, true);
    for t = 1:tmax
        newp(sub, t, model) = tempnewp(mysub, t, model);
        newdeltaq(sub, t, model) = tempdeltaq(mysub, t, model);
        newout(sub, t, model) = tempout(mysub, t, model);
    end
poolsub = poolsub + nagent;
end
%% ------------------------ Plots ----------------------------- %
figure('Renderer', 'painters', 'Position', [10 10 1600 1200]);
suptitle('21 subjects randomly drawn (risk neutral)');
subjectlist = randsample(1:105, 21, false);
for i = 1:length(subjectlist)
subplot(7, 3, i);
reversalplot(...
    newp(subjectlist(i), :, model)',... %data
    [],... %time when reversal occurs
    [],... %time when cond changes
    ones(3) * 0.5,... % chance lvl
    'k',... %curve color
    0.7,... %linewidth
    0.7,... % alpha
    0, 1,... % ylims
    9,... %fontsize,
    sprintf('subject %d / \\alpha=%.2f, \\beta=%.2f',...
        subjectlist(i),...
        round(fit_params(subjectlist(i), 2, model), 2),...
        round(fit_params(subjectlist(i), 1, model), 2)),...
    'trials',... %xlabel
    'p of option 2'... % ylabel
    );
    %plot(newout(subjectlist(i), :, model), 'Color', [0.4660    0.6740    0.1880]);
end
figure('Renderer', 'painters', 'Position', [10 10 1600 1200]);
suptitle('21 subjects randomly drawn (risk neutral)');
for i = 1:length(subjectlist)
subplot(7, 3, i);
reversalplot(...
    newdeltaq(subjectlist(i), :, model)',... %data
    [],... %time when reversal occurs
    [],... %time when cond changes
    ones(3) * 0.5,... % chance lvl
    'k',... %curve color
    0.7,... %linewidth
    0.7,... % alpha
    0, 1,... % ylims
    9,... %fontsize,
    sprintf('subject %d / \\alpha=%.2f, \\beta=%.2f',...
        subjectlist(i),...
        round(fit_params(subjectlist(i), 2, model), 2),...
        round(fit_params(subjectlist(i), 1, model), 2)),...
    'trials',... %xlabel
    '\Delta Q'... % ylabel
    );
end


% figure('Renderer', 'painters', 'Position', [10 10 1600 1200]);
% suptitle('risk');
% subplot(3, 2, model);
% reversalplot(...
%     newqvalues(:, :, model)',... %data
%     [],... %time when reversal occurs
%     [0, 4*48],... %time when cond changes
%     zeros(3),... % chance lvl
%     'k',... %curve color
%     0.7,... %linewidth
%     0.7,... % alpha
%     -1, 1,... % ylims
%     9,... %fontsize,
%     models{model},... %title,
%     'trials',... %xlabel
%     'Q of option 2'... % ylabel
%     );
% 
% figure('Renderer', 'painters', 'Position', [10 10 1600 1200]);
% suptitle('risk');
% subplot(3, 2, model);
% reversalplot(...
%     newout(:, :, model)',... %data
%     [],... %time when reversal occurs
%     [0, 4*48],... %time when cond changes
%     zeros(3),... % chance lvl
%     'k',... %curve color
%     0.7,... %linewidth
%     0.7,... % alpha
%     0, 5,... % ylims
%     9,... %fontsize,
%     models{model},... %title,
%     'trials',... %xlabel
%     'Outcome'... % ylabel
%     );
