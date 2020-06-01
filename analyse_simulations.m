close all
clear all

addpath './simulation'
addpath './fit'
addpath './utils'

%% ----------------------- Variables to modify -------------------------- % 
whichmodel = [1, 2, 3, 4, 5, 8];
nmodel = length(whichmodel);
nagent = 10;
%% ---------------------- Define variables --------------------- ---------%

condlabels = {'risk', 'statusquo1', 'statusquo2'};
models = {'QLearning', 'Asymmetric',...
    'AsymmetricPessimistic', 'Priors', 'Impulsive Perseveration',...
    'Gradual Perseveration', 'Full'};
    
% when reversal occurs
t_reversal = [];
t_cond = [0, 24, 2*24, 3*24, 4*24, 5*24, 6*24, 7*24, (8*24+1)];

%% ----------------------Loop over conditions----------------------------% 

[con, con2, cho, out, nsub, p, correct, qvalues] = load_data('sim', 'risk');
tmax = length(p{1}(:, 1));

try
    %% ------------------- Try to load formatted data --------------- % 
    data = load('data/analyse/risk');
    newdata = data.data;
    newp = newdata('newp');
    newcorrect = newdata('newcorrect');
    newqvalues = newdata('newqvalues');
    
    %disp(std(ndata(nagent, :, 1)));
catch
    %% --------------------- Reformat data in order to analyse ------ %
    % Format probabilitites and correct choices
    % into matrices
    tempnewp = zeros(nsub, tmax, nmodel);
    tempnewcorrect = zeros(nsub, tmax, nmodel);
    tempqvalues = zeros(nsub, tmax, nmodel);
    for sub = 1:nsub
        tempnewp(sub, :, :) = p{sub}(:, whichmodel);    
        tempnewcorrect(sub, :,  :) = correct{sub}(:, whichmodel);
        tempqvalues(sub, :, :) = qvalues{sub}(:, whichmodel);
    end
    % Average agents with identical cognitive parameters
    newp = zeros(nsub/nagent, tmax, nmodel);
    newcorrect = zeros(nsub/nagent, tmax, nmodel);
    newqvalues = zeros(nsub/nagent, tmax, nmodel);
    poolsub = [1, nagent];
    for sub = 1:nsub/nagent          
        for model = 1:nmodel
            for t = 1:tmax
                newp(sub, t, model) = mean(...
                    tempnewp(poolsub(1):poolsub(2), t, model));
                newcorrect(sub, t, model) = mean(...
                    tempnewcorrect(poolsub(1):poolsub(2), t, model));
                newqvalues(sub, t, model) = mean(...
                    tempqvalues(poolsub(1):poolsub(2), t, model));
            end
        end
       poolsub = poolsub + nagent;
    end
    % save
    data = containers.Map(...
        {'newp', 'newcorrect', 'newqvalues'},...
        {newp, newcorrect, newqvalues});
    save('data/analyse/risk', 'data');
end

%% ------------------------ Plots ----------------------------- %
% probabilities
% --------------------------------------------------------------
f = figure('Renderer', 'painters', 'Position', [10 10 1800 1000]);
suptitle('risk');
for model = 1:nmodel
    subplot(round(nmodel/2), 2, model);
    reversalplot(...
        newp(:, :, model)',... %data
        t_reversal,... %time when reversal occurs
        t_cond,... %time when cond changes
        ones(3) * 0.5,... % chance lvl
        'k',... %curve color
        0.7,... %linewidth
        0.7,... % alpha
        0, 1,... % ylims
        9,... %fontsize,
        models{whichmodel(model)},... %title,
        'trials',... %xlabel
        'p of option 2'... % ylabel
    ); 
     
end
saveas(f, 'fig/2sessions_analysis_proba_risk.png');

% --------------------------------------------------------------
% correct choices
% --------------------------------------------------------------
f = figure('Renderer', 'painters', 'Position', [10 10 1800 1000]);
suptitle('risk');
for model = 1:nmodel
    subplot(round(nmodel/2), 2, model);
    reversalplot(...
        newcorrect(:, :, model)',... %data
        t_reversal,... %time when reversal occurs
        t_cond,... %time when cond changes
        ones(3) * 0.5,... % chance lvl
        'k',... %curve color
        0.7,... %linewidth
        0.7,... % alpha
        0, 1,... % ylims
        9,... %fontsize,
        models{whichmodel(model)},... %title,
        'trials',... %xlabel
        '% of choices maximizing utility'... % ylabel
    ); 
end
saveas(f, 'fig/2sessions_analysis_choice_risk.png');


% --------------------------------------------------------------
% Q values
% --------------------------------------------------------------
f = figure('Renderer', 'painters', 'Position', [10 10 1800 1000]);
suptitle('risk');
for model = 1:nmodel
    subplot(round(nmodel/2), 2, model);
    reversalplot(...
        newqvalues(:, :, model)',... %data
        t_reversal,... %time when reversal occurs
        t_cond,... %time when cond changes
        zeros(3),... % chance lvl
        'k',... %curve color
        0.7,... %linewidth
        0.7,... % alpha
        -1, 1,... % ylims
        9,... %fontsize,
        models{whichmodel(model)},... %title,
        'trials',... %xlabel
        'Q of option 2'... % ylabel
    );    
end
saveas(f, 'fig/2sessions_analysis_qvalue_risk.png');
