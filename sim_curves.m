close all
clear all

addpath './simulation'
addpath './fit'
addpath './utils'

%% ----------------------- Variables to modify -------------------------- % 
whichmodel = [1, 2, 3, 4, 5];
nmodel = length(whichmodel);
nagent = 10;
%% ---------------------- Define variables --------------------- ---------%

condlabels = {'risk', 'statusquo1', 'statusquo2'};
models = {'RW', 'RW\pm',...
    'RW\pm_{\omega^-}', 'RW_\omega', 'RW_\phi',...
    'Bayesian'};
    
% when reversal occurs
t_reversal = [];
t_cond = [0, 20, 2*20, 3*20, 4*20, 5*20, 6*20, 7*20, 8*20, 9*(20+1)];

%% ----------------------Loop over conditions----------------------------% 

[con, con2, cho, out, nsub, nsim, p, correct, qvalues] = load_data('sim', 'sim_online_exp');
tmax = length(p{1}(:, 1));

try
    %% ------------------- Try to load formatted data --------------- % 
    data = load('data/analyse/online_exp');
    newdata = data.data;
    newp = newdata('newp');
    newcorrect = newdata('newcorrect');
    newqvalues = newdata('newqvalues');
    newcho = newdata('newcho');
    
    %disp(std(ndata(nagent, :, 1)));
catch
    %% --------------------- Reformat data in order to analyse ------ %
    % Format probabilitites and correct choices
    % into matrices
    tempnewp = zeros(nsub, tmax, nmodel);
    tempnewcorrect = zeros(nsub, tmax, nmodel);
    tempqvalues = zeros(nsub, tmax, nmodel);
    tempnewcho = zeros(nsub, tmax, nmodel);
    for sub = 1:nsub
        tempnewp(sub, :, :) = p{sub}(:, whichmodel);    
        tempnewcorrect(sub, :,  :) = correct{sub}(:, whichmodel);
        tempqvalues(sub, :, :) = qvalues{sub}(:, whichmodel);
        tempnewcho(sub, :, :) = cho{sub}(:, whichmodel);
    end
    % Average agents with identical cognitive parameters
    newp = zeros(nsub/nagent, tmax, nmodel);
    newcorrect = zeros(nsub/nagent, tmax, nmodel);
    newqvalues = zeros(nsub/nagent, tmax, nmodel);
    newcho = zeros(nsub/nagent, tmax, nmodel);
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
                newcho(poolsub(1):poolsub(2), t, model) = ...
                    tempnewcho(poolsub(1):poolsub(2), t, model);
            end
        end
       poolsub = poolsub + nagent;
    end
    % save
    data = containers.Map(...
        {'newp', 'newcorrect', 'newqvalues', 'newcho'},...
        {newp, newcorrect, newqvalues, newcho});
    save('data/analyse/risk_online_cond', 'data');
end

%% ------------------------ Plots ----------------------------- %
% -----------------------------------------
% probabilities
% -----------------------------------------
f = figure('Renderer', 'painters', 'Position', [10 10 1500 1200]);

bounds = [1:20; 21:40; 41:60 ; 61:80; 81:100; 101:120; 121:140; 141:160; 161:180];
co = {
        [0.8500    0.3250    0.0980],...
        [107/255 196/255 103/255], ...
        [0.9290    0.6940    0.1250],...
        [0.4940    0.1840    0.5560],...
        [0    0.4470    0.7410],...
        [0.4660    0.6740    0.1880],...
        [0.3010    0.7450    0.9330],...
        [0.6350    0.0780    0.1840]...
    };

% 1: P(r=0|a=1)=1 ; P(r=1|a=2)=0.5, P(r=-1|a=2)=0.5
% 2: P(r=0|a=1)=1 ; P(r=1|a=2)=0.5, P(r=-1|a=2)=0.5', ...
% 3: P(r=1|a=1)=0.75, P(r=-1|a=1)=0.25 ; P(r=1|a=2)=0.25, P(r=-1|a=2)=0.75',...
% 4: P(r=0|a=1)=0.5, P(r=1|a=1)=0.5 ; P(r=0|a=2)=0.5, P(r=-1|a=2)=0.5'...
imodel = 0;
conds = [9 8 7 3 2 1 6 5 4];
j = 0;
for i = conds
    j = j + 1;
    subplot(3, 3, j)
    for model = 1:nmodel
        [L(model), sem] = reversalplot(...
            (newcho(:, bounds(i, :), model) == 2)',... %data
            [],... %time when reversal occurs
            [],... %time when cond changes
            ones(3) * 0.5,... % chance lvl
            co{model},... %curve color
            0.9,... %linewidth
            0.3,... % alpha
            0, 1,... % ylims
            15,... %fontsize,
            '',... %title,
            'trials',... %xlabel
            'P(risky)'... % ylabel
        );     
        hold on
    end
    %title(titles(i));
    %l = legend(L, models{1:nmodel});
    %l.FontSize = 40;
end