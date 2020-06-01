close all
clear all

addpath './simulation'
addpath './fit'
addpath './utils'

%% ----------------------- Variables to modify -------------------------- % 
whichmodel = 1:6;
nmodel = length(whichmodel);
nagent = 20;
%% ---------------------- Define variables --------------------- ---------%

condlabels = {'risk', 'statusquo1', 'statusquo2'};
models = {'QLearning', 'Asymmetric', 'AsymmetricPessimistic',...
        'Perseveration', 'Priors', 'Full'};
    
% when reversal occurs
t_reversal = containers.Map(...
    {'risk', 'statusquo1', 'statusquo2'},...
    {[], [36, 56, 64, 78, 84, 90, 132, 154, 162, 176, 182, 188],...
    [30, 60, 90, 126, 156, 186]});

t_cond = containers.Map(...
    condlabels,...
    repelem({[0, 24, 2*24, 3*24, 4*24, 5*24, 6*24, 7*24, (8*24+1)]}, 3));

%% ----------------------Loop over conditions----------------------------% 

for i = 1:length(condlabels)
    if strcmp(condlabels{i}, 'risk')
    [con, con2, cho, out, nsub, p, correct, qvalues] = load_data('sim2', sprintf('%s_new', condlabels{i}));
    tmax = length(p{1}(:, 1));
    %disp(cho{1}(3*48:tmax, 1));

    try
        %% ------------------- Try to load formatted data --------------- % 
        data = load(sprintf('data/analyse_2sessions/%s_new', condlabels{i}));
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
        for sub = i:nsub
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
        save(sprintf('data/analyse_2sessions/%s', condlabels{i}), 'data');
    end
    
    
    %% ------------------------ Plots ----------------------------- %
    
    % probabilities
    % --------------------------------------------------------------
    f = figure('Renderer', 'painters', 'Position', [10 10 1800 1000]);
    suptitle(condlabels{i});
    for model = 1:nmodel
        subplot(nmodel/2, 2, model);
        reversalplot(...
            newp(:, :, model)',... %data
            t_reversal(condlabels{i}),... %time when reversal occurs
            t_cond(condlabels{i}),... %time when cond changes
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
    saveas(f, sprintf('fig/2sessions_analysis_proba_%s.png', condlabels{i}));
  
    % correct choices
    % --------------------------------------------------------------
    f = figure('Renderer', 'painters', 'Position', [10 10 1800 1000]);
    suptitle(condlabels{i});
    for model = 1:nmodel
        subplot(nmodel/2, 2, model);
        reversalplot(...
            newcorrect(:, :, model)',... %data
            t_reversal(condlabels{i}),... %time when reversal occurs
            t_cond(condlabels{i}),... %time when cond changes
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
    saveas(f, sprintf('fig/2sessions_analysis_choice_%s.png', condlabels{i}));

    
    % Q values
    % --------------------------------------------------------------
    f = figure('Renderer', 'painters', 'Position', [10 10 1800 1000]);
    suptitle(condlabels{i});
    for model = 1:nmodel
        subplot(nmodel/2, 2, model);
        reversalplot(...
            newqvalues(:, :, model)',... %data
            t_reversal(condlabels{i}),... %time when reversal occurs
            t_cond(condlabels{i}),... %time when cond changes
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
    saveas(f, sprintf('fig/2sessions_analysis_qvalue_%s.png', condlabels{i}));
    end
    
end        
