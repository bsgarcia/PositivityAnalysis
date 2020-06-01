close all
clear all

addpath './simulation'
addpath './fit'
addpath './utils'

%% ----------------------- Variables to modify -------------------------- % 
whichmodel = 1:6;
nmodel = length(whichmodel);
%% ---------------------- Define variables --------------------- ---------%

condlabels = {'risk', 'statusquo1', 'statusquo2'};
models = {'QLearning', 'Asymmetric', 'AsymmetricPessimistic',...
        'Perseveration', 'Priors', 'Full'};
    
t_reversal = containers.Map(...
    condlabels,...
    {[], [72, 112, 128, 156, 168, 180], [60, 120, 180]});
t_cond = containers.Map(...
    condlabels,...
    repelem({[0, 48, 2*48, 3*48, (4*48)+1]}, 3));


%% ----------------------Loop over conditions----------------------------% 

for i = 1:length(condlabels)
    [con, con2, cho, out, nsub, p, correct, qvalues] = load_data('sim', condlabels{i});
    tmax = length(p{1}(:, 1));
    %disp(cho{1}(3*48:tmax, 1));

    try
        %% ------------------- Try to load formatted data --------------- % 
        data = load(sprintf('data/analyse/%s', condlabels{i}));
        newdata = data.data;
        newp = newdata('newp');
        newcorrect = newdata('newcorrect');
        newqvalues = newdata('newqvalues');
        
        %disp(std(ndata(50, :, 1)));
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
        newp = zeros(nsub/50, tmax, nmodel);
        newcorrect = zeros(nsub/50, tmax, nmodel);
        newqvalues = zeros(nsub/50, tmax, nmodel);
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
                end
            end
           poolsub = poolsub + 50;
        end
        % save
        data = containers.Map(...
            {'newp', 'newcorrect', 'newqvalues'},...
            {newp, newcorrect, newqvalues});
        save(sprintf('data/analyse/%s', condlabels{i}), 'data');
    end
    
    
    %% ------------------------ Plots ----------------------------- %
    
    % probabilities
    % --------------------------------------------------------------
    f = figure('Renderer', 'painters', 'Position', [10 10 1600 1200]);
    suptitle(condlabels{i});
    imodel = 0;
    for model = whichmodel
        imodel = imodel + 1;
        subplot(nmodel/2, 2, imodel);
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
            models{model},... %title,
            'trials',... %xlabel
            'p of option 2'... % ylabel
        ); 
         
    end
    saveas(f, sprintf('analysis_proba_%s.png', condlabels{i}));

    % correct choices
    % --------------------------------------------------------------
    f = figure('Renderer', 'painters', 'Position', [10 10 1600 1200]);
    suptitle(condlabels{i});
    imodel = 0;
    for model = whichmodel
        imodel = imodel + 1;
        subplot(nmodel/2, 2, imodel);
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
            models{model},... %title,
            'trials',... %xlabel
            '% of choices maximizing utility'... % ylabel
        ); 
    end
    
    saveas(f, sprintf('analysis_choice_%s.png', condlabels{i}));

    
    % Q values
    % --------------------------------------------------------------
    f = figure('Renderer', 'painters', 'Position', [10 10 1600 1200]);
    suptitle(condlabels{i});
    imodel = 0;
    for model = whichmodel
        imodel = imodel + 1;
        subplot(nmodel/2, 2, imodel);
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
            models{model},... %title,
            'trials',... %xlabel
            'Q of option 2'... % ylabel
        );    
    end
    saveas(f, sprintf('analysis_qvalue_%s.png', condlabels{i}));

end        
