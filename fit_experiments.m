% This function find the best fitting model/parameters_lpp
close all
clear all

addpath './simulation'
addpath './fit'
addpath './utils'

%% Load experiment data
% --------------------------------------------------------------------
[con, con2, cho, out, nsubs] = load_data('exp');


%% Define variables
% --------------------------------------------------------------------
whichmodel = [1, 2, 5, 6, 7, 8];
nmodel = length(whichmodel);
subjecttot = nsubs;
% --------------------------------------------------------------------

% 1: basic df=2
% 2: asymmetric neutral df=3
% 3: asymmetric pessimistic df=3
% 4: priors df=3
% 5: impulsive perseveration df=3
% 6: gradual perseveration df=4
% 7: Semi-full df=4
% 8: full df=5
% 9: Bayesian df=3

models = {'RW', 'RW\pm',...
    'RW\pm_{\omega^-}', 'RW_\omega', 'RW_\phi',...
    'RW_\tau', 'Semi-Full', 'Full', 'Bayesian'};

paramlabels = {
    '\beta', '\alpha+', '\alpha-', '\omega', '\phi',...
    '\tau', '\sigma_{xi}', '\sigma_{\epsilon}'}; 
nparam = length(paramlabels);

try
    data = load('data/fit_exp/data_usingfit');
    lpp = data.data('lpp');
    parameters_lpp = data.data('parameters_lpp');  %% Optimization parameters_lpp 
    parameters_ll = data.data('parameters_ll');
    ll = data.data('ll');
    hessian_ll = data.data('hessian_ll');
    hessian_lpp = data.data('hessian_lpp');
    whichmodel = data.data('whichmodel');
catch
    parameters_lpp = zeros(subjecttot, nparam, nmodel);
    parameters_ll = zeros(subjecttot, nparam, nmodel);
    
    ll = zeros(subjecttot, nmodel);
    lpp = zeros(subjecttot, nmodel);
    report = zeros(subjecttot, nmodel);
    gradient = cell(subjecttot, nmodel);
    hessian_lpp = cell(subjecttot, nmodel);
    hessian_ll = cell(subjecttot, nmodel);
    options = optimset(...
        'Algorithm',...
        'interior-point',...
        'Display', 'off',...
        'MaxIter', 5000,...
        'MaxFunEval', 5000);

    w = waitbar(0, 'Fitting subject');
    
    for nsub = 1:subjecttot
        
        waitbar(...
            nsub/subjecttot,...  % Compute progression
            w,...
            sprintf('%s%d', 'Fitting subject ', nsub)...
            );
        
        for model = whichmodel
            [
                p,...
                l,...
                rep,...
                output,...
                lmbda,...
                grad,...
                hess,...
                ] = fmincon(...
                @(x) getlpp_usingfit(x, con{nsub}, cho{nsub}, out{nsub}, model),...
                [1, .5, .5, 0, 0, .5, .15, .15],...
                [], [], [], [],...
                [0, 0, 0, -1, -2, 0, 0, 0],...
                [Inf, 1, 1, 1, 2, 1, 1, 1],...
                [],...
                options...
                );
            parameters_lpp(nsub, :, model) = p;
            lpp(nsub, model) = l;
            report(nsub, model) = rep;
            gradient{nsub, model} = grad;
            hessian_lpp{nsub, model}= hess;
            
            [
                p1,...
                l1,...
                rep1,...
                output,...
                lmbda,...
                grad1,...
                hess1,...
                ] = fmincon(...
                @(x) getll(x, con{nsub}, cho{nsub}, out{nsub}, model),...
                [1, .5, .5, 0, 0, .5, .15, .15],...
                [], [], [], [],...
                [0, 0, 0, -1, -2, 0, 0, 0],...
                [Inf, 1, 1, 1, 2, 1, 1, 1],...
                [],...
                options...
                );
            parameters_ll(nsub, :, model) = p1;
            ll(nsub, model) = l1;
            hessian_ll{nsub, model} = hess1;
        end
    end
    %%Save the data
    data = containers.Map({'whichmodel', 'parameters_ll', 'll', 'hessian_ll',...
        'parameters_lpp', 'lpp', 'hessian_lpp'},...
        {whichmodel, parameters_ll, ll, hessian_ll, parameters_lpp, lpp, hessian_lpp});
    save('data/fit_exp/data', 'data_usingfit');
    close(w);
    
end

%parameters_lpp = parameters_lpp(best_ids, :, :);
figure
subplot(1, 2, 1)
barplot_param_comparison(parameters_lpp, [2, 3, 5], 7, paramlabels, models, 0.5);
subplot(1, 2, 2)
barplot_param_comparison(parameters_lpp, [2, 3, 5, 6], 8, paramlabels, models, 0.5);

figure
subplot(1, 2, 1)
barplot_param_comparison(parameters_ll, [2, 3, 5], 7, paramlabels, models, 0.5);
subplot(1, 2, 2)
barplot_param_comparison(parameters_ll, [2, 3, 5, 6], 8, paramlabels, models, 0.5);

% Correlations
% ----------------------------------------------------------------------------
% bias1 = (parameters_lpp(:, 2, 2) - parameters_lpp(:, 3, 2))./...
%     (parameters_lpp(:, 2, 2) + parameters_lpp(:, 3, 2)) ;
% bias2 = (parameters_lpp(:, 2, 3) - parameters_lpp(:, 3, 3)) ./...
%     (parameters_lpp(:, 2, 3) + parameters_lpp(:, 3, 3));
% perse = parameters_lpp(:, 5, 5);
% prior = parameters_lpp(:, 4, 4);
% color = [107/255 196/255 103/255];
% color = [107/255 196/255 103/255];
% 
% subplot(2, 3, 4)
% scatterCorr(bias1, perse, color, 0.5, 2, 1);
% xlabel('\alpha+ > \alpha-', 'FontSize', 20);
% ylabel('\phi', 'FontSize', 20);
% title('')
% set(gca, 'Fontsize', 30);
% subplot(2, 3, 5)
% scatterCorr(bias1, prior, color, 0.5, 2, 1);
% xlabel('\alpha+ > \alpha-', 'FontSize', 20);
% ylabel('\omega', 'FontSize', 20);
% title('')
% set(gca, 'Fontsize', 30);
% subplot(2, 3, 6)
% scatterCorr(bias1, bias2, color, 0.5, 2, 1);
% xlabel( '\alpha+ > \alpha-', 'FontSize', 20);
% ylabel('\alpha+ > \alpha- ; \omega = -1', 'FontSize', 20);
% title('')
% set(gca, 'Fontsize', 30);
%figure
%violinplot([parameters_lpp(:, 2, 2), parameters_lpp(:, 2, 3)],...
%    {'\alpha+', '\alpha-'});


% figure
% llasy = abs(ll(:, 2) - ll(:, 1)); %./ (ll(:, 2) + ll(:, 1)); 
% llper = abs(ll(:, 5) - ll(:, 1)); %./ (ll(:, 5) + ll(:, 1)); 
% scatterCorr(llasy, llper, color, 0.5, 2, 1);

%error('ssssssssss');

% --------------------------------------------------------------------

% 1: basic df=2
% 2: asymmetric neutral df=3
% 3: asymmetric pessimistic df=3
% 4: priors df=3
% 5: impulsive perseveration df=3
% 6: gradual perseveration df=4
% 7: Semi-full df=4
% 8: full df=5
% 9: Bayesian df=3
i = 0;
nfpm = [2, 3, 3, 3, 3, 4, 4, 5, 3];

for n = whichmodel
    i = i + 1;
    bic(1:85, i) = -2 * -ll(1:85, n) + nfpm(n) * log(96);
    bic(86:105, i) = -2 * -ll(86:105, n) + nfpm(n) * log(96*2);
    aic(:, i)= -2 * -ll(:, n) + 2*nfpm(n);
    me(:, i) = -lpp(:, n) + (nfpm(n)/2)*log(2*pi) - (1/2)*log(...hess
       arrayfun(@(x) det(cell2mat(x)), {hessian_lpp{:, n}})');
    %
end

VBA_groupBMC(-me');
%VBA_groupBMC(-aic');
VBA_groupBMC(-bic');


% --------------------------------------------------------------------
function parametersfitbarplot(parameters_lpp, nmodel, whichmodel, models) 
    % %% to correct
    nfpm = [2, 3, 3, 3, 3, 4, 6, 3];
    parameters_lpp(:, 1, :) = 1/parameters_lpp(:, 1, :);
    y = zeros(nmodel, 5);
    for model = whichmodel
        switch model
            case {1, 2, 3, 7}
                y(model, 1:nfpm(model), 1) = mean(parameters_lpp(:, 1:nfpm(model), model), 1);
                s(model, 1:nfpm(model), 1) = sem(parameters_lpp(:, 1:nfpm(model), model));
            case 4
                y(model, 1:nfpm(model), 1) = mean(parameters_lpp(:, [1, 2, 4], model), 1);
                s(model, 1:nfpm(model), 1) = sem(parameters_lpp(:, [1, 2, 4], model));
            case 5
                y(model, 1:nfpm(model), 1) = mean(parameters_lpp(:, [1, 2, 5], model), 1);
                s(model, 1:nfpm(model), 1) = sem(parameters_lpp(:, [1, 2, 5], model));
            case 6
                y(model, 1:nfpm(model), 1) = mean(parameters_lpp(:, [1, 2, 5, 6], model), 1);
                s(model, 1:nfpm(model), 1) = sem(parameters_lpp(:, [1, 2, 5, 6], model));
            case 8 
                y(model, 1:nfpm(model), 1) = mean(parameters_lpp(:, [1, 7, 8], model), 1);
                s(model, 1:nfpm(model), 1) = sem(parameters_lpp(:, [1, 7, 8], model));
        end
    end
    figure
    hBar = bar(y, 'FaceAlpha', 0.7);
    set(gca, 'XTickLabels', models);
    hold on

    % Finding the number of groups and the number of bars in each group
    ngroups = size(y, 1);
    nbars = size(y, 1);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    try
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, y(:, i), s(:, i), 'k', 'linestyle', 'none', 'linewidth', 1);
    end
    catch
    end
    box off
end


function matrix = computeposterior(criterion, nmodel, models, whichmodel)
    %set options
    options.modelNames = {models{whichmodel}};
    options.DisplayWin = true;

    [posterior, outcome] = VBA_groupBMC(-criterion'./2, options);
    for fittedmodel = 1:nmodel
        matrix(fittedmodel, 1, 1) = mean(posterior.r(fittedmodel, :));
    end
end

function barplot_param_comparison(parameters, param_idx, model_idx, labels, models, ymax)
    %hold on
    nsub = size(parameters, 1);
    for i = 1:length(param_idx)
        y(i, :) = reshape(parameters(:, param_idx(i), model_idx), [], 1);
        means(i) = mean(y(i, :));
        errors(i) = sem(y(i, :));
    end
    b = bar(means, 'EdgeColor', 'black');
    box off
    hold on
    box off
    e = errorbar(means, errors, 'Color', 'black', 'LineWidth', 3, 'LineStyle', 'none');
    %hold off
    box off
    b.FaceColor = 'flat';
    b.CData(1, :) = [107/255 196/255 103/255];
    b.CData(2, :) = [149/255 230/255 146/255];
    if length(param_idx) == 3
        b.CData(3, :) = [149/255 240/255 146/255];
        set(gca, 'XTickLabel',{labels{param_idx(1)}, labels{param_idx(2)}, labels{param_idx(3)}});
        set(gca, 'FontSize', 30);
    elseif length(param_idx) == 4
        b.CData(3, :) = [149/255 240/255 146/255];
        b.CData(4, :) = [60/255 240/255 146/255];
        set(gca, 'XTickLabel',{labels{param_idx(1)},...
            labels{param_idx(2)}, labels{param_idx(3)}, labels{param_idx(4)}});
        set(gca, 'FontSize', 30);
    else
        set(gca, 'XTickLabel',{labels{param_idx(1)}, labels{param_idx(2)}});
        set(gca, 'FontSize', 30);
    end

    title(models{model_idx});
    %ymin = ylim;
    %disp(ymin(1));
    %ylim([ymin(1), 1]); 
    box off
    ax1 = gca;
    hold(ax1, 'all');

    for i = 1:length(param_idx)      
        ax(i) = axes('Position',get(ax1,'Position'),'XAxisLocation','top',...
         'YAxisLocation','right','Color','none','XColor','k','YColor','k');
          
        hold(ax(i), 'all');
        
        X = ones(1, nsub)-Shuffle(linspace(-0.15, 0.15, nsub));
        scatter(...
            X + (i-1),...
            y(i, :),...
             'filled', 'Parent', ax1, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1,...
             'MarkerFaceColor', [107/255 220/255 103/255],...
             'MarkerEdgeColor', [107/255 220/255 103/255]);
        set(gca, 'xtick', []);
        set(gca, 'box', 'off');
        set(ax(i), 'box', 'off');
        box off
    end
    uistack(e, 'top');
end

