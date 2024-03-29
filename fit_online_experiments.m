% This function find the best fitting model/parameters

% 1: basic df=2
% 2: asymmetric neutral df=3
% 3: asymmetric pessimistic df=3
% 4: priors df=3
% 5: impulsive perseveration df=3
% 6: gradual perseveration df=3
% 7: full df=5
% 8: Bayesian df=3

close all
clear all

addpath './simulation'
addpath './fit'
addpath './utils'

%% Load experiment data
% --------------------------------------------------------------------
[cho, out, con] = getdata('data/data_online_exp/learningdata');


%% Modifiable variables
% --------------------------------------------------------------------
whichmodel = [1, 2, 5, 6, 7, 8];
% --------------------------------------------------------------------

%% Run
% --------------------------------------------------------------------
nmodel = length(whichmodel);
subjecttot = length(cho);

models = {'RW', 'RW\pm',...
    'RW\pm_{\omega^-}', 'RW_\omega', 'RW_\phi',...
    'RW_{\tau}', 'Full', 'Bayesian'};

paramlabels = {
    '\beta', '\alpha+', '\alpha-', '\alpha_{rel}', '\omega', '\phi',...
    '\tau', '\sigma_{xi}', '\sigma_{\epsilon}'}; 
nparam = length(paramlabels);

try
    data = load('data/fit_exp/online_exp');
    lpp = data.data('lpp');
    parameters = data.data('parameters');  %% Optimization parameters 
    ll = data.data('ll');
    hessian = data.data('hessian');
catch
    runfit(subjecttot, nparam, nmodel, whichmodel, con, cho, out);
end

%parametersfitbarplot(parameters, nmodel, whichmodel, models) 

%% Plots
% % --------------------------------------------------------------------
figure
params = {[2], [2, 3], [2, 3], [2, 4], [2, 5], [2, 5, 6], [2, 3, 5], [2, 3, 5, 6]};
alternatives = whichmodel(2:end);
j = 0;
for i = alternatives
    j = j +1; 
    subplot(2, 4, j)
    barplot_param_comparison(parameters, params{i}, i, paramlabels, models, 0.5);
end

% -------------- Params ----------------------------------------------- % 
bias1 = (parameters(:, 2, 2) - parameters(:, 3, 2))./...
    (parameters(:, 2, 2) + parameters(:, 3, 2)) ;
%bias2 = (parameters(:, 2, 3) - parameters(:, 3, 3)) ./...
%    (parameters(:, 2, 3) + parameters(:, 3, 3));
perse = parameters(:, 6, 5);
tau = (parameters(:, 7, 6) .* parameters(:, 6, 6)) ./...
(parameters(:, 7, 6) ./ parameters(:, 6, 6));
prior = parameters(:, 5, 4);
color = [107/255 196/255 103/255];

subplot(2, 3, 4)
scatterCorr(bias1, perse, color, 0.5, 2, 1);
xlabel('\alpha+ > \alpha-', 'FontSize', 20);
ylabel('\phi', 'FontSize', 20);
title('')
set(gca, 'Fontsize', 30);
subplot(2, 3, 5)
scatterCorr(bias1, prior, color, 0.5, 2, 1);
xlabel('\alpha+ > \alpha-', 'FontSize', 20);
ylabel('\omega', 'FontSize', 20);
title('')
set(gca, 'Fontsize', 30);
subplot(2, 3, 6)
scatterCorr(bias1, tau, color, 0.5, 2, 1);
xlabel( '\alpha+ > \alpha-', 'FontSize', 20);
ylabel('\phi * \tau', 'FontSize', 20);
title('')
set(gca, 'Fontsize', 30);
%figure
%violinplot([parameters(:, 2, 2), parameters(:, 2, 3)],...
%    {'\alpha+', '\alpha-'});


% figure
% llasy = abs(ll(:, 2) - ll(:, 1)); %./ (ll(:, 2) + ll(:, 1)); 
% llper = abs(ll(:, 5) - ll(:, 1)); %./ (ll(:, 5) + ll(:, 1)); 
% scatterCorr(llasy, llper, color, 0.5, 2, 1);

%error('ssssssssss');

% --------------------------------------------------------------------

i = 0;
nfpm = [3, 4, 4, 4, 4, 5, 7, 4];

for n = whichmodel
    i = i + 1;
   % bic(1:85, i) = -2 * -ll(:, n) + nfpm(n) * log(96);
    bic(:, i) = -2 * -ll(:, n) + nfpm(n) * log(180);
    aic(:, i)= -2 * -ll(:, n)...
            + 2*nfpm(n);
    me(:, i) = -lpp(:, n) + (nfpm(n)/2)*log(2*pi) - (1/2)*log(...
       arrayfun(@(x) det(cell2mat(x)), {hessian{:, n}})');
    %
end
% --------------------------------------------------------------------
% figure
% bar(mean(aic, 1));
% ylabel('AIC');
%VBA_groupBMC(-aic');
VBA_groupBMC(-bic');

%VBA_groupBMC(-bic'./2);



%% Functions
% --------------------------------------------------------------------

function runfit(subjecttot, nparam, nmodel, whichmodel, con, cho, out)
    parameters = zeros(subjecttot, nparam, nmodel);
    ll = zeros(subjecttot, nmodel);
    lpp = zeros(subjecttot, nmodel);
    report = zeros(subjecttot, nmodel);
    gradient = cell(subjecttot, nmodel);
    hessian = cell(subjecttot, nmodel);
    options = optimset(...
        'Algorithm',...
        'interior-point',...
        'Display', 'off',...
        'MaxIter', 10000,...
        'MaxFunEval', 10000);

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
                @(x) getlpp(x, con{nsub}, cho{nsub}, out{nsub}, model),...
                [1, .5, .5, .5, 0, 0, .5, .15, .15],...
                [], [], [], [],...
                [0, 0, 0, 0, -2, -2, 0, 0, 0],...
                [Inf, 1, 1, 1, 2, 2, 1, 1, 1],...
                [],...
                options...
                );
            parameters(nsub, :, model) = p;
            lpp(nsub, model) = l;
            report(nsub, model) = rep;
            gradient{nsub, model} = grad;
            hessian{nsub, model}= hess;
            
            [
                p1,...
                l1,...
                rep1,...
                grad1,...
                hess1,...
                ] = fmincon(...
                @(x) getll(x, con{nsub}, cho{nsub}, out{nsub}, model),...
                [1, .5, .5, .5, 0, 0, .5, .15, .15],...
                [], [], [], [],...
                [0, 0, 0, 0, -2, -2, 0, 0, 0],...
                [Inf, 1, 1, 1, 2, 2, 1, 1, 1],...
                [],...
                options...
                );
            ll(nsub, model) = l1;
        end
    end
    %%Save the data
    data = containers.Map({'parameters', 'lpp' 'll', 'hessian'}, {parameters, lpp, ll, hessian});
    save('data/fit_exp/online_exp', 'data');
    close(w);
   
end
% --------------------------------------------------------------------
function parametersfitbarplot(parameters, nmodel, whichmodel, models) 
    % %% to correct
    nfpm = [3, 4, 4, 4, 4, 5, 7, 4];
    parameters(:, 1, :) = 1/parameters(:, 1, :);
    y = zeros(nmodel, max(nfpm(whichmodel)));
    i = 0;
    for model = whichmodel
        i = i + 1;
        switch model
            case 1
                y(i, 1:nfpm(model), 1) = mean(parameters(:, [1, 2, 4], model), 1);
                s(i, 1:nfpm(model), 1) = sem(parameters(:, [1, 2, 4], model));
            case {2, 3, 7}
                y(i, 1:nfpm(model), 1) = mean(parameters(:, 1:nfpm(model), model), 1);
                s(i, 1:nfpm(model), 1) = sem(parameters(:, 1:nfpm(model), model));
            case 4
                y(i, 1:nfpm(model), 1) = mean(parameters(:, [1, 2, 4, 5], model), 1);
                s(model, 1:nfpm(model), 1) = sem(parameters(:, [1, 2, 4, 5], model));
            case 5
                y(i, 1:nfpm(model), 1) = mean(parameters(:, [1, 2, 4, 6], model), 1);
                s(i, 1:nfpm(model), 1) = sem(parameters(:, [1, 2, 4, 6], model));
            case 6
                y(i, 1:nfpm(model), 1) = mean(parameters(:, [1, 2, 4, 6, 7], model), 1);
                s(i, 1:nfpm(model), 1) = sem(parameters(:, [1, 2, 4, 6, 7], model));
            case 8 
                y(i, 1:nfpm(model), 1) = mean(parameters(:, [1, 8, 9], model), 1);
                s(i, 1:nfpm(model), 1) = sem(parameters(:, [1, 8, 9], model));
        end
    end
    figure
    hBar = bar(y, 'FaceAlpha', 0.7);
    set(gca, 'XTickLabels', {models{whichmodel}});
    hold on

    % Finding the number of groups and the number of bars in each group
    ngroups = size(y, 2);
    nbars = size(y, 1);
    % Calculating the width for each bar group
    groupwidth = min(0.8, nbars/(nbars + 1.5));
    % Set the position of each error bar in the centre of the main bar
    % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
    
    for i = 1:nbars
        % Calculate center of each bar
        x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
        errorbar(x, y(:, i, 1), s(:, i, 1), 'k', 'linestyle', 'none', 'linewidth', 1);
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
        set(gca, 'FontSize', 25);
    elseif length(param_idx) == 4
        b.CData(3, :) = [149/255 240/255 146/255];
        b.CData(4, :) = [60/255 240/255 146/255];
        set(gca, 'XTickLabel',{labels{param_idx(1)},...
            labels{param_idx(2)}, labels{param_idx(3)}, labels{param_idx(4)}});
        set(gca, 'FontSize', 25);
    else
        set(gca, 'XTickLabel',{labels{param_idx(1)}, labels{param_idx(2)}});
        set(gca, 'FontSize', 25);
    end

    title(models{model_idx});
    %ymin = ylim;
    %disp(ymin(1));
    %ylim([ymin(1), 1]); 
    box off
    ax1 = gca;
    hold(ax1, 'all');
    set(ax1, 'box', 'off');

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
        %set(ax, 'bof', 'off');
        %set(ax2, 'ytick', []);
        box off
%         ax1 = gca;
%         %set(ax1, 'XColor', 'r', 'YColor', 'r'); 
%         hold(ax1, 'all');  %   <--------------------------------
%         ax2 = axes('Position',get(ax1,'Position'),'XAxisLocation','top',...
%         'YAxisLocation','right','Color','none','XColor','k','YColor','k');
%         hold(ax2, 'all');  %   <--------------------------------
%         X = ones(1, nsub)-Shuffle(linspace(-0.15, 0.15, nsub));
%         scatter(...
%             X,...
%             y1,...
%             'filled', 'Parent', ax1, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1,...
%             'MarkerFaceColor', [107/255 220/255 103/255], 'MarkerEdgeColor', [107/255 220/255 103/255]);
%         s = scatter(...
%             X + 1,...
%             y2,...
%             'filled', 'Parent', ax1, 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1,...
%             'MarkerFaceColor', [107/255 214/255 103/255], 'MarkerEdgeColor', [107/255 214/255 103/255]);
%         set(gca, 'xtick', []);
%         set(gca, 'box', 'off');
%         set(ax1, 'box', 'off');
%         set(ax2, 'ytick', []);
%         box off
    end
    uistack(e, 'top');
end

function [cho, out, con] = getdata(file)
    data = load(file);
    data = data.learningdata(:, 1:18);
    ncond = max(data(:, 13));
    nsession = max(data(:, 18)) -1 ;
    sub_ids = unique(data(:, 2));
    i = 1;
    k = 1;
    tmaxsession = 60;
    for id = 1:length(sub_ids)
        sub = sub_ids(id);
        mask_sub = data(:, 2) == sub;
        if ismember(sum(mask_sub), [213, 228])
            t = 1;
            for sess = 0:nsession
                mask_sess = data(:, 18) == sess;
                mask = logical(mask_sub .* mask_sess);
                [noneed, trialorder] = sort(data(mask, 12));

                tempcho = data(mask, 9); 
                tempcho = tempcho(trialorder);
                tempout = data(mask, 7); 
                tempout = tempout(trialorder);
                tempcon = data(mask, 13);
                tempcon = tempcon(trialorder);

                for j = 1:tmaxsession

                    cho{i}(t) = tempcho(j);
                    out{i}(t) = tempout(j);
                    con{i}(t) = tempcon(j) + 1;

                    t = t + 1;
                end
            end
            if length(cho{i}) ~= 180
                error('No good length');
            end
%             if (sum(out{i}(:)) < 4)
%                 best_ids(k) = i;
%                 k = k + 1;
%             end
            i = i + 1;
        end 
    end
end



