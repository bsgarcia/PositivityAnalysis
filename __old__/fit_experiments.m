% This function find the best fitting model/parameters
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
whichmodel = 1:6;
nmodel = length(whichmodel);
models = {'QLearning', 'Asymmetric',...
    'AsymmetricPessimistic', 'Perseveration'...
     'Priors', 'Full'};

try
    data = load('data/fit_exp/data');
    lpp = data.data('lpp');
    parameters = data.data('parameters');
catch 
    %% Optimization
    subjecttot = nsubs;
    nparam = 5;
    parameters = zeros(subjecttot, nparam, nmodel);
    lpp = zeros(subjecttot, nmodel);
    report = zeros(subjecttot, nmodel);
    gradient = cell(subjecttot, nmodel);
    hessian = cell(subjecttot, nmodel);


    % 1: basic df=2
    % 2: asymmetric neutral df=3
    % 3: asymmetric pessimistic df=3
    % 4: perseveration df=3
    % 5: priors df=3
    % 6: full df=5

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
                grad,...
                hess,...
             ] = fmincon(...
                    @(x) getlpp(x, con{nsub}, cho{nsub}, out{nsub}, model),...
                    [1, .5, .5, 0, 0, .5],...
                    [], [], [], [],...
                    [0, 0, 0, -2, -1, 0],...
                    [Inf, 1, 1, 2, 1, 1],...
                    [],...
                    options...
                );
            parameters(nsub, :, model) = p;
            lpp(nsub, model) = l;
            report(nsub, model) = rep;
            gradient{nsub, model} = grad;
            hessian{nsub, model}= hess;

        end
    end
    data = containers.Map({'parameters', 'lpp'}, {parameters, lpp});
    save('data/fit_exp/data', 'data');
    close(w);
end
% nfpm=[2 3]; % number of free parameters

%% to correct
nfpm = [2, 3, 3, 3, 3, 5];
i = 0;
for n = whichmodel
    i = i + 1;
    bic(1:85, i) = -2 * -lpp(1:85, n) + nfpm(n) * log(96);
    bic(86:105, i) = -2 * -lpp(86:105, n) + nfpm(n) * log(96*2);
    % l2 is already positive
    % aic(:,n)=lpp(:,n)+2*nfpm(n);
    %
end

% results_bic = computeposterior(bic', nmodel, models, whichmodel);
%results_aic = computeposterior(aic, nmodel);
%
%BarsAndErrorPlot(results_bic, 'k', 'k', 1, 1, 0, 1, 1, 14, '', '', '');
% 1/beta
parameters(:, 1, :) = 1/parameters(:, 1, :);
y = zeros(nmodel, 5);
for model = whichmodel
    switch model
        case {1, 2, 6}
            y(model, 1:nfpm(model), 1) = mean(parameters(:, 1:nfpm(model), model), 1);
            s(model, 1:nfpm(model), 1) = sem(parameters(:, 1:nfpm(model), model));
        case 3
            y(model, 1:nfpm(model), 1) = mean(parameters(:, 1:nfpm(model), model), 1);
            s(model, 1:nfpm(model), 1) = sem(parameters(:, 1:nfpm(model), model));

            %y(model, nfpm+1, 1) = -1;
        case 4
            y(model, 1:nfpm(model), 1) = mean(parameters(:, [1, 2, 4], model), 1);
            s(model, 1:nfpm(model), 1) = sem(parameters(:, [1, 2, 4], model));

        case 5
            y(model, 1:nfpm(model), 1) = mean(parameters(:, [1, 2, 5], model), 1);
            s(model, 1:nfpm(model), 1) = sem(parameters(:, [1, 2, 5], model));

    end
end
% x = 1:size(y, 1) * size(y, 2);
% x = reshape(x, size(y));
figure
hBar = bar(y, 'FaceAlpha', 0.7);
set(gca, 'XTickLabels', models);
hold on

% Finding the number of groups and the number of bars in each group
ngroups = size(y, 1);
nbars = size(y, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, y(:, i), s(:, i), 'k', 'linestyle', 'none', 'linewidth', 1);
end
box off
% for i = 1:6
%     errorbar([0.1, 0.2, 0.3, 0.4, 0.5], y(i, :), [0.5, 0.5, 0.5, 0.5, 0.5], 'LineStyle', 'None', 'Color', 'k');
% end
%er.Color = [0 0 0];                            
%er.LineStyle = 'none';  
%for k1 = 1:size(y, 1)-1
%    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
%    ydt(k1,:) = hBar(k1).YData;
%end
%hold on
%errorbar(ctr, ydt, s, '.r')
%hold off

figure;
%% What is it?
bias1 = (parameters(:, 2, 2) - parameters(:, 3, 2)) ./ (parameters(:, 2, 2) ...
    + parameters(:, 3, 2));
bias2 = (parameters(:, 2, 3) - parameters(:, 3, 3)) ./ (parameters(:, 2, 3) ...
    + parameters(:, 3, 3));
perse = parameters(:, 4, 4);
prior = parameters(:, 5, 5);

subplot(1, 2, 1)
scatterCorr(bias1, perse, 2, 1);
subplot(1, 2, 2)
scatterCorr(bias1, prior, 2, 1);

function matrix = computeposterior(criterion, nmodel, models, whichmodel)
    %set options
    options.modelNames = {models{whichmodel}};
    options.DisplayWin = true;

    [posterior, outcome] = VBA_groupBMC(-criterion, options);
    for fittedmodel = 1:nmodel
        matrix(fittedmodel, 1, 1) = mean(posterior.r(fittedmodel, :));
    end
end



% 
% 
% %% ----------------------- Plots ----------------------------------

% 
% %%
% aaaa = [
%     (parameters(:, 2, 6) - parameters(:, 3, 6)) ./ (parameters(:, 2, 6) +...
%     parameters(:, 3, 6)), parameters(:, 4:5, 6)
% ];
% 
% %
% figure;
% BarsAndErrorPlot(...
%     aaaa',...
%     [0, 0, 0],...
%     [0.5, 0.5, 0.5],...
%     2,...
%     0.5,...
%     -1,...
%     1,...
%     -0.2,...
%     14,...
%     '',...
%     '',...
%     ''...
% );
% 
% %%
% figure;
% 
% subplot(1, 2, 1)
% 
% scatterCorr(...
%     (parameters(:, 2, 6) - parameters(:, 3, 6))./(parameters(:, 2, 6)...
%     + parameters(:, 3, 6)), parameters(:, 4, 6), 1, 1);
% 
% subplot(1, 2, 2)
% 
% scatterCorr(...
%     (parameters(:, 2, 6) - parameters(:, 3, 6))./(parameters(:, 2, 6)...
%     + parameters(:, 3, 6)), parameters(:, 5, 6), 1, 1);
% 
% 
% %bestmodel=find(bicgroup==min(bicgroup));
% 
% %bestparameters(:,:)=parameters(:,:,bestmodel);
% 
% %% nex is the one to download anyxay
% 
% %parameters2(:,1,2)=1./parameters2(:,1,2);
% parameters(:, 1, 2) = 1 ./ parameters(:, 1, 2);
% 
% figure;
% BarsAndErrorPlotThree(...
%     parameters(:, :, 2)',...
%     [0, 0, 0],...
%     [0.5, 0.5, 0.5],...
%     2,...
%     0.5,...
%     0,...
%     1,...
%     -0.2,...
%     14,...
%     '', '', '', '', '', '');
% hold on
% %BarsAndErrorPlotThree(parameters2(:,:,2)',[0.5 0.5 0.5],[0 0 0],2,0.5,0,1,+0.2,14,'','','','','','');
% 
% %%
% %figure;
% %BarsAndErrorPlotTwo(ll',[0 0 0],[0.5 0.5 0.5],2,0.5,0,1,-0.2,14,'','','','','','');
% %hold on
% %BarsAndErrorPlotTwo(ll2',[0.5 0.5 0.5],[0 0 0],2,0.5,0,1,+0.2,14,'','','','','','');
% 
% [outcome, post] = VBA_groupBMC((-bic(:, [1, 2, 4, 5, 6]))');
% 
% %%
% 
% %figure;
% %scatterCorr(parameters(:,2,2)-parameters(:,3,2),parameters2(:,2,2)-parameters2(:,3,2),1,[0 0 1]);
