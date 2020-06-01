clear all
close all

addpath './simulation'
addpath './fit'
addpath './utils'

% Variables to modify
%% ---------------------------------------------------------------------- % 
whichmodel = [1, 2, 5, 6, 7, 8];
%% ---------------------------------------------------------------------- %
%% Set  variables
nmodel = length(whichmodel);
paramlabels = {
    '\beta', '\alpha+', '\alpha-', '\omega', '\phi',...
    '\tau', '\sigma_{xi}', '\sigma_{\epsilon}'}; 
nparam = length(paramlabels);
tmax = 4*48;
models = {'QLearning', 'Asymmetric', 'AsymmetricPessimistic', 'Priors', ...
    'Impulsive Perseveration', 'Gradual Perseveration', 'Semi-Full', 'Full', 'Kalman'};
nsim = 100;
% 1: basic df=2
% 2: asymmetric neutral df=3
% 3: asymmetric pessimistic df=3
% 4: priors df=3
% 5: impulsive perseveration df=3
% 6: gradual perseveration df=4
% 7: Semi-full df=4
% 8: full df=5
% 9: Bayesian df=3
options = optimset( ...
    'Algorithm', ...
    'interior-point', ...
    'Display', 'off', ...
    'MaxIter', 2000, ...
    'MaxFunEval', 2000);
%'UseParallel', true);

w = waitbar(0, 'Get data');

% try
%     [ll, parameters] = getdata('35');
% catch
%     [ll, parameters] = runfit('35', nmodel, whichmodel, nparam,...
%         options, w);
% end
conds = {'35_opt'};
for j = 1:length(conds) 
    name = conds{j};
    try
        data = getdata(name);
    catch
        data = runfit(name, nmodel, whichmodel, nparam,...
            options, w);
    end
        
    dd = data{1}('ll');
    subjecttot = length(dd(1, 1, :));

    loweraicmatrix = zeros(nmodel, nmodel);
    ppmatrix = zeros(nmodel, nmodel);
    for i = 1:nsim
        d = data{i}('ll');
        [bic, aic] = computebicaic(d, tmax, nmodel, subjecttot, whichmodel);
        loweraicmatrix = computeloweraic(loweraicmatrix, aic, whichmodel, nmodel, subjecttot);
        ppmatrix = computehigherpp(ppmatrix, aic, whichmodel, nmodel, subjecttot);
    end

    xlabels = {models{whichmodel}};
    ylabels = {models{whichmodel}};
    %loweraicmatrix = loweraicmatrix ./ subjecttot;
    %ppmatrix = ppmatrix ./ (subjecttot .* 100);

    plotheatmap(loweraicmatrix, xlabels, ylabels,...
        'Fitted model', 'Simulated data using model',...
        sprintf('AIC %s', name),...
        '% model wins');

    plotheatmap(ppmatrix, xlabels, ylabels,...
        'Fitted model', 'Simulated data using model',...
        sprintf('maxPP cond=%s', name),...
        '% mean of maxPP');
end

%% ---------------------- Plot functions -------------------------------- % 

function plotheatmap(data, xlabels, ylabels, xaxislabel, yaxislabel, titlelabel,...
    cname)
    f = figure('Renderer', 'painters', 'Position', [10 10 800 700]);
    h = heatmap(xlabels, ylabels, data);    
    title(titlelabel);
    ax = gca;
    axs = struct(ax);   %and ignore the warning
    try
        warning('off', 'last');
    catch     
    end
    xlabel(xaxislabel);   
    ylabel(yaxislabel);
    %then ax becomes the handle to the heatmap
    c = axs.Colorbar;    %now you have a handle to the colorbar object    
    c.Label.String = cname;
    %set(c.Label, 'Rotation', -360+-90);
    %saveas(f, sprintf('fig/2sessions_recovery_%s_%d.png', titlelabel, length(xlabels)));
end

%% ---------------------- computation functions -------------------------------- % 


function [bic, aic, me] = computebicaic(ll, tmax, nmodel, subjecttot, whichmodel)
    % 1: basic df=2
    % 2: asymmetric neutral df=3
    % 3: asymmetric pessimistic df=3
    % 4: priors df=3
    % 5: impulsive perseveration df=3
    % 6: gradual perseveration df=4
    % 7: Semi-full df=4
    % 8: full df=5
    % 9: Bayesian df=3
    nfpm = [2, 3, 3, 3, 3, 4, 4, 5, 3];

    bic = zeros(nmodel, nmodel, subjecttot);
    aic = zeros(nmodel, nmodel, subjecttot);

    for fittedmodel = whichmodel
        bic(fittedmodel, whichmodel, :) = -2 * -ll(fittedmodel, whichmodel, :)...
            + nfpm(fittedmodel) * log(tmax);
        
        aic(fittedmodel, whichmodel, :) = -2 * -ll(fittedmodel, whichmodel, :)...
            + 2*nfpm(fittedmodel);
       
    end
end

function loweraicmatrix = computeloweraic(loweraicmatrix, aic, whichmodel, nmodel, subjecttot)

    for datamodel = whichmodel
        for sub = 1:subjecttot
            [useless, argmin] = min(aic(:, datamodel, sub));
             loweraicmatrix(datamodel, argmin) = loweraicmatrix(datamodel, argmin) + 1;
        end
    end
    %loweraicmatrix = (loweraicmatrix ./ subjecttot);
end


function matrix = computehigherpp(matrix, criterion, whichmodel, nmodel, subjecttot)
    i = 0;
    for datamodel = whichmodel
        i = i + 1;
        %set options
%         %options.modelNames = models(whichmodel);
%         options.figName = sprintf(...
%             '%s condition, data generated using %s', cond, models{datamodel});
        options.DisplayWin = false;
        
        formatedmatrix(1:nmodel, 1:subjecttot) = -criterion(whichmodel, datamodel, :);
        
        [posterior, outcome] = VBA_groupBMC(formatedmatrix./2, options);
        for fittedmodel = 1:nmodel
                matrix(i, fittedmodel) = matrix(i, fittedmodel)...
                    + (argmax(posterior.r(fittedmodel, :)) == fittedmodel);
        end
        %matrix = matrix ./ subjecttot;
    end
end

function matrix = computeposterior(...
    criterion, whichmodel, nmodel, models, subjecttot, cond)
    i = 0;
    for datamodel = whichmodel
        i = i + 1;
        %set options
        %options.modelNames = models(whichmodel);
        options.figName = sprintf(...
            '%s condition, data generated using %s', cond, models{datamodel});
        options.DisplayWin = false;
        
        formatedmatrix(1:nmodel, 1:subjecttot) = -criterion(whichmodel, datamodel, :);
        
        [posterior, outcome] = VBA_groupBMC(formatedmatrix./2, options);
        for fittedmodel = 1:nmodel
            matrix(i, fittedmodel) = mean(posterior.r(fittedmodel, :));
        end
    end
end


function bicmatrix = computepercentagewinning(...
    bic, whichmodel, nmodel, subjecttot)

    bicmatrix = zeros(nmodel, nmodel);
    for datamodel = whichmodel       
        for sub = 1:subjecttot
            [useless, argmin] = min(bic(:, datamodel, sub));
             bicmatrix(datamodel, argmin) = bicmatrix(datamodel, argmin) + 1;
        end
    end
    bicmatrix = bicmatrix ./ (subjecttot/100);
end

%% --------------------------- Get fit data functions ------------ % 

function data = getdata(file)
        data = load(sprintf('data/fit_sim/%s', file));
        data = data.data;
end

function data = runfit(file, nmodel, whichmodel, nparam,...
    options, w)
    [
        con, ...
        con2,....
        cho, ...
        out, ...
        nsubs, ...
        nsim, ...
    ] = load_data('sim',  file);
    
    data = cell(1);

    for n = 1:nsim
        waitbar( ...
                n/nsim, ... % Compute progression
                w, ...
                sprintf('Fitting sim %d for cond %s ', n, file) ...
        );

        subjecttot = nsubs;

        for nsub = 1:subjecttot
            
            for fittedmodel = whichmodel
                    templl = [];
                    tempparam = [];
                    temphess = [];
                    for datamodel = whichmodel
                        [
                            p, ...
                            l, ...
                            rep, ...
                            output,...
                            lambda,...
                            grad, ...
                            hess, ...
                        ] = fmincon( ...
                                    @(x) ...
                                    getll( ...
                                        x, ...
                                        con{n}{nsub}(:, :, datamodel), ...
                                        cho{n}{nsub}(:, :, datamodel), ...
                                        out{n}{nsub}(:, :, datamodel), ...
                                        fittedmodel ...
                                    ), ...
                                [1, .5, .5, 0, 0, .5, .15, .15],...
                                [], [], [], [],...
                                [0, 0, 0, -1, -2, 0, 0, 0],...
                                [Inf, 1, 1, 1, 2, 1, 1, 1],...
                                [], ...
                                options ...
                            );
                        tempparam(datamodel, :) = p;
                        templl(datamodel) = l;
                        temphess{datamodel} = hess;

                    end
                    ll(fittedmodel, :, nsub) = templl;
                    parameters(fittedmodel, :, nsub, :) = tempparam;
                    hessian{fittedmodel, :, nsub} = temphess;
            end
        end
        data{n} = containers.Map({'parameters', 'll', 'hessian',},....
                {parameters, ll, hessian}...
        );
    end
    save(sprintf('data/fit_sim/%s', file), 'data');
end
