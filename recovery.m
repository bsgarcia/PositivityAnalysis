cclear all
close all

addpath './simulation'
addpath './fit'
addpath './utils'

% Variables to modify
%% ---------------------------------------------------------------------- % 
allmodel = 1:8;
whichmodel = [1, 2, 3, 4, 5];
%% ---------------------------------------------------------------------- %
%% Set  variables
nallmodel = length(allmodel);
nmodel = length(whichmodel);
nparam = 8;
tmax = 196;
models = {'QLearning', 'Asymmetric', 'AsymmetricPessimistic', 'Priors', ...
    'Impulsive Perseveration', 'Gradual Perseveration', 'Full'};

% 1: basic df=2
% 2: asymmetric neutral df=3
% 3: asymmetric pessimistic df=3
% 4: perseveration df=3
% 5: priors df=3
% 6: full df=5

options = optimset( ...
    'Algorithm', ...
    'interior-point', ...
    'Display', 'off', ...
    'MaxIter', 2000, ...
    'MaxFunEval', 2000);
%'UseParallel', true);

w = waitbar(0, 'Get data');

try
    [ll, parameters] = getdata('conf');
catch
    [ll, parameters] = runfit('conf', nmodel, allmodel, whichmodel, nparam,...
        options, w);
end

subjecttot = length(ll(1, 1, :));

[bic, aic] = computebicaic(ll, tmax, nmodel, subjecttot, whichmodel);

% Compute the posterior probabilities
bicmatrix = computeposterior(...
    bic, whichmodel, nmodel, models, subjecttot, 'online');

aicmatrix = computeposterior(...
    aic, whichmodel, nmodel, models, subjecttot, 'online');

xlabels = {models{whichmodel}};
ylabels = {models{whichmodel}};

% plot bic recovery for one condition
% plotheatmap(bicmatrix, xlabels, ylabels,...
%     'Fitted model', 'Simulated data using model x',...
%     sprintf('BIC'),...
%     'mean of posterior probabilities');

% plot aic recovery for one condition
plotheatmap(aicmatrix, xlabels, ylabels,...
    'Fitted model', 'Simulated data',...
    sprintf('AIC'),...
    'mean of posterior probabilities');
% plot aic recovery for one condition
plotheatmap(bicmatrix, xlabels, ylabels,...
    'Fitted model', 'Simulated data',...
    sprintf('BIC'),...
    'mean of posterior probabilities');
    
close(w);


%% ---------------------- Plot functions -------------------------------- % 

function plotheatmap(data, xlabels, ylabels, xaxislabel, yaxislabel, titlelabel,...
    cname)
    f = figure('Renderer', 'painters', 'Position', [10 10 800 700]);
    h = heatmap(xlabels, ylabels, data);
    %set(gca,'defaulttextinterpreter','latex');
    %set(gca, 'defaultAxesTickLabelInterpreter','latex');  
    %load('greenish');
    %colormap(flip(mymap));
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
    saveas(f, sprintf('fig/2sessions_recovery_%s_%d.png', titlelabel, length(xlabels)));
end

%% ---------------------- computation functions -------------------------------- % 


function [bic, aic] = computebicaic(ll, lpp, tmax, nmodel, subjecttot, whichmodel)
    % 1: basic df=2
    % 2: asymmetric neutral df=3
    % 3: asymmetric pessimistic df=3
    % 4: perseveration df=3
    % 5: priors df=3
    % 6: full df=5
    nfpm = [2, 3, 3, 3, 3, 4, 6, 3];

    bic = zeros(nmodel, nmodel, subjecttot);
    aic = zeros(nmodel, nmodel, subjecttot);
    %lme = zeros(nmodel, nmodel, subjecttot);

    for fittedmodel = whichmodel
        bic(fittedmodel, whichmodel, :) = -2 * -ll(fittedmodel, whichmodel, :)...
            + nfpm(fittedmodel) * log(tmax);
        
        aic(fittedmodel, whichmodel, :) = -2 * -ll(fittedmodel, whichmodel, :)...
            + 2*nfpm(fittedmodel);
        
        me(:, i) = -lpp(:, n) + (nfpm(n)/2)*log(2*pi) - (1/2)*log(...
            arrayfun(@(x) det(cell2mat(x)), {hessian{:, n}})');
       
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

function [ll, parameters] = getdata(file)
        data = load(sprintf('data/fit_sim/%s', file));
        ll = data.data('lpp');
        parameters = data.data('parameters');
end

function [ll, parameters] = runfit(file, nmodel, allmodel, whichmodel, nparam,...
    options, w)
    [
        con, ...
        con2,....
        cho, ...
        out, ...
        nsubs, ...
    ] = load_data('sim',  file);

    subjecttot = nsubs;
    parameters = repelem({zeros(subjecttot, nparam, nmodel)}, nmodel);
    ll = zeros(nmodel, nmodel, subjecttot);
    report = repelem({zeros(subjecttot, nmodel)}, nmodel);
    gradient = repelem({cell(subjecttot, nmodel)}, nmodel);
    hessian = repelem({cell(subjecttot, nmodel)}, nmodel);

    for nsub = 1:subjecttot
        if ~(mod(nsub, 5))
            waitbar( ...
                nsub/subjecttot, ... % Compute progression
                w, ...
                sprintf('Fitting subject %d for cond %s ', nsub, file) ...
            );
        end
        parfor fittedmodel = allmodel
            if ismember(fittedmodel, whichmodel)
                templl = [];
                for datamodel = whichmodel
                    [
                        p, ...
                        l, ...
                        rep, ...
                        grad, ...
                        hess, ...
                    ] = fmincon( ...
                                @(x) ...
                                getllnorel( ...
                                    x, ...
                                    con{nsub}(:, :, datamodel), ...
                                    cho{nsub}(:, :, datamodel), ...
                                    out{nsub}(:, :, datamodel), ...
                                    fittedmodel ...
                                ), ...
                            [1, .5, .5, 0, 0, .5, .15, .15],...
                            [], [], [], [],...
                            [0, 0, 0, -1, -2, 0, 0, 0],...
                            [Inf, 1, 1, 1, 2, 1, 1, 1],...
                            [], ...
                            options ...
                        );

                    parameters{fittedmodel}(nsub, :, datamodel) = p;
                    templl(datamodel) = l;
                    report{fittedmodel}(nsub, datamodel) = rep;
                    gradient{fittedmodel}{nsub, datamodel} = grad;
                    hessian{fittedmodel}{nsub, datamodel} = hess;

                end
                ll(fittedmodel, :, nsub) = templl;
            end
        end
    end
    data = containers.Map({'parameters', 'll'},....
            {parameters, ll}...
    );
    save(sprintf('data/fit_sim/%s', file), 'data');
end
