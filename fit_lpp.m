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
% 6: gradual perseveration df=3
% 7: fulpp df=5
% 8: Bayesian df=3

models = {'RW', 'RW\pm',...
    'RW\pm_{\omega^-}', 'RW_\omega', 'RW_\phi',...
    'RW_\tau', 'Semi-Full', 'Full', 'Bayesian'};

paramlabels = {
    '\beta', '\alpha+', '\alpha-', '\omega', '\phi',...
    '\tau', '\sigma_{xi}', '\sigma_{\epsilon}'}; 
nparam = length(paramlabels);

try
    data = load('data/fit_exp/data_new_ms_lpp');
    parameters_lpp = data.data('parameters_lpp');
    lpp = data.data('lpp');
    hessian_lpp = data.data('hessian_lpp');
catch
    parameters_lpp = zeros(subjecttot, nparam, nmodel);
    
    lpp = zeros(subjecttot, nmodel);
    report = zeros(subjecttot, nmodel);
    gradient = cell(subjecttot, nmodel);
    hessian_lpp = cell(subjecttot, nmodel);
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
                p1,...
                l1,...
                rep1,...
                output,...
                lmbda,...
                grad1,...
                hess1,...
                ] = fmincon(...
                @(x) getlpp(x, con{nsub}, cho{nsub}, out{nsub}, model),...
                [1, .5, .5, 0, 0, .5, .15, .15],...
                [], [], [], [],...
                [0, 0, 0, -1, -2, 0, 0, 0],...
                [300, 1, 1, 1, 2, 1, 1, 1],...
                [],...
                options...
                );
            parameters_lpp(nsub, :, model) = p1;
            lpp(nsub, model) = l1;
            hessian_lpp{nsub, model} = hess1;
        end
    end
    %%Save the data
    data = containers.Map({'parameters_lpp','lpp', 'hessian_lpp'},...
        {parameters_lpp, lpp, hessian_lpp});
    save('data/fit_exp/data_new_ms_lpp', 'data');
    close(w);
    
end

params = {
    [1, 2],...
    [1, 2, 3],...
    [1, 2, 3],...
    [1, 2, 4],...
    [1, 2, 5],...
    [1, 2, 5, 6],...
    [1, 2, 3, 5],...
    [1, 2, 3, 5, 6]
};

% for model = whichmodel  
%     figure;
%     suptitle(models{model});
%     i = 0;
%     for p = params{model}
%         
%         i = i +1;
%         subplot(1, length(params{model}), i);
%         
%         if ismember(p, [2, 3])
%             dist = 'beta';
%         elseif ismember(p, [1])
%             dist = 'gamma';
%         elseif ismember(p, [6])
%             dist = 'gamma';
%         else
%              parameters_lpp(:, p, model) = ...
%                  rescale(parameters_lpp(:, p, model), 0, 1);
%             dist = 'normal';
%         end
%         disp(model);
%         histfit(parameters_lpp(:, p, model), 10, dist);
%         ylim([0, 100]);
%         title(gca, paramlabels{p});
%     end
% end
% clear dist;
dist = zeros(1, 2);
for p = [1, 2, 3, 5, 6]
    if ismember(p, [2, 3])
        dist(p, :) = betafit(parameters_lpp(:, p, 7));
    elseif ismember(p, [1, 6])
        if p == 6
            dist(p, :) = gamfit(parameters_lpp(:, p, 6));
        else
            dist(p, :) = gamfit(parameters_lpp(:, p, 7));
        end
    else
        dist(p, :) = normfit(parameters_lpp(:, p, 7));
    end
end
    
  
    
    
    