close all
clear all
format longG
addpath './simulation'
addpath './fit'
addpath './utils'

ncond = 9;
rewards = cell(ncond, 1, 1);
probs = cell(ncond, 1, 1);

% j = 1;
% for i = [0.75, 0.5, 0.25]
%     rewards{j} = {[0, 0], [-1, 1]};
%     probs{j} = {[0.5, 0.5], [1-i, i]};
%     j = j + 1;
% end
% for i = [0.75, 0.5, 0.25]
%     rewards{j} = {[-1, -1], [-2, 0]};
%     probs{j} = {[0.5, 0.5], [1-i, i]};
%     j = j + 1;
% end
% for i = [0.75, 0.5, 0.25]
%     rewards{j} = {[1, 1], [0, 2]};
%     probs{j} = {[0.5, 0.5], [1-i, i]};
%     j = j + 1;
% end

data = load('learningdata');
data = data.learningdata(:, 1:18);
ncond = max(data(:, 13));
nsession = max(data(:, 18));
sub_ids = unique(data(:, 2));
i = 0;
for id = 1:length(sub_ids)
    sub = sub_ids(id);
    mask_sub = data(:, 2) == sub;
%     if ~(ismember(sum(mask_sub), [213, 228]))
%         disp(sprintf('%d : %d', sum(mask_sub), sub));
%     end
    if ismember(sum(mask_sub), [213, 228])
        i = i + 1;
        for cond = 0:ncond 
            mask_sub = data(:, 2) == sub;
            mask_cond = data(:, 13) == cond;
            mask_sess1 = ismember(data(:, 18), [0, 1, 2]);
            mask_sess2 = data(:, 18) == 3;
            mask_exp = logical(mask_sub .* mask_cond .* mask_sess1);
            mask_lot = logical(mask_sub .* mask_cond .* mask_sess2);
            %[noneed, trialorder] = sort(data(mask, 12));
            tempcho = data(:, 9); 
            cho_exp(i, :, cond+1) = tempcho(mask_exp);
            cho_lot(i, :, cond+1) = tempcho(mask_lot);
        end
    end
end
disp(sprintf('N subjects %d', i));
i = 0;
j = 1;
conds = [9 8 7 3 2 1 6 5 4];
p_risky_exp = zeros(3, 3); 
p_risky_lot = zeros(3, 3); 
for cond = conds
    
    i = i + 1;
    disp(i);
    disp(j);
    p_risky_exp(j, i) = sum(reshape(cho_exp(:, :, cond) == 2, [], 1)) / length(reshape(...
        cho_exp(:, :, cond), [], 1));
    p_risky_lot(j, i) = sum(reshape(cho_lot(:, :, cond) == 2, [], 1)) / length(reshape(...
       cho_lot(:, :, cond), [], 1));
   

   if mod(i, 3) == 0
       j = j + 1;
       i = 0;
   end
end
%disp(sprintf('disp subject %d', i));
%disp(sprintf('disp subject %d', length(sub_ids)));
ylabels = {'R', 'N', 'P'};
xlabels = {.25, .5, .75};
yaxislabel = 'Framing';
xaxislabel = 'P win risky';
cname = 'p(risky)';

plotheatmap(p_risky_exp, xlabels, ylabels, xaxislabel, yaxislabel, 'Experiment',...
    cname);
plotheatmap(p_risky_lot, xlabels, ylabels, xaxislabel, yaxislabel, 'Lotteries',...
    cname);
plotheatmap([p_risky_lot-p_risky_exp], xlabels, ylabels, xaxislabel, yaxislabel, 'Gap',...
    cname);

function plotheatmap(data, xlabels, ylabels, xaxislabel, yaxislabel, titlelabel,...
    cname)
    colormap default
    f = figure('Renderer', 'painters', 'Position', [10 10 800 700]);
    h = heatmap(xlabels, ylabels, data);
        colormap default

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
    %ylim([0, 1]);
    %then ax becomes the handle to the heatmap
    c = axs.Colorbar;    %now you have a handle to the colorbar object    
    c.Label.String = cname;
    %set(c.Label, 'Rotation', -360+-90);
    %saveas(f, sprintf('fig/2sessions_recovery_%s_%d.png', titlelabel, length(xlabels)));
        colormap default

end

