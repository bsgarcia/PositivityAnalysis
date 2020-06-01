close all
clear all
format longG
addpath './simulation'
addpath './fit'
addpath './utils'

ncond = 9;
rewards = cell(ncond, 1, 1);
probs = cell(ncond, 1, 1);

j = 1;
for i = [0.75, 0.5, 0.25]
    rewards{j} = {[0, 0], [-1, 1]};
    probs{j} = {[0.5, 0.5], [1-i, i]};
    j = j + 1;
end
for i = [0.75, 0.5, 0.25]
    rewards{j} = {[-1, -1], [-2, 0]};
    probs{j} = {[0.5, 0.5], [1-i, i]};
    j = j + 1;
end
for i = [0.75, 0.5, 0.25]
    rewards{j} = {[1, 1], [0, 2]};
    probs{j} = {[0.5, 0.5], [1-i, i]};
    j = j + 1;
end

data = load('learningdata');
data = data.learningdata(:, 1:18);
ncond = max(data(:, 13));
nsession = max(data(:, 18));
sub_ids = unique(data(:, 2));
i = 1;
for id = 1:length(sub_ids)
    sub = sub_ids(id);
    mask_sub = data(:, 2) == sub;
    if ~(ismember(sum(mask_sub), [213, 228]))
        disp(sprintf('%d : %d', sum(mask_sub), sub));
    end
    if ismember(sum(mask_sub), [213, 228])
        for cond = 0:ncond 
            mask_sub = data(:, 2) == sub;
            mask_cond = data(:, 13) == cond;
            mask_sess = ismember(data(:, 18), [0, 1, 2]);
            mask = logical(mask_sub .* mask_cond .* mask_sess);
            [noneed, trialorder] = sort(data(mask, 12));
            tempcho = data(mask, 9); 
            cho(i, :, cond+1) = tempcho(trialorder);
            tempout = data(mask, 7); 
            out(i, :, cond+1) = tempout(trialorder);
            tempcorr = data(mask, 10);
            corr(i, :, cond+1) = tempcorr(trialorder);   
            temprew = data(mask, 17);
            rew(i, :, cond+1) = temprew(trialorder);
        end
         if sum(out(i, :)) > 4
            i = i + 1;
         end
    end
end
error('');

disp(sprintf('disp subject %d', i));
disp(sprintf('disp subject %d', length(sub_ids)));

color = get(gca,'ColorOrder');
colors = repelem({color(1, :), color(7, :), color(5, :)}, 3);

i = 1;
conds = [9 8 7 3 2 1 6 5 4];
for cond = conds
    subplot(3, 3, i)
    if ismember(cond, [2, 5, 8])
        reversalplot(...
            (cho(:, :, cond)==2)',... %data
            [],... %time when reversal occurs
            [],... %time when cond changes
            ones(3) * 0.5,... % chance lvl
            colors{cond},... %curve color
            0.9,... %linewidth
            0.3,... % alpha
            0, 1,... % ylims
            15,... %fontsizemat,
            probs{cond}{2}(2),... %title,
            'trials',... %xlabel
            'risky choice rate'... % ylabel
    );
    else
    reversalplot(...
        corr(:, :, cond)',... %data
        [],... %time when reversal occurs
        [],... %time when cond changes
        ones(3) * 0.5,... % chance lvl
        colors{cond},... %curve color
        0.9,... %linewidth
        0.3,... % alpha
        0, 1,... % ylims
        15,... %fontsizemat,
        probs{cond}{2}(2),... %title,
        'trials',... %xlabel
        'correct choice'... % ylabel
    );
    end
    i = i + 1;
end


% i = 1;
% for cond = conds
%     subplot(3, 3, i)
%         reversalplot(...
%             (cho(:, :, cond)==2)',... %data
%             [],... %time when reversal occurs
%             [],... %time when cond changes
%             ones(3) * 0.5,... % chance lvl
%             colors{cond},... %curve color
%             0.9,... %linewidth
%             0.3,... % alpha
%             0, 1,... % ylims
%             15,... %fontsizemat,
%             probs{cond}{2}(2),... %title,
%             'trials',... %xlabel
%             'risky choice rate'... % ylabel
%     );
%     i = i + 1;
% end

% data = load('learningdata');
% data = data.learningdata(:, 1:18);
% ncond = max(data(:, 13));
% nsession = max(data(:, 18));
% sub_ids = unique(data(:, 2));
% i = 1;
% clear cho;
% clear out;
% clear rew;
% clear corr;
% for id = 1:length(sub_ids)
%     sub = sub_ids(id);
%     mask_sub = data(:, 2) == sub;
%     disp(sum(mask_sub));
%     if ismember(sum(mask_sub), [213, 228])
%         for cond = 0:ncond 
%             mask_sub = data(:, 2) == sub;
%             mask_cond = data(:, 13) == cond;
%             mask_sess = data(:, 18) == 3;
%             mask = logical(mask_sub .* mask_cond .* mask_sess);
%             [noneed, trialorder] = sort(data(mask, 12));
%             tempcho = data(mask, 9); 
%             cho(i, :, cond+1) = tempcho(trialorder);
%             tempout = data(mask, 7); 
%             out(i, :, cond+1) = tempout(trialorder);
%             tempcorr = data(mask, 10);
%             corr(i, :, cond+1) = tempcorr(trialorder);   
%             temprew = data(mask, 17);
%             rew(i, :, cond+1) = temprew(trialorder);
%         end
%         i = i + 1;
%     end
% end
% 
% figure
% 
% i = 1;
% conds = [9 8 7 3 2 1 6 5 4];
% for cond = conds
%     subplot(3, 3, i)
%     
%         reversalplot(...
%             (cho(:, :, cond)==2)',... %data
%             [],... %time when reversal occurs
%             [],... %time when cond changes
%             ones(3) * 0.5,... % chance lvl
%             colors{cond},... %curve color
%             0.9,... %linewidth
%             0.3,... % alpha
%             0, 1,... % ylims
%             15,... %fontsizemat,
%             probs{cond}{2}(2),... %title,
%             'trials',... %xlabel
%             'risky choice rate'... % ylabel
%     );
%     
%     i = i + 1;
% end
% 
