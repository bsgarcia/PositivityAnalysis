clear all
close all
format longg
data = load('data/data_online_exp/learningdata');
%ids = loadcsv('learning_data.csv')
format longg

data = data.learningdata(:, 1:18);
ncond = max(data(:, 13));
nsession = max(data(:, 18));
sub_ids = unique(data(:, 2));
i = 1;
%matrix = zeros(1, 5);
for id = 1:length(sub_ids)
    sub = sub_ids(id);
    mask_sub = data(:, 2) == sub;
    if sum(mask_sub) > 100
        mask_sub = data(:, 2) == sub;
        mask_sess = ismember(data(:, 18), [0, 1, 2, 3]);
        mask = logical(mask_sub .* mask_sess);
        out = data(mask, 7);
        matrix(i, 1) = sub;
        matrix(i, 2) = (sum(out) * 7.5) / 100;
        %matrix{i, 5} = length(out);
        i = i + 1;
    end
end
format longg
save('matrix');

writetable(array2table(matrix));
