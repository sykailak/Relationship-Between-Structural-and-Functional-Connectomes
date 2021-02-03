
%% Load structural data

str_data_dir = dir('connectomes-data/Task2Data/*_WFA_68.csv');

for t = 1:size(str_data_dir,1) %1~19
    str_data{t} = readmatrix(str_data_dir(t).name);
end

clearvars str_data_dir
%% Load functional data
f_data_dir = dir('connectomes-data/Task2Data/*_rsfMRI_68.csv');

for t = 1:size(f_data_dir,1) %1~19
    f_data{t} = readmatrix(f_data_dir(t).name);
end

clearvars f_data_dir t



%% Estimate the structural connectivity density for each vertex (the sum across rows or
% columns of each structural connectivity matrix). Do the same for the functional connectivity
% matrices, and then model functional connectivity density based on structural connectivity
% density, independently for each vertex. Do you find any strong associations between the
% modalities using this approach?









%% preprocess data for easier data handling

% Concatenate all connectivity matrices into (68,68,19)
str_all = zeros(68,68,19);
f_all = zeros(68,68,19);

for i = 1:68
    for j = 1:68
        for s = 1:19
            str_all(i,j,s) = str_data{s}(i,j);
            f_all(i,j,s) = f_data{s}(i,j);
        end
    end
end

clearvars i j s


%% Get all UNIQUE connections, not including the diagonal connections
% Modify into shape (19, 2278)

str = zeros(19,2278); % 2278 = (1/2)*68*67
f = zeros(19,2278);
k = 1;
for i = 1:68
    for j = 1:68
        if i>j
            str(:,k) = squeeze(str_all(i,j,:));
            f(:,k) = squeeze(f_all(i,j,:));
            k=k+1
        end
    end
end

clearvars k i j
%% 
X = str;
Y = f(:,3);
rng default % For reproducibility 
[b,fitinfo] = lasso(X,Y,'CV',5);
lassoPlot(b,fitinfo,'PlotType','Lambda','XScale','log');

lassoPlot(b,fitinfo,'PlotType','CV');
legend('show') % Show legend

