


% Manually threshold the FA map at different values between 0.1 and 0.8,
% setting subthreshold voxels to zero. Discuss the spatial distribution of 
% voxels remaining above threshold as the latter increases.

FAmap = load_nii('connectomes-data/Task1Data/tractor/diffusion/dti_FA.nii.gz') %(96,96,60)


% load parcellation data
d_parcellation = load_nii('connectomes-data/Task1Data/tractor/diffusion/parcellation.nii.gz')
%(96,96,60)



%%

%Because FAmap is 3D
disp('range in FAmap img:')
disp([min(min(min(FAmap.img))), max(max(max(FAmap.img)))]);
disp('size of FAmap img:');
disp(size(FAmap.img));

% set threshold range
threshold = [0.1:0.1:0.8];

for t = 1:length(threshold)
    % create copy
    FAmap_store = FAmap.img;
    parcel_store = d_parcellation.img;
    % voxels with values lower than threshould equal to 0
    %FAmap_store(FAmap_store < threshold(t)) = 0;
    tf_matrix{t} = FAmap_store > threshold(t);
    parcel_store = parcel_store(tf_matrix{t});
    % store values in ROI in a cell
    %fa_maps{t} = FAmap_store;
    % store INDICES of ROI values in a list (536349,1)
    %non_zeros{t} = find(fa_maps{t});
    %roi_indices = non_zeros{t};
    
    %Take regions affected by thresholds
    %for y = 1:length(roi_indices)
    %    roi_inc{t}(y) = d_parcellation.img(roi_indices(y)); %which rois are present in each map
    %end 
end
%%
for t = 1:8
    subplot(2,4,t)
    mat = tf_matrix{t}(:,:,30)
    [r, c] = size(mat);                          % Get the matrix size
    imagesc((1:c)+0.5, (1:r)+0.5, mat);          % Plot the image
    colormap(gray);                              % Use a gray colormap
    axis equal                                   % Make axes grid sizes equal
    title(sprintf('FA threshold of %.1f', threshold(t)))
    set(gca, 'XTick', [], 'YTick', [], ...  % Change some axes properties
             'XLim', [1 c+1], 'YLim', [1 r+1], ...
             'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
end
%% Create cell to store number of elements for each roi in each thresholded map 

% Load region information
fileId = fopen('connectomes-data/Task1Data/tractor/functional/parcellation.lut');
region_info = textscan(fileId, '%d %s %s %s %s %s %s', 'Headerlines', 3); % skip the first unknown region
region_ids = region_info{1}; %unique regions
region_type = region_info{5}; %types of regions

for x = 1:size(roi_inc,2) %1~8
    % select roi for threshold
    roi_loc = roi_inc{x};
    for y = 1:length(region_ids)
        %stores number of elements for each roi in each thresholded map
        spat_loc{x}{y,1} = (region_ids(y));
        spat_loc{x}{y,2} = sum(roi_loc(:) == region_ids(y)); 
        spat_loc{x}{y,3} = convertCharsToStrings(region_type(y));
    end
end
%% Spatial distribution of regions for different thresholds 

for t = 1:length(threshold)
    region_freq=[];
    region_freq=[region_freq,spat_loc{t}{:,2}];
    region_name=[];
    region_name=[region_name,spat_loc{t}{:,3}];
    
    region_list = unique(region_name);
    
    table=[region_list;zeros(1,6)]';
    for i = 1:length(region_list)
        this_reg = region_list(i);
        idx = region_name==this_reg;
        table(i,2) = sum(region_freq(idx));
    end 
    
    final_tab{t}=table;
    subplot(4,2,t)
    bar(str2double(table(:,2)));
    region_list(1) = "cerebellar\newlinecortex";
    region_list(2) = "cerebral\newlinecortex";
    region_list(3) = "cerebrospinal\newlinefluid";
    region_list(4) = "deep grey\newlinematter";
    region_list(6) = "white matter";

    set(gca,'xticklabel', region_list);
    ylim([0 10^4])
    title(sprintf('Threshold %0.1f' ,threshold(t)));


end

sgtitle('Number of streamlines for different FA thresholds');



%%
for t = 1:8
    subplot(2,4,t)
    mat = FA_graphs{t}
    [r, c] = size(mat);                          % Get the matrix size
    imagesc((1:c)+0.5, (1:r)+0.5, mat);          % Plot the image
    colormap(gray);                              % Use a gray colormap
    axis equal                                   % Make axes grid sizes equal
    title(sprintf('FA threshold of %.1f', threshold(t)))
    set(gca, 'XTick', [], 'YTick', [], ...  % Change some axes properties
             'XLim', [1 c+1], 'YLim', [1 r+1], ...
             'GridLineStyle', '-', 'XGrid', 'on', 'YGrid', 'on');
end


%% reading the csv files

filelist = {'FA.1_graph.csv'; 'FA.2_graph.csv'; 'FA.3_graph.csv';...
            'FA.4_graph.csv'; 'FA.5_graph.csv'; 'FA.6_graph.csv'; 'FA.7_graph.csv';...
            'FA.8_graph.csv'};
        
for file = 1:size(filelist,1) %1~8
    
    %first three lines are headers
    FA_graphs{file} = readmatrix([pwd '/connectomes-data/Task1Data/connectomes/' filelist{file}], 'HeaderLines',3);
    
    %edge density in index 1, number of vertices in 2, number of edges in 3
    %[density{file}(1), density{file}(2), density{file}(3)] = ...
    density{file}=density_und(FA_graphs{file}); 
     %density_und(FA_graphs{file}); 
 
    %Average shortest path
    %distance matrices and characteristic path length
    dist_mat{file} = distance_bin(FA_graphs{file}); 
    %network characteristic path length in (1), global efficiency in (2)
    [char_path{file}(1), char_path{file}(2)]  = charpath(dist_mat{file},0,0); 
    
    
    %network characteristic path length in (1), global efficiency in (2)
    %[Not inverse]
    
    %global efficiency
    %efficiency{file} = 1 / efficiency_bin(FA_graphs{file}); 
    
    %mean clustering coefficient
    mean_cluster_coeff{file} = mean(clustering_coef_bu(FA_graphs{file})); 
    
end

% As FA increases...
% Density decreases
% average shortest path increase except for the last
% efficiency increases
% mean cluster coneff non monotonic



%% Plotting

figure;
%[ha, pos] = tight_subplot(3,1,[.01 .03],[.1 .01],[.01 .01])
[ha, pos] = tight_subplot(3,1,[.1 .1],[.1 .05],[.05 .01])


%for ii = 1:3;
axes(ha(1)); 
plot([density{:}])
title({'','','Density'},'FontSize',15)
set(gca,'XTickLabel',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8] );
axes(ha(2));
plot([char_path{1}(1),char_path{2}(1),char_path{3}(1),char_path{4}(1),char_path{5}(1),char_path{6}(1),char_path{7}(1),char_path{8}(1)])
%xlabel('FA thresholds') 
title({'','','Average shortest path'},'FontSize',15)
set(gca,'XTickLabel',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8] );

axes(ha(3));
plot([mean_cluster_coeff{:}])
xlabel('FA thresholds') 
title({'','','Mean clustering coefficient'},'FontSize',15)
set(gca,'XTickLabel',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8] );



%% Calculate a weighted FUNCTIONAL connectome, by first extracting an average time series in
% each cortical region, and then deriving a correlation matrix using the shrinkage approach
% of Schäfer & Strimmer (2005).

% load parcellation data
parcellation = load_nii('connectomes-data/Task1Data/tractor/functional/parcellation.nii.gz');

% load functional data
data_func = load_nii('connectomes-data/Task1Data/tractor/functional/data.nii.gz');
%(64,64,30,15)

% Load region information
fileId = fopen('connectomes-data/Task1Data/tractor/functional/parcellation.lut');
region_info = textscan(fileId, '%d %s %s %s %s %s %s', 'Headerlines', 3); % skip the first unknown region
region_ids = region_info{1}; %unique regions
num_regions = size(region_ids);


% get the average time series for each region:
avg_time_series = zeros(num_regions(1), 15);
for i = 1:num_regions(1)
    id = region_ids(i);
    % one time-series for each region (averaged over the voxels in that region)
    % (we have 15 time steps)
    
    % Create True/False mask for the region over (64,64,30)
    region = (parcellation.img == id); 
    
    % match the size of data_func.img and region to (64,64,30,15) to
    % successfully select voxels of the region in question
    voxel_list = data_func.img(region(:,:,:,ones(1,15)));
    
    % The masked result is a vector of 15's time series. We reshape.
    size_voxel_list = size(voxel_list); %(n,1)
    region_voxels = reshape(voxel_list, [size_voxel_list(1)/15, 15]); 
    
    % Take mean and store the results
    region_avg_time_series = mean(region_voxels, 1);
    avg_time_series(i,:) = region_avg_time_series;
end


%% Correlation matrix (function definition at the end of the script)

[n, p] = size(avg_time_series'); %15,108
sx = makeMeanZero(avg_time_series'); sx = makeStdOne(sx); % convert S to R
[r, vr] = varcov(sx);

for i = 1:length(threshold)
% Eqn on p4 of Schafer and Strimmer 2005
% lambda denotes shrinkage intensity
lambda = threshold(i);
lambda = min(lambda, 1); lambda = max(lambda, 0); 
Rhat = (1-lambda)*r;
Rhat(logical(eye(p))) = 1; %diagonal values are 1
cor_matrix{i} = Rhat;
end



%% threshold resulting matrix at a correlation 0.1 & binarise it

for t = 1:length(threshold)
    % create copy
    cor_store = cor_matrix{t};
    % voxels with values lower than threshould equal to 0
    cor_store(cor_store < 0.1) = 0;
    cor_store(cor_store >= 0.1) = 1;
    binary_cor{t}=cor_store;

end


%% calculate graph metrics 
       
for file = 1:size(binary_cor,2) %1~8
    
    %edge density in index 1, number of vertices in 2, number of edges in 3
    %[cor_density{file}(1), cor_density{file}(2), cor_density{file}(3)] = ...
    % density_und(binary_cor{file}); 
 
     cor_density{file} = density_und(binary_cor{file}); 

    %Average shortest path
    %distance matrices and characteristic path length
    cor_dist_mat{file} = distance_bin(binary_cor{file}); 
    %network characteristic path length in (1), global efficiency in (2)
    [cor_char_path{file}(1), cor_char_path{file}(2)]  = charpath(cor_dist_mat{file},0,0); 
    %network characteristic path length in (1), global efficiency in (2)
    %[Not inverse]
    
    %global efficiency
    cor_efficiency{file} = 1 / efficiency_bin(binary_cor{file}); 
    
    %mean clustering coefficient
    letssee{file} = clustering_coef_bu(binary_cor{file});
    cor_mean_cluster_coeff{file} = mean(clustering_coef_bu(binary_cor{file})); 
    
end

% As shrinkage intensity increases...
% edge density decreases
% Average shortest path increases
% global efficiency increases
% mean clustering coefficient decreases


%% plotting



figure;

subplot(4,1,1)
plot([cor_density{:}])
xlabel('FA thresholds') 
title({'','','Density'},'FontSize',15)

subplot(4,1,2)
plot([cor_char_path{1}(1),cor_char_path{2}(1),cor_char_path{3}(1),cor_char_path{4}(1),cor_char_path{5}(1),cor_char_path{6}(1),cor_char_path{7}(1),cor_char_path{8}(1)])
xlabel('FA thresholds') 
title({'','','Average shortest path'},'FontSize',15)

subplot(4,1,3)
plot([cor_char_path{1}(2),cor_char_path{2}(2),cor_char_path{3}(2),cor_char_path{4}(2),cor_char_path{5}(2),cor_char_path{6}(2),cor_char_path{7}(2),cor_char_path{8}(2)])
xlabel('FA thresholds') 
title({'','','Efficiency'},'FontSize',15)


subplot(4,1,4)
plot([cor_mean_cluster_coeff{:}])
xlabel('FA thresholds') 
title({'','','Mean clustering coefficient'},'FontSize',15)

%% Plotting

figure;
%[ha, pos] = tight_subplot(3,1,[.01 .03],[.1 .01],[.01 .01])
[ha, pos] = tight_subplot(3,1,[.1 .1],[.1 .05],[.05 .01])


%for ii = 1:3;
axes(ha(1)); 
plot([cor_density{:}])
title({'','','Density'},'FontSize',15)
set(gca,'XTickLabel',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8] );

axes(ha(2));
plot([cor_char_path{1}(1),cor_char_path{2}(1),cor_char_path{3}(1),cor_char_path{4}(1),cor_char_path{5}(1),cor_char_path{6}(1),cor_char_path{7}(1),cor_char_path{8}(1)])
%xlabel('FA thresholds') 
title({'','','Average shortest path'},'FontSize',15)
set(gca,'XTickLabel',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8] );

axes(ha(3));
plot([cor_mean_cluster_coeff{:}])
xlabel('FA thresholds') 
title({'','','Mean clustering coefficient'},'FontSize',15)
set(gca,'XTickLabel',[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8] );



%%
test1 = (clustering_coef_bu(binary_cor{7}))'
%%
test = (clustering_coef_bu(binary_cor{8}))'

%% Consider the effects of retaining or discarding negative correlations with
% absolute value greater than 0.1, and what significance these edges may have 
% in the context of brain activity.


% Discarding negative correlations greater than |0.1|
for t = 1:length(threshold)
    % create copy
    cor_store = cor_matrix{t};
    cor_store=(cor_store > -0.1);
    discard_negcor{t} = cor_store;
end

for file = 1:size(discard_negcor,2) %1~8
    
    %edge density in index 1, number of vertices in 2, number of edges in 3
    %[neg_density{file}(1), neg_density{file}(2), neg_density{file}(3)] = ...
    % density_und(discard_negcor{file}); 
 
    neg_density{file} = density_und(discard_negcor{file}); 

    %Average shortest path
    %distance matrices and characteristic path length
    neg_dist_mat{file} = distance_bin(discard_negcor{file}); 
    %network characteristic path length in (1), global efficiency in (2)
    [neg_char_path{file}(1), neg_char_path{file}(2)]  = charpath(neg_dist_mat{file},0,0); 
    %network characteristic path length in (1), global efficiency in (2)
    %[Not inverse]
    
    %global efficiency
    neg_efficiency{file} = 1 / efficiency_bin(discard_negcor{file}); 
    
    %mean clustering coefficient
    neg_mean_cluster_coeff{file} = mean(clustering_coef_bu(discard_negcor{file})); 
    
end


%% plotting


figure;

subplot(4,1,1)
plot([neg_density{:}])
xlabel('FA thresholds') 
title({'','','Density'},'FontSize',15)

subplot(4,1,2)
plot([neg_char_path{1}(1),neg_char_path{2}(1),neg_char_path{3}(1),neg_char_path{4}(1),neg_char_path{5}(1),neg_char_path{6}(1),neg_char_path{7}(1),neg_char_path{8}(1)])
xlabel('FA thresholds') 
title({'','','Average shortest path'},'FontSize',15)

subplot(4,1,3)
plot([neg_char_path{1}(2),neg_char_path{2}(2),neg_char_path{3}(2),neg_char_path{4}(2),neg_char_path{5}(2),neg_char_path{6}(2),neg_char_path{7}(2),neg_char_path{8}(2)])
xlabel('FA thresholds') 
title({'','','Efficiency'},'FontSize',15)


subplot(4,1,4)
plot([neg_mean_cluster_coeff{:}])
xlabel('FA thresholds') 
title({'','','Mean clustering coefficient'},'FontSize',15)




%%






%%


% NOT Discarding negative correlations greater than |0.1|

for file = 1:size(discard_negcor,2) %1~8
    
    %edge density in index 1, number of vertices in 2, number of edges in 3
    %[neg_density{file}(1), neg_density{file}(2), neg_density{file}(3)] = ...
    % density_und(discard_negcor{file}); 
 
    not_neg_density{file} = density_und(cor_matrix{file}); 

    %Average shortest path
    %distance matrices and characteristic path length
    not_neg_dist_mat{file} = distance_bin(cor_matrix{file}); 
    %network characteristic path length in (1), global efficiency in (2)
    [not_neg_char_path{file}(1), not_neg_char_path{file}(2)]  = charpath(not_neg_dist_mat{file},0,0); 
    %network characteristic path length in (1), global efficiency in (2)
    %[Not inverse]
    
    %global efficiency
    not_neg_efficiency{file} = 1 / efficiency_bin(cor_matrix{file}); 
    
    %mean clustering coefficient
    not_neg_mean_cluster_coeff{file} = mean(clustering_coef_bu(cor_matrix{file})); 
    
end


%% plotting


figure;

subplot(4,1,1)
plot([not_neg_density{:}])
xlabel('FA thresholds') 
title({'','','Density'},'FontSize',15)

subplot(4,1,2)
plot([not_neg_char_path{1}(1),not_neg_char_path{2}(1),not_neg_char_path{3}(1),not_neg_char_path{4}(1),not_neg_char_path{5}(1),not_neg_char_path{6}(1),not_neg_char_path{7}(1),not_neg_char_path{8}(1)])
xlabel('FA thresholds') 
title({'','','Average shortest path'},'FontSize',15)

subplot(4,1,3)
plot([not_neg_char_path{1}(2),not_neg_char_path{2}(2),not_neg_char_path{3}(2),not_neg_char_path{4}(2),not_neg_char_path{5}(2),not_neg_char_path{6}(2),not_neg_char_path{7}(2),not_neg_char_path{8}(2)])
xlabel('FA thresholds') 
title({'','','Efficiency'},'FontSize',15)


subplot(4,1,4)
plot([not_neg_mean_cluster_coeff{:}])
xlabel('FA thresholds') 
title({'','','Mean clustering coefficient'},'FontSize',15)

%% Function Definitions
function [S, VS] = varcov(x)
% s(i,j) = cov X(i,j)
% vs(i,j) = est var s(i,j)
[n,p] = size(x);
xc = makeMeanZero(x); 
S = cov(xc);
XC1 = repmat(reshape(xc', [p 1 n]), [1 p 1]); % size p*p*n !
XC2 = repmat(reshape(xc', [1 p n]),  [p 1 1]); % size p*p*n !
VS = var(XC1 .* XC2, 0,  3) * n/((n-1)^2);
end

function xc = makeMeanZero(x)
% make column means zero
[n,p] = size(x);
m = mean(x);
xc = x - ones(n, 1)*m; 
end

function xc = makeStdOne(x)
% make column variances one
[n,p] = size(x);
sd = ones(n, 1)*std(x);
xc = x ./ sd; 
end


