clearvars -except data_path AgeFolderList dataFolderList
close all

%data_path = "Z:\Jui-Yen Huang\Derived Files\Data\In vivo Calcium Imaging\20200827-0918_jRGECO1a_FGF9-Ncre-GAD30-7\20200827_jRGECO1a_FGF9-Ncre-GAD30-7_P11\P11_30-7ALBL_200";
%data_path = "Z:\Jui-Yen Huang\Derived Files\Data\In vivo Calcium Imaging\20210914_0924_jRGECO1a_C57-C2-1\20210924_jRGECO1a_C57-C2-1_P21\P21_C2-1_BR_L3_180_stim"
load_path = fullfile(data_path,'.\Modol_outputs\all.mat');
%%load data and construct similarity matrix
%here we use the covariance matrix generated in modol paper
% load ('.\Modol_outputs\all.mat'); % this loads the covariance matrix from modol output
load(load_path);

clearvars -except CovM Race start_time end_time detected_events AgeFolderList dataFolderList data_path
%%
%%identify current wowrking folder for naming file
[root_path,animal_path,~] = fileparts(data_path); %pwd: display the current folder
animal_ID = split(animal_path, '_');

%make the folder in the current working directory to save all outpus
store_folder_path = fullfile(data_path,'DBSCAN_output');
if ~exist(store_folder_path, 'dir')
    mkdir(store_folder_path)
end

% construct save path
save_path = fullfile(store_folder_path,'output_mat.mat');
%% 
% In DBSCAN (Density-Based Spatial Clustering of Applications with Noise), choosing 
% the right values for epsilon (ε) and minimum points (minPts) is crucial for 
% effective clustering. 
% 
% In our case /goal: 
%% 
% * *Temporal Coherence*: If neuronal patterns tend to persist for a few frames, 
% cluster sizes should reflect this. For example, |minPts| could be chosen as 
% the average duration of such patterns in frames. --in our case, we do not expect 
% a pattern need to persist for how many frames. While it must be at least 3 point. 
% So we set *Minimum Points (minPts)* as 3.
% * *Epsilon (ε)*: Determines the threshold for two frames to belong to the 
% same cluster. Use the *k-distance graph* to estimate find this value for specific 
% dataset
%% 
% In addition, the input matrix must larger than 0. As other clustering method 
% all use Covariance Matrix, and the range is -1 to 1. Have to conver it become 
% positive linerly.

% Convert Correlation to Distance
distMatrix = 1 - CovM; % Convert correlation to distance
distMatrix = distMatrix / max(distMatrix(:)); % Normalize distances to [0, 1]

% Compute k-Distance Graph for Epsilon Selection
k = 5; % Choose k
kDist = zeros(size(distMatrix, 1), 1);
for i = 1:size(distMatrix, 1)
    sortedD = sort(distMatrix(i, :));
    kDist(i) = sortedD(k); % Distance to the k-th nearest neighbor
end
sortedDist = sort(kDist);

% Find knee point for epsilon
startPoint = [1, sortedDist(1)];
endPoint = [length(sortedDist), sortedDist(end)];
distances = zeros(size(sortedDist));
for i = 1:length(sortedDist)
    point = [i, sortedDist(i)];
    distances(i) = abs(det([endPoint - startPoint; point - startPoint])) / norm(endPoint - startPoint);
end

% Identify the best epsilon with a refinement condition
% find the appropriate epilison which is below 0.25
remainingPeaks = distances; % Copy distances to track remaining peaks
validEpsilon = false;
%lastCheckedEpsilon = NaN; % Initialize to store the last checked epsilon

while ~validEpsilon
    [~, kneeIdx] = max(remainingPeaks); % Find the current maximum peak
    epsilon = sortedDist(kneeIdx);     % Corresponding epsilon value
    lastCheckedEpsilon = epsilon;      % Store the last checked epsilon

    if epsilon <= 0.25
        validEpsilon = true; % Stop if epsilon is valid
    else
        remainingPeaks(kneeIdx) = -inf; % Remove this peak and find the next one
        if all(remainingPeaks == -inf)
            disp('No valid epsilon found below 0.25.'); % Fail-safe if no valid epsilon exists
            epsilon = 0.1; % Use the last checked epsilon
            validEpsilon = true; % Exit the loop after setting epsilon
        end
    end
end

% Plot k-Distance Graph as needed
clf
figure;
hold on
plot(sortedDist, 'b-', 'LineWidth', 1.5);
plot(kneeIdx, epsilon, 'ro', 'MarkerSize', 10, 'LineWidth', 2); % Highlight the knee point
hold off
xlabel('Points sorted by distance');
ylabel(sprintf('%d-th Nearest Neighbor Distance', k));
title('k-Distance Graph');
grid on;

fig_save_path_k_Distance_graph = fullfile(store_folder_path,'k_Distance_graph.fig');
saveas (gcf,fig_save_path_k_Distance_graph);
close all

% Apply DBSCAN by using the epsilon identifed with K distance graph
minPts = 3; % using the minimum points
labels = dbscan(distMatrix, epsilon, minPts, 'Distance', 'precomputed');

if all(isnan(labels), 'all')
    save(save_path);
    disp("Matrix cannot be clustered");
    return; % This will end the execution of the function or script 
end

% Analyze Clusters as needed
numClusters = numel(unique(labels(labels > 0))); % Exclude noise
%disp(['Number of clusters: ', num2str(numClusters)]);
%disp(['Noise points: ', num2str(sum(labels == -1))]);

% Visualize Clustering
clf
scatter(1:length(labels), labels, 15, labels, 'filled');
colormap('jet');
xlabel('Frame Index');
ylabel('Cluster Label');
title('DBSCAN Clustering Results');
grid on;

fig_save_path_DBSCAN_clustering_results = fullfile(store_folder_path,'DBSCAN_clustering_results.fig');
saveas (gcf,fig_save_path_DBSCAN_clustering_results);
close all
%%
clustered_idx = find(labels ==1);
if isempty(clustered_idx)
    disp("no cluster")
    C0 = {};
    save (save_path); 
    return
end
%%
%% Below start using the output from DBSCAN and analyzed by the rest of Model code
[NCell,NRace] = size(Race); %total cell/neuron number; total frame number
IDX2 = labels; % Cluster labels from DBSCAN
NCl = max(IDX2, [], 'all'); % Determine the number of clusters (including noise)

% Handle noise points (-1): Assign them to a new cluster index (NCl + 1)
noise_idx = find(IDX2 == -1);
if ~isempty(noise_idx)
    IDX2(IDX2 == -1) = NCl + 1; % Reassign noise points to a new cluster
    NCl = NCl + 1; % Update the total number of clusters
end

% Calculate silhouette values for all points
silhV_all = silh(CovM, IDX2); % Silhouette values for all clusters

% Calculate median silhouette values for each cluster
sCl_all = zeros(1, NCl); % Preallocate array for silhouette values
for i = 1:NCl
    sCl_all(i) = median(silhV_all(IDX2 == i));
end

% Remove the silhouette value of the noise cluster (last cluster)
if ~isempty(noise_idx)
    sCl_all = sCl_all(1:end-1);
end

% Sort silhouette values in descending order
[sCl_all, xCl] = sort(sCl_all, 'descend');
sCl = rmmissing(sCl_all); % Remove NaN values, if any

% Determine the number of valid clusters after sorting
NCl_beforeStat = length(sCl);
%%
%Taking the index output and use that follow the rest of Modol code
[~,x2] = sort(IDX2);
% Sorted normalized covariance matrix
MSort = CovM(x2,x2);
%Msort = distMatrix(x2,x2);

%detected events clusters and their scores
R = cell(0);
CellScore = zeros(NCell,NCl);
CellScoreN = zeros(NCell,NCl);
for i = 1:NCl
    R{i} = find(IDX2==i);
    CellScore(:,i) = sum(Race(:,R{i}),2);
    CellScoreN(:,i) = CellScore(:,i)/length(R{i});
end
%Assign cells to cluster with which it most likely spikes
[~,CellCl] = max(CellScoreN,[],2);
%Remove cells with less than 2 spikes in a given cluster
CellCl(max(CellScore,[],2)<2) = 0;
[X1,x1] = sort(CellCl);

%plot the initial clustering result
clf
f = figure ('visible','on');
subplot(2,1,1)
imagesc(MSort)
colormap jet
axis image
xlabel('RACE #')
ylabel('RACE #')
title('Covariance between GCEs')
colorbar

subplot(2,1,2)
imagesc(Race(x1,x2),[-1 1.2])
xlabel('RACE #')
ylabel('Cell #')
title('Detected Events sorted by GCE clusters and neurons')

fig_save_path_initial_cluster = fullfile(store_folder_path,'Initial_cluster.fig');
saveas (gcf,fig_save_path_initial_cluster);
close all
%%
%%For statistical significance
NShuf = 5000; %for real analysis
%NShuf = 100; %for test
disp ('NShuf')
tic

%Count number of participation to each cluster
CellP = zeros(NCell,NCl); CellR = zeros(NCell,NCl);
for i = 1:NCl
    CellP(:,i) = sum(Race(:,IDX2 == i),2);
    CellR(:,i) = CellP(:,i)/sum(IDX2 == i);
end

%Test for statistical significance
CellCl = zeros(NCl,NCell); %Binary matrix of cell associated to clusters
for j = 1:NCell
    %Random distribution among Clusters
    RClr = zeros(NCl,NShuf);
    Nrnd = sum(Race(j,:) ~= 0);
    if Nrnd == 0
        continue
    end
    for l = 1:NShuf
        Random = randperm(NRace);
        Random = Random(1:Nrnd);
        Racer = zeros(1,NRace);
        Racer(Random) = 1;
        for i = 1:NCl
            RClr(i,l) = sum(Racer(:,IDX2 == i),2);
        end
    end
    RClr = sort(RClr,2);
    %ThMin = mean(Random) - 2*std(Random);
    %Proba above 95th percentile
    ThMax = RClr(:,floor(NShuf*(1-0.05/NCl))); 
    for i = 1:NCl
        CellCl(i,j) = double(CellP(j,i)>ThMax(i));% - double(RCl(:,j)<ThMin);
    end
end

toc

A0 = find(sum(CellCl) == 0); %Cells not in any cluster
A1 = find(sum(CellCl) == 1); %Cells in one cluster
A2 = find(sum(CellCl) >= 2); %Cells in several clusters
%%
for i = A2
    [~,idx] = max(CellR(i,:));
    CellCl(:,i) = 0;
    CellCl(idx,i) = 1;
end
C0 = cell(0);
k = 0;
inds = [];
for i = 1:NCl
    if length(find(CellCl(i,:)))>2
        k = k+1;
        C0{k} = find(CellCl(i,:));
        inds = [inds; k];
    end
end

%Participation rate to its own cluster
CellParticip = max(CellR([A1 A2],:),[],2);

%C0 cell assary here is the neurons corresponding to each cluster
%%
NCl = length(C0); %the NCl here is the good cluster # that passed statistical test
if ~NCl
    NCl = 0;   
    save (save_path);    
    disp ('There were no significant clusters found!! Cannot run this cell...')
    return; % This will end the execution of the function or script 
end

%[NCell,NRace] = size(Race);
%Cell count in each cluster
RCl = zeros(NCl,NRace);
PCl = zeros(NCl,NRace);
for i = 1:NCl
    RCl(i,:) = sum(Race(C0{i},:));
end

RCln = zeros(NCl,NRace);
for j = 1:NRace
    %Random distribution among Clusters
    RClr = zeros(NCl,NShuf);
    Nrnd = sum(Race(:,j) ~= 0); % Changed since we don't have binary data
    if ~Nrnd % neuron doesn't fire in time period
        continue
    end
    for l = 1:NShuf
        Random = randperm(NCell);
        Random = Random(1:floor(Nrnd));
        Racer = zeros(NCell,1);
        Racer(Random) = 1;
        for i = 1:NCl
            RClr(i,l) = sum(Racer(C0{i}));
        end
    end
    % ThMin = mean(Random) - 2*std(Random);
    RClr = sort(RClr,2);
    % ThMin = mean(Random) - 2*std(Random);
    %Proba above 95th percentile
    ThMax = RClr(:,round(NShuf*(1-0.05/NCl)));
    for i = 1:NCl
        PCl(i,j) = double(RCl(i,j)>ThMax(i));% - double(RCl(:,j)<ThMin);
    end
    %Normalize (probability)
    RCln(:,j) = RCl(:,j)/sum(Race(:,j));
end

%%
if ~NCl
   save (save_path);
   disp('There were no significant clusters found!! Cannot run this cell...')
   return; % This will end the execution of the function or script 
end
%%
%% Sorting ROI corresponding to the assigned/significant GCE
% Times that significantly recruit 0 cell assemblies; will not plot this
Cl0 = find(sum(PCl,1) == 0);
% Times that significantly recruit 1 cell assembly
Cl1 = find(sum(PCl,1) == 1);
% etc.
Cl2 = find(sum(PCl,1) == 2);
Cl3 = find(sum(PCl,1) == 3);
Cl4 = find(sum(PCl,1) == 4);

Bin = 2.^(0:NCl-1);

%Sort Cl1
[~,x01] = sort(Bin*PCl(:,Cl1));
Cl1 = Cl1(x01);

%Sort Cl2
[~,x02] = sort(Bin*PCl(:,Cl2));
Cl2 = Cl2(x02);

%Sort Cl3
[~,x03] = sort(Bin*PCl(:,Cl3));
Cl3 = Cl3(x03);

RList = [Cl1 Cl2 Cl3 Cl4];
%x1 from DetectRace

[X1,x1] = sort(Bin*CellCl(inds, :));
%%
if ~NCl
    error('There were no significant clusters found!! Cannot run this cell...')
    
end
%%
%Plot the post-sorted Rastermap (with significant GCEs)
clf
figure ("Visible","on");
imagesc(Race(x1,RList) ~= 0);
colormap(flipud(gray))
title('Sorted Rastermap (with significant GCEs)')

fig_save_path_Sorted_Rastermap = fullfile(store_folder_path,'Sorted_Rastermap.fig');
saveas(gcf, fig_save_path_Sorted_Rastermap);
close all
%%
%%Below are for making figures for representing figure
%Define a common colormap:
%non-cluster-no fill; clustered color fill
newCmap = [255 255 255; 227 26 28; 31 120 180; 178 223 138; 51 160 44; 251 154 153; 166 206 227; 253 191 111; 255 127 0; 202 178 214; 106 61 154; 255 255 153; 177 89 40; 0 0 0;];
newCmap = newCmap ./ 255;

%sorting GCE for making figure
sortInds = [];
for ndx = 1:length(C0)
    sortInds = [sortInds, C0{ndx}];
end
for ndx = 1:size(Race, 1)
    if ~ismember(ndx, sortInds)
    sortInds = [sortInds, ndx];
    end
end

Race_sort = double(Race(sortInds,RList));
if size(RList, 2) == 0
          
end

PCl_sort = PCl(:, RList);
     
% Assign cluster number on rastermap (so colors show up in imagesc)
for jdx = 1:size(PCl_sort, 2)
    for kdx = 1:length(C0)
        neurs = C0{kdx};
        sortedNeurs = [];
        for neur = neurs
           sortedNeurs = [sortedNeurs, find(sortInds == neur)];
        end
        if PCl_sort(kdx, jdx) == 1
           for mdx = sortedNeurs
               if Race_sort(mdx, jdx)
                  Race_sort(mdx, jdx) = kdx + 1;
               end
           end
        end
    end
end

%Plotting figure with sorted cell and event
clf
Figure = figure('visible','on');
imagesc(Race_sort);
scMap = [0 0 0; newCmap(1:max(Race_sort(:)), :)];
scMap(2, :) = 0.3;
colormap(scMap)
xlabel("Sorted GCE")
ylabel("Contour #")

title (animal_path,'Interpreter','none');
Figure.Position = [400 400 400 300];

fig_save_path_cluster_on_detected_events = fullfile(store_folder_path,'cluster_on_detected_events.fig');
saveas(gcf, fig_save_path_cluster_on_detected_events);

fig_save_path_cluster_on_detected_events_jpeg = fullfile(store_folder_path,'cluster_on_detected_events.jpeg');
saveas(gcf, fig_save_path_cluster_on_detected_events_jpeg);
close all
%%
%Plot ROI corresponding to its color-coded cluster
%prepare ROI coordinate
load_suite2p_path = fullfile(data_path,"suite2p/plane0/Fall.mat");
load(load_suite2p_path);%load data from suite2p for ROI coordinate
%load ("./suite2p/plane0/Fall.mat");%load data from suite2p for ROI coordinate

ROI_idx = find(iscell(:,1) == 1);
ROI_ID = ROI_idx-1;
stats_for_ROIs = stat(ROI_idx)';

yx = NaN(length(stats_for_ROIs),2);

  for z = 1:length(yx)
      yx(z,1:2) = stats_for_ROIs{z,1}.med;
  end

xy = zeros(size(yx, 1), size(yx, 2));
xy(:,1) = yx(:,2);
xy(:,2) = yx(:,1);
xy = round(xy);

%Plot ROI figure
clf
Figure = figure('visible','on');
scatter(xy(A0, 1), xy(A0, 2), 'w', 'o', 'filled', 'MarkerEdgeColor', 'k')
  hold on;
  for color = 1:length(C0)
     scatter(xy(C0{color},1), xy(C0{color},2), 40, newCmap(color + 1, :), 'o', 'filled', 'MarkerEdgeColor', 'k');
  end
  labels = {};
  for jdx = 1:NCl+1
      if jdx == 1
         labels{jdx} = "No Cluster";
      else
         labels{jdx} = string(jdx - 1);
      end
  end
 
xlabel("X position (um)")
ylabel("Y position (um)")     
title (animal_path,'Interpreter','none');
Figure.Position = [400 400 400 300];

fig_save_path_cluster_on_ROIs = fullfile(store_folder_path,'cluster_on_ROIs.fig');
saveas (gcf,fig_save_path_cluster_on_ROIs);

fig_save_path_cluster_on_ROIs_jpeg = fullfile(store_folder_path,'cluster_on_ROIs.jpeg');
saveas (gcf,fig_save_path_cluster_on_ROIs_jpeg);
close all 
save (save_path);