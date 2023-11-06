function Modol_script_uniform(sim_type, path)
%%
addpath ('Cossart Lab');
addpath("DependentFurntion_BCT_Modol");
%%
parpool(48)
%%
path_char = char(path);
if path_char(end) ~= '\' & path_char(end) ~= '/'
    path = strcat(path_char, "/");
end
%%
detected_events = readmatrix(strcat(path, 'detected_events.xlsx')) ~= 0;
%%
try
    start_time = readmatrix(strcat(path, "start_time.xlsx"));
    end_time = readmatrix (strcat(path, "end_time.xlsx"));
catch ME
    error("There is no start/end time data")
end
%%
Race = detected_events(:,floor(start_time/0.0656):floor(end_time/0.0656) - 1);

%% New Workflow

% Assigning Type of Similarity Matrix
if strcmp (sim_type, 'Cosine') 
    %Similarity Matrix for community detection
    rho =  CosineM(Race);
    nreps = 100;
    %make the folder in the current working directory to save all outpus
    folder_name= 'Cosine_CommDetectUniform_Modol';
    output_path = strcat(path, folder_name, "/");
    if ~exist(output_path, 'dir')
       mkdir(output_path)
    end         
else
    rho = CovarM(Race);
    nreps = 500;
    %make the folder in the current working directory to save all outpus
    folder_name= 'Cov_CommDetectUniform_Modol';
    output_path = strcat(path, folder_name, "/");
    if ~exist(output_path, 'dir')
       mkdir(output_path)
    end 
end

%%identify current wowrking folder for naming file
[root_path,animal_path,~] = fileparts(pwd); %pwd: display the current folder
animal_ID = split(animal_path, '_');

%Community detection--uniform null model
%Ref: compares observed matrix with matrix where all cells are considered to be correlated at magnitude equal to gamma
% (see Bazzi et al 2015, Community detection in temporal ...)

%gamma is a parameter that one can vary--for community size
%In our case, it is difficult to decide which gamma should, so we run different gamma value (gammavals)
gammavals = linspace(0,1,11); %use different gamma 
%gammavals = linspace(0,0.2,3); %shorter for testing 
<<<<<<< HEAD
% nreps = 500; the repeat # for each gamma value used for community_louvain function, should run several iterations
=======
nreps = 100; % the repeat # for each gamma value used for community_louvain function, should run several iterations
>>>>>>> upstream/main


%Pre-allocate matrix / array for data storage 
% store the ci from each given gamma value (each vector) and repeat (colume in vector) into the array.
ci_all = cell(1,length(gammavals)); 
% store the AgreementMatrix from each given gamma (each vector) into the array.
Agreement_all = cell(1,length(gammavals));
% store the q from each given gamma value (row) and repeat (colume) into the array.
q_all = zeros(length(gammavals),nreps);
% store the agreement matrix from all given gamma value and their repeat.
d = zeros(length(rho));
% store the agreement matrix from a given gamma value for its repeat, temp for each given gamma value.
Agreement_temp = zeros(length(rho));
disp(nreps)
gcp; %Activate the parallel pool
for igamma = 1:length(gammavals) % use different gamma for community_louvain
    disp(igamma)    
    gamma = gammavals(igamma);
    b = (rho - gamma).*~eye(length(rho));
    
    %tic
    
    ci = zeros(length(rho),nreps);
    q = zeros(1,nreps);
    parfor irep = 1:nreps
           [ci(:,irep),q(:,irep)] = community_louvain(ones(size(rho)),[],[],b);
    end
    
    % calculate the agreement matrix for the repeat of given gamma value
    Agreement_temp = Agreement_temp + agreement(ci);  
    
    %toc    
    % generate agreement matrix for all repeat from all given gamma value. 
    d = d + agreement(ci); 
    
    %save data from the given gamma value into array/vector            
    ci_all{igamma} = ci;
    q_all(igamma,:) = q;
    Agreement_all{igamma} = Agreement_temp;                      
end

%% using the concensus matrix
%normalize the concensus matrix into probability
d_sum = sum (d,"all");
d_Pro = d/d_sum;
d_Pro_max = max (d_Pro,[], 'all');
% The value of d_Pro is small, thus using tau = 0, not setting threshold
% Using the ciu from consensus_und function to assign community
ciu = consensus_und(d_Pro,0,10);

%Plot the grid community by using the function in BCT 
clf
figure('visible','on');
[gx,gy,idxorder] = grid_communities(ciu); % reorder matrix by community
imagesc(rho(idxorder,idxorder));
hold on;
plot(gx,gy,'k');
colormap(fcn_cmapjet);
colorbar;

saveas (gcf,fullfile(output_path,'Grid_community.fig'));
saveas (gcf,fullfile(output_path,'Grid_community.jpeg'));
save (fullfile(output_path, 'CommDetectUniform.mat'));
close all
%%
%% Below start using the output from community detection and analyzed by the rest of Model code
%assign the variable for further step
[NCell,NRace] = size(Race); %total cell/neuron number; total frame number
IDX2 = ciu;
IDX2 = transpose (IDX2);

% communit/cluster number
NCl= max(IDX2); % this is the initial community / cluster number
%%
%Calculate silhValue of each cluster/community after community detection before statistic test
%the silhV of all, using the silh function from Modol code
silhV_all = silh(rho,IDX2); % Replaced CovM with rho AB - 03/24/2022  

sCl_all = zeros(1,NCl); % to store the silV of each cluster
for i = 1:NCl
    sCl_all(i) = median(silhV_all(IDX2==i));
end

[sCl_all,xCl] = sort(sCl_all,'descend'); %sort the silhValue of all cluster/community
sCl = rmmissing (sCl_all); % remove NaN
NCl_beforeStat = length(sCl); % This is the cluster that has silhouette value
%%
%Taking the index output from community detection and use that follow the rest of Modol code
[~,x2] = sort(IDX2);
% Sorted normalized covariance matrix
MSort = rho(x2,x2);% Replaced CovM with rho AB - 03/24/2022  

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

saveas (gcf,fullfile(output_path,'Initial_cluster.fig'));
close all
%%
%%For statistical significance
NShuf = 5000; %for real analysis

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
    save (fullfile(output_path,'all.mat'))
    exit()
    % error('There were no significant clusters found!! Cannot run this cell...')
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
   save (fullfile(output_path,'all.mat'));% Edited AB - 03/24/2022
   disp('There were no significant clusters found!! Cannot run this cell...')
    
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

saveas(gcf, fullfile(output_path, 'Sorted_Rastermap.fig'));% Edited AB - 03/24/2022
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

saveas(gcf, fullfile(output_path, 'cluster_on_detected_events.fig'));% Edited AB - 03/24/2022
saveas(gcf, fullfile(output_path, 'cluster_on_detected_events.jpeg'));% Edited AB - 03/24/2022
close all
save (fullfile(output_path,'all.mat')); % Edited AB - 03/24/2022
