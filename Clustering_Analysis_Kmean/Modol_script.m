function Modol_script(path)
%%
addpath 'Cossart Lab'
%%
%c = parcluster('local');
%slurm_id = getenv('SLURM_JOB_ID');
%store_loc = strcat("~/.matlab/local_cluster_jobs/R2020a/", slurm_id);
%mkdir(store_loc);
%c.JobStorageLocation = store_loc;
parpool(32)
%%
path_char = char(path);
if path_char(end) ~= '\' & path_char(end) ~= '/'
    path = strcat(path_char, "/");
end
output_path = strcat(path, "Modol_outputs/");
if ~exist(output_path, 'dir')
   mkdir(output_path)
end
%%
detected_events = readmatrix(strcat(path, 'detected_events.xlsx')) ~= 0;
%%
try
    start_time = readmatrix(strcat(path, "start_time.xlsx"));
    end_time = readmatrix (strcat(path, "end_time.xlsx"));
catch ME
    error("There is no start/end time data")
    %start_time = 1*0.0656;
    %end_time = length(detected_events)*0.0656;
end
%%
Race = detected_events(:,floor(start_time/0.0656):floor(end_time/0.0656) - 1);
%%
[NCell,NRace] = size(Race);
CovM = CovarM(Race);
[IDX2,sCl,~,~] = kmeansopt(CovM,10,'precomp');
% Number clusters
NCl = max(IDX2);

[~,x2] = sort(IDX2);
% Sorted normalized covariance matrix
MSort = CovM(x2,x2);

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
%%
figure
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
SSS='Clusters.fig';
St=strcat(output_path,SSS);

savefig(St)
%%
SSS='Clusters.mat';
St=strcat(output_path,SSS);
save(St,'IDX2')
%%
profile on
NTrials = 1000;
sClrnd = zeros(1,NTrials);
for i = 1:NTrials
    disp(i)
%     rnd_ind = randperm(size(CovM,1));
%     M_swap = CovM(rnd_ind,:);
%     M_swap = M_swap(:,rnd_ind);
    sClrnd(i) = kmeansoptrnd(Race,10,NCl,'var');
end

NClOK = sum(sCl>prctile(sClrnd(1:NTrials-1),95));
sClOK = sCl(1:NClOK)';

save(St,'NClOK', 'sClOK', '-append')

RaceOK = Race(:,IDX2<=NClOK);
NRaceOK = size(RaceOK,2);
profile off
%%
NShuf = 5000;
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
    %         ThMin = mean(Random) - 2*std(Random);
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
save(St,'CellParticip', 'C0', '-append')
%%
NCl = length(C0);
if ~NCl
    NCl = 0;
    SSS='all.mat';
    SSS=strcat(output_path,SSS);
    save(SSS,'-v7.3');
    exit()
    % error('There were no significant clusters found!! Cannot run this cell...')
end
[NCell,NRace] = size(Race);
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
    %         ThMin = mean(Random) - 2*std(Random);
    RClr = sort(RClr,2);
    %         ThMin = mean(Random) - 2*std(Random);
    %Proba above 95th percentile
    ThMax = RClr(:,round(NShuf*(1-0.05/NCl)));
    for i = 1:NCl
        PCl(i,j) = double(RCl(i,j)>ThMax(i));% - double(RCl(:,j)<ThMin);
    end
    %Normalize (probability)
    RCln(:,j) = RCl(:,j)/sum(Race(:,j));
end
save(St, 'PCl', '-append');
%%
if ~NCl
    disp('There were no significant clusters found!! Cannot run this cell...')
end
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
figure
imagesc(Race(x1,RList) ~= 0)
% imagesc(Race(x1,RList))
% colormap hot
colormap(flipud(gray))
title('Sorted Rastermap (with significant GCEs)')
% colorbar

SSS='CellCluster.fig';
St=strcat(output_path,SSS);
savefig(St)
close all;

SSS='all.mat';
SSS=strcat(output_path,SSS);
save(SSS,'-v7.3');
