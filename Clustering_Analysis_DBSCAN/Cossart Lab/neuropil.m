d1=length(Cn(1,:));
d2=length(Cn(:,1));

BW1=zeros(d1,d2);
BW2=zeros(d1,d2);


for i=1:length(TraceCellsOld(:,1))
 cont = medfilt1(Coor{i}')';
 for j=1:size(cont,2)
 BW1(cont(2,j),cont(1,j))=1;
 BW2=bwperim(BW1);
 end
 BW3{i}=imfill(BW2,'holes');
        BW1=zeros(d1,d2);
        BW2=zeros(d1,d2);
end

BWupdated=zeros(d1,d2);
BWtemp2={};
for i=1:length(TraceCellsOld(:,1))
    se1 = strel('line',35,90);
    se2 = strel('line',35,0);
    BWtemp=imdilate(BW3{i},[se1,se2]);
    BWtemp2{i}=BWtemp-BW3{i};
    BWupdated=[BW3{i}+BWupdated];
end

BWupdated(BWupdated>0)=1;

BWLoc={};
for i=1:length(TraceCellsOld(:,1))
BWLoc{i}=BWupdated.*BWtemp2{i};
end


for i=1:length(TraceCellsOld(:,1))
    wts{i}=find(BWLoc{i}==1);
end

nonz=zeros(1,length(Y(1,1,:)));

for i=1:length(TraceCellsOld(:,1))
t=wts{i};
if sum(t)==0
    neupil(i,:)=nonz;
else
tcel=[];
[I1,J1] = ind2sub([d1,d2],t);
tcel=zeros(length(I1),length(Y(1,1,:)));
for j=length(I1)
tcel(j,:)=Y(I1(j),J1(j),:);
end
    neupil(i,:)=sum(tcel)/length(tcel(:,1));
end
end
    

att2={};
parfor i = 1:length(TraceCellsOld(:,1))
if sum(neupil(i,:))>1
SAMPLES = cont_ca_sampler(neupil(i,:));    %% MCMC  
att2{i}=SAMPLES.params.spiketimes_;
else
    att2{i}=NaN;
end
end
    
    
d1=length(TraceCellsOld(:,1));
d2=2*length(neupil(1,:));
spikenums2=zeros(d1,d2);
for i = 1:length(TraceCellsOld(:,1))
    ac=att2{i};
    if length(ac)>1
    acc = 0.5;
    result = round(ac/acc)*acc;
    result=result*2;
    for j=1:length(result)
        spikenums2(i,result(j))=1;
    end
    end
end

for i=1:length(locs)
    for j=1:length(att2)
        axxx=spikenums2(j,locs(i)-3:locs(i)+3);
        ax=min(find(axxx==1));
        if ax>=1
            t1_2(i,j)=ax;
        end
    end
end

t11_2=t1-4;
t666=t11(1:length(TraceCellsOld(:,1)),:)-t11_2;

