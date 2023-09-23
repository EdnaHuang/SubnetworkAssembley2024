numClus=max(IDX2);
for i=1:numClus
    clu=find(IDX2==i);
    bla=[];
    for j=1:length(clu-1)
       bla(:,:,end+1:end+4)=Y(:,:,round(locs(clu(j))/2):round(locs(clu(j))/2)+3);
       if j==length(clu-1)
           bla=bla(:,:,2:end);
           avgClus{i}=max(bla,[],3);
       end
    end
end


celltype=[];
d1=length(Cn(1,:));
d2=length(Cn(:,1));
figure('units','normalized','outerposition',[0 0 1 1]);
imshow(mat2gray((avgClus{2})))
hold on;

BW4=zeros(d1,d2);
for i=1:length(CellCl(1,:))
a=sum(CellCl(:,i));
%no participation in cluster
if a==0
    celltype(1,i)=0;
else
     bb=find(CellCl(:,i)==1);
     imshow(mat2gray((avgClus{bb})))
 cont = medfilt1(Coor{i}')';
BW1=zeros(d1,d2);
 BW2=zeros(d1,d2);
for j=1:size(cont,2)
 BW1(cont(2,j),cont(1,j))=1;
    se1 = strel('line',5,90);
    se2 = strel('line',5,0);
    BWtemp=imdilate(BW1,[se1,se2]);
        BWtemp= imfill(BWtemp,'holes');
 BW2=bwperim(BWtemp);
end
B = bwboundaries(BW2);
visboundaries(B,'Color','r','LineWidth', 0.0005)

B2=bwboundaries(BW4);
visboundaries(B2,'LineWidth',0.2,'Color','g')
BW4=BW1+BW4;

[x,y,z] = ginput(1);

celltype(1,i)=z;

end
end
