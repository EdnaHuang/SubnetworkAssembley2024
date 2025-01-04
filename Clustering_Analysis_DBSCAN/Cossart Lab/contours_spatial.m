figure;
imshow(Cn)
hold on;

for i=1:length(aaa)
a=sum(CellCl(:,i));
%no participation in cluster
if a==0
end

if a==1
 bb=find(CellCl(:,i)==1);
 cont = medfilt1(Coor{aaa(i),1}')';
 
 
 if bb==1
 BW1=zeros(d1,d2);
 BW2=zeros(d1,d2);
for j=1:size(cont,2)
 BW1(cont(2,j),cont(1,j))=1;
 BW2=bwperim(BW1);
end
if celltype(1,i)==1
B = bwboundaries(BW2);
visboundaries(B,'Color','r','LineWidth', 0.2,'LineStyle','-.')
else
B = bwboundaries(BW2);
visboundaries(B,'Color','r','LineWidth', 0.2)
end
end

 end
 
  if bb==2
 BW1=zeros(d1,d2);
 BW2=zeros(d1,d2);
for j=1:size(cont,2)
 BW1(cont(2,j),cont(1,j))=1;
 BW2=bwperim(BW1);
end
if celltype(1,i)==1
B = bwboundaries(BW2);
visboundaries(B,'Color','g','LineWidth', 0.2,'LineStyle','-.')
else
B = bwboundaries(BW2);
visboundaries(B,'Color','g','LineWidth', 0.2)
end
 end
 
 
 if bb==3
 BW1=zeros(d1,d2);
 BW2=zeros(d1,d2);
for j=1:size(cont,2)
 BW1(cont(2,j),cont(1,j))=1;
 BW2=bwperim(BW1);
end
if celltype(1,i)==1
B = bwboundaries(BW2);
visboundaries(B,'Color','y','LineWidth', 0.2,'LineStyle','-.')
else
B = bwboundaries(BW2);
visboundaries(B,'Color','y','LineWidth', 0.2)
end
 end
 
end