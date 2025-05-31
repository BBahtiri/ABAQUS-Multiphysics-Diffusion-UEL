%% Script developed to read and modify Abaqus input files %%

clear variables
clc

% Open the file
fid = fopen('CT.inp','r');

file=textscan(fid,'%s','Delimiter','\n');
%fclose(fid);
flines=(file{1});

lif=length(flines); 
%fclose(fid);

% Search the lines containing *ELEMENT
FFLINES=upper(flines);
Lin_element_cells=strfind(FFLINES,['*' upper('Element')]); 
Lin_element_zeros=cellfun(@isempty,Lin_element_cells); 
start_elements=find(Lin_element_zeros==0); %numbers of the lines that contains *Element

Lin_order_cells=strfind(FFLINES,'*');
Lin_order_zeros=cellfun(@isempty,Lin_order_cells); 
Lin_orders=find(Lin_order_zeros==0); %numbers of the lines that contains *Element

for i=1:length(start_elements)
 end_elements(i)=Lin_orders(find(Lin_orders==start_elements(i))+1);
 p=1;
 str_a=FFLINES{start_elements(i)}(~isspace(FFLINES{start_elements(i)}));%eliminate spaces
 flag=0;
 if length(str_a)==(length('Element')+1)
  flag=1;    
 else
  if str_a(length('Element')+2)==','
   flag=1;
  end
 end

 if flag==1
  for b=(start_elements(i)+1):(end_elements(i)-2) 
   first=strread(FFLINES{b}, '%f', 'delimiter', ',');
   b=b+1;
   second=strread(FFLINES{b}, '%f', 'delimiter', ',');
   Element{i}(p,:)=vertcat(first,second);
   p=p+1;  
  end
  Element{i}(2:2:end,:) = [];
 end
end

LE=Element{1}; % Element numbering and connectivity
[NUME, ~] = size(LE);
Order=floor(log10(NUME));
offset=10*10^Order;
LEnew=zeros(NUME,21);
for i=1:NUME
 LEnew(i,:)= [(offset+i) LE(i,2) LE(i,3) LE(i,4) LE(i,5) LE(i,6) LE(i,7)...
     LE(i,8) LE(i,9) LE(i,10) LE(i,11) LE(i,12) LE(i,13) LE(i,14) LE(i,15) LE(i,16)...
     LE(i,17) LE(i,18) LE(i,19) LE(i,20) LE(i,21)];   
end    
for i=1:length(Element)
 % A new file is created so as to avoid overwriting the previous one
 NameFile=strcat('VisualMesh.inp');
 % We open the new input file and start making changes
 dlmwrite(NameFile,LEnew,'precision', 8);
end

fclose all;