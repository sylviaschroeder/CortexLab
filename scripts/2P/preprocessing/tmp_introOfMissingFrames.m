% do this at line 108 of blockReg2P.m
p1=IMG(:,:,4,273:end);
p2=IMG(:,:,1,273:end);
p3=IMG(:,:,2,273:end);
p4=IMG(:,:,3,273:end);
IMG(:,:,1,273)=(IMG(:,:,1,272)+IMG(:,:,4,273))./2;
IMG(:,:,1,274:end)=p1(:,:,1,1:end-1);
IMG(:,:,2,273:end)=p2;
IMG(:,:,3,273:end)=p3;
IMG(:,:,4,273:end)=p4;


% in 16th tiff file, frame 172 in plane1 is missing, which is frame
% 171*4+1=685 in the whole stack (data)
% do this at line 179 of blockReg2P.m
d1=data(:,:,1:4:end);
d4=data(:,:,4:4:end);
insert=(d1(:,:,171)+d4(:,:,172))./2;
drest=data(:,:,685:end);
d=data;
d(:,:,685)=insert;
d(:,:,686:end)=drest(:,:,1:end-1);
rememberR=drest(:,:,end);
rememberG=drest(:,:,end);
data=d;

% for each following tiff file:
% do this at line 226 of blockReg2P.m
data=cat(3,rememberR,data(:,:,1:end));
rememberR=data(:,:,end);
data(:,:,end)=[];

data=cat(3,rememberG,data(:,:,1:end));
rememberG=data(:,:,end);
data(:,:,end)=[];