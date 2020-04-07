% run kilosort.m

obj.P.chanMap.connected(obj.P.chanMap.ycoords<600) = false;
obj.P.chanMap.connected(obj.P.chanMap.ycoords>1300) = false;
obj.ops.Th=[8 4];
% obj.ops.mergeThreshold=.25; %.1;
% obj.ops.ccsplit=.97; %.99;
% obj.ops.lam = 100; %120;