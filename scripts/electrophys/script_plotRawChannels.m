
% parameter information was taken from params.m
nChInFile  = 384;
filename = 'SS061_20160511_g0_t0.imec_AP_CAR.bin';
dataType = 'int16';

d = dir(filename);
nSamp = d.bytes/2/nChInFile ;

mmf = memmapfile(filename, 'Format', {dataType, [nChInFile nSamp], 'x'});

data = mmf.Data.x(195:203,1:1000*30000);

figure
hold on
for k=1:2:9
    plot((1:size(data,2))/30000, double(data(k,:))./10+1950+(k-1)*10)
end
xlabel('Time (s)')
ylabel('Voltage')
title('Channels 195:2:203')
figure
hold on
for k=2:2:9
    plot((1:size(data,2))/30000,double(data(k,:))./10+1950+(k-1)*10)
end
xlabel('Time (s)')
ylabel('Voltage')
title('Channels 196:2:202')