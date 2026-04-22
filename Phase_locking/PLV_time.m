%Ricardo Bruña et al 2018 J. Neural Eng. 15 056011
function [plv] = PLV_time(eegdata)
% - eegdata: channels x samples x trials(time windows)
% disp(size(eegdata))
[nc, ns, nt] = size(eegdata); 
ndat = eegdata./ abs(eegdata); 
plv = zeros(nc, nc, nt); 
    for t = 1: nt 
        plv(:,:, t) = abs(ndat(:,:, t) * ndat(:,:, t)') / ns; 
    end
end

%another method
% [nc, ns, nt] = size(eeg_data_timewindow); 
% phs = angle(eeg_data_timewindow);  
% plv = zeros(nc, nc, nt); 
% for t = 1: nt 
%     tplv = complex(zeros(nc)); 
%     for s = 1: ns 
%         dphs = bsxfun(@minus, phs(:, s, t), phs(:, s, t)'); 
%         tplv = tplv + exp(1i * dphs); 
%     end
%     plv(:,:, t) = abs(tplv / ns);
% end