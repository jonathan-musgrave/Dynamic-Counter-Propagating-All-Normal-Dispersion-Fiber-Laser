% Calculate the Carrier frequency of a [nw,nz] array usingthe COM of the
% intensity profile. Returns an [1,nz] array of the center frequencies
function [w0] = calc_CarrierFreq(UU,w)
    [nw,nz] = size(UU);
    w0 = zeros(1,nz);
    for i = 1:nz
       w0(i) = calc_COM(w',abs(UU(:,i)).^2);
    end
end