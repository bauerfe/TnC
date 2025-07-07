
function [SLAdyn] = DynamicSLA(Sl, dflotm1,curve,boost)
    %curve = 15; % to set the slope 
    %boost = 3; % boost when leaf age is 0
    %dfloall = 0:300;
    SLAdyn = Sl .* (1 + boost * exp(-dflotm1 / curve));
end


