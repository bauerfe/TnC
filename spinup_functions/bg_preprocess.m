%! Pre-process some variables

function [Se, Se_fc, Psi_s, Ts, V, VT] = preprocess(Ta, Tdp, O, V, Soil_Param, Phy, SPAR, Bio_Zs)
    jj=0;
    for i=1:24:length(Ta)
        jj=jj+1;
        pdind = [max(1,i-24):i-1];
        [Se_bio(jj),Se_fc,Psi_bio(jj),Tdp_bio(jj),VSUM(jj),VTSUM(jj)]=Biogeo_environment(Tdp(pdind,:),O(pdind,:),V(pdind,:),Soil_Param,Phy,SPAR,Bio_Zs);
        Aew(jj) = 0.000008575*exp(11.67*Se_bio(jj)./Se_fc)*(Se_bio(jj)>0.2);
    end
    Aew(Aew>1)=1;
    Se = (log(nanmean(Aew)/0.000008575)/11.67)*Se_fc;
    Psi_s = nanmean(Psi_bio);
    Ts = nanmean(Tdp_bio);
    V = nanmean(VSUM);
    VT = nanmean(VTSUM);
end