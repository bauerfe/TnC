function plot_results(B, R_litter, R_litter_sur, R_microbe, R_bacteria, R_ew, VOL, BfixN, Min_N, Min_P, RmycAM, RmycEM, ...
    N2flx, NH4_Uptake, NO3_Uptake, P_Uptake, K_Uptake, LEAK_NH4, LEAK_NO3, LEAK_P, LEAK_K, LEAK_DOC, ...
    LEAK_DON, LEAK_DOP, Lk, Zbio, rsd)

    rsd = mean(rsd);
    ZBIOG = 0.001 * Zbio;
    
    figure(106)
    subplot(3,1,1)
    plot(B(:,1),'b','LineWidth', 1.5);
    hold on; grid on;  ylabel('gC/m^2')
    plot(B(:,2),'g','LineWidth', 1.5);
    plot(B(:,3),'r','LineWidth', 1.5);
    legend('Ab. Met','Ab Str','Ab Str Lig')
    subplot(3,1,2)
    plot(B(:,6),'b','LineWidth', 1.5);
    hold on; grid on;  ylabel('gC/m^2')
    plot(B(:,7),'g','LineWidth', 1.5);
    plot(B(:,8),'r','LineWidth', 1.5);
    legend('Be. Met','Be Str','Be Str Lig')
    subplot(3,1,3)
    plot(B(:,4),'b','LineWidth', 1.5);
    hold on; grid on;  ylabel('gC/m^2')
    plot(B(:,5),'g','LineWidth', 1.5);
    legend('Ab. Wood','Ab Wood Lign')
    
    figure(107)
    subplot(2,2,1)
    plot(B(:,9),'g','LineWidth', 1.5);
    hold on; grid on;
    plot(B(:,10),'k','LineWidth', 1.5);
    plot(B(:,11),'y','LineWidth', 1.5);
    legend('SOM POC Lign','SOM POC - Cell','SOM MOC')
    subplot(2,2,2)
    plot(B(:,12),'k','LineWidth', 1.5);
    hold on; grid on;
    plot(B(:,13),'r','LineWidth', 1.5);
    legend('DOC-B','DOC-F')
    subplot(2,2,3)
    plot(B(:,18),'g','LineWidth', 1.5);
    hold on; grid on;
    plot(B(:,19),'k','LineWidth', 1.5);
    plot(B(:,20),'r','LineWidth', 1.5);
    plot(B(:,21),'b','LineWidth', 1.5);
    plot(B(:,22),'y','LineWidth', 1.5);
    title('Carbon Pool')
    xlabel('Days'); ylabel('gC/m^2')
    legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Earthworms')
    subplot(2,2,4)
    plot(B(:,14),'g','LineWidth', 1.5);
    hold on; grid on;
    plot(B(:,15),'k','LineWidth', 1.5);
    plot(B(:,16),'r','LineWidth', 1.5);
    plot(B(:,17),'b','LineWidth', 1.5);
    title('Carbon Pool')
    xlabel('Days'); ylabel('gC/m^2')
    legend('EM-POC-B','EM-POC-F','EM-MOC-B','EM-MOC-F')
    
    figure(108)
    subplot(2,2,1)
    plot(B(:,23),'b','LineWidth', 1.5);
    hold on; grid on;  ylabel('gN/m^2')
    plot(B(:,24),'g','LineWidth', 1.5);
    plot(B(:,25),'r','LineWidth', 1.5);
    legend('Ab. Lit','Ab Wod','Be Lit')
    subplot(2,2,2)
    plot(B(:,26),'g','LineWidth', 1.5);
    hold on; grid on;
    legend('SOM')
    subplot(2,2,3)
    plot(B(:,27),'g','LineWidth', 1.5);
    hold on; grid on;
    plot(B(:,28),'k','LineWidth', 1.5);
    plot(B(:,29),'r','LineWidth', 1.5);
    plot(B(:,30),'b','LineWidth', 1.5);
    title('Nitrogen Pool')
    xlabel('Days'); ylabel('gN/m^2')
    legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza')
    subplot(2,2,4)
    plot(B(:,31),'g','LineWidth', 1.5);
    hold on; grid on;
    plot(B(:,32),'k','LineWidth', 1.5);
    plot(B(:,33),'b','LineWidth', 1.5);
    title('Nitrogen Pool')
    xlabel('Days'); ylabel('gN/m^2')
    legend('NH4+ ','NO3-','DON')
    
    figure(109)
    subplot(3,2,1)
    plot(B(:,35),'b','LineWidth', 1.5);
    hold on; grid on;  ylabel('gP/m^2')
    plot(B(:,36),'g','LineWidth', 1.5);
    plot(B(:,37),'r','LineWidth', 1.5);
    legend('Ab. Lit','Ab Wod','Be Lit')
    subplot(3,2,2)
    plot(B(:,38),'g','LineWidth', 1.5);
    hold on; grid on;
    legend('SOM')
    subplot(3,2,3)
    plot(B(:,39),'g','LineWidth', 1.5);
    hold on; grid on;
    plot(B(:,40),'k','LineWidth', 1.5);
    plot(B(:,41),'r','LineWidth', 1.5);
    plot(B(:,42),'b','LineWidth', 1.5);
    title('Phosporus Pool')
    xlabel('Days'); ylabel('gP/m^2')
    legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza')
    subplot(3,2,4)
    plot(B(:,43),'g','LineWidth', 1.5);
    hold on; grid on;
    plot(B(:,47),'r','LineWidth', 1.5);
    title('Phosporus Pool')
    xlabel('Days'); ylabel('gP/m^2')
    legend('Mineral','DOP')
    subplot(3,2,5)
    plot(B(:,44),'k','LineWidth', 1.5);
    hold on; grid on;
    title('Phosporus Pool')
    xlabel('Days'); ylabel('gP/m^2')
    legend('Primary Material')
    subplot(3,2,6)
    plot(B(:,46),'k','LineWidth', 1.5);
    hold on; grid on;
    plot(B(:,45),'g','LineWidth', 1.5);
    xlabel('Days'); ylabel('gP/m^2')
    legend('Occluded','Secondary')
    
    
    figure(110)
    subplot(3,2,1)
    plot(B(:,48),'b','LineWidth', 1.5);
    hold on; grid on;  ylabel('gK/m^2')
    plot(B(:,49),'g','LineWidth', 1.5);
    plot(B(:,50),'r','LineWidth', 1.5);
    legend('Ab. Lit','Ab Wod','Be Lit')
    subplot(3,2,2)
    plot(B(:,51),'g','LineWidth', 1.5);
    hold on; grid on;
    legend('SOM')
    subplot(3,2,3)
    plot(B(:,52),'g','LineWidth', 1.5);
    hold on; grid on;
    title('Potassium Pool')
    xlabel('Days'); ylabel('gK/m^2')
    legend('Mineral Solution ')
    subplot(3,2,4)
    plot(B(:,53),'g','LineWidth', 1.5);
    hold on; grid on;
    plot(B(:,54),'k','LineWidth', 1.5);
    title('Potassium Pool')
    xlabel('Days'); ylabel('gK/m^2')
    legend('Excheangeable ','Non-Excheangeable')
    subplot(3,2,5)
    plot(B(:,55),'g','LineWidth', 1.5);
    hold on; grid on;
    legend('Primary Minerals')
    
    
    
    figure(111)
    plot((B(:,2)+B(:,3)+B(:,1))./B(:,23),'g','LineWidth', 1.5);
    hold on; grid on;
    plot((B(:,5)+B(:,4))./B(:,24),'m','LineWidth', 1.5);
    plot((B(:,6)+B(:,7)+B(:,8))./B(:,25),'k','LineWidth', 1.5);
    plot((B(:,9)+B(:,10)+B(:,11))./B(:,26),'b','LineWidth', 1.5);
    plot(B(:,18)./B(:,27),'r','LineWidth', 1.5);
    plot(B(:,19)./B(:,28),'y','LineWidth', 1.5);
    plot(B(:,20)./B(:,29),'c','LineWidth', 1.5);
    plot(B(:,21)./B(:,30),'Color',[0.168 0.50586 0.3372],'LineWidth', 1.5);
    plot(B(:,22)./B(:,34),'Color',[0.06 0.7 0.6],'LineWidth', 1.5);
    title('C:N Ratio')
    xlabel('Days'); ylabel('C:N')
    legend('AG Litter','AG Wood','BG Litter','SOM','Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Earthworms')
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    figure(112)
    plot((B(:,2)+B(:,3)+B(:,1))./B(:,35),'g','LineWidth', 1.5);
    hold on; grid on;
    plot((B(:,5)+B(:,4))./B(:,36),'m','LineWidth', 1.5);
    plot((B(:,6)+B(:,7)+B(:,8))./B(:,37),'k','LineWidth', 1.5);
    plot((B(:,9)+B(:,10)+B(:,11))./B(:,38),'b','LineWidth', 1.5);
    plot(B(:,18)./B(:,39),'r','LineWidth', 1.5);
    plot(B(:,19)./B(:,40),'y','LineWidth', 1.5);
    plot(B(:,20)./B(:,41),'c','LineWidth', 1.5);
    plot(B(:,21)./B(:,42),'Color',[0.168 0.50586 0.3372],'LineWidth', 1.5);
    title('C:P Ratio')
    xlabel('Days'); ylabel('C:P')
    legend('AG Litter','AG Wood','BG Litter','SOM','Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza')
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    figure(113)
    plot((B(:,2)+B(:,3)+B(:,1))./B(:,48),'g','LineWidth', 1.5);
    hold on; grid on;
    plot((B(:,5)+B(:,4))./B(:,50),'m','LineWidth', 1.5);
    plot((B(:,6)+B(:,7)+B(:,8))./B(:,37),'k','LineWidth', 1.5);
    plot((B(:,9)+B(:,10)+B(:,11))./B(:,51),'b','LineWidth', 1.5);
    title('C:K Ratio')
    xlabel('Days'); ylabel('C:K')
    legend('AG Litter','AG Wood','BG Litter','SOM')
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(114)
    subplot(2,2,1)
    plot(R_litter,'r','LineWidth', 1.5);
    hold on; grid on;
    plot(R_litter-R_litter_sur,'m','LineWidth', 1.5);
    plot(R_microbe,'b','LineWidth', 1.5);
    title('Respiration Het.')
    plot(R_ew,'g','LineWidth', 1.5);
    xlabel('Days'); ylabel('[gC/m2 day]')
    legend('Litter','Litter below','Microbe','Earthworms')
    subplot(2,2,2)
    plot(VOL,'r','LineWidth', 1.5);
    hold on; grid on;
    plot(N2flx,'b','LineWidth', 1.5);
    title('N - Fluxes')
    xlabel('Days'); ylabel('[gN/m2 day]')
    legend('NH_4 Vol.','N_2')
    subplot(2,2,3)
    plot(Min_N,'k','LineWidth', 1.5);
    hold on ;  grid on 
    title('N - Fluxes')
    xlabel('Days'); ylabel('[gN/m2 day]')
    legend('Min-N')
    subplot(2,2,4)
    plot(Min_P,'k','LineWidth', 1.5);
    hold on ;  grid on 
    title('P - Fluxes')
    xlabel('Days'); ylabel('[gP/m2 day]')
    legend('Min-P')
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    figure(115)
    subplot(3,1,1)
    plot(NH4_Uptake,'r','LineWidth', 1.5);
    hold on; grid on;
    plot(NO3_Uptake,'b','LineWidth', 1.5);
    title('N Uptake.')
    xlabel('Days'); ylabel('[gN/m2 day]')
    legend('NH_4','NO_3')
    subplot(3,1,2)
    plot(P_Uptake,'g','LineWidth', 1.5);
    hold on; grid on;
    title('P Uptake')
    xlabel('Days'); ylabel('[gP/m2 day]')
    subplot(3,1,3)
    plot(K_Uptake,'m','LineWidth', 1.5);
    hold on; grid on;
    title('K Uptake')
    xlabel('Days'); ylabel('[gK/m2 day]')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure(116)
    subplot(2,2,1)
    plot(LEAK_DOC,'b','LineWidth', 1.5);
    hold on; grid on;
    title('DOC Leaching')
    xlabel('Days'); ylabel('[gC/m2 day]')
    subplot(2,2,2)
    plot(LEAK_NH4,'r','LineWidth', 1.5);
    hold on; grid on;
    plot(LEAK_NO3,'b','LineWidth', 1.5);
    hold on; grid on;
    plot(LEAK_DON,'m','LineWidth', 1.5);
    title('N Leaching')
    xlabel('Days'); ylabel('[gN/m2 day]')
    legend('NH_4','NO_3','DON')
    subplot(2,2,3)
    plot(LEAK_P,'g','LineWidth', 1.5);
    hold on; grid on;
    plot(LEAK_DOP,'m','LineWidth', 1.5);
    title('P Leaching')
    legend('P','PON')
    xlabel('Days'); ylabel('[gP/m2 day]')
    subplot(2,2,4)
    plot(LEAK_K,'m','LineWidth', 1.5);
    hold on; grid on;
    title('K Leaching')
    xlabel('Days'); ylabel('[gK/m2 day]')
    %%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    Lkday=Lk;
    figure(120)
    subplot(2,2,1)
    plot(LEAK_DOC./(Lkday)*1000,'b','LineWidth', 1.5);
    hold on; grid on;
    title('DOC Conc.')
    xlabel('Days'); ylabel('[mg/l]')
    subplot(2,2,2)
    plot(LEAK_NH4./(Lkday)*1000,'r','LineWidth', 1.5);
    hold on; grid on;
    plot(LEAK_NO3./(Lkday)*1000,'b','LineWidth', 1.5);
    hold on; grid on;
    plot(LEAK_DON./(Lkday)*1000,'m','LineWidth', 1.5);
    title('N Conc.')
    xlabel('Days'); ylabel('[mg/l]')
    legend('NH_4','NO_3','DON')
    subplot(2,2,3)
    plot(LEAK_P./(Lkday)*1000*1000,'g','LineWidth', 1.5);
    hold on; grid on;
    plot(LEAK_DOP./(Lkday)*1000*1000,'m','LineWidth', 1.5);
    title('P Conc.')
    legend('P','DOP')
    xlabel('Days');  ylabel('[ug/ l]')
    subplot(2,2,4)
    plot(LEAK_K./(Lkday)*1000,'m','LineWidth', 1.5);
    hold on; grid on;
    title('K Conc.')
    xlabel('Days'); ylabel('[mg/ l]')
    
    
    figure(121)
    subplot(3,2,1)
    plot((B(:,18)+B(:,19)+B(:,20) + B(:,21))./(B(:,9)+B(:,10)+B(:,11)+B(:,12)+B(:,13)),'b','LineWidth', 1.5);
    hold on; grid on;  ylabel('[-]')
    legend('Microbial/Substrate')
    subplot(3,2,2)
    plot((B(:,22))./(B(:,18)+B(:,19)+B(:,20) + B(:,21)),'b','LineWidth', 1.5);
    hold on; grid on;  ylabel('[-]')
    legend('Earthworms/Microbial')
    subplot(3,2,3)
    plot(sum(B(:,6:21),2)./(ZBIOG*rsd),'b','LineWidth', 1.5);
    hold on; grid on;  ylabel('SOC [gC /kg soil]')
    subplot(3,2,4)
    plot(sum(B(:,6:21),2),'b','LineWidth', 1.5);
    hold on; grid on;  ylabel('SOC [gC /m^2]')
    subplot(3,2,5)
    plot(((B(:,19) + B(:,20) + B(:,21))./B(:,18)),'b','LineWidth', 1.5);
    hold on; grid on;  ylabel('[-]')
    legend('Fungi/Bacteria')
    subplot(3,2,6)
    plot( B(:,19)./(B(:,20)+B(:,21)),'b','LineWidth', 1.5);
    hold on; grid on;  ylabel('[-]')
    legend('Saprotrophic/Mycorrhiza')
    
    Rlitter_sub= R_litter-R_litter_sur; 
    rb = R_bacteria./R_microbe; 
    rf = 1 -rb; 
    TSR = Rlitter_sub + R_microbe + R_ew; 
    figure(122)
    subplot(1,1,1)
    plot((R_bacteria + rb.*Rlitter_sub)./TSR,'r','LineWidth', 1.5);
    hold on; grid on;
    plot((rf.*R_microbe + rf.*Rlitter_sub)./TSR,'b','LineWidth', 1.5);
    title('Respiration Het.')
    plot(R_ew./TSR,'g','LineWidth', 1.5);
    xlabel('Days'); ylabel('[%]')
    legend('Bacteria','Fungi','Earthworms')
end