
fc=9*10^8; %carrier frequence

[Ae_Ru4,Ae_Ru6,Ae_Ru8]  =pass_loss(fc);

aAe=mean(Ae_Ru4,2); %the average of Ae after running 10000 times
aAeb=mean(Ae_Ru6,2); %the average of the best case after running 10000 times
aAew=mean(Ae_Ru8,2); %the average of the worst case after running 10000 times

fc2g = 2*10^9;

[Ae_Ru4_2g,Ae_Ru6_2g,Ae_Ru8_2g]  =pass_loss(fc2g);

aAe_2g=mean(Ae_Ru4_2g,2); %the average of Ae after running 10000 times
aAeb_2g=mean(Ae_Ru6_2g,2); %the average of the best case after running 10000 times
aAew_2g=mean(Ae_Ru8_2g,2); %the average of the worst case after running 10000 times


Ru = 0.1:0.05:1;
figure()

semilogy(Ru,aAe,'k',Ru,aAe_2g,'b--',Ru,aAeb,'g',Ru,aAew,'r',Ru,aAeb_2g,'b--o',Ru,aAew_2g,'c--'); %set the best case is green,the worst line is red

legend('fc =900MHZ ','fc = 2GHZ')
xlabel('Normalized Reuse Distance Ru');grid;
ylabel('ASE[Bits/Sec/Hz/Km^2]');
title('ASE in Uplink Transmission(a=2,b=2,R=200m),');



function [Ae_Ru4,Ae_Ru6,Ae_Ru8] = pass_loss(frequency)
            
    %function [Ae] = pass_loss(frequency)
    R=0.2; %the cell radius(KM
    R0=0.02; %the closest distance the mobile can be from the the BS antenna
    hB=0.01; % BS antenna height
    hm=0.002; % mobile antenna height
    
    Ru_iteration=0; 
    a=2;
    b=2;
    K=1;
    St=1;
    NI=6;
    R_iteration = 0;
    
    %lamdac=(3*10^5)/frequency; %calculate lamdac
    lamdac=(3*10^5)/frequency; %calculate lamdac
    g=(4*hB*hm)/lamdac; %calculate the break point of the pass-loss curve
    
    
    for R = 0.1:0.05:1;
        Ru_iteration=Ru_iteration+1;
        
        
        %for Ru=4:2:8;  % this is the normalized Ru
        Ru_4 = 4;
        Ru_6 = 6;
        Ru_8 = 8;
        D_Ru_4=Ru_4.*R; %reuse distance
        D_Ru_6=Ru_6.*R; %reuse distance
        D_Ru_8=Ru_8.*R; %reuse distance
        
        
        for i=1:1:10000; % repeat the process for 10,000 times
            %%%%%%%%%% STEP 1%%%%%
            u=rand(1,1);
            r= R0+(R - R0)*sqrt(u); % Step 1, the Random Distance of the user
            
            %%%%%%%%%%STEP 2%%%%%%%%%%%
            ui = rand(6,1);  %Step 2, ,6 means NI = 6
            vi = rand(6,1);  % Step 2,
            xi = R0+(R - R0)*sqrt(ui);  %the distance of the of the interferer UI
            thetai=2*pi*vi; %%the distance of the of the interferer UI
            
            %%%%%%%%%%STEP 3%%%%%%%%%%%
            %ri=(D*D+(xi).^2+(2*xi*D).*sin(thetai)).^(1/2); %the actual distance of the interferer
            ri_Ru4=(D_Ru_4.^2+(xi).^2+(2*xi*D_Ru_4).*sin(thetai)).^(1/2);
            ri_Ru6=(D_Ru_6.^2+(xi).^2+(2*xi*D_Ru_6).*sin(thetai)).^(1/2);
            ri_Ru8=(D_Ru_8.^2+(xi).^2+(2*xi*D_Ru_8).*sin(thetai)).^(1/2);
            %%%%%%%%%%STEP 4 & Step 5%%%%%%%%%%%
            Sd=(K/((r^a)*(1+r/g)^b))*St;
            
            S_R_interfered_Ru4 =(K./((ri_Ru4.^a).*(1+ri_Ru4./g).^b))*St; 
            S_R_interfered_Ru6 =(K./((ri_Ru6.^a).*(1+ri_Ru6./g).^b))*St; 
            S_R_interfered_Ru8 =(K./((ri_Ru8.^a).*(1+ri_Ru8./g).^b))*St; 
            
            gammad_Ru4  = Sd./(mean(S_R_interfered_Ru4));
            gammad_Ru6  = Sd./(mean(S_R_interfered_Ru6));
            gammad_Ru8  = Sd./(mean(S_R_interfered_Ru8));
            %%%%%%%%%% STETP 6 %%%%%%%%%%%
            
            %Ae=(4./(pi*(Ru.^2)*(R*R))).*(log2(1+gammad)); %the ASE of genera simulations
            Ae_Ru4(Ru_iteration,i)=(4./(pi*(Ru_4.^2)*(R*R))).*(log2(1+gammad_Ru4)); %the ASE of genera simulations
            Ae_Ru6(Ru_iteration,i)=(4./(pi*(Ru_6.^2)*(R*R))).*(log2(1+gammad_Ru6)); %the ASE of genera simulations
            Ae_Ru8(Ru_iteration,i)=(4./(pi*(Ru_8.^2)*(R*R))).*(log2(1+gammad_Ru8)); %the ASE of genera simulations
            
            
            end
     end
end
