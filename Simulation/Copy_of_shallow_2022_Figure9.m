R=200; %the cell radius(m)
R0=20; %the closest distance the mobile can be from the the BS antenna
hB=10; % BS antenna height
hm=2; % mobile antenna height

%R=0.2; %the cell radius(KM)
%R0=0.02; %the closest distance the mobile can be from the the BS antenna
%hB=0.01; % BS antenna height
%hm=0.002; % mobile antenna height

fc=9*10^8; %carrier frequence
lamdac=(3*10^8)/fc; %calculate lamdac
g=(4*hB*hm)/lamdac; %calculate the break point of the pass-loss curve

c=0;
a=2;
b=2;
K=1;
St=1;
NI=6;

Xi=10/log(10); %the constant number of xi
sigmad=2.5119^4; %sigama is measured in dB
sigmaI=2.5119^4;
W0=1;

for Ru=2:0.2:10
    D=Ru.*R;
    c=c+1;
    for i=1:1:10000;
        u=rand(1,1);
        r=R0+(R-R0)*(u^(1/2)); %user's position to the BS
        ui=rand(6,1); %interferer's random distance to the BSi
        vi=rand(6,1); %interferer's random angle to the BSi
        xi=R0+(R-R0)*(ui.^(1/2)); %interferer's position to the BSi
        thetai=2*pi*vi; %interferer'angle to the BSi
        ri=(D*D+(xi).^2+(2*D*xi).*sin(thetai)).^(1/2); %the distance from the ithinterferer to theBS

        Sd=(K/((r^a)*(1+r/g)^b))*St; % equation 1!(St)(K)
        Si=(K./((ri.^a).*(1+ri./g).^b))*St ; %6 of received power level from ithinterfering mobile
        SI=sum(Si); %total interfering power
        Sib=(K/(((D+R)^a)*(1+(D+R)/g)^b))*St ; %6 of received power level fromith interfering mobile
        SIb=6*(Sib); %total interfering power
        Siw=(K/(((D-R)^a)*(1+(D-R)/g)^b))*St ; %6 of received power level fromith interfering mobile
        SiW=6*(Siw); %total interfering power


        %Calculate log parameters simulation
        udd=1/(r.^a.*(1+(r/g)).^b)*St;
        uii=1./(ri.^a.*(1+(ri/g)).^b)*St;
        mdd=Xi*log(udd);
        mii=Xi*log(uii);
        
        %log mean function for desired user simulation
        mudd = log((mdd^2)/sqrt(sigmad.^(1/2)+mdd^2));
        muii = log((mii.^2)./sqrt(sigmad.^(1/2)+mii.^2));
        %log variance function for interferer simulation
        sigmadd = sqrt(log(sigmad.^(1/2)/(mdd^2)+1));
        sigmai = sqrt(log(sigmad.^(1/2)./(mii.^2)+1));
        
        
        log_ud=lognrnd(mudd,sigmadd,[1,1]);
        log_ui=lognrnd(muii,sigmai,[6,1]);
        SII=sum(Si.*log_ui);
        gammadi=(Sd*log_ud)/SII;
        
        %Calculate log parameters best-case
        udd_b=1/(r.^a.*(1+(r/g)).^b)*St;
        uii_b=1./((D+R)^a*(1+((D+R)/g))^b)*St;
        mdd_b=Xi*log(udd_b);
        mii_b=Xi*log(uii_b);
        %log function for desired user simulation
        
        mudd_b = log((mdd_b^2)/sqrt(sigmad.^(1/2)+mdd_b^2));
        muii_b = log((mii_b^2)/sqrt(sigmad^(1/2)+mii_b^2));
        %log function for interferer simulation
        sigmadd_b = sqrt(log(sigmad.^(1/2)/(mdd_b^2)+1));
        sigmai_b = sqrt(log(sigmad^(1/2)./(mii_b^2)+1));
        
        log_udb=lognrnd(mudd_b,sigmadd_b,[1,1]);
        log_uib=lognrnd(muii_b,sigmai_b,[1,1]);
        SIIb=6*(Sib*log_uib);
        gammadib=(Sd*log_udb)/SIIb;
                                
        
        
        %Calculate log parameters best-case
        udd_w=1/(r.^a.*(1+(r/g)).^b)*St;
        uii_w=1./((D-R)^a*(1+((D-R)/g))^b)*St;
        mdd_w=Xi*log(udd_w);
        mii_w=Xi*log(uii_w);
        %log function for desired user simulation
        mudd_w = log((mdd_w^2)/sqrt(sigmad.^(1/2)+mdd_w^2));
        muii_w = log((mii_w^2)/sqrt(sigmad^(1/2)+mii_w^2));
        %log function for interferer simulation
        sigmadd_w = sqrt(log(sigmad^(1/2)/(mdd_w^2)+1));
        sigmai_w = sqrt(log(sigmad^(1/2)/(mii_w^2)+1));
        log_udw=lognrnd(mudd_w,sigmadd_w,[1,1]);
        log_uiw=lognrnd(muii_w,sigmai_w,[1,1]);
        SIIw=6*(Siw*log_uiw);
        
        gammadiw=(Sd*log_udw)/SIIw; 
        gammad=Sd/SI; %CIR power ratio
        gammab=((((R*(Ru+1))/r)^a)*((g+(Ru+1)*R)/(g+r))^b)/NI;%the best case ofgammad
        gammaw=((((R*(Ru-1))/r)^a)*((g+(Ru-1)*R)/(g+r))^b)/NI;%the worst case ofgammad
        
        
        Ae(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammad)); %the generalsimulations
        Aeb(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammab)); %the best caseinterference
        Aew(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammaw));%the worst caseinterference
        Ael(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammadi)); %the generalsimulations
        Aeblb(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammadib)); %the best caseinterference
        Aewlw(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammadiw));%the worst caseinterference
    end
end
aAe=mean(Ae,2)*10^6;%the average of Ae after running 10000 times
aAeb=mean(Aeb,2)*10^6; %the average of the best of Ae after running 10000times
aAew=mean(Aew,2)*10^6; %the average of the worst of Ae after running10000 times
aAeld=mean(Ael,2)*10^6;%the average of Ae after running 10000 times
aAebldb=mean(Aeblb,2)*10^6; %the average of the best of Ae after running10000 times
aAewldw=mean(Aewlw,2)*10^6; %the average of the worst of Ae afterrunning 10000 times
Ru=2:0.2:10;
plot(Ru,aAe,'k',Ru,aAeb,'g',Ru,aAew,'r',Ru,aAeld,'--b',Ru,aAebldb,'--c',Ru,aAewldw,'--m');%set the best case is green,the worst line is red
legend('Simulation No Shadowing ','Worst-Case No Shadowing','Best-CaseNoShadowing','Simulation Shadowing','Best-Case Shadowing','Worst-Case Shadowing')
xlabel('Normalized Reuse Distance Ru');grid;
ylabel('ASE[Bits/Sec/Hz/Km^2]');
title('Figure 9 Effect of Shadowing,');