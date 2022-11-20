R=200; %the cell radius(m)
R0=20; %the closest distance the mobile can be from the the BS antenna
hB=10; % BS antenna height
hm=2; % mobile antenna height

fc=9*10^8; %carrier frequence
lamdac=(3*10^8)/fc; %calculate lamdac
g=(4*hB*hm)/lamdac; %calculate the break point of the pass-loss curve

c=0;
a=2;
b=2;
K=1;
St=1;
NI=6;
W0=1;
md = 1;
mi = 1.0;

Xi=10/log(10); %the constant number of xi
sigmad_shadow = 6; %sigama is measured in dB
sigmaI_shadow = 6;

%syms j x
%sum_d_sysm = symsum(1/(md+j)^2,j,0,Inf);
%sum_i_sysm = symsum(1/(mi+j)^2,j,0,Inf);
%sum_d = double(sum_d_sysm);
%sum_i  =double(sum_i_sysm);
%sigmad = Xi^2*sum_d+sigmad_shadow;
%sigmaI = Xi^2*sum_i+sigmaI_shadow;

for Ru=2.0:0.2:10;
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
        mdd= log(udd);
        mii= log(uii);
        log_ud=lognrnd(mdd,sigmad_shadow/Xi,[1,1]);
        log_ui=lognrnd(mii,sigmaI_shadow/Xi,[6,1]);
        gamma_b_d = log_ud/md;
        gamma_b_i = log_ui/mi;
        PDF_SD = gamrnd(md,gamma_b_d,1,1);
        PDF_SI = gamrnd(mi,gamma_b_i,6,1);
        gammad_multi = PDF_SD/sum(PDF_SI);
        

        %Calculate log parameters best-case
        udd_b=1/(r.^a.*(1+(r/g)).^b)*St;
        uii_b=1./((D+R)^a*(1+((D+R)/g))^b)*St;
        mdd_b=log(udd_b);
        mii_b=log(uii_b);
        log_udb=lognrnd(mdd_b,sigmad_shadow/Xi,[1,1]);
        log_uib=lognrnd(mii_b,sigmaI_shadow/Xi,[6,1]);
        gamma_b_best = log_udb/md;
        gamma_b_best_i = log_uib/mi;
        PDF_SD_best = gamrnd(md,gamma_b_best,1,1);
        PDF_SI_best = gamrnd(mi,gamma_b_best_i,6,1);
        gammadib = PDF_SD_best/sum(PDF_SI_best);
        
        
        %Calculate log parameters best-case
        udd_w=1/(r.^a.*(1+(r/g)).^b)*St;
        uii_w=1./((D-R)^a*(1+((D-R)/g))^b)*St;
        mdd_w=log(udd_w);
        mii_w=log(uii_w);
        log_udw=lognrnd(mdd_w,sigmad_shadow/Xi,[1,1]);
        log_uiw=lognrnd(mii_w,sigmaI_shadow/Xi,[6,1]);
        gamma_b_worst = log_udw/md;
        gamma_b_worst_i = log_uiw/mi;
        PDF_SD_worst = gamrnd(md,gamma_b_worst,1,1);
        PDF_SI_worst = gamrnd(mi,gamma_b_worst_i,6,1);
        gammadiw=PDF_SD_worst/sum(PDF_SI_worst); 
        

        gammad=Sd/SI; %CIR power ratio
        gammab=((((R*(Ru+1))/r)^a)*((g+(Ru+1)*R)/(g+r))^b)/NI;%the best case ofgammad
        gammaw=((((R*(Ru-1))/r)^a)*((g+(Ru-1)*R)/(g+r))^b)/NI;%the worst case ofgammad
        
        
        Ae(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammad)); %the generalsimulations
        Aeb(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammab)); %the best caseinterference
        Aew(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammaw));%the worst caseinterference
        
        Ael(c,i)=(4/(pi*(Ru.^2)*(R*R)))*(log2(1+gammad_multi)); %the generalsimulations
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

Ru=2.0:0.2:10;

plot(Ru,aAe,'k',Ru,aAeb,'g',Ru,aAew,'r',Ru,aAeld,'--b',Ru,aAebldb,'--k',Ru,aAewldw,'--m');%set the best case is green,the worst line is red
legend('Simulation ','Worst-Case','Best-Case','Simulation Combination','Best-Case Combination','Worst-Case Combination')
xlabel('Normalized Reuse Distance Ru');grid;
ylabel('ASE[Bits/Sec/Hz/Km^2]');
title('Figure Effect of Combination,');