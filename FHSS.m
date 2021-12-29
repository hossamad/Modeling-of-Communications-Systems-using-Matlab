clear;clc;close all;
%% Single user part
SP=[16 32 64];
for kk3=1:length(SP)
%% Testing the FHSS system with spreading factor =SP(m),2 symbol per hopping and without jamming
L = SP(kk3);        %no of frequency bands which will be used for hopping
n = 1e4;            %no of data symbols in the simulation
Lh = 1;             %no of hops per symbol, therefore it also is the symbol duration.
m = 1;              %no of users
%generating information bits
s_data=round(rand(n,m));
% note that exp(jx)=cos(x)+jsin(x).
xbase1 = [exp(j*2*pi*rand(Lh*n,1))];
xbase0 = [exp(j*2*pi*rand(Lh*n,1))];
% modulate using two orthogonal frequencies(xbase1 and xbase0):
%"kron(x1,x2)" function multiply the x2 vector by each element in the x1 ..
% .. vector except zeros in the x1 vector.
%-----------------------------------
%"kron(s_data,ones(Lh,1)).*xbase1" convert the '1' data symbols into ones..
%.. have a "Lh" duration and then modulates them by multiplying them by the..
%.. 2nd orthogonal freq. for BFSK.
%--------------------------------------
%"kron((1-s_data), ones(Lh,1)).*xbase0" convert the '0' data symbols into ones..
%.. have a "Lh" duration and then modulates them by multiplying them by the..
%.. 1st orthogonal freq. for BFSK.
xmodsig = [kron(s_data,ones(Lh,1)).*xbase1 kron((1-s_data), ones(Lh,1)).*xbase0];
clear xbase0 xbase1;
%% genrating a random hopping sequence nLh long
Phop = round(rand(Lh*n, 1)*(L-1))+1;%PN hopping pattern (random numbers from from 1 to L)
% repeat each number in the generated random sequence, because we need each..
%.. band to contain 2 data bits (symbols) of BFSK.
   Phop_repeated =zeros(Lh*n,1);
    repeat=0;
    k=1;
    for i= 1:Lh*n  
        if repeat == 0
            Phop_repeated (i)=Phop(k);
            RepeatedNumber=Phop(k);
            repeat=1;
            k=1+k;
        else
            Phop_repeated(i)=RepeatedNumber;
            repeat=0;
        end
    end
%divde the '1' modulated signals into 'L' bands.
Xsiga = sparse(1:Lh*n, Phop_repeated, xmodsig(:,1));
%divde the '0' modulated signals into 'L' bands.
Xsigb = sparse(1:Lh*n, Phop_repeated, xmodsig(:,2));
% note that length(Xsiga)=length(Xsiga)=Lh*n.
% Generating noise sequences for both frequency channels
noise1=randn(Lh*n,1)+j*randn(Lh*n,1);
noise2=randn(Lh*n,1)+j*randn(Lh*n,1);
% divide the generated noise into L bands.
Nsiga=sparse(1:Lh*n,Phop_repeated,noise1);
Nsigb=sparse(1:Lh*n,Phop_repeated,noise2);
clear noisel noise2 ;
BER_WithoutJamming=[];
BER_az=[];
for i=1:10
    Eb2N(i)=i;                   % ( Eb/N in dB ) 
    Eb2N_num=10^(Eb2N(i)/10);    % Eb/N in numeral 
    Var_n=1/(2*Eb2N_num);        %1/SNR is the noise variance
    signois=sqrt(Var_n);          % standard deviation 
    ych1=Xsiga+signois*Nsiga;    % AWGN complex channels 
    ych2=Xsigb+signois*Nsigb;    % AWGN channels 
    for kk=0:n-1
        Yvec1=[]; Yvec2=[];
        for kk2=1:Lh 
            Yvec1=[Yvec1 ych1(kk*Lh+kk2, Phop_repeated(kk*Lh+kk2))];
            Yvec2=[Yvec2 ych2(kk*Lh+kk2, Phop_repeated(kk*Lh+kk2))];
        end
        Q=full(Yvec1);
        I=full(Yvec2);
        ydim1=Q*Q';
        ydim2=I*I';
        dec(kk+1) = (ydim1>ydim2);%detected values
        
    end

% Compute BER from simulation 
decTransposed=transpose(dec);
BER_WithoutJamming(i) = sum( (decTransposed ~= s_data) );
% Compare against analytical BER .
BER_az(i)=0.5*exp(-Eb2N_num/2);
end
BER_WithoutJamming=BER_WithoutJamming./(n);


%% Testing the FHSS system with spreading factor =SP(m),2 symbol per hopping and jamming
BER_WithJamming=[];
Xsiga_JammingCase=Xsiga;
Xsigb_JammingCase=Xsigb;
% Add 2 jammed channels (randomly picked) 
nch1=round(rand*(L-1))+1;
if nch1==1
    nch2=nch1+1;
else
    nch2=nch1-1;
end
%compute number of number of element in these 2 bands
t1=0;t2=0;
for hh=1:length(Phop)
    if Phop(hh)==nch1
        t1=t1+1;
    end
    if Phop(hh)==nch2
        t2=t2+1;
    end
end
jammingTone1=300*cos(2*pi*t1);
jammingTone2=300*cos(2*pi*t2);
Xsiga_JammingCase(:,nch1)=Xsiga_JammingCase(:,nch1)+jammingTone1;
Xsigb_JammingCase(:,nch1)=Xsigb_JammingCase(:,nch1)+jammingTone2;
Xsiga_JammingCase(:,nch2)=Xsiga_JammingCase(:,nch2)+jammingTone1;
Xsigb_JammingCase(:,nch2)=Xsigb_JammingCase(:,nch2)+jammingTone2;
% Generating the channel noise (AWGN) 
for i=1:10
    Eb2N(i)=i;                   % ( Eb/N in dB ) 
    Eb2N_num=10^(Eb2N(i)/10);    % Eb/N in numeral 
    Var_n=1/(2*Eb2N_num);        %1/SNR is the noise variance
    signois=sqrt(Var_n);          % standard deviation 
    ych1=Xsiga_JammingCase+signois*Nsigb;    % AWGN complex channels 
    ych2=Xsigb_JammingCase+signois*Nsigb;    % AWGN channels 
    for kk=0:n-1
        Yvec1=[]; Yvec2=[];
        for kk2=1:Lh  
            Yvec1=[Yvec1 ych1(kk*Lh+kk2, Phop_repeated(kk*Lh+kk2))];
            Yvec2=[Yvec2 ych2(kk*Lh+kk2, Phop_repeated(kk*Lh+kk2))];
        end
        Q_jammedCase=full(Yvec1);
        I_jammedCase=full(Yvec2);
        ydim1=Q_jammedCase*Q_jammedCase';
        ydim2=I_jammedCase*I_jammedCase';
        dec_jammedCase(kk+1) = (ydim1>ydim2);%detected values
        
    end

% Compute BER from simulation 
decTransposed_jammedCase=transpose(dec_jammedCase);
BER_WithJamming(i) = sum( (decTransposed_jammedCase ~= s_data) );
end
BER_WithJamming=BER_WithJamming./(n);
%store results 
    if kk3==1
        BER_SP16_SU_NJ=BER_WithoutJamming;
        BER_SP16_SU_J=BER_WithJamming;
    elseif kk3==2
        BER_SP32_SU_NJ=BER_WithoutJamming;
        BER_SP32_SU_J=BER_WithJamming; 
    else
        BER_SP64_SU_NJ=BER_WithoutJamming;
        BER_SP64_SU_J=BER_WithJamming;  
    end
end
figure
figber=semilogy(Eb2N,BER_SP16_SU_NJ,Eb2N,BER_SP32_SU_NJ,Eb2N,BER_SP64_SU_NJ,...
    Eb2N,BER_SP16_SU_J,Eb2N,BER_SP32_SU_J,Eb2N,BER_SP64_SU_J);
set(figber, 'Linewidth', 2);
legend('without jamming, Spread factor=16', 'without jamming, Spread factor=32',...
    'without jamming, Spread factor=64','with jamming, Spread factor=16', ...
    'with jamming, Spread factor=32','with jamming, Spread factor=64');
fx=xlabel('E_b/N (dB)');
fy=ylabel('Bit error rate');grid on; 
title('BER of the FHSS system with and without jamming (single user)')


%% 2 USERS
for kk3=1:length(SP)
%% Testing the FHSS system with spreading factor =SP(m),2 symbol per hopping and without jamming
L = SP(kk3);        %no of frequency bands which will be used for hopping
n = 1e4;            %no of data symbols in the simulation
Lh = 1;             %no of hops per symbol, therefore it also is the symbol duration.
m =2;              %no of users
%generating information bits
s_data=round(rand(n,m));
% note that exp(jx)=cos(x)+jsin(x).
xbase1 = [exp(j*2*pi*rand(Lh*n,1))];
xbase0 = [exp(j*2*pi*rand(Lh*n,1))];
% modulate using two orthogonal frequencies(xbase1 and xbase0):
xmodsig1 = [kron(s_data(:,1),ones(Lh,1)).*xbase1 kron((1-s_data(:,1)), ones(Lh,1)).*xbase0];
xmodsig2 = [kron(s_data(:,2),ones(Lh,1)).*xbase1 kron((1-s_data(:,2)), ones(Lh,1)).*xbase0];
%% genrating a random hopping sequence nLh long
Phop = round(rand(Lh*n/2, 1)*(L-1))+1;%PN hopping pattern (random numbers from from 1 to L)
Phop_repeated =zeros(Lh*n,1);
repeat=0;
k=1;
for i= 1:Lh*n  
	if repeat == 0
    	Phop_repeated (i)=Phop(k);
    	RepeatedNumber=Phop(k);
        repeat=1;
    	k=1+k;
    else
        Phop_repeated(i)=RepeatedNumber;
        repeat=0;
	end
end
%divde the '1' modulated signals into 'L' bands.
Xsiga = sparse(1:Lh*n, Phop_repeated, xmodsig1(:,1));
%divde the '0' modulated signals into 'L' bands.
Xsigb = sparse(1:Lh*n, Phop_repeated, xmodsig1(:,2));
%divde the '1' modulated signals into 'L' bands.
Xsiga2 = sparse(1:Lh*n, Phop_repeated, xmodsig2(:,1));
%divde the '0' modulated signals into 'L' bands.
Xsigb2 = sparse(1:Lh*n, Phop_repeated, xmodsig2(:,2));
% Generating noise sequences for both frequency channels
noise1=randn(Lh*n,1)+j*randn(Lh*n,1);
noise2=randn(Lh*n,1)+j*randn(Lh*n,1);
% divide the generated noise into L bands.
Nsiga=sparse(1:Lh*n,Phop_repeated,noise1);
Nsigb=sparse(1:Lh*n,Phop_repeated,noise2);
clear noisel noise2 ;
BER_WithoutJamming=[];
BER_az=[];
for i=1:10
    Eb2N(i)=i;                   % ( Eb/N in dB ) 
    Eb2N_num=10^(Eb2N(i)/10);    % Eb/N in numeral 
    Var_n=1/(2*Eb2N_num);        %1/SNR is the noise variance
    signois=sqrt(Var_n);         % standard deviation 
    ych1=Xsiga+signois*Nsiga;    % AWGN complex channels 
    ych2=Xsigb+signois*Nsigb;    % AWGN channels 
    for kk=0:n-1
        Yvec1=[]; Yvec2=[];
    for kk2=1:Lh 
    	Yvec1=[Yvec1 ych1(kk*Lh+kk2, Phop_repeated(kk*Lh+kk2))];
    	Yvec2=[Yvec2 ych2(kk*Lh+kk2, Phop_repeated(kk*Lh+kk2))];
    end
    Q=full(Yvec1);
    I=full(Yvec2);
    ydim1=Q*Q';
    ydim2=I*I';
    dec(kk+1) = (ydim1>ydim2);%detected values
    end 
    %user 2 same channel
    ych12=Xsiga2+signois*Nsiga;    % AWGN complex channels 
    ych22=Xsigb2+signois*Nsigb;    % AWGN channels 
    for kk=0:n-1
        Yvec12=[]; Yvec22=[];
    for kk2=1:Lh 
    	Yvec12=[Yvec12 ych12(kk*Lh+kk2, Phop_repeated(kk*Lh+kk2))];
    	Yvec22=[Yvec22 ych22(kk*Lh+kk2, Phop_repeated(kk*Lh+kk2))];
    end
        Q2=full(Yvec12);
        I2=full(Yvec22);
        ydim12=Q2*Q2';
        ydim22=I2*I2';
        dec2(kk+1) = (ydim12>ydim22);%detected values 
    end

% Compute BER from simulation 
decTransposed=transpose(dec);
BER_WithoutJamming1(i) = sum( (decTransposed ~= s_data(:,1)) );
% Compute BER from simulation 
decTransposed2=transpose(dec2);
BER_WithoutJamming2(i) = sum( (decTransposed2 ~= s_data(:,2)) );
end
BER_WithoutJamming1=BER_WithoutJamming1./(n);
BER_WithoutJamming2=BER_WithoutJamming2./(n);

%% Testing the FHSS system with spreading factor =SP(m),2 symbol per hopping and jamming
BER_WithJamming1=[];
BER_WithJamming2=[];
% Add a jammed channel (randomly picked) 
nch1=round(rand*(L-1))+1;
Xsiga_JammingCase1=Xsiga;
Xsigb_JammingCase1=Xsigb;
Xsiga_JammingCase2=Xsiga2;
Xsigb_JammingCase2=Xsigb2;
% Generating noise sequences for both frequency channels 
% Add 2 jammed channels (randomly picked) 
nch1=round(rand*(L-1))+1;
if nch1==1
    nch2=nch1+1;
else
    nch2=nch1-1;
end
%compute number of number of element in these 2 bands
t1=0;t2=0;
for hh=1:length(Phop)
    if Phop(hh)==nch1
        t1=t1+1;
    end
    if Phop(hh)==nch2
        t2=t2+1;
    end
end
jammingTone1=300*cos(2*pi*t1);
jammingTone2=300*cos(2*pi*t2);
%u1
Xsiga_JammingCase1(:,nch1)=Xsiga_JammingCase1(:,nch1)+jammingTone1;
Xsigb_JammingCase1(:,nch1)=Xsigb_JammingCase1(:,nch1)+jammingTone2;
Xsiga_JammingCase1(:,nch2)=Xsiga_JammingCase1(:,nch2)+jammingTone1;
Xsigb_JammingCase1(:,nch2)=Xsigb_JammingCase1(:,nch2)+jammingTone2;
%u2
Xsiga_JammingCase2(:,nch1)=Xsiga_JammingCase2(:,nch1)+jammingTone1;
Xsigb_JammingCase2(:,nch1)=Xsigb_JammingCase2(:,nch1)+jammingTone2;
Xsiga_JammingCase2(:,nch2)=Xsiga_JammingCase2(:,nch2)+jammingTone1;
Xsigb_JammingCase2(:,nch2)=Xsigb_JammingCase2(:,nch2)+jammingTone2;
% Generating the channel noise (AWGN) 
for i=1:10
    Eb2N(i)=i;                   % ( Eb/N in dB ) 
    Eb2N_num=10^(Eb2N(i)/10);    % Eb/N in numeral 
    Var_n=1/(2*Eb2N_num);        %1/SNR is the noise variance
    signois=sqrt(Var_n);          % standard deviation 
    %user1
    ych1=Xsiga_JammingCase1+signois*Nsigb;    % AWGN complex channels 
    ych2=Xsigb_JammingCase1+signois*Nsigb;    % AWGN channels 
    %USER1
    for kk=0:n-1
    Yvec1=[]; Yvec2=[];
    for kk2=1:Lh  
    	Yvec1=[Yvec1 ych1(kk*Lh+kk2, Phop_repeated(kk*Lh+kk2))];
    	Yvec2=[Yvec2 ych2(kk*Lh+kk2, Phop_repeated(kk*Lh+kk2))];
    end
    Q_jammedCase=full(Yvec1);
    I_jammedCase=full(Yvec2);
    ydim1=Q_jammedCase*Q_jammedCase';
    ydim2=I_jammedCase*I_jammedCase';
    dec_jammedCase1(kk+1) = (ydim1>ydim2);%detected values  
    end
    
    %user 2
    ych12=Xsiga_JammingCase2+signois*Nsigb;    % AWGN complex channels 
    ych22=Xsigb_JammingCase2+signois*Nsigb;    % AWGN channels 
    for kk=0:n-1
    Yvec12=[]; Yvec22=[];
    for kk2=1:Lh  
        	Yvec12=[Yvec12 ych12(kk*Lh+kk2, Phop_repeated(kk*Lh+kk2))];
            Yvec22=[Yvec22 ych22(kk*Lh+kk2, Phop_repeated(kk*Lh+kk2))];
    end
    Q_jammedCase2=full(Yvec12);
    I_jammedCase2=full(Yvec22);
    ydim12=Q_jammedCase2*Q_jammedCase2';
    ydim22=I_jammedCase2*I_jammedCase2';
    dec_jammedCase2(kk+1) = (ydim12>ydim22);%detected values
    end
% Compute BER from simulation 
decTransposed_jammedCase1=transpose(dec_jammedCase1);
BER_WithJamming1(i) = sum( (decTransposed_jammedCase1 ~= s_data(:,1)) );
% Compute BER from simulation 
decTransposed_jammedCase2=transpose(dec_jammedCase2);
BER_WithJamming2(i) = sum( (decTransposed_jammedCase2 ~= s_data(:,2)) );
end
BER_WithJamming1=BER_WithJamming1./(n);
BER_WithJamming2=BER_WithJamming2./(n);
%store results 
    if kk3==1
        BER_SP16_2U_NJ1=BER_WithoutJamming1;
        BER_SP16_2U_J1=BER_WithJamming1;
        BER_SP16_2U_NJ2=BER_WithoutJamming2;
        BER_SP16_2U_J2=BER_WithJamming2;
    elseif kk3==2
        BER_SP32_2U_NJ1=BER_WithoutJamming1;
        BER_SP32_2U_J1=BER_WithJamming1; 
        BER_SP32_2U_NJ2=BER_WithoutJamming2;
        BER_SP32_2U_J2=BER_WithJamming2; 
    else
        BER_SP64_2U_NJ1=BER_WithoutJamming1;
        BER_SP64_2U_J1=BER_WithJamming1; 
        BER_SP64_2U_NJ2=BER_WithoutJamming2;
        BER_SP64_2U_J2=BER_WithJamming2; 
    end
end
%PLOT USER 1
figure
figber=semilogy(Eb2N,BER_SP16_2U_NJ1,Eb2N,BER_SP32_2U_NJ1,Eb2N,BER_SP64_2U_NJ1,...
    Eb2N,BER_SP16_2U_J1,Eb2N,BER_SP32_2U_J1,Eb2N,BER_SP64_2U_J1);
set(figber, 'Linewidth', 2);
legend('NoJam, SF=16', 'NoJam, SF=32',...
    'U1,NoJam, SF=64','U1,withJam, SF=16', ...
    'U1,withJam, SF=32','U1,withJam, SF=64');
fx=xlabel('E_b/N (dB)');
fy=ylabel('Bit error rate');grid on; 
title('BER of the FHSS system with and without jamming (the first user)')
%PLOT USER 2
figure
figber=semilogy(Eb2N,BER_SP16_2U_NJ2,Eb2N,BER_SP32_2U_NJ2,Eb2N,BER_SP64_2U_NJ2,...
    Eb2N,BER_SP16_2U_J2,Eb2N,BER_SP32_2U_J2,Eb2N,BER_SP64_2U_J2);
set(figber, 'Linewidth', 2);
legend('NoJam, SF=16', 'NoJam, SF=32',...
    'NoJam, SF=64','withJam, SF=16', ...
    'withJam, SF=32','withJam, SF=64');
fx=xlabel('E_b/N (dB)');
fy=ylabel('Bit error rate');grid on; 
title('BER of the FHSS system with and without jamming (the second user)')
figure
figber=semilogy(Eb2N,BER_SP16_2U_NJ1,Eb2N,BER_SP32_2U_NJ1,Eb2N,BER_SP64_2U_NJ1,...
    Eb2N,BER_SP16_2U_NJ2,Eb2N,BER_SP32_2U_NJ2,Eb2N,BER_SP64_2U_NJ2);
set(figber, 'Linewidth', 2);
legend('user1, NoJam, SF=16', 'user1, NoJam, SF=32',...
    'user1, NoJam, SF=64','user2, NoJam, SF=16', 'user2, NoJam, SF=32',...
    'user2, NoJam, SF=64');
fx=xlabel('E_b/N (dB)');
fy=ylabel('Bit error rate');grid on; 
title('BER of the FHSS system user1 vs user2 (without jamming)')
figure
figber=semilogy(Eb2N,BER_SP16_2U_J1,Eb2N,BER_SP32_2U_J1,Eb2N,BER_SP64_2U_J1,...
    Eb2N,BER_SP16_2U_J2,Eb2N,BER_SP32_2U_J2,Eb2N,BER_SP64_2U_J2);
set(figber, 'Linewidth', 2);
legend('user1, withJam, SF=16', ...
    'user1, withJam, SF=32','user1, withJam, SF=64', ...
    'user2, withJam, SF=16', ...
    'user2, withJam, SF=32','user2, withJam, SF=64');
fx=xlabel('E_b/N (dB)');
fy=ylabel('Bit error rate');grid on; 
title('BER of the FHSS system user1 vs user2 (with jamming)')
