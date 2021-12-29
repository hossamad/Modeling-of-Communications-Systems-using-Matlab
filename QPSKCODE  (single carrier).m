clc;   clear all; 
%%%%%%%%%%%%%%%%%   Single Carrier   %%%%%%%%%%%%%%%%%%%%%
NumOfBitsInBitStream=36e5;% A number that is divisible by 3 (for repetition code) and 2 (for QPSK symbols)
NumOfSymbols=NumOfBitsInBitStream/2;
SNR=-30:1.5:20; %(IN dB)
SNRInLinear=10.^((SNR/10));
%% Uncoded
BER=zeros(1,length(SNR));
%generate Rayleigh for uncoded QPSK
v_1=randn(NumOfSymbols,1) ; 
v_2=randn(NumOfSymbols,1) ; 
R=sqrt(((v_1).^2+(v_2).^2)/2);
% define symbol energy
Eo=sum(R.^2)/length(R); Eb=Eo/2;
for i= 1:length(SNR)
    QPSKConstellation=zeros(NumOfSymbols,1);
    % generate randam data bits   
    data=(randn(NumOfBitsInBitStream,1)>0);
    % mapping
    NameOfSymbol=0;
    count=1;
    for k=1:2:NumOfBitsInBitStream
        %data(K)>>b_1   data(K+1)>>b_0
       if data(k)==0 && data(k+1)==0 
            NameOfSymbol=2;
       elseif data(k)==0 && data(k+1)==1
            NameOfSymbol=3;
       elseif data(k)==1 && data(k+1)==0
            NameOfSymbol=1;
       elseif data(k)==1 && data(k+1)==1
            NameOfSymbol=4;
       end
    QPSKConstellation(count,1)=sqrt(Eo)* (cos((2*NameOfSymbol-1)*(pi/4)) - 1i*sin((2*NameOfSymbol-1)*(pi/4)));
    count=count+1;
    end
    % add noise (%as know all channel info)
    n=sqrt(Eb/SNRInLinear(i)).*( randn(size(QPSKConstellation)) + 1i*randn(size(QPSKConstellation)) );
    received_sig= (R.*QPSKConstellation+n)./R;
    % demapping
    FinalBits=zeros(NumOfBitsInBitStream,1);
    v=1;
    for k=1:NumOfSymbols
       if  real(received_sig(k,1)) > 0
           if imag(received_sig(k,1)) > 0
                FinalBits(v)=1; v=v+1;    FinalBits(v)=1; v=v+1;
           else
                FinalBits(v)=1; v=v+1;    FinalBits(v)=0; v=v+1;  
           end
       else
           if imag(received_sig(k,1)) > 0
                FinalBits(v)=0; v=v+1;    FinalBits(v)=1; v=v+1;
           else
                FinalBits(v)=0; v=v+1;    FinalBits(v)=0; v=v+1; 
           end
       end
    end
    BER(i)= sum((FinalBits ~= data));
end
BER=BER/(2*NumOfSymbols);
%% Coded
SERofCodedQPSK=zeros(length(SNR),1);
%generate Rayleigh for coded QPSK
v_1_ForCodedQPSK=randn(NumOfSymbols*3,1);   
v_2_ForCodedQPSK=randn(NumOfSymbols*3,1); 
R_ForCodedQPSK=sqrt(((v_1_ForCodedQPSK).^2+(v_2_ForCodedQPSK).^2)/2);
% use the same symbol energy, that used in uncoded, divided by 3 to make..
% .. the comparision fair in the power wise.
E_forCoded=Eo/3; Eb_forCoded=E_forCoded/2;
for i= 1:length(SNR)
    QPSKConstellation_beforeRepetition=zeros(NumOfSymbols,1);
    % generate randam data bits   
    data=(randn(NumOfBitsInBitStream,1)>0);
    % mapping
    NameOfSymbol=0;
    count=1;
    for k=1:2:NumOfBitsInBitStream
        %data(K)>>b_1   data(K+1)>>b_0
       if data(k)==0 && data(k+1)==0 
            NameOfSymbol=2;
       elseif data(k)==0 && data(k+1)==1
            NameOfSymbol=3;
       elseif data(k)==1 && data(k+1)==0
            NameOfSymbol=1;
       elseif data(k)==1 && data(k+1)==1
            NameOfSymbol=4;
       end
    QPSKConstellation_beforeRepetition(count,1)=sqrt(E_forCoded)* (cos((2*NameOfSymbol-1)*(pi/4)) - 1i*sin((2*NameOfSymbol-1)*(pi/4)));
    count=count+1;
    end
    % repeat symbols
    QPSKConstellation_AfterRepetition=zeros(NumOfSymbols*3,1);
    repeat=0;
    k=1;
    for h= 1:NumOfSymbols*3    
        if repeat == 0
            QPSKConstellation_AfterRepetition(h)=QPSKConstellation_beforeRepetition(k);
            RepeatedSymbol=QPSKConstellation_beforeRepetition(k);
            repeat=1+repeat;
            k=1+k;
        else
            QPSKConstellation_AfterRepetition(h)=RepeatedSymbol;
            repeat=1+repeat;
        end
        if repeat== 3
            repeat=0;
        end
    end
    % Add noise which is generated using the same symbol energy divided by 3
    n2=sqrt(Eb_forCoded/SNRInLinear(i)).*( randn(size(QPSKConstellation_AfterRepetition)) + 1i*randn(size(QPSKConstellation_AfterRepetition)) );
    received_sigofCoded= (R_ForCodedQPSK.*QPSKConstellation_AfterRepetition+n2)./R_ForCodedQPSK;
    % demapping
    DemappedQPSKConstellation=zeros(NumOfSymbols*3,1);
    count=1;
    for k=1:NumOfSymbols*3 
        NameOfSymbol=0;
        if  real(received_sigofCoded(k,1)) > 0
           if imag(received_sigofCoded(k,1)) > 0
                NameOfSymbol=4;
           else
                NameOfSymbol=1;  
           end
       else
           if imag(received_sigofCoded(k,1)) > 0
                NameOfSymbol=3;
           else
                NameOfSymbol=2;  
           end
        end
    DemappedQPSKConstellation(count,1)=sqrt(E_forCoded) * (cos((2*NameOfSymbol-1)*(pi/4)) - 1i*sin((2*NameOfSymbol-1)*(pi/4)));
    count=count+1;
    end
    %compute SER
    k=1;
    while k <= (NumOfSymbols*3-2)
        err=0;
        if DemappedQPSKConstellation(k) ~= QPSKConstellation_AfterRepetition(k) 
            err=err+1;
        end
        k=k+1;
        if DemappedQPSKConstellation(k) ~= QPSKConstellation_AfterRepetition(k)  
            err=err+1;
        end
        k=k+1;
        if DemappedQPSKConstellation(k) ~= QPSKConstellation_AfterRepetition(k)  
            err=err+1;
        end
        k=k+1;
        if err==2 || err==3
            SERofCodedQPSK(i)=SERofCodedQPSK(i)+1;
        end 
    end  
end
BERofCodedQPSK=SERofCodedQPSK./(2*NumOfSymbols);
figure
semilogy(SNR,BER,'b',SNR,BERofCodedQPSK,'r', 'LineWidth', 2)
xlabel('Eb/No (dB)');ylabel('Bit error rate');title('QPSK (single carrier)');
legend('Uncoded','Coded');grid on;  

