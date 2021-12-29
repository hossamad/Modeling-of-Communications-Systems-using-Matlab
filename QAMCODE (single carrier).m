clear all; clc;
%%%%%%%%%%%%%%%%%   Single Carrier   %%%%%%%%%%%%%%%%%%%%%
NumOfBitsInBitStream=8*1e5;% A number that is divisible by 4 
NumOfSymbols=NumOfBitsInBitStream/4;
SNR=-10:1:25; %(IN dB)
SNRInLinear=10.^((SNR/10));
%% Uncoded
SER=zeros(1,length(SNR));
%generate Rayleigh for uncoded QPSK
v_1=randn(NumOfSymbols,1) ; 
v_2=randn(NumOfSymbols,1) ; 
R=sqrt(((v_1).^2+(v_2).^2)/2);
Eo=sum(R.^2)/length(R); Eb=2.5*Eo;
for h= 1:length(SNR)
    % generate randam data bits   
    data=(randn(1,NumOfBitsInBitStream)>0);
    % mapping
    QAMConstellation  = MapperFunction_16QAM( data,NumOfSymbols,Eo );
    % add noise
    n=sqrt(Eb/SNRInLinear(h)) .*( randn(size(QAMConstellation)) + ...
        1i*randn(size(QAMConstellation)) );
    receivedSig = (R.*QAMConstellation + n)./R;
    ReceivedSymbols(:,1)=real(receivedSig);
    ReceivedSymbols(:,2)=imag(receivedSig);
    % demapping
    a_demapped=zeros(NumOfSymbols,1);
    b_demapped=zeros(NumOfSymbols,1);
    Data_bits_reshaped = zeros(NumOfSymbols,4);
    temp = zeros(1 , 4);%temprorary vector to store each symbol(4-bits)
    for i = 1:NumOfSymbols
        %first row 
        if(ReceivedSymbols(i,1)>2*sqrt(Eo) && ReceivedSymbols(i,2)>2*sqrt(Eo))
            a_demapped(i,1)= 3;     b_demapped(i,1)=3 ;
        elseif(ReceivedSymbols(i,1)>2*sqrt(Eo) && ReceivedSymbols(i,2)<2*sqrt(Eo)...
                && ReceivedSymbols(i,2)>0)
            a_demapped(i,1)=3 ;     b_demapped(i,1)= 1;
        elseif(ReceivedSymbols(i,1)>2*sqrt(Eo) && ReceivedSymbols(i,2)<0 &&...
                ReceivedSymbols(i,2)>-2*sqrt(Eo) )
            a_demapped(i,1)=3 ;     b_demapped(i,1)= -1;
        elseif(ReceivedSymbols(i,1)>2*sqrt(Eo) && ReceivedSymbols(i,2)<-2*sqrt(Eo) )
            a_demapped(i,1)=3 ;     b_demapped(i,1)=-3 ;
        %second row
        elseif(ReceivedSymbols(i,1)<2*sqrt(Eo)  && ReceivedSymbols(i,1)>0 ...
                && ReceivedSymbols(i,2)>2*sqrt(Eo))
            a_demapped(i,1)= 1;     b_demapped(i,1)= 3;
        elseif(ReceivedSymbols(i,1)<2*sqrt(Eo) && ReceivedSymbols(i,1)>0 ...
                && ReceivedSymbols(i,2)<2*sqrt(Eo) && ReceivedSymbols(i,2)>0)
            a_demapped(i,1)= 1;     b_demapped(i,1)= 1;
        elseif(ReceivedSymbols(i,1)<2*sqrt(Eo) && ReceivedSymbols(i,1)>0...
            && ReceivedSymbols(i,2)<0 && ReceivedSymbols(i,2)>-2*sqrt(Eo) )
            a_demapped(i,1)= 1;     b_demapped(i,1)= -1;
        elseif(ReceivedSymbols(i,1)<2*sqrt(Eo) && ReceivedSymbols(i,1)>0 ...
                && ReceivedSymbols(i,2)<-2*sqrt(Eo) )
            a_demapped(i,1)= 1;     b_demapped(i,1)= -3; 
        %third row
        elseif(ReceivedSymbols(i,1)<0 && ReceivedSymbols(i,1)>-2*sqrt(Eo) ...
                && ReceivedSymbols(i,2)>2*sqrt(Eo))
            a_demapped(i,1)= -1;     b_demapped(i,1)= 3;
        elseif(ReceivedSymbols(i,1)<0 && ReceivedSymbols(i,1)>-2*sqrt(Eo) ...
                && ReceivedSymbols(i,2)<2*sqrt(Eo) && ReceivedSymbols(i,2)>0)
            a_demapped(i,1)= -1;     b_demapped(i,1)= 1;
        elseif(ReceivedSymbols(i,1)<0 && ReceivedSymbols(i,1)>-2*sqrt(Eo) ...
                && ReceivedSymbols(i,2)<0 && ReceivedSymbols(i,2)>-2*sqrt(Eo) )
            a_demapped(i,1)= -1;     b_demapped(i,1)= -1;
        elseif(ReceivedSymbols(i,1)<0 && ReceivedSymbols(i,1)>-2*sqrt(Eo) ...
                && ReceivedSymbols(i,2)<-2*sqrt(Eo) )
            a_demapped(i,1)= -1;     b_demapped(i,1)=-3;
        %fourh row
        elseif(ReceivedSymbols(i,1)<-2*sqrt(Eo) && ReceivedSymbols(i,2)>2*sqrt(Eo))
            a_demapped(i,1)=-3 ;     b_demapped(i,1)=3 ;
        elseif(ReceivedSymbols(i,1)<-2*sqrt(Eo) && ReceivedSymbols(i,2)<2*sqrt(Eo) ...
                && ReceivedSymbols(i,2)>0)
            a_demapped(i,1)= -3;     b_demapped(i,1)= 1;
        elseif(ReceivedSymbols(i,1)<-2*sqrt(Eo) && ReceivedSymbols(i,2)>-2*sqrt(Eo) ...
                && ReceivedSymbols(i,2)<0)
            a_demapped(i,1)= -3;     b_demapped(i,1)= -1;
        elseif(ReceivedSymbols(i,1)<-2*sqrt(Eo) && ReceivedSymbols(i,2)<-2*sqrt(Eo) )
            a_demapped(i,1)=-3;     b_demapped(i,1)=-3;
        end
    
    end
    Demapped_QAM = sqrt(Eo) * (a_demapped+1i*b_demapped);
    % compute SER
    SER(h) = sum((real(Demapped_QAM) ~= real(QAMConstellation)) + (imag(Demapped_QAM) ...
        ~= imag(QAMConstellation)));
end
BER=SER/(4*NumOfSymbols);


%% Coded
SERofCodedQAM=zeros(length(SNR),1);
%generate Rayleigh for coded QPSK
v_1_ForCodedQAM=randn(NumOfSymbols*3,1);   
v_2_ForCodedQAM=randn(NumOfSymbols*3,1); 
R_ForCodedQAM=sqrt(((v_1_ForCodedQAM).^2+(v_2_ForCodedQAM).^2)/2);
% use the same symbol energy, that used in uncoded, divided by 3 to make..
% .. the comparision fair in the power wise.
E_forCoded=Eo/3; Eb_forCoded=2.5*E_forCoded;
for h= 1:length(SNR)
    % generate randam data bits   
    dataR=(randn(1,NumOfBitsInBitStream)>0);
    % mapping
    QAMConstellation_beforeRepetition = MapperFunction_16QAM( dataR,NumOfSymbols,E_forCoded );
    % repeat symbols
    QAMConstellation_AfterRepetition=zeros(NumOfSymbols*3,1);
    repeat=0;
    k=1;
    for i= 1:NumOfSymbols*3    
        if repeat == 0
            QAMConstellation_AfterRepetition(i)=QAMConstellation_beforeRepetition(k);
            RepeatedSymbol=QAMConstellation_beforeRepetition(k);
            repeat=1+repeat;
            k=1+k;
        else
            QAMConstellation_AfterRepetition(i)=RepeatedSymbol;
            repeat=1+repeat;
        end
        if repeat== 3
            repeat=0;
        end
    end
    % add noise     
    n_R=sqrt(Eb_forCoded/SNRInLinear(h)) .*( randn(size(QAMConstellation_AfterRepetition)) ...
        + 1i*randn(size(QAMConstellation_AfterRepetition)) );
    receivedSigofCoded = (R_ForCodedQAM.*QAMConstellation_AfterRepetition + n_R)./R_ForCodedQAM;
    ReceivedSymbols_coded(:,1)=real(receivedSigofCoded);
    ReceivedSymbols_coded(:,2)=imag(receivedSigofCoded);
    % demapping
    a_demappedR=zeros(NumOfSymbols*3,1); b_demappedR=zeros(NumOfSymbols*3,1);
    for i = 1:NumOfSymbols*3
        %first row 
        if(ReceivedSymbols_coded(i,1)>2*sqrt(E_forCoded) && ...
                ReceivedSymbols_coded(i,2)>2*sqrt(E_forCoded))
            a_demappedR(i,1)= 3;     b_demappedR(i,1)=3 ;
        elseif(ReceivedSymbols_coded(i,1)>2*sqrt(E_forCoded) && ...
                ReceivedSymbols_coded(i,2)<2*sqrt(E_forCoded) && ReceivedSymbols_coded(i,2)>0)
            a_demappedR(i,1)=3 ;     b_demappedR(i,1)= 1;
        elseif(ReceivedSymbols_coded(i,1)>2*sqrt(E_forCoded) && ...
                ReceivedSymbols_coded(i,2)<0 && ReceivedSymbols_coded(i,2)>-2*sqrt(E_forCoded) )
            a_demappedR(i,1)=3 ;     b_demappedR(i,1)= -1;
        elseif(ReceivedSymbols_coded(i,1)>2*sqrt(E_forCoded) && ...
                ReceivedSymbols_coded(i,2)<-2*sqrt(E_forCoded) )
            a_demappedR(i,1)=3 ;     b_demappedR(i,1)=-3 ;
        %second row
        elseif(ReceivedSymbols_coded(i,1)<2*sqrt(E_forCoded)  && ...
                ReceivedSymbols_coded(i,1)>0 && ReceivedSymbols_coded(i,2)>2*sqrt(E_forCoded))
            a_demappedR(i,1)= 1;     b_demappedR(i,1)= 3;
        elseif(ReceivedSymbols_coded(i,1)<2*sqrt(E_forCoded) && ...
                ReceivedSymbols_coded(i,1)>0 && ReceivedSymbols_coded(i,2)<2*sqrt(E_forCoded) &&...
            ReceivedSymbols_coded(i,2)>0)
            a_demappedR(i,1)= 1;     b_demappedR(i,1)= 1;
        elseif(ReceivedSymbols_coded(i,1)<2*sqrt(E_forCoded) && ...
                ReceivedSymbols_coded(i,1)>0&& ReceivedSymbols_coded(i,2)<0 &&...
                ReceivedSymbols_coded(i,2)>-2*sqrt(E_forCoded) )
            a_demappedR(i,1)= 1;     b_demappedR(i,1)= -1;
        elseif(ReceivedSymbols_coded(i,1)<2*sqrt(E_forCoded) &&...
                ReceivedSymbols_coded(i,1)>0 && ReceivedSymbols_coded(i,2)<-2*sqrt(E_forCoded) )
            a_demappedR(i,1)= 1;     b_demappedR(i,1)= -3; 
        %third row
        elseif(ReceivedSymbols_coded(i,1)<0 && ReceivedSymbols_coded(i,1)>-2*sqrt(E_forCoded)  && ...
                ReceivedSymbols_coded(i,2)>2*sqrt(E_forCoded))
            a_demappedR(i,1)= -1;     b_demappedR(i,1)= 3;
        elseif(ReceivedSymbols_coded(i,1)<0 && ReceivedSymbols_coded(i,1)>-2*sqrt(E_forCoded) && ...
                ReceivedSymbols_coded(i,2)<2*sqrt(E_forCoded) && ReceivedSymbols_coded(i,2)>0)
            a_demappedR(i,1)= -1;     b_demappedR(i,1)= 1;
        elseif(ReceivedSymbols_coded(i,1)<0 && ReceivedSymbols_coded(i,1)>-2*sqrt(E_forCoded) &&...
                ReceivedSymbols_coded(i,2)<0 && ReceivedSymbols_coded(i,2)>-2*sqrt(E_forCoded) )
            a_demappedR(i,1)= -1;     b_demappedR(i,1)= -1;
        elseif(ReceivedSymbols_coded(i,1)<0 && ReceivedSymbols_coded(i,1)>-2*sqrt(E_forCoded) && ...
                ReceivedSymbols_coded(i,2)<-2*sqrt(E_forCoded) )
            a_demappedR(i,1)= -1;     b_demappedR(i,1)=-3;
        %fourh row
        elseif(ReceivedSymbols_coded(i,1)<-2*sqrt(E_forCoded) && ...
                ReceivedSymbols_coded(i,2)>2*sqrt(E_forCoded))
            a_demappedR(i,1)=-3 ;     b_demappedR(i,1)=3 ;
        elseif(ReceivedSymbols_coded(i,1)<-2*sqrt(E_forCoded) &&...
                ReceivedSymbols_coded(i,2)<2*sqrt(E_forCoded) && ReceivedSymbols_coded(i,2)>0)
            a_demappedR(i,1)= -3;     b_demappedR(i,1)= 1;
        elseif(ReceivedSymbols_coded(i,1)<-2*sqrt(E_forCoded) && ...
                ReceivedSymbols_coded(i,2)>-2*sqrt(E_forCoded) && ...
                ReceivedSymbols_coded(i,2)<0)
            a_demappedR(i,1)= -3;     b_demappedR(i,1)= -1;
        elseif(ReceivedSymbols_coded(i,1)<-2*sqrt(E_forCoded) && ...
                ReceivedSymbols_coded(i,2)<-2*sqrt(E_forCoded) )
            a_demappedR(i,1)=-3;     b_demappedR(i,1)=-3;
        end
    end
    Demapped_QAMR = sqrt(E_forCoded) .* (a_demappedR+1i*b_demappedR);
    %compute SER
    k=1;
    while k <= (NumOfSymbols*3-2)
        err=0;
        if Demapped_QAMR(k) ~= QAMConstellation_AfterRepetition(k) 
            err=err+1;
        end
        k=k+1;
        if Demapped_QAMR(k) ~= QAMConstellation_AfterRepetition(k)  
            err=err+1;
        end
        k=k+1;
        if Demapped_QAMR(k) ~= QAMConstellation_AfterRepetition(k)  
            err=err+1;
        end
        k=k+1;
        if err==2 || err==3
            SERofCodedQAM(h)=SERofCodedQAM(h)+1;
        end 
    end
end
BERofCodedQAM=SERofCodedQAM./NumOfBitsInBitStream;
figure
semilogy(SNR, BER, 'b',SNR,BERofCodedQAM,'r', 'LineWidth', 2)
xlabel('Eb/No (dB)');ylabel('Bit error rate');title('QAM (single carrier)');
legend('Uncoded','coded');grid on;