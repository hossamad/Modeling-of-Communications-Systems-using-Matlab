 clear all; clc;
%%%%%%%%%%%%%%%%%   Multicarriers QAM  %%%%%%%%%%%%%%%%%%%%%
NumOfBitsInBitStream=400*256;% A number that is divisible by 256 (for ofdm symbols)
NumOfOFDMSymbols=NumOfBitsInBitStream/(16*16);
NumOfSymbols=NumOfBitsInBitStream/4;
SNR=-5:1.5:15; %(IN dB)
SNRInLinear=10.^((SNR/10));
N=64;


%% Uncoded (flat channel)
BER_QAMUncodedflatchann=zeros(1,length(SNR));
Eo=2;  Eb=2.5*Eo;
flatFadingChannel=(randn()+1i*randn());
for i= 1:length(SNR)
    % generate randam data bits   
    data=(randn(NumOfBitsInBitStream,1)>0);   
    %interleaver
    InterleavedData=zeros(1,NumOfBitsInBitStream);
    colsNum=16; rowsNum=16;  NumOfElements=16*16;
    for m=0:NumOfOFDMSymbols-1
        % convert the 16*16 elements from the data vector into matrix 16*16
        ReshapedMatrix=zeros(rowsNum,colsNum);
        for h=0:rowsNum-1
            ReshapedMatrix(h+1,:)=data( h*colsNum+1+NumOfElements*m : ...
                h*colsNum+colsNum+NumOfElements*m ,1);
        end
        % read column by column
        for h=0:colsNum-1
            InterleavedData(1, h*rowsNum+1+NumOfElements*m : h*rowsNum+rowsNum+NumOfElements*m )...
                =ReshapedMatrix(:,h+1);
        end                     
    end
     % mapping
     QAMConstellation = MapperFunction_16QAM( InterleavedData,NumOfSymbols,Eo );
    %IFFT
    %take 64 symbols and get 64 IFFT of them
    IFFToutput=zeros(NumOfSymbols,1);
    for k=0:(NumOfSymbols/64)-1
        IFFToutput((k*64)+1:(k*64)+64,1)=ifft(QAMConstellation((k*64)+1 : (k*64)+64 ,1),N);
    end
    %add cyclic extension
    CyclicExtensionLen=N/4;
    NewLenAfterCyclicExtension=NumOfSymbols+ (CyclicExtensionLen)*NumOfOFDMSymbols;
    TXOutput=zeros(1,NewLenAfterCyclicExtension);
    for k=0:( NumOfOFDMSymbols - 1 )
       buffer=zeros(1, N+CyclicExtensionLen );
       buffer( 1,1:CyclicExtensionLen  )=IFFToutput( (k*N)+(N-CyclicExtensionLen)+1 : k*N+N ,1);
       buffer( 1,CyclicExtensionLen+1:N+CyclicExtensionLen  )=IFFToutput( (k*N)+1 : k*N+N ,1);
       TXOutput(1,k*(N+CyclicExtensionLen)+1 : ...
           k*(N+CyclicExtensionLen)+N+CyclicExtensionLen )=buffer;  
    end
     %% flat fading channel
    n=sqrt(Eb/SNRInLinear(i)).*( randn(size(TXOutput)) + 1i*randn(size(TXOutput)) );
    rxedSymbolsFromfFlatChannel=(flatFadingChannel.*TXOutput+n)./flatFadingChannel;
    % get the FFT of OFDM symbols  
    SymbolsAfterFFT=zeros(1,NumOfSymbols);
    for k=0:( NewLenAfterCyclicExtension/(N+CyclicExtensionLen) - 1 )
       NewBuffer=zeros(1,N);
       %update indeces
       SymBegin=k*(N+CyclicExtensionLen)+CyclicExtensionLen+1;
       symEnd=k*(N+CyclicExtensionLen)+(N+CyclicExtensionLen);
       %get FFT FOR one OFDM symbol
       NewBuffer(1,:)=rxedSymbolsFromfFlatChannel(1, SymBegin :symEnd );
       SymbolsAfterFFT(1, k*N+1:k*N+N )=fft(NewBuffer,N);  
    end
    % demapping
    DemapperOutput  = DemapperFunction_16QAM( SymbolsAfterFFT,NumOfSymbols,NumOfBitsInBitStream,Eo );
    DemapperOutput2=transpose(DemapperOutput);
    % deinterleaver
    for m=0:NumOfOFDMSymbols-1
    % convert the 256 elements from the FinalBits vector into matrix 8*16
        colsNum=16;  rowsNum=16;  NumofElems=16*16;
        ReshapedMatrix=zeros(rowsNum,colsNum);
        for h=0:colsNum-1
            ReshapedMatrix(:,h+1)=DemapperOutput2( h*rowsNum+1+NumofElems*m :...
                h*rowsNum+rowsNum+NumofElems*m ,1);
        end
        % read row by row
        for h=0:rowsNum-1
            DeinterleavedData(  h*colsNum+1+NumofElems*m : h*colsNum+colsNum+NumofElems*m ,1)...
                =ReshapedMatrix(h+1,:);
        end                     
    end
    BER_QAMUncodedflatchann(i)=sum( (DeinterleavedData~=data) );
end
BER_QAMUncodedflatchann=BER_QAMUncodedflatchann./NumOfBitsInBitStream;



%% Coded (flat channel)
%update NumOfOFDMSymbols
NumOfOFDMSymbols=(NumOfBitsInBitStream*3)/256;
NumOfSymbols=(NumOfBitsInBitStream*3)/4;
E_forCoded=Eo/3; Eb_forCoded=2.5*E_forCoded;
BERofMC_QAM_FC=zeros(1,length(SNR));
for i= 1:length(SNR)
    % generate randam data bits   
    data_FC_QAM=(randn(NumOfBitsInBitStream,1)>0);   
    % repeat bits
    RepeatedData_FC_QAM=zeros(NumOfBitsInBitStream*3,1);
    repeat=0;
    k=1;
    for h= 1:NumOfBitsInBitStream*3    
        if repeat == 0
            RepeatedData_FC_QAM(h)=data_FC_QAM(k);
            RepeatedBit=data_FC_QAM(k);
            repeat=1+repeat;
            k=1+k;
        else
            RepeatedData_FC_QAM(h)=RepeatedBit;
            repeat=1+repeat;
        end
        if repeat== 3
            repeat=0;
        end
    end   
    %interleaver 
    InterleavedData_FC_QAM=zeros(1,NumOfBitsInBitStream*3);
    colsNum=16; rowsNum=16;  NumOfElements=16*16;
    for m=0:NumOfOFDMSymbols-1
        % convert the 16*16 elements from the data vector into matrix 16*16
        ReshapedMatrix_FC_QAM=zeros(rowsNum,colsNum);
        for h=0:rowsNum-1
            ReshapedMatrix_FC_QAM(h+1,:)=RepeatedData_FC_QAM( h*colsNum+1+NumOfElements*m :...
                h*colsNum+colsNum+NumOfElements*m ,1);
        end
        % read column by column
        for h=0:colsNum-1
            InterleavedData_FC_QAM(1, h*rowsNum+1+NumOfElements*m :...
                h*rowsNum+rowsNum+NumOfElements*m )...
                =ReshapedMatrix_FC_QAM(:,h+1);
        end                     
    end
    % mapping
    QAMConstellation_FC_QAM  = MapperFunction_16QAM( InterleavedData_FC_QAM,NumOfSymbols,E_forCoded );
    %IFFT
    %take 64 symbols and get 64 IFFT of them
    IFFToutput_FC_QAM=zeros(NumOfSymbols,1);
    for k=0:(NumOfSymbols/64)-1
        IFFToutput_FC_QAM((k*64)+1:(k*64)+64,1)=...
            ifft(QAMConstellation_FC_QAM((k*64)+1 : (k*64)+64 ,1),N);
    end
    %add cyclic extension
    CyclicExtensionLen=N/4;
    NewLenAfterCyclicExtension = NumOfSymbols + (CyclicExtensionLen)*NumOfOFDMSymbols;
    TXOutput_FC_QAM=zeros(1,NewLenAfterCyclicExtension);
    for k=0:( NumOfOFDMSymbols - 1 )
       buffer=zeros(1, N+CyclicExtensionLen );
       buffer( 1,1:CyclicExtensionLen  )=IFFToutput_FC_QAM(...
           (k*N)+(N-CyclicExtensionLen)+1 : k*N+N ,1);
       buffer( 1,CyclicExtensionLen+1:N+CyclicExtensionLen  )=...
           IFFToutput_FC_QAM( (k*N)+1 : k*N+N ,1);
       TXOutput_FC_QAM(1, k*(N+CyclicExtensionLen)+1 : k*(N+CyclicExtensionLen)+...
           N+CyclicExtensionLen )=buffer;  
    end
    %% flat fading channel
    n_FC_QAM=sqrt(Eb_forCoded/SNRInLinear(i)).*( randn(size(TXOutput_FC_QAM)) +...
        1i*randn(size(TXOutput_FC_QAM)) );
    OutputFlatChannel_FC_QAM=(flatFadingChannel.*TXOutput_FC_QAM+n_FC_QAM)./flatFadingChannel; 
    % get the FFT of OFDM symbols  
    SymbolsAfterFFT_FC_QAM=zeros(NumOfSymbols,1);
    for k=0:( NewLenAfterCyclicExtension/(N+CyclicExtensionLen) - 1 )
       %update indeces
       SymBegin=k*(N+CyclicExtensionLen)+CyclicExtensionLen+1;
       symEnd=k*(N+CyclicExtensionLen)+(N+CyclicExtensionLen);
       %get FFT FOR one OFDM symbol
       SymbolsAfterFFT_FC_QAM( k*N+1:k*N+N ,1)=fft( ...
           OutputFlatChannel_FC_QAM(1, SymBegin :symEnd ) ,N);  
    end
       % demapping
 DemapperOutput_FC_QAM  = DemapperFunction_16QAM( SymbolsAfterFFT_FC_QAM,NumOfSymbols,...
     NumOfBitsInBitStream*3,E_forCoded );
    DemapperOutput2_FC_QAM=transpose(DemapperOutput_FC_QAM);
    % deinterleaver
    for m=0:NumOfOFDMSymbols-1
    % convert the 256 elements from the FinalBits vector into matrix 8*16
        colsNum=16;  rowsNum=16;  NumofElems=16*16;
        ReshapedMatrix_FC_QAM=zeros(rowsNum,colsNum);
        for h=0:colsNum-1
            ReshapedMatrix_FC_QAM(:,h+1)=DemapperOutput2_FC_QAM( h*rowsNum+1+NumofElems*m : ...
                h*rowsNum+rowsNum+NumofElems*m ,1);
        end
        % read row by row
        for h=0:rowsNum-1
            DeinterleavedData_FC_QAM(  h*colsNum+1+NumofElems*m : h*colsNum+colsNum+NumofElems*m ,1)=...
                ReshapedMatrix_FC_QAM(h+1,:);
        end                     
    end
        k=1;
    while k <= (NumOfBitsInBitStream*3-2)
        err=0;
        if DeinterleavedData_FC_QAM(k) ~= RepeatedData_FC_QAM(k) 
            err=err+1;
        end
        k=k+1;
        if DeinterleavedData_FC_QAM(k) ~= RepeatedData_FC_QAM(k)  
            err=err+1;
        end
        k=k+1;
        if DeinterleavedData_FC_QAM(k) ~= RepeatedData_FC_QAM(k)  
            err=err+1;
        end
        k=k+1;
        if err==2 || err==3
           BERofMC_QAM_FC(i)=BERofMC_QAM_FC(i)+1;
        end 
    end
end
BERofMC_QAM_FC=BERofMC_QAM_FC./NumOfBitsInBitStream;


%% Uncoded (Frequency selective Channel)
NumOfOFDMSymbols=NumOfBitsInBitStream/256;
NumOfSymbols=NumOfBitsInBitStream/4;
BERofUncodedQAM_FreqSelec=zeros(length(SNR),1);
freqSelectiveChannel=[0.8 0 0 0 0 0 0 0 0 0 0.6];
for i= 1:length(SNR)
    % generate randam data bits   
    data_FS=(randn(NumOfBitsInBitStream,1)>0);   
    %interleaver
    InterleavedData_FS=zeros(1,NumOfBitsInBitStream);
    colsNum=16; rowsNum=16;  NumOfElements=16*16;
    for m=0:NumOfOFDMSymbols-1
        % convert the 16*16 elements from the data vector into matrix 16*16
        ReshapedMatrix_FS=zeros(rowsNum,colsNum);
        for h=0:rowsNum-1
            ReshapedMatrix_FS(h+1,:)=data_FS( h*colsNum+1+NumOfElements*m : ...
                h*colsNum+colsNum+NumOfElements*m ,1);
        end
        % read column by column
        for h=0:colsNum-1
            InterleavedData_FS(1, h*rowsNum+1+NumOfElements*m : h*rowsNum+rowsNum+NumOfElements*m )=...
                ReshapedMatrix_FS(:,h+1);
        end                     
    end
     % mapping
     QAMConstellation_FS = MapperFunction_16QAM( InterleavedData_FS,NumOfSymbols,Eo );
    %IFFT
    %take 64 symbols and get 64 IFFT of them
    IFFToutput_FS=zeros(NumOfSymbols,1);
    for k=0:(NumOfSymbols/64)-1
        IFFToutput_FS((k*64)+1:(k*64)+64,1)=ifft(QAMConstellation_FS((k*64)+1 : (k*64)+64 ,1),N);
    end
    %add cyclic extension
    CyclicExtensionLen=N/4;
    NewLenAfterCyclicExtension=NumOfSymbols+ (CyclicExtensionLen)*NumOfOFDMSymbols;
    TXOutput_FS=zeros(1,NewLenAfterCyclicExtension);
    for k=0:( NumOfOFDMSymbols - 1 )
       buffer=zeros(1, N+CyclicExtensionLen );
       buffer( 1,1:CyclicExtensionLen  )=IFFToutput_FS( (k*N)+(N-CyclicExtensionLen)+1 : k*N+N ,1);
       buffer( 1,CyclicExtensionLen+1:N+CyclicExtensionLen  )=IFFToutput_FS( (k*N)+1 : k*N+N ,1);
       TXOutput_FS(1,k*(N+CyclicExtensionLen)+1 : k*(N+CyclicExtensionLen)+...
           N+CyclicExtensionLen )=buffer;  
    end
    % Frequency selective channel
    n=sqrt(Eb/SNRInLinear(i)).*( randn(size(TXOutput_FS)) + 1i*randn(size(TXOutput_FS)) );
    ChannelinputFS=TXOutput_FS+n;
    for k=0:NumOfOFDMSymbols-1
        buffer2FS=ChannelinputFS(1,k*(N+CyclicExtensionLen)+1 : ...
            k*(N+CyclicExtensionLen)+N+CyclicExtensionLen );
        buffer3FS=conv(freqSelectiveChannel,buffer2FS);
        buffer4FS=deconv(buffer3FS,freqSelectiveChannel);
        FFTInputFS(1,k*(N+CyclicExtensionLen)+1 : k*(N+CyclicExtensionLen)+...
            N+CyclicExtensionLen)=buffer4FS;
    end
    % get the FFT of OFDM symbols  
    SymbolsAfterFFT_FS=zeros(1,NumOfSymbols);
    for k=0:( NewLenAfterCyclicExtension/(N+CyclicExtensionLen) - 1 )
       NewBuffer=zeros(1,N);
       %update indeces
       SymBegin=k*(N+CyclicExtensionLen)+CyclicExtensionLen+1;
       symEnd=k*(N+CyclicExtensionLen)+(N+CyclicExtensionLen);
       %get FFT FOR one OFDM symbol
       NewBuffer(1,:)=FFTInputFS(1, SymBegin :symEnd );
       SymbolsAfterFFT_FS(1, k*N+1:k*N+N )=fft(NewBuffer,N);  
    end
    % demapping
    DemapperOutput_FS  = DemapperFunction_16QAM( SymbolsAfterFFT_FS,NumOfSymbols,...
        NumOfBitsInBitStream,Eo );
    DemapperOutput2_FS=transpose(DemapperOutput_FS);
    % deinterleaver
    for m=0:NumOfOFDMSymbols-1
    % convert the 256 elements from the FinalBits vector into matrix 8*16
        colsNum=16;  rowsNum=16;  NumofElems=16*16;
        ReshapedMatrix_FS=zeros(rowsNum,colsNum);
        for h=0:colsNum-1
            ReshapedMatrix_FS(:,h+1)=DemapperOutput2_FS( h*rowsNum+1+NumofElems*m : ...
                h*rowsNum+rowsNum+NumofElems*m ,1);
        end
        % read row by row
        for h=0:rowsNum-1
            DeinterleavedData_FS(  h*colsNum+1+NumofElems*m : h*colsNum+colsNum+NumofElems*m ,1)=...
                ReshapedMatrix_FS(h+1,:);
        end                     
    end
    BERofUncodedQAM_FreqSelec(i)=sum( (DeinterleavedData_FS~=data_FS) );
end
BERofUncodedQAM_FreqSelec=BERofUncodedQAM_FreqSelec./NumOfBitsInBitStream;




%% Coded (Frequency selective Channel)
NumOfOFDMSymbols=(NumOfBitsInBitStream*3)/256;
NumOfSymbols=(NumOfBitsInBitStream*3)/4;
BERofCodedQAM_FreqSelec=zeros(length(SNR),1);
for i= 1:length(SNR)
    % generate randam data bits   
    data_FSC_QAM=(randn(NumOfBitsInBitStream,1)>0);   
    % repeat bits
    RepeatedData_FSC_QAM=zeros(NumOfBitsInBitStream*3,1);
    repeat=0;
    k=1;
    for h= 1:NumOfBitsInBitStream*3    
        if repeat == 0
            RepeatedData_FSC_QAM(h)=data_FSC_QAM(k);
            RepeatedBit=data_FSC_QAM(k);
            repeat=1+repeat;
            k=1+k;
        else
            RepeatedData_FSC_QAM(h)=RepeatedBit;
            repeat=1+repeat;
        end
        if repeat== 3
            repeat=0;
        end
    end   
    %interleaver 
    InterleavedData_FSC_QAM=zeros(1,NumOfBitsInBitStream*3);
    colsNum=16; rowsNum=16;  NumOfElements=16*16;
    for m=0:NumOfOFDMSymbols-1
        % convert the 16*16 elements from the data vector into matrix 16*16
        ReshapedMatrix_FSC_QAM=zeros(rowsNum,colsNum);
        for h=0:rowsNum-1
            ReshapedMatrix_FSC_QAM(h+1,:)=RepeatedData_FSC_QAM( h*colsNum+1+NumOfElements*m :...
                h*colsNum+colsNum+NumOfElements*m ,1);
        end
        % read column by column
        for h=0:colsNum-1
            InterleavedData_FSC_QAM(1, h*rowsNum+1+NumOfElements*m :...
                h*rowsNum+rowsNum+NumOfElements*m )...
                =ReshapedMatrix_FSC_QAM(:,h+1);
        end                     
    end
    % mapping
    QAMConstellation_FSC_QAM  = MapperFunction_16QAM( InterleavedData_FSC_QAM,...
        NumOfSymbols,E_forCoded );
    %IFFT
    %take 64 symbols and get 64 IFFT of them
    IFFToutput_FSC_QAM=zeros(NumOfSymbols,1);
    for k=0:(NumOfSymbols/64)-1
        IFFToutput_FSC_QAM((k*64)+1:(k*64)+64,1)=ifft...
            (QAMConstellation_FSC_QAM((k*64)+1 : (k*64)+64 ,1),N);
    end
    %add cyclic extension
    CyclicExtensionLen=N/4;
    NewLenAfterCyclicExtension = NumOfSymbols + (CyclicExtensionLen)*NumOfOFDMSymbols;
    TXOutput_FSC_QAM=zeros(1,NewLenAfterCyclicExtension);
    for k=0:( NumOfOFDMSymbols - 1 )
       buffer=zeros(1, N+CyclicExtensionLen );
       buffer( 1,1:CyclicExtensionLen  )=IFFToutput_FSC_QAM( ...
           (k*N)+(N-CyclicExtensionLen)+1 : k*N+N ,1);
       buffer( 1,CyclicExtensionLen+1:N+CyclicExtensionLen  )=...
           IFFToutput_FSC_QAM( (k*N)+1 : k*N+N ,1);
       TXOutput_FSC_QAM(1, k*(N+CyclicExtensionLen)+1 :...
           k*(N+CyclicExtensionLen)+N+CyclicExtensionLen )=buffer;  
    end
    % Frequency selective channel
    n=sqrt(Eb_forCoded/SNRInLinear(i)).*( randn(size(TXOutput_FSC_QAM)) + ...
        1i*randn(size(TXOutput_FSC_QAM)) );
    ChannelinputFSC=TXOutput_FSC_QAM+n;
    for k=0:NumOfOFDMSymbols-1
        buffer2FSC=ChannelinputFSC(1,k*(N+CyclicExtensionLen)+1 : ...
            k*(N+CyclicExtensionLen)+N+CyclicExtensionLen );
        buffer3FSC=conv(freqSelectiveChannel,buffer2FSC);
        buffer4FSC=deconv(buffer3FSC,freqSelectiveChannel);
        FFTInputFSC(1,k*(N+CyclicExtensionLen)+1 : k*(N+CyclicExtensionLen)+...
            N+CyclicExtensionLen)=buffer4FSC;
    end
    % get the FFT of OFDM symbols  
    SymbolsAfterFFT_FSC_QAM=zeros(NumOfSymbols,1);
    for k=0:( NewLenAfterCyclicExtension/(N+CyclicExtensionLen) - 1 )
       %update indeces
       SymBegin=k*(N+CyclicExtensionLen)+CyclicExtensionLen+1;
       symEnd=k*(N+CyclicExtensionLen)+(N+CyclicExtensionLen);
       %get FFT FOR one OFDM symbol
       SymbolsAfterFFT_FSC_QAM( k*N+1:k*N+N ,1)=fft( FFTInputFSC(1, SymBegin :symEnd ) ,N);  
    end
       % demapping
 DemapperOutput_FSC_QAM  = DemapperFunction_16QAM( SymbolsAfterFFT_FSC_QAM,NumOfSymbols,...
     NumOfBitsInBitStream*3,E_forCoded );
    DemapperOutput2_FSC_QAM=transpose(DemapperOutput_FSC_QAM);
    % deinterleaver
    for m=0:NumOfOFDMSymbols-1
    % convert the 256 elements from the FinalBits vector into matrix 8*16
        colsNum=16;  rowsNum=16;  NumofElems=16*16;
        ReshapedMatrix_FSC_QAM=zeros(rowsNum,colsNum);
        for h=0:colsNum-1
            ReshapedMatrix_FSC_QAM(:,h+1)=DemapperOutput2_FSC_QAM...
                ( h*rowsNum+1+NumofElems*m : h*rowsNum+rowsNum+NumofElems*m ,1);
        end
        % read row by row
        for h=0:rowsNum-1
            DeinterleavedData_FSC_QAM(  h*colsNum+1+NumofElems*m :...
                h*colsNum+colsNum+NumofElems*m ,1)...
                =ReshapedMatrix_FSC_QAM(h+1,:);
        end                     
    end
        k=1;
    while k <= (NumOfBitsInBitStream*3-2)
        err=0;
        if DeinterleavedData_FSC_QAM(k) ~= RepeatedData_FSC_QAM(k) 
            err=err+1;
        end
        k=k+1;
        if DeinterleavedData_FSC_QAM(k) ~= RepeatedData_FSC_QAM(k)  
            err=err+1;
        end
        k=k+1;
        if DeinterleavedData_FSC_QAM(k) ~= RepeatedData_FSC_QAM(k)  
            err=err+1;
        end
        k=k+1;
        if err==2 || err==3
           BERofCodedQAM_FreqSelec(i)=BERofCodedQAM_FreqSelec(i)+1;
        end 
    end
end
BERofCodedQAM_FreqSelec=BERofCodedQAM_FreqSelec./NumOfBitsInBitStream;

% plot results
figure
semilogy(SNR,BER_QAMUncodedflatchann,SNR,BERofMC_QAM_FC,'LineWidth', 2)
xlabel('Eb/No (dB)');ylabel('Bit error rate');
title('QAM (Multicarriers, flat channel)');
legend('Uncoded','coded');grid on;
figure
semilogy(SNR,BERofUncodedQAM_FreqSelec,SNR,BERofCodedQAM_FreqSelec,'LineWidth', 2)
xlabel('Eb/No (dB)');ylabel('Bit error rate');
title('QAM (Multicarriers, frequency selective channel)');
legend('Uncoded','coded');grid on;
figure
semilogy(SNR,BER_QAMUncodedflatchann,SNR,BERofMC_QAM_FC,SNR,...
    BERofUncodedQAM_FreqSelec,SNR,BERofCodedQAM_FreqSelec,'LineWidth', 2)
xlabel('Eb/No (dB)');ylabel('Bit error rate');title('QAM (Multicarriers)');
legend('Uncoded flat channel','coded flat channel',...
    'Uncoded freq selective channel','coded freq selective channel');grid on;