 clear all; clc;
%%%%%%%%%%%%%%%%%   Multicarriers QPSK  %%%%%%%%%%%%%%%%%%%%%
NumOfBitsInBitStream=998400;% A number that is divisible 128 (for ofdm symbols)
NumOfOFDMSymbols=NumOfBitsInBitStream/128;
SNR=-4:1:15; %(IN dB)
SNRInLinear=10.^((SNR/10));
N=64;
%% Uncoded (Flat Channel)
NumOfSymbols=NumOfBitsInBitStream/2;
NumOfOFDMSymbols=NumOfBitsInBitStream/128;
BERofUncodedQPSK_flatChan=zeros(1,length(SNR));
Eo=2; Eb=Eo/2;
% declare flat channel
flatFadingChannel=(randn()+1i*randn());
for i= 1:length(SNR)
    % generate randam data bits   
    data=(randn(NumOfBitsInBitStream,1)>0);   
    %interleaver
    InterleavedData=zeros(1,NumOfBitsInBitStream);
    for m=0:NumOfOFDMSymbols-1
        % convert the 128 elements from the data vector into matrix 8*16
        ReshapedMatrix=zeros(8,16);
        for h=0:7
            ReshapedMatrix(h+1,:)=data( h*16+1+128*m : h*16+16+128*m ,1);
        end
        % read column by column
        for h=0:15
            InterleavedData(1, h*8+1+128*m : h*8+8+128*m )=ReshapedMatrix(:,h+1);
        end                     
    end
    InputToMapper=InterleavedData;
    %mapping
    QPSKConstellation=zeros(NumOfSymbols,1);
    count=1;
    for k=1:2:NumOfBitsInBitStream
       if InputToMapper(k)==0 && InputToMapper(k+1)==0 
            NameOfSymbol=2;
       elseif InputToMapper(k)==0 && InputToMapper(k+1)==1
            NameOfSymbol=3;
       elseif InputToMapper(k)==1 && InputToMapper(k+1)==0
            NameOfSymbol=1;
       elseif InputToMapper(k)==1 && InputToMapper(k+1)==1
            NameOfSymbol=4;
       end
    QPSKConstellation(count,1)=sqrt(Eo)* (cos((2*NameOfSymbol-1)*(pi/4)) - ...
        1i*sin((2*NameOfSymbol-1)*(pi/4)));
    count=count+1;
    end
    %IFFT by taking each 64 symbols and get 64 IFFT of them
    IFFToutput=zeros(NumOfSymbols,1);
    for k=0:(NumOfSymbols/64)-1
        IFFToutput((k*64)+1:(k*64)+64,1)=ifft(QPSKConstellation((k*64)+1 : (k*64)+64 ,1),N);
    end
    %add cyclic extension
    CyclicExtensionLen=N/4;
    NewLenAfterCyclicExtension=NumOfSymbols+ (CyclicExtensionLen)*NumOfOFDMSymbols;
    TXOutput=zeros(1,NewLenAfterCyclicExtension);
    for k=0:( NumOfOFDMSymbols - 1 )
       buffer=zeros(1, N+CyclicExtensionLen );
       buffer( 1,1:CyclicExtensionLen  )=IFFToutput( (k*N)+(N-CyclicExtensionLen)+1 : k*N+N ,1);
       buffer( 1,CyclicExtensionLen+1:N+CyclicExtensionLen  )=IFFToutput( (k*N)+1 : k*N+N ,1);
       TXOutput(1,k*(N+CyclicExtensionLen)+1 : k*(N+CyclicExtensionLen)+...
           N+CyclicExtensionLen )=buffer;  
    end
    %% flat fading channel
    n=sqrt(Eb/SNRInLinear(i)).*( randn(size(TXOutput)) + 1i*randn(size(TXOutput)) );
    OutputFlatChannel=(flatFadingChannel.*TXOutput+n)./flatFadingChannel;
    % get the FFT of OFDM symbols  
    SymbolsAfterFFT=zeros(NumOfSymbols,1);
    for k=0:( NewLenAfterCyclicExtension/(N+CyclicExtensionLen) - 1 )
       %update indeces
       SymBegin=k*(N+CyclicExtensionLen)+CyclicExtensionLen+1;
       symEnd=k*(N+CyclicExtensionLen)+(N+CyclicExtensionLen);
       %get FFT FOR one OFDM symbol
       SymbolsAfterFFT( k*N+1:k*N+N ,1)=fft( OutputFlatChannel(1, SymBegin :symEnd ) ,N);  
    end
    %demapping
    received_sig=SymbolsAfterFFT;
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
    % deinterleaver
    k=0;
    for m=0:NumOfOFDMSymbols-1
        % convert the 128 elements from the FinalBits vector into matrix 8*16
        colsNum=16;
        rowsNum=8;
        NumofElems=128;
        ReshapedMatrix=zeros(rowsNum,colsNum);
        for h=0:colsNum-1
            ReshapedMatrix(:,h+1)=FinalBits( h*rowsNum+1+NumofElems*m : ...
                h*rowsNum+rowsNum+NumofElems*m ,1);
        end
        % read row by row
        for h=0:rowsNum-1
            DeinterleavedData(  h*colsNum+1+NumofElems*m :...
                h*colsNum+colsNum+NumofElems*m ,1)=ReshapedMatrix(h+1,:);
        end                     
    end
    BERofUncodedQPSK_flatChan(i)=sum( (DeinterleavedData~=data) );
end
BERofUncodedQPSK_flatChan=BERofUncodedQPSK_flatChan./(NumOfBitsInBitStream);
 



%% Coded (Flat Channel)
E_forCoded=Eo/3;  Eb_forCoded=E_forCoded/2;
NumOfOFDMSymbols=(NumOfBitsInBitStream*3)/128;
NumOfSymbols=(NumOfBitsInBitStream*3)/2;
BERofCodedQPSK_flatChann=zeros(length(SNR),1);
for i= 1:length(SNR)
    % generate randam data bits   
    data=(randn(NumOfBitsInBitStream,1)>0);   
    % repeat bits
    RepeatedData=zeros(NumOfBitsInBitStream*3,1);
    repeat=0;
    k=1;
    for h= 1:NumOfBitsInBitStream*3    
        if repeat == 0
            RepeatedData(h)=data(k);
            RepeatedBit=data(k);
            repeat=1+repeat;
            k=1+k;
        else
            RepeatedData(h)=RepeatedBit;
            repeat=1+repeat;
        end
        if repeat== 3
            repeat=0;
        end
    end
    %interleaver    
    InterleavedData=zeros(1,NumOfBitsInBitStream*3);
    for m=0:NumOfOFDMSymbols-1
        % convert the 128 elements from the RepeatedData vector into matrix 8*16
        ReshapedMatrix=zeros(8,16);
        for h=0:7
            ReshapedMatrix(h+1,:)=RepeatedData( h*16+1+128*m : h*16+16+128*m ,1);
        end
        % read column by column
        for h=0:15
            InterleavedData(1, h*8+1+128*m : h*8+8+128*m )=ReshapedMatrix(:,h+1);
        end                     
    end
    InputToMapper=InterleavedData;
    %mapping
    QPSKConstellation=zeros(NumOfSymbols,1);
    count=1;
    for k=1:2:NumOfBitsInBitStream*3
       if InputToMapper(k)==0 && InputToMapper(k+1)==0 
            NameOfSymbol=2;
       elseif InputToMapper(k)==0 && InputToMapper(k+1)==1
            NameOfSymbol=3;
       elseif InputToMapper(k)==1 && InputToMapper(k+1)==0
            NameOfSymbol=1;
       elseif InputToMapper(k)==1 && InputToMapper(k+1)==1
            NameOfSymbol=4;
       end
    QPSKConstellation(count,1)=sqrt(E_forCoded)* (cos((2*NameOfSymbol-1)*(pi/4)) - 1i*sin((2*NameOfSymbol-1)*(pi/4)));
    count=count+1;
    end
    %IFFT
    %take 64 symbols and get 64 IFFT of them
    IFFToutput=zeros(NumOfSymbols,1);
    for k=0:(NumOfSymbols/64)-1
        IFFToutput((k*64)+1:(k*64)+64,1)=ifft(QPSKConstellation((k*64)+1 : (k*64)+64 ,1),N);
    end
    %add cyclic extension
    CyclicExtensionLen=N/4;
    NewLenAfterCyclicExtension = NumOfSymbols + (CyclicExtensionLen)*NumOfOFDMSymbols;
    TXOutput=zeros(1,NewLenAfterCyclicExtension);
    for k=0:( NumOfOFDMSymbols - 1 )
       buffer=zeros(1, N+CyclicExtensionLen );
       buffer( 1,1:CyclicExtensionLen  )=IFFToutput( (k*N)+(N-CyclicExtensionLen)+1 : k*N+N ,1);
       buffer( 1,CyclicExtensionLen+1:N+CyclicExtensionLen  )=IFFToutput( (k*N)+1 : k*N+N ,1);
       TXOutput(1, k*(N+CyclicExtensionLen)+1 : k*(N+CyclicExtensionLen)+N+CyclicExtensionLen )=buffer;  
    end
    %% flat fading channel
    n=sqrt(Eb_forCoded/SNRInLinear(i)).*( randn(size(TXOutput)) + 1i*randn(size(TXOutput)) );
    OutputFlatChannel=(flatFadingChannel.*TXOutput+n)./flatFadingChannel; 
    % get the FFT of OFDM symbols  
    SymbolsAfterFFT=zeros(NumOfSymbols,1);
    for k=0:( NewLenAfterCyclicExtension/(N+CyclicExtensionLen) - 1 )
       %update indeces
       SymBegin=k*(N+CyclicExtensionLen)+CyclicExtensionLen+1;
       symEnd=k*(N+CyclicExtensionLen)+(N+CyclicExtensionLen);
       %get FFT FOR one OFDM symbol
       SymbolsAfterFFT( k*N+1:k*N+N ,1)=fft( OutputFlatChannel(1, SymBegin :symEnd ) ,N);  
    end
    %demapping
    received_sig=SymbolsAfterFFT;
    FinalBits=zeros(NumOfBitsInBitStream*3,1);
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
    % deinterleaver
    for m=0:NumOfOFDMSymbols-1
        % convert the 128 elements from the FinalBits vector into matrix 8*16
        colsNum=16; rowsNum=8; NumofElems=128;
        ReshapedMatrix=zeros(rowsNum,colsNum);
        for h=0:colsNum-1
            ReshapedMatrix(:,h+1)=FinalBits( h*rowsNum+1+NumofElems*m : h*rowsNum+rowsNum+NumofElems*m ,1);
        end
        % read row by row
        for h=0:rowsNum-1
            DeinterleavedData(  h*colsNum+1+NumofElems*m : h*colsNum+colsNum+NumofElems*m ,1)=ReshapedMatrix(h+1,:);
        end                     
    end
    %compute BER
    k=1;
    while k <= (NumOfBitsInBitStream*3-2)
        err=0;
        if DeinterleavedData(k) ~= RepeatedData(k) 
            err=err+1;
        end
        k=k+1;
        if DeinterleavedData(k) ~= RepeatedData(k)  
            err=err+1;
        end
        k=k+1;
        if DeinterleavedData(k) ~= RepeatedData(k)  
            err=err+1;
        end
        k=k+1;
        if err==2 || err==3
            BERofCodedQPSK_flatChann(i)=BERofCodedQPSK_flatChann(i)+1;
        end 
    end  
end
BERofCodedQPSK_flatChann=BERofCodedQPSK_flatChann./(NumOfBitsInBitStream);





%% Uncoded (Frequency selective Channel)
NumOfOFDMSymbols=NumOfBitsInBitStream/128;
NumOfSymbols=NumOfBitsInBitStream/2;
BERofUncodedQPSK_FreqSelec=zeros(length(SNR),1);
freqSelectiveChannel=[0.8 0 0 0 0 0 0 0 0 0 0.6];
for i= 1:length(SNR)
    % generate randam data bits   
    dataFS=(randn(NumOfBitsInBitStream,1)>0);   
    %interleaver
    InterleavedDataFS=zeros(1,NumOfBitsInBitStream);
    for m=0:NumOfOFDMSymbols-1
        % convert the 128 elements from the data vector into matrix 8*16
        ReshapedMatrixFS=zeros(8,16);
        for h=0:7
            ReshapedMatrixFS(h+1,:)=dataFS( h*16+1+128*m : h*16+16+128*m ,1);
        end
        % read column by column
        for h=0:15
            InterleavedDataFS(1, h*8+1+128*m : h*8+8+128*m )=ReshapedMatrixFS(:,h+1);
        end                     
    end
    InputToMapperFS=InterleavedDataFS;
    %mapping
    QPSKConstellationFS=zeros(NumOfSymbols,1);
    count=1;
    for k=1:2:NumOfBitsInBitStream
       if InputToMapperFS(k)==0 && InputToMapperFS(k+1)==0 
            NameOfSymbol=2;
       elseif InputToMapperFS(k)==0 && InputToMapperFS(k+1)==1
            NameOfSymbol=3;
       elseif InputToMapperFS(k)==1 && InputToMapperFS(k+1)==0
            NameOfSymbol=1;
       elseif InputToMapperFS(k)==1 && InputToMapperFS(k+1)==1
            NameOfSymbol=4;
       end
    QPSKConstellationFS(count,1)=sqrt(Eo)* (cos((2*NameOfSymbol-1)*(pi/4)) - 1i*sin((2*NameOfSymbol-1)*(pi/4)));
    count=count+1;
    end
    %IFFT by taking each 64 symbols and get 64 IFFT of them
    IFFToutputFS=zeros(NumOfSymbols,1);
    for k=0:(NumOfSymbols/64)-1
        IFFToutputFS((k*64)+1:(k*64)+64,1)=ifft(QPSKConstellationFS((k*64)+1 : (k*64)+64 ,1),N);
    end
    %add cyclic extension
    CyclicExtensionLen=N/4;
    NewLenAfterCyclicExtension=NumOfSymbols+ (CyclicExtensionLen)*NumOfOFDMSymbols;
    TXOutputFS=zeros(1,NewLenAfterCyclicExtension);
    for k=0:( NumOfOFDMSymbols - 1 )
       bufferFS=zeros(1, N+CyclicExtensionLen );
       bufferFS( 1,1:CyclicExtensionLen  )=IFFToutputFS( (k*N)+(N-CyclicExtensionLen)+1 : k*N+N ,1);
       bufferFS( 1,CyclicExtensionLen+1:N+CyclicExtensionLen  )=IFFToutputFS( (k*N)+1 : k*N+N ,1);
       TXOutputFS(1,k*(N+CyclicExtensionLen)+1 : k*(N+CyclicExtensionLen)+N+CyclicExtensionLen )=bufferFS;  
    end
    % Frequency selective channel
    n=sqrt((Eo)/SNRInLinear(i)).*( randn(size(TXOutputFS)) + 1i*randn(size(TXOutputFS)) );
    ChannelinputFS=TXOutputFS+n;
    for k=0:NumOfOFDMSymbols-1
        buffer2FS=ChannelinputFS(1,k*(N+CyclicExtensionLen)+1 : k*(N+CyclicExtensionLen)+N+CyclicExtensionLen );
        buffer3FS=conv(freqSelectiveChannel,buffer2FS);
        buffer4FS=deconv(buffer3FS,freqSelectiveChannel);
        FFTInputFS(1,k*(N+CyclicExtensionLen)+1 : k*(N+CyclicExtensionLen)+N+CyclicExtensionLen)=buffer4FS;
    end
    % get the FFT of OFDM symbols  
    SymbolsAfterFFTFS=zeros(NumOfSymbols,1);
    for k=0:( NewLenAfterCyclicExtension/(N+CyclicExtensionLen) - 1 )
       %update indeces
       SymBegin=k*(N+CyclicExtensionLen)+CyclicExtensionLen+1;
       symEnd=k*(N+CyclicExtensionLen)+(N+CyclicExtensionLen);
       %get FFT FOR one OFDM symbol
       SymbolsAfterFFTFS( k*N+1:k*N+N ,1)=fft( FFTInputFS(1, SymBegin :symEnd ) ,N);  
    end
    %demapping
    received_sig=SymbolsAfterFFTFS;
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
    % deinterleaver
    k=0;
    for m=0:NumOfOFDMSymbols-1
        % convert the 128 elements from the FinalBits vector into matrix 8*16
        colsNum=16;  rowsNum=8;  NumofElems=128;
        ReshapedMatrixFS=zeros(rowsNum,colsNum);
        for h=0:colsNum-1
            ReshapedMatrixFS(:,h+1)=FinalBits( h*rowsNum+1+NumofElems*m : h*rowsNum+rowsNum+NumofElems*m ,1);
        end
        % read row by row
        for h=0:rowsNum-1
            DeinterleavedDataFS(  h*colsNum+1+NumofElems*m : h*colsNum+colsNum+NumofElems*m ,1)=ReshapedMatrixFS(h+1,:);
        end                     
    end
    BERofUncodedQPSK_FreqSelec(i)=sum( (DeinterleavedDataFS~=dataFS) );
end
BERofUncodedQPSK_FreqSelec=BERofUncodedQPSK_FreqSelec./(NumOfBitsInBitStream);

%% Coded (Frequency selective Channel)
E_forCoded=Eo/3;  Eb_forCoded=E_forCoded/2;
NumOfOFDMSymbols=(NumOfBitsInBitStream*3)/128;
NumOfSymbols=(NumOfBitsInBitStream*3)/2;
BERofCodedQPSK_FreqSelec=zeros(length(SNR),1);
for i= 1:length(SNR)
    % generate randam data bits   
    data_CFS=(randn(NumOfBitsInBitStream,1)>0);   
    % repeat bits
    RepeatedData_CFS=zeros(NumOfBitsInBitStream*3,1);
    repeat=0;
    k=1;
    for h= 1:NumOfBitsInBitStream*3    
        if repeat == 0
            RepeatedData_CFS(h)=data_CFS(k);
            RepeatedBit=data_CFS(k);
            repeat=1+repeat;
            k=1+k;
        else
            RepeatedData_CFS(h)=RepeatedBit;
            repeat=1+repeat;
        end
        if repeat== 3
            repeat=0;
        end
    end
    %interleaver    
    InterleavedData_CFS=zeros(1,NumOfBitsInBitStream*3);
    for m=0:NumOfOFDMSymbols-1
        % convert the 128 elements from the RepeatedData vector into matrix 8*16
        ReshapedMatrix_CFS=zeros(8,16);
        for h=0:7
            ReshapedMatrix_CFS(h+1,:)=RepeatedData_CFS( h*16+1+128*m : h*16+16+128*m ,1);
        end
        % read column by column
        for h=0:15
            InterleavedData_CFS(1, h*8+1+128*m : h*8+8+128*m )=ReshapedMatrix_CFS(:,h+1);
        end                     
    end
    InputToMapper_CFS=InterleavedData_CFS;
    %mapping
    QPSKConstellation_CFS=zeros(NumOfSymbols,1);
    count=1;
    for k=1:2:NumOfBitsInBitStream*3
       if InputToMapper_CFS(k)==0 && InputToMapper_CFS(k+1)==0 
            NameOfSymbol=2;
       elseif InputToMapper_CFS(k)==0 && InputToMapper_CFS(k+1)==1
            NameOfSymbol=3;
       elseif InputToMapper_CFS(k)==1 && InputToMapper_CFS(k+1)==0
            NameOfSymbol=1;
       elseif InputToMapper_CFS(k)==1 && InputToMapper_CFS(k+1)==1
            NameOfSymbol=4;
       end
    QPSKConstellation_CFS(count,1)=sqrt(E_forCoded)* (cos((2*NameOfSymbol-1)*(pi/4)) - 1i*sin((2*NameOfSymbol-1)*(pi/4)));
    count=count+1;
    end
    %IFFT
    %take 64 symbols and get 64 IFFT of them
    IFFToutput_CFS=zeros(NumOfSymbols,1);
    for k=0:(NumOfSymbols/64)-1
        IFFToutput_CFS((k*64)+1:(k*64)+64,1)=ifft(QPSKConstellation_CFS((k*64)+1 : (k*64)+64 ,1),N);
    end
    %add cyclic extension
    CyclicExtensionLen=N/4;
    NewLenAfterCyclicExtension = NumOfSymbols + (CyclicExtensionLen)*NumOfOFDMSymbols;
    TXOutput_CFS=zeros(1,NewLenAfterCyclicExtension);
    for k=0:( NumOfOFDMSymbols - 1 )
       buffer=zeros(1, N+CyclicExtensionLen );
       buffer( 1,1:CyclicExtensionLen  )=IFFToutput_CFS( (k*N)+(N-CyclicExtensionLen)+1 : k*N+N ,1);
       buffer( 1,CyclicExtensionLen+1:N+CyclicExtensionLen  )=IFFToutput_CFS( (k*N)+1 : k*N+N ,1);
       TXOutput_CFS(1, k*(N+CyclicExtensionLen)+1 : k*(N+CyclicExtensionLen)+N+CyclicExtensionLen )=buffer;  
    end
    %% Frequency selective channel
    n_CFS=sqrt(Eb_forCoded/SNRInLinear(i)).*( randn(size(TXOutput_CFS)) + 1i*randn(size(TXOutput_CFS)) );
    Channelinput_CFS=TXOutput_CFS+n_CFS;
    for k=0:NumOfOFDMSymbols-1
        buffer2_CFS=Channelinput_CFS(1,k*(N+CyclicExtensionLen)+1 : k*(N+CyclicExtensionLen)+N+CyclicExtensionLen );
        buffer3_CFS=conv(freqSelectiveChannel,buffer2_CFS);
        buffer4_CFS=deconv(buffer3_CFS,freqSelectiveChannel);
        FFTInput_CFS(1,k*(N+CyclicExtensionLen)+1 : k*(N+CyclicExtensionLen)+N+CyclicExtensionLen)=buffer4_CFS;
    end 
    % get the FFT of OFDM symbols  
    SymbolsAfterFFT_CFS=zeros(NumOfSymbols,1);
    for k=0:( NewLenAfterCyclicExtension/(N+CyclicExtensionLen) - 1 )
       %update indeces
       SymBegin=k*(N+CyclicExtensionLen)+CyclicExtensionLen+1;
       symEnd=k*(N+CyclicExtensionLen)+(N+CyclicExtensionLen);
       %get FFT FOR one OFDM symbol
       SymbolsAfterFFT_CFS( k*N+1:k*N+N ,1)=fft( FFTInput_CFS(1, SymBegin :symEnd ) ,N);  
    end
    %demapping
    received_sig=SymbolsAfterFFT_CFS;
    FinalBits=zeros(NumOfBitsInBitStream*3,1);
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
    % deinterleaver
    for m=0:NumOfOFDMSymbols-1
        % convert the 128 elements from the FinalBits vector into matrix 8*16
        colsNum=16; rowsNum=8; NumofElems=128;
        ReshapedMatrix_CFS=zeros(rowsNum,colsNum);
        for h=0:colsNum-1
            ReshapedMatrix_CFS(:,h+1)=FinalBits( h*rowsNum+1+NumofElems*m : h*rowsNum+rowsNum+NumofElems*m ,1);
        end
        % read row by row
        for h=0:rowsNum-1
            DeinterleavedData_CFS(  h*colsNum+1+NumofElems*m : h*colsNum+colsNum+NumofElems*m ,1)=ReshapedMatrix_CFS(h+1,:);
        end                     
    end
    %compute BER
    k=1;
    while k <= (NumOfBitsInBitStream*3-2)
        err=0;
        if DeinterleavedData_CFS(k) ~= RepeatedData_CFS(k) 
            err=err+1;
        end
        k=k+1;
        if DeinterleavedData_CFS(k) ~= RepeatedData_CFS(k)  
            err=err+1;
        end
        k=k+1;
        if DeinterleavedData_CFS(k) ~= RepeatedData_CFS(k)  
            err=err+1;
        end
        k=k+1;
        if err==2 || err==3
            BERofCodedQPSK_FreqSelec(i)=BERofCodedQPSK_FreqSelec(i)+1;
        end 
    end  
end
BERofCodedQPSK_FreqSelec=BERofCodedQPSK_FreqSelec./(NumOfBitsInBitStream);



% plot results
figure
semilogy(SNR,BERofUncodedQPSK_flatChan,SNR,BERofCodedQPSK_flatChann,'LineWidth', 2)
xlabel('Eb/No (dB)');ylabel('Bit error rate');
title('QPSK (Multicarriers,flat channel)');
legend('Uncoded','coded');grid on;
figure
semilogy(SNR,BERofUncodedQPSK_FreqSelec,SNR,BERofCodedQPSK_FreqSelec,'LineWidth', 2)
xlabel('Eb/No (dB)');ylabel('Bit error rate');
title('QPSK (Multicarriers, frequency selective channel)');
legend('Uncoded','coded');grid on;
figure
semilogy(SNR,BERofUncodedQPSK_flatChan,SNR,BERofCodedQPSK_flatChann,SNR,...
    BERofUncodedQPSK_FreqSelec,SNR,BERofCodedQPSK_FreqSelec,'LineWidth', 2)
xlabel('Eb/No (dB)');ylabel('Bit error rate');
title('QPSK (Multicarriers, frequency selective channel)');
legend('Uncoded flat channel','coded flat channel','Uncoded freq selective channel',...
    'coded freq selective channel');grid on;