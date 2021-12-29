function [ DemapperOutput ] = DemapperFunction_16QAM( SymbolsAfterFFT,NumOfSymbols,NumOfBitsInBitStream,Eo )
   % demapping
    DemapperInput(:,1)=real( SymbolsAfterFFT);
    DemapperInput(:,2)=imag( SymbolsAfterFFT);
    DemapperOutput = zeros(1,NumOfBitsInBitStream);
    v=1;
    for u = 1:NumOfSymbols
        %first row 
        if(DemapperInput(u,1)>2*sqrt(Eo) && DemapperInput(u,2)>2*sqrt(Eo))
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=1; v=v+1; 
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=1;v=v+1;
        elseif(DemapperInput(u,1)>2*sqrt(Eo) && DemapperInput(u,2)<2*sqrt(Eo)...
                && DemapperInput(u,2)>0)
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=1; v=v+1;  
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=1;v=v+1;
        elseif(DemapperInput(u,1)>2*sqrt(Eo) && DemapperInput(u,2)<0 ...
                && DemapperInput(u,2)>-2*sqrt(Eo) )
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=1; v=v+1; 
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=0;v=v+1;
        elseif(DemapperInput(u,1)>2*sqrt(Eo) && DemapperInput(u,2)<-2*sqrt(Eo) )
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=1; v=v+1; 
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=0;v=v+1;
        %second row
        elseif(DemapperInput(u,1)<2*sqrt(Eo)  && DemapperInput(u,1)>0 ...
                && DemapperInput(u,2)>2*sqrt(Eo))
           DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=1; v=v+1; 
           DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=1;v=v+1;
        elseif(DemapperInput(u,1)<2*sqrt(Eo) && DemapperInput(u,1)>0 ...
                && DemapperInput(u,2)<2*sqrt(Eo) && DemapperInput(u,2)>0)
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=1; v=v+1; 
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=1;v=v+1;
        elseif(DemapperInput(u,1)<2*sqrt(Eo) && DemapperInput(u,1)>0 ...
            && DemapperInput(u,2)<0 && DemapperInput(u,2)>-2*sqrt(Eo) )
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=1; v=v+1; 
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=0;v=v+1;
        elseif(DemapperInput(u,1)<2*sqrt(Eo) && DemapperInput(u,1)>0 ...
                && DemapperInput(u,2)<-2*sqrt(Eo) )
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=1; v=v+1; 
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=0;v=v+1;
        %third row
        elseif(DemapperInput(u,1)<0 && DemapperInput(u,1)>-2*sqrt(Eo)  ...
                && DemapperInput(u,2)>2*sqrt(Eo))
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=0; v=v+1; 
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=1;v=v+1;
        elseif(DemapperInput(u,1)<0 && DemapperInput(u,1)>-2*sqrt(Eo) ...
                && DemapperInput(u,2)<2*sqrt(Eo) && DemapperInput(u,2)>0)
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=0; v=v+1; 
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=1;v=v+1;
        elseif(DemapperInput(u,1)<0 && DemapperInput(u,1)>-2*sqrt(Eo) ...
                && DemapperInput(u,2)<0 && DemapperInput(u,2)>-2*sqrt(Eo) )
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=0; v=v+1; 
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=0;v=v+1;
        elseif(DemapperInput(u,1)<0 && DemapperInput(u,1)>-2*sqrt(Eo) ...
                && DemapperInput(u,2)<-2*sqrt(Eo) )
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=0; v=v+1; 
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=0;v=v+1;
        %fourh row
        elseif(DemapperInput(u,1)<-2*sqrt(Eo) && DemapperInput(u,2)>2*sqrt(Eo))
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=0; v=v+1; 
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=1;v=v+1;
        elseif(DemapperInput(u,1)<-2*sqrt(Eo) && DemapperInput(u,2)<2*sqrt(Eo) ...
                && DemapperInput(u,2)>0)
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=0; v=v+1;
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=1;v=v+1;
        elseif(DemapperInput(u,1)<-2*sqrt(Eo) && DemapperInput(u,2)>-2*sqrt(Eo) ...
                && DemapperInput(u,2)<0)
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=0; v=v+1; 
            DemapperOutput(1,v)=1; v=v+1; DemapperOutput(1,v)=0;v=v+1;
        elseif(DemapperInput(u,1)<-2*sqrt(Eo) && DemapperInput(u,2)<-2*sqrt(Eo) )
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=0; v=v+1; 
            DemapperOutput(1,v)=0; v=v+1; DemapperOutput(1,v)=0;v=v+1;
        end
    end


end

