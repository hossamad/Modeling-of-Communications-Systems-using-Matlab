function [ QAMConstellation ] = MapperFunction_16QAM( InputToMapper,NumOfSymbols,Eo )
    a=zeros(NumOfSymbols,1);  b=zeros(NumOfSymbols,1); 
    % reshaping
    for r = 0 : 1 :(NumOfSymbols-1)
        Data_bits_reshaped(r+1,:) = InputToMapper(1,(4*r)+1 :(4*r)+4);
    end
    for r=1:NumOfSymbols
        if Data_bits_reshaped(r,1)==0 && Data_bits_reshaped(r,2)==0 &&...
                Data_bits_reshaped(r,3)==0 && Data_bits_reshaped(r,4)==0
            a(r,1)=-3;     b(r,1)=-3;
        elseif Data_bits_reshaped(r,1)==0 && Data_bits_reshaped(r,2)==0 &&...
                Data_bits_reshaped(r,3)==0 && Data_bits_reshaped(r,4)==1
            a(r,1)=-3 ;     b(r,1)=3 ;
        elseif Data_bits_reshaped(r,1)==0 && Data_bits_reshaped(r,2)==0 &&...
                Data_bits_reshaped(r,3)==1 && Data_bits_reshaped(r,4)==0
            a(r,1)=-3 ;     b(r,1)=-1 ;
        elseif Data_bits_reshaped(r,1)==0 && Data_bits_reshaped(r,2)==0 &&...
                Data_bits_reshaped(r,3)==1 && Data_bits_reshaped(r,4)==1
            a(r,1)=-3 ;     b(r,1)=1 ;
        elseif Data_bits_reshaped(r,1)==0 && Data_bits_reshaped(r,2)==1 &&...
                Data_bits_reshaped(r,3)==0 && Data_bits_reshaped(r,4)==0
            a(r,1)=3 ;     b(r,1)=-3 ;
        elseif Data_bits_reshaped(r,1)==0 && Data_bits_reshaped(r,2)==1 &&...
                Data_bits_reshaped(r,3)==0 && Data_bits_reshaped(r,4)==1
            a(r,1)=3 ;     b(r,1)=3 ;
        elseif Data_bits_reshaped(r,1)==0 && Data_bits_reshaped(r,2)==1 &&...
                Data_bits_reshaped(r,3)==1 && Data_bits_reshaped(r,4)==0
            a(r,1)=3 ;     b(r,1)=-1 ;
        elseif Data_bits_reshaped(r,1)==0 && Data_bits_reshaped(r,2)==1 &&...
                Data_bits_reshaped(r,3)==1 && Data_bits_reshaped(r,4)==1
            a(r,1)=3 ;     b(r,1)=1 ;
        elseif Data_bits_reshaped(r,1)==1 && Data_bits_reshaped(r,2)==0 &&...
                Data_bits_reshaped(r,3)==0 && Data_bits_reshaped(r,4)==0
            a(r,1)=-1 ;     b(r,1)=-3 ;
        elseif Data_bits_reshaped(r,1)==1 && Data_bits_reshaped(r,2)==0 &&...
                Data_bits_reshaped(r,3)==0 && Data_bits_reshaped(r,4)==1
            a(r,1)=-1 ;     b(r,1)=3 ;
        elseif Data_bits_reshaped(r,1)==1 && Data_bits_reshaped(r,2)==0 &&...
                Data_bits_reshaped(r,3)==1 && Data_bits_reshaped(r,4)==0
            a(r,1)=-1 ;     b(r,1)=-1 ;
        elseif Data_bits_reshaped(r,1)==1 && Data_bits_reshaped(r,2)==0 &&...
                Data_bits_reshaped(r,3)==1 && Data_bits_reshaped(r,4)==1
            a(r,1)=-1 ;     b(r,1)=1 ;
        elseif Data_bits_reshaped(r,1)==1 && Data_bits_reshaped(r,2)==1 &&...
                Data_bits_reshaped(r,3)==0 && Data_bits_reshaped(r,4)==0
            a(r,1)=1 ;     b(r,1)=-3 ;
        elseif Data_bits_reshaped(r,1)==1 && Data_bits_reshaped(r,2)==1 &&...
            Data_bits_reshaped(r,3)==0 && Data_bits_reshaped(r,4)==1
            a(r,1)=1 ;     b(r,1)=3 ;
        elseif Data_bits_reshaped(r,1)==1 && Data_bits_reshaped(r,2)==1 &&...
                Data_bits_reshaped(r,3)==1 && Data_bits_reshaped(r,4)==0
            a(r,1)=1 ;     b(r,1)=-1 ;
        elseif Data_bits_reshaped(r,1)==1 && Data_bits_reshaped(r,2)==1 &&...
                Data_bits_reshaped(r,3)==1 && Data_bits_reshaped(r,4)==1
            a(r,1)=1 ;     b(r,1)=1 ;               
        end
    end
    QAMConstellation = sqrt(Eo) * (a+1i*b);

end

