function [outMatrix]=GenerateBandDiagMatrix(band, matrixRowSize, typeStart, typeEnd)

    shift = floor(length(band)/2);
    extension = 2*shift;
    outMatrix = zeros(matrixRowSize, matrixRowSize);
    if(size(band,1) == 1 || size(band,2) == 1)
        if(typeStart == 'Identity')
            for i=1:shift
                outMatrix(i,i)=1;
            end
        elseif(typeStart == 'dudxZero')            
            centralElement = (length(band)+1)/2;
            for i = 1:shift       
                pos = length(band(centralElement:-1:1)) - 1;
                outMatrix(i,1:i+pos) =  band(centralElement+1-i:end);
                pos = length(band(centralElement-i:-1:1)) - 1;
                outMatrix(i,2:2+pos) = outMatrix(i,i+1:i+1+pos) + band(centralElement-i:-1:1);
            end
        end
        for i=1 + shift:matrixRowSize - shift
            outMatrix(i,i-shift:i+extension-shift) = band;
        end
        if(typeEnd == 'Identity')
            for i=matrixRowSize - shift + 1:matrixRowSize
                outMatrix(i,i)=1;
            end
        elseif(typeEnd == 'dudxZero')
            centralElement = (length(band)+1)/2;
            for i = 1:shift       
                pos = length(band(centralElement:-1:1)) - 1;
                outMatrix(i,end-pos-i+1:end) =  band(centralElement+1-i:end);
                pos = length(band(centralElement-i:-1:1)) - 1;
                outMatrix(i,end-pos-1:end) = outMatrix(i,end-pos-1:end) + band(centralElement-i:-1:1);
            end
        end
    else    
        fdSize = size(band);
        centralFiniteDiffPosition = (fdSize(1)+1)/2;
        midPoint = centralFiniteDiffPosition; %fdSize(2)/2;
        shift = floor(fdSize(1)/2);
        extension = 2*shift;
        for i=1:shift
            %outMatrix(i,1:1+fdSize(2)-1) = band(i,:);
            outMatrix(i,i)=1;
        end
        for i=1 + shift:matrixRowSize - shift
            outMatrix(i,i-shift:i+extension-shift) = band(centralFiniteDiffPosition,1:fdSize(1));
        end
        %iBand = centralFiniteDiffPosition;
        for i=matrixRowSize - shift + 1:matrixRowSize
            %iBand = iBand+1;
            %outMatrix(i,end-fdSize(2)+1:end) = band(iBand,:);
            outMatrix(i,i)=1;
        end
    end
end