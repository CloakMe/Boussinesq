function [outMatrix]=GenerateBandDiagMatrix(band, matrixRowSize, typeStart, typeEnd)

    if(nargin == 2)
        typeStart = 'Identity';
        typeEnd = 'Identity';
    end
    if(nargin == 3)
        typeEnd = 'Identity';
    end
    shift = floor(length(band)/2);
    extension = 2*shift;
    outMatrix = zeros(matrixRowSize, matrixRowSize);
    if(size(band,1) == 1 || size(band,2) == 1)
        if(strcmp(typeStart, 'Identity'))
            for i=1:shift
                outMatrix(i,i)=1;
            end
        elseif(strcmp(typeStart, 'dudxZero'))
            centralElement = (length(band)+1)/2;
            for i = 1:shift      
                pos = length(band(centralElement:-1:1)) - 1;
                outMatrix(i,1:i+pos) = band(centralElement+1-i:end);
                pos = length(band(centralElement-i:-1:1)) - 1;
                outMatrix(i,2:2+pos) = outMatrix(i,2:2+pos) + band(centralElement-i:-1:1);
            end
        end
        for i=1 + shift:matrixRowSize - shift
            outMatrix(i,i-shift:i+extension-shift) = band;
        end
        if(strcmp(typeEnd, 'Identity')) 
            for i=matrixRowSize - shift + 1:matrixRowSize
                outMatrix(i,i)=1;
            end
        elseif(strcmp(typeEnd, 'dudxZero') )
            centralElement = (length(band)+1)/2;
            endEl = length(band) - 1;
            for i=matrixRowSize - shift + 1:matrixRowSize
                pos = length(band(1:endEl)) - 1;
                outMatrix(i,end-pos:end) =band(1:endEl);
                startEl = matrixRowSize + 1 - i;
                pos = length(band(endEl+1:end)) - 1;
                outMatrix(i,end-1-pos:end-1) = outMatrix(i,end-1-pos:end-1) + band(end:-1:endEl+1);
                endEl = endEl - 1;
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