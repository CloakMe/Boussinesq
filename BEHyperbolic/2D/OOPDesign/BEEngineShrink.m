classdef (ConstructOnLoad) BEEngineShrink < BEEngine
   % Class help goes here
  properties ( SetAccess = protected, GetAccess = public)
    dhb
  end 
  
  methods
    
    function this = BEEngineShrink( dscrtParams, eqParams, ic )
        this = this@BEEngine( dscrtParams, eqParams, ic );
        [ ff, this.dhb ] = BEUtilities.GetFinDiffMat( this.sx, this.h );
    end
    
    function [ this ] = UpdateMatricesSizes( this, newSizeX, newSizeY )
        %tic
        numAuxPnts = this.order / 2;
        this.dhb = this.dhb( numAuxPnts+1:end - numAuxPnts, numAuxPnts+1:end - numAuxPnts );
        [ this.eigenFinDiffMat, w ] = eig( -this.dhb );

        [ this.IminusDHdiag, this.minusDHdiag ] =...
            this.DiagonalizeAndGetDiagOfFinDiffMat( newSizeX );
        this.vdah = zeros( newSizeX, newSizeY );
        this.sx = newSizeX;
        this.sy = newSizeY;
        this.x = this.x( numAuxPnts+1:end - numAuxPnts );
        this.y = this.y( numAuxPnts+1:end - numAuxPnts );
        %toc
    end
    
    function [ mat1, mat2, mat3, mat4 ] = UpdateMatrices( this, coeff, mat1, mat2, mat3, mat4 )
        
        numAuxPnts = this.order / 2;
        mat1 = mat1( coeff*numAuxPnts+1:end - coeff*numAuxPnts,...
                      coeff*numAuxPnts+1:end - coeff*numAuxPnts );
        mat2 = mat2( coeff*numAuxPnts+1:end - coeff*numAuxPnts,...
                      coeff*numAuxPnts+1:end - coeff*numAuxPnts );
        
        if( nargin == 4 )
            mat3 = 0; mat4 = 0;
            return;
        end
        
        mat3 = mat3( coeff*numAuxPnts+1:end - coeff*numAuxPnts,...
                      coeff*numAuxPnts+1:end - coeff*numAuxPnts );
        if( nargin == 5 )
            mat4 = 0;
            return;
        end
        mat4 = mat4( coeff*numAuxPnts+1:end - coeff*numAuxPnts,...
                      coeff*numAuxPnts+1:end - coeff*numAuxPnts );
    end
    
  end
  
  methods ( Abstract = true )
    % virtual
    BESolver( this )
  end
  
end
    