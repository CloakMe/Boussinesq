classdef TestBEDomainUtilsP2Spec
    methods(Static)          
        function [ result ] = SimpleBEDomainTestYder1(st,en)
            [ domainUtilsP2Spec ] = TestBEDomainUtilsP2Spec.GetDomainUtilsP2SpecInst(st,en,2);         
            
            M = [  1    1    1;...
                   2    2    2;...
                   1    1    2.5];   
            [ res ] = domainUtilsP2Spec.YDer( M, [1 -2 1] );
            expectedxRes = [  0     0     0;...
                              0     0     0;...
                              0     1.5   0];                        
            result = min(min(res == expectedxRes) );
            
            if( result == 1 )
                disp( 'SimpleBEDomainTestYder1 is OK.' );
            else
                disp( 'SimpleBEDomainTestYder1 Failed!' );
            end
        end 
        
        function [ result ] = SimpleBEDomainTestYder2(st,en)
            
            [ domainUtilsP2Spec ] = TestBEDomainUtilsP2Spec.GetDomainUtilsP2SpecInst(st,en,2);        
            
            M = [  1    1    1    2    1;...
                   2    2    2    2    2;...
                   1    1    1    1    1;...
                   2    2    2    1    1;...
                   2    2    2    1    1];   
            [ res ] = domainUtilsP2Spec.YDer( M, [1 -2 1] );
            expectedxRes = [   0    0    1   -2    0;...
                               0    0    0    0    0;...
                               0    0    0    0    0;...
                               0    0   -1    1    0;...
                               0    0   -1    1    0];                         
            result = min(min(res == expectedxRes) );
            
            if( result == 1 )
                disp( 'SimpleBEDomainTestYder2 is OK.' );
            else
                disp( 'SimpleBEDomainTestYder2 Failed!' );
            end
        end 
        
        function [ result ] = SimpleBEDomainTestYder3(st,en)
            
            [ domainUtilsP2Spec ] = TestBEDomainUtilsP2Spec.GetDomainUtilsP2SpecInst(st,en,4);        
            
            M = [  1    1    1    2    1;...
                   2    2    2    2    2;...
                   1    1    1    1    1;...
                   2    2    2    1    1;...
                   2    2    2    1    1];   
            [ res ] = domainUtilsP2Spec.YDer( M, [-1/12    4.0/3   -2.5000    4.0/3   -1/12] );
            expectedxRes = [   0    0    4/3   0    0;...
                               0    0    0     0    0;...
                               0    0    0     0    0;...
                               0    0   -1.25  0    0;...
                               0    0   -1.25  0    0];                         
            result = min(min( abs( res - expectedxRes ) < 10^(-10) ) );
            
            if( result == 1 )
                disp( 'SimpleBEDomainTestYderiv1 is OK.' );
            else
                disp( 'SimpleBEDomainTestYderiv1 Failed!' );
            end
        end 
        
        function [ result ] = SimpleBEDomainTestYder4(st,en)
            
            [ domainUtilsP2Spec ] = TestBEDomainUtilsP2Spec.GetDomainUtilsP2SpecInst(st,en,2);        
            
            M = [  1    1    1    2    1;...
                   2    2    2    2    2;...
                   1    1    1    1    1;...
                   2    2    2    1    1;...
                   2    2    2    1    1];   
               
            domainUtilsP2Spec = domainUtilsP2Spec.SetDeltaTimeDerivatives( M, 1.0 );
            [ res ] = domainUtilsP2Spec.YDeriv( M, [1 -2 1], 0 );
            
            expectedxRes = [   0    0    1   -2    0;...
                               0    0    0    0    0;...
                               0    0    0    0    0;...
                               0    0   -1    1    0;...
                               0    0   -1    1    0];                         
            result = min(min(res == expectedxRes) );
            
            if( result == 1 )
                disp( 'SimpleBEDomainTestYderiv2 is OK.' );
            else
                disp( 'SimpleBEDomainTestYderiv2 Failed!' );
            end
        end 
        
        function [ result ] = SimpleBEDomainTestYderiv1(st,en)
            
            [ domainUtilsP2Spec ] = TestBEDomainUtilsP2Spec.GetDomainUtilsP2SpecInst(st,en,4);        
            
            M = [  1    1    1    2    1;...
                   2    2    2    2    2;...
                   1    1    1    1    1;...
                   2    2    2    1    1;...
                   2    2    2    1    1];   
               
            domainUtilsP2Spec = domainUtilsP2Spec.SetDeltaTimeDerivatives( M, 1.0 );
            [ res ] = domainUtilsP2Spec.YDeriv( M, [-1/12    4.0/3   -2.5000    4.0/3   -1/12], 0 );
            
            expectedxRes = [   -0.3049    0.2028    1.3333   -0.0560    0.7010;...
                               -0.0005   -1.1354   -0.0000    2.1086   -0.2493;...
                                0.1021    0.3195   -0.0000    1.1020    0.1869;...
                               -0.0005   -1.1354   -1.2500    1.0543   -0.1246;...
                               -0.6099    0.4056   -1.2500   -0.0280    0.7010];                         
            result = min(min(  abs( res - expectedxRes ) < 10^(-3)  ) );
            
            if( result == 1 )
                disp( 'SimpleBEDomainTestYderiv2 is OK.' );
            else
                disp( 'SimpleBEDomainTestYderiv2 Failed!' );
            end
        end 
        
    end
    
    methods(Static, Access = private)       

        function [ domainUtilsP2Spec ] = GetDomainUtilsP2SpecInst(st,en,ord)
            x = st:en;
            y = st:en;
            domainUtilsP2Spec = BEDomainUtilsP2Spec( x, y, ord, 3, 0.3, 0.1, 0.2 );
        end 
    end
end