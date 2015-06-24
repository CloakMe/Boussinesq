classdef TestBEDomainUtilsP2edges
    methods(Static)          
        function [ result ] = SimpleBEDomainTestYder1(st,en)
            [ domainUtilsP2edges ] = TestBEDomainUtilsP2edges.GetDomainUtilsP2edgesInst(st,en,2);         
            
            M = [  1    1    1;...
                   2    2    2;...
                   1    1    2.5];   
               
            domainUtilsP2edges = domainUtilsP2edges.SetDeltaTimeDerivatives( M, 1.0 );
            [ res1 ] = domainUtilsP2edges.YDeriv( M, [1 -2 1], 0 );
            expectedxRes1 = [  -0.5677         0    1.0543;...
                                2.1978         0    0.5495;...
                               -0.5677    1.5000    2.6357];  
                           
            [ res2 ] = domainUtilsP2edges.XDeriv( M, [1 -2 1], 0 );
            expectedxRes2 = [ 0           0           0
                              -2.0000   -2.0000   -0.5000
                              0         -1.3650       0];  

                          
            result = min(min( abs( res1 - expectedxRes1 ) < 10^(-3) ) );
            result = result + min(min( abs( res2 - expectedxRes2 ) < 10^(-3) ) );
            if( result == 2 )
                disp( 'SimpleBEDomainTestYder1 is OK.' );
            else
                disp( 'SimpleBEDomainTestYder1 Failed!' );
            end
        end 
        
        function [ result ] = SimpleBEDomainTestYder3(st,en)
            
            [ domainUtilsP2Edges ] = TestBEDomainUtilsP2edges.GetDomainUtilsP2edgesInst(st,en,4);        
            
            M = [  1    1    1    2    1;...
                   2    2    2    2    2;...
                   1    1    1    1    1;...
                   2    2    2    1    1;...
                   2    2    2    1    1];   
            domainUtilsP2Edges = domainUtilsP2Edges.SetDeltaTimeDerivatives( M, 1 );
            [ res ] = domainUtilsP2Edges.YDeriv( M, [-1/12    4.0/3   -2.5000    4.0/3   -1/12], 0 );
            expectedxRes = [   -0.3049    0.2028    4/3   -0.0560    0.7010;...
                               -0.0005   -1.1354   -0.0000    2.1086   -0.2493;...
                               -2.8388   -2.8388   -0.0000   -1.3736   -1.4652;...
                               -0.0005   -1.1354   -1.25    1.0543   -0.1246;...
                               -0.6099    0.4056   -1.25   -0.0280    0.7010];                         
            result = min(min( abs( res - expectedxRes ) < 10^(-3) ) );
            
            if( result == 1 )
                disp( 'SimpleBEDomainTestYderiv1 is OK.' );
            else
                disp( 'SimpleBEDomainTestYderiv1 Failed!' );
            end
        end 
        
        function [ result ] = SimpleBEDomainTestYder4(st,en)
            
            [ domainUtilsP2Edges ] = TestBEdomainUtilsP2Edges.GetdomainUtilsP2EdgesInst(st,en,2);        
            
            M = [  1    1    1    2    1;...
                   2    2    2    2    2;...
                   1    1    1    1    1;...
                   2    2    2    1    1;...
                   2    2    2    1    1];   
               
            domainUtilsP2Edges = domainUtilsP2Edges.SetDeltaTimeDerivatives( M, 1.0 );
            [ res ] = domainUtilsP2Edges.YDeriv( M, [1 -2 1], 0 );
            
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
            
            [ domainUtilsP2Edges ] = TestBEdomainUtilsP2Edges.GetdomainUtilsP2EdgesInst(st,en,4);        
            
            M = [  1    1    1    2    1;...
                   2    2    2    2    2;...
                   1    1    1    1    1;...
                   2    2    2    1    1;...
                   2    2    2    1    1];   
               
            domainUtilsP2Edges = domainUtilsP2Edges.SetDeltaTimeDerivatives( M, 1.0 );
            [ res ] = domainUtilsP2Edges.YDeriv( M, [-1/12    4.0/3   -2.5000    4.0/3   -1/12], 0 );
            
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
        
        function [ result ] = SimpleBEDomainTestYderiv2(st,en)
            
            [ domainUtilsP2Edges ] = TestBEDomainUtilsP2edges.GetDomainUtilsP2edgesInst(st,en,2);        
            
            M = [  1    1    1    2    1;...
                   2    2    2    2    2;...
                   1    1    1    1    1;...
                   2    2    2    1    1;...
                   2    2    2    1    1];   
            domainUtilsP2Edges = domainUtilsP2Edges.SetDeltaTimeDerivatives( M, 1 );
            [ res ] = domainUtilsP2Edges.YDeriv( M, [1 -2 1], 0 );
            expectedxRes = [   -0.3049         0    1.0000   -2.0000    0.7010;...
                                2.1978         0         0         0    2.1978;...
                               -2.1978         0         0         0   -1.0989;...
                                1.0989         0   -1.0000    1.0000         0;...
                               -0.6099         0   -1.0000    1.0000    0.7010];                         
            result = min(min( abs( res - expectedxRes ) < 10^(-3) ) );
            
            if( result == 1 )
                disp( 'SimpleBEDomainTestYder2 is OK.' );
            else
                disp( 'SimpleBEDomainTestYder2 Failed!' );
            end
        end 
        
        
    end
    
    methods(Static, Access = private)       

        function [ domainUtilsP2edges ] = GetDomainUtilsP2edgesInst(st,en,ord)
            x = st:en;
            y = st:en;
            domainUtilsP2edges = BEDomainUtilsP2Edges( x, y, ord, 3, 0.3, 0.1, 0.2 );
        end 
    end
end