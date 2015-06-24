classdef TestBEDomainUtilsP2   
    methods(Static)          
        function [ result ] = SimpleBEDomainTest(st,en)
            [ domainUtilsP2 ] = TestBEDomainUtilsP2.GetDomainUtilsP2inst(st,en,2);        
            [ res ] = domainUtilsP2.TestFun( domainUtilsP2.X, domainUtilsP2.Y );
            expectedxRes = [  2     1     2;...
                              1     0     1;...
                              2     1     2];                        
            result = min(min(res == expectedxRes) );
            
            [ res ] = domainUtilsP2.GetTestFunMid();
            result = result + min(min(res == expectedxRes) );
            if( result == 2 )
                disp( 'SimpleBEDomainTest is OK.' );
            else
                disp( 'SimpleBEDomainTest Failed!' );
            end
        end 
        
        function [ result ] = DerivBEDomainTest(st,en)
            [ domainUtilsP2 ] = TestBEDomainUtilsP2.GetDomainUtilsP2inst(st,en,2);            
            [ res ] = domainUtilsP2.TestFun( domainUtilsP2.X, domainUtilsP2.Y );

            res = domainUtilsP2.YDerivative( res,...
                                       domainUtilsP2.GetTestFunLeft(),...
                                       domainUtilsP2.GetTestFunRight(),...
                                       [ 1 -2 1] );         
            expectedxRes = [  2     2     2;...
                              2     2     2;...
                              2     2     2];
            result = min(min(res == expectedxRes) );
            
            if( result == 1 )
                disp( 'DerivBEDomainTest is OK.' );
            else
                disp( 'DerivBEDomainTest Failed!' );
            end
        end 
        
        function [ result ] = BEBndDerFunTest(st,en)
            
            [ domainUtilsP2 ] = TestBEDomainUtilsP2.GetDomainUtilsP2inst(st,en,2);        
            t = 0;
            [ res ] = domainUtilsP2.DersBnd( domainUtilsP2.X, domainUtilsP2.Y, t );
                                      
            t = -.1;
            [ res_1 ] = domainUtilsP2.DersBnd( domainUtilsP2.X, domainUtilsP2.Y, t );
            t = .1;
            [ res1 ] = domainUtilsP2.DersBnd( domainUtilsP2.X, domainUtilsP2.Y, t );
            secTDer = (1/.1^2)*( res1(:,:,1) + res_1(:,:,1) - 2*res(:,:,1) );
            
            Diff = secTDer - res(:,:,3);
                      
            if( max( max (Diff ) ) < .00000001)
                disp( 'BEBndFunTest is OK.' );
            else
                disp( 'BEBndFunTest Failed!' );
            end
        end
        
        function [ result ] = BEBndDerFunTest2(st,en)
            
            [ domainUtilsP2 ] = TestBEDomainUtilsP2.GetDomainUtilsP2inst(st,en,4);   
            
            fd = BEUtilities.GetFinDiffCoeff([ -2,-1,0,1,2], 4);
            
            t = 0;
            [ res ] = domainUtilsP2.GetDersBndMid( t );
                                      
            t = -.2;
            [ res_2 ] = domainUtilsP2.GetDersBndMid( t );
            t = -.1;
            [ res_1 ] = domainUtilsP2.GetDersBndMid( t );
            t = .1;
            [ res1 ] = domainUtilsP2.GetDersBndMid( t );
            t = .2;
            [ res2 ] = domainUtilsP2.GetDersBndMid( t );
            fourthTDer = (1/.1^4)*( fd(1)*res2(:,:,1) + fd(2)*res1(:,:,1) +...
                fd(3)*res(:,:,1) + fd(4)*res_1(:,:,1) + fd(5)*res_2(:,:,1) );
            
            Diff = fourthTDer - res(:,:,5);
                      
            if( max( max (Diff ) ) < .00000001)
                disp( 'BEBndDerFunTest2 is OK.' );
            else
                disp( 'BEBndDerFunTest2 Failed!' );
            end
        end
        
    end
    
    methods(Static, Access = private)          
        function [ domainUtilsP2 ] = GetDomainUtilsP2inst(st,en,ord)
            x = st:en;
            y = st:en;
            domainUtilsP2 = BEDomainUtilsP2( x, y, ord, 3, 0.3, 0.1 );
        end 
    end
    
    %TestBEDomainUtilsP2.SimpleBEDomainTest()
end

