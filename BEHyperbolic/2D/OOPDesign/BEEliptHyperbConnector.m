%return
clear; clc;

bndCutSize = 0;
% ChristovIC_40_80_bt1_c090_h010_O(h^6)
% partialPath = 'BEEliptic\Boussinesq2D\SavedWorkspaces\';
partialPath = 'BEEliptic\Boussinesq2D\ZeroBoundary\ChristovIC_40_bt1_c090\Oh4\';
waveFactory = WaveFactory( partialPath, 'ChristovIC_40_ZB1_bt1_c090_h010_O(h^4)', bndCutSize, 0 ); %
%waveFactory = WaveFactory( 'BestFitIC' );

tEnd=10.0;
SavingTheSolution = 1;
fprintf('SavingTheSolution = %.1d\n', SavingTheSolution);

tau = waveFactory.h/10;
if( waveFactory.beta == 3 && waveFactory.order == 2 )
    tau = waveFactory.h/40;
elseif( waveFactory.beta == 1 && waveFactory.order == 2 )
    tau = waveFactory.h/200;
end
fprintf('tau = %.6f, h = %.6f, tau/h = %.6f\n', tau, waveFactory.h, waveFactory.h/tau);

%turnOnCaxis = 0;
%waveFactory.PlotSingleWave( turnOnCaxis );

estep = max(floor((1/tau)/10),1); %zapazwat se 20 stypki za edinitsa vreme
dscrtParams = BEDiscretizationParameters( waveFactory.x, waveFactory.y ,waveFactory.h, waveFactory.order,...
                                         tau, tEnd, estep );
eqParams = BEEquationParameters( waveFactory.alpha, waveFactory.beta1, waveFactory.beta2, waveFactory.c );
ic = BEInitialCondition( waveFactory.u_t0 , waveFactory.dudt_t0, waveFactory.mu, waveFactory.theta );   
%engine = BEEngineEnergySaveZeroBnd( dscrtParams, eqParams, ic ); %BEEngineTaylorSoftBnd %BEEngineEnergySaveSoftBnd
engine = BEEngineTaylor( dscrtParams, eqParams, ic );
% _____________________________________
tic

%(VS) vector scheme: -->  O(tau + h^2)
%(VS) vector scheme: -->  O(tau + h^2)
%[tt, max_vv, t ,v1l, v2l]  = BE2D_v6(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0); ver = 6;
%(VC) Explicit method with variable change applied -->  O(tau^2 + h^2)  tau<function(h,beta)<h ..
%[tt, max_v, t, vl]  = BE2D_v4(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0);  ver = 4;
%(NVC) Explicit method NO variable change -->  O(tau^2 + h^2) 
%[tt, max_v, t, vl]  = BE2D_v3(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0);  ver = 3;
%Taylor method variable change applied --> O(tau^4 + h^2)  tau<function(h,beta)<h ..
%[tt, max_v, t, vl]  = BE2D_t1(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0);  ver = 1;
% Taylor method variable change applied --> O(tau^4 + h^4)  tau<function(h,beta)<h ..
[engine, tt, max_v, t, EN, II, vl, dvl] = engine.BESolver( );
% Taylor method variable change applied --> O(tau^4/tau^4 + h^8)  tau<function(h,beta)<h ..
%[tt, max_v, t, EN, II, vl, dvl]  =  BE2D_t8(x,y,h,tau,t_end,beta1,beta2,al,c,c1,ord,0,estep,u_t0,dudt_t0);  ver = 2;
% Taylor method variable change applied --> O(tau^4 + h^4)  tau<function(h,beta)<h ..
%[tt, max_v, t, vl, dvl]  =  CH2D_t2(x,y,h,tau,t_end,beta1,beta2,al,c^2,ord,estep,u_t0,dudt_t0);  ver = 222;
%(NVC sit) Explicit method NO variable change -->  O(tau^2 + h^2) 
%[tt, max_v, t, vl]  = BE2D_v3_sit(x,y,h,tau,t_end,beta1,beta2,al,estep,u_t0,dudt_t0);  ver = 33;
%(VC sit) Explicit method with variable change applied -->  O(tau^2 + h

 %{
  xInd(1) = 1;
  sx = length(engine.x)
  sy = length(engine.y)
  xInd(2) = sx;
  yInd(1) = 1;
  yInd(2) = sy;
  %}
elapsedTime = toc;
fprintf('Elapsed Time = %.2f min \n', elapsedTime/60);

curSz = size( vl );

x = engine.x;
y = engine.y;

bt = waveFactory.beta1/waveFactory.beta2;
%DrawEnergyForHyperbolicBE( engine, tt );

useZeroBoundary = 0;
if( waveFactory.mu == 0 )
    useZeroBoundary = 1;
end
fprintf('Use ZeroBoundary = %d\n', useZeroBoundary);
if( SavingTheSolution == 1)
	try
        save (['SavedWorkspaces\' engine.GetName() '_' num2str(floor(waveFactory.compBox.x_end2)) ...
            '_ZB' num2str(useZeroBoundary) '_bt' ...
            num2str(bt) '_c0' num2str(floor(waveFactory.c*100)) ...
          '_h0' num2str(waveFactory.h*100,'%.02d') '_O(h^' num2str(  waveFactory.order  ) ')']);

    catch ex
        fprintf('Trying to save only the most important solution parameters, ...\n');
        save (['SavedWorkspaces\' engine.GetName() '_' num2str(floor(waveFactory.compBox.x_end2)) ...
             '_ZB'  num2str(useZeroBoundary)  '_bt' ...
             num2str(bt) '_c0' num2str(floor(waveFactory.c*100)) '_h0' ...
             num2str(waveFactory.h*100,'%.02d') '_O(h^' num2str(  waveFactory.order  ) ')'], ...
             'tau', 'x', 'y', 'tt', 'max_v', 't', 'EN', 'II', 'vl', 'dvl');
        fprintf('parameters saved!\n');
    end
end

topView = 1;
viewTypeX = 0;
viewTypeY = 90;
%MovieForBEHyperbolic( viewTypeX, viewTypeY, tt, x, y );

Q = 41;
if (false)
    figure(11)
    mesh(x, y(1:Q), vl(:,1:Q)');
    view( viewTypeX, viewTypeY );
    title('Bottom domain boundary near y=-40');
    xlabel('x');            ylabel('y');
    figure(12)
    mesh(x(1:Q), y, vl(1:Q,:)');
    view( viewTypeX, viewTypeY );    
    title('Left domain boundary near x=-40');
    xlabel('x');            ylabel('y');

    figure(13)
    mesh(x, y(end-Q:end), vl(:,end-Q:end)');
    view( viewTypeX, viewTypeY );    
    title('Top domain boundary near y=40');
    xlabel('x');            ylabel('y');    
    figure(14)
    mesh(x(end-Q:end), y, vl(end-Q:end,:)');
    view( viewTypeX, viewTypeY );
    title('Right domain boundary near x=40');
    xlabel('x');            ylabel('y');
end

figure(15)
mesh(x,y,vl')
view( viewTypeX, viewTypeY );
title('solution');
xlabel('x');            ylabel('y');
colorbar;
figure(16)
%hold on;
plot(t(1:end-1),max_v,'k', t(1), max_v+0.005, 'k', t(1), max_v-0.005, 'k' );
%hold off;
title('Evolution of the maximum');
xlabel('time "t"');  ylabel('max(u_h)');
figure(17)
%hold on;
plot(t(1:end-1),EN,'k',t(1),EN(2)+EN(2)/1000.0,t(end),EN(end)-EN(end)/1000.0 )
%hold off;
title('Energy functional');
xlabel('time "t"');  ylabel('EN');

figure(18)
%hold on;
[indeces, shift] = BEUtilities.GetCommonIndexArray( t, II );
indeces(1) = [];
plot(t(indeces+shift),II(indeces),'k',t(1),II(2)+II(2)/1000.0,t(end),II(end)-II(end)/1000.0 )
%hold off;
title('Integral');
xlabel('time "t"');  ylabel('Integral');

return;
%----------------------------------------------------------------------------------------
CreatePictures( viewTypeX, viewTypeY, tt, x, y ); 

waveFactory = WaveFactory( 'ChristovIC_80_bt1_c090_h020_O(h^4)', 5 );

    tau = 0.1;
    tEnd=44.0;
    tt = tau:tau:tEnd;
        
    compBoxToShow = waveFactory.compBox;
    compBoxToShow.x_st = waveFactory.x( 1 );
    compBoxToShow.x_end = waveFactory.x( end );
    compBoxToShow.y_st = waveFactory.y( 1 );
    compBoxToShow.y_end = waveFactory.y( end );

    topView = 1;
    viewTypeX = 81;
    viewTypeY = 18;
    MovieForBEHyperbolic( compBoxToShow, viewTypeX, viewTypeY, tt, waveFactory.x, waveFactory.y );

    
    flipU = fliplr(this.u_t0);
    flipU = this.u_t0 + flipU;
    this.dudt_t0 = this.dudt_t0 + fliplr(this.dudt_t0);
    this.mu = 2*c1;

    figure(14)
    mesh(x,y,flipU')
    title('solution');
    xlabel('x');            ylabel('y');

    figure(15)
    mesh(x,y,this.u_t0')
    title('solution');
    xlabel('x');            ylabel('y');

   
    figure(17)
%hold on;
plot(tt,EN,'k',tt(1),EN(2)+EN(2)/1000.0,tt(end),EN(end)-EN(end)/1000.0 )
%hold off;
title('Energy functional');
xlabel('time "t"');  ylabel('EN');