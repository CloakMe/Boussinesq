clear; clc;
%return

bndCutSize = 0; 
% ChristovIC_80_bt1_c090_h020_O(h^6)
% ChristovIC_40_bt1_c090_h020_O(h^4)
% ChristovIC_80_bt3_c052_h020_O(h^6)
% ChristovIC_36_bt5_c040_h020_O(h^6)
% ChristovIC_40_bt069_c090_h020_O(h^6) 
% ChristovIC_40_bt05_c085_h020_O(h^6)
waveFactory = WaveFactory( 'ChristovIC_40_bt3_c045_h010_O(h^6)', bndCutSize );
%waveFactory = WaveFactory( 'BestFitIC' );

    tau = 0.05;
    tEnd=10.0;
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
  elapsedTime = toc
  curSz = size( vl );

  x = engine.x;
  y = engine.y;
      
  topView = 1;
  viewTypeX = 90;
  viewTypeY = 90;
  MovieForBEHyperbolic( viewTypeX, viewTypeY, tt, x, y );
       
   Q = 41;
    figure(11)
    mesh(x, y(1:Q), vl(:,1:Q)');
    title('Bottom domain boundary');
    xlabel('x');            ylabel('y (y_{start})');
    figure(12)
    mesh(x(1:Q), y, vl(1:Q,:)');
    title('Left domain boundary');
    xlabel('x (x_{start})');            ylabel('y');

    figure(13)
    mesh(x, y(end-Q:end), vl(:,end-Q:end)');
    title('Top domain boundary');
    xlabel('x');            ylabel('y (y_{end})');
    title('Right domain boundary');
    figure(14)
    mesh(x(end-Q:end), y, vl(end-Q:end,:)');
    xlabel('x (x_{end})');            ylabel('y');
    
    figure(15)
    mesh(x,y,vl')
    title('solution');
    xlabel('x');            ylabel('y');
    figure(16)
    %hold on;
    plot(t(1:end-1),max_v,'k')
    %hold off;
    title('Evolution of the maximum');
    xlabel('time "t"');  ylabel('max(v)');
    figure(17)
    %hold on;
    plot(tt,EN,'k',tt(1),EN(1)+EN(1)/1000.0,tt(end),EN(end)-EN(end)/1000.0 )
    %hold off;
    title('Energy functional');
    xlabel('time "t"');  ylabel('EN');
        
    bt = waveFactory.beta1/waveFactory.beta2;
    %DrawEnergyForHyperbolicBE( engine, tt );
    save (['SavedWorkspaces\Hyperb_' num2str(floor(waveFactory.compBox.x_end2)) '_bt' ...
        num2str(bt) '_c0' num2str(floor(waveFactory.c*100)) ...
      '_h0' num2str(waveFactory.h*100) '_O(h^' num2str(  waveFactory.order  ) ')']);
    return;
%----------------------------------------------------------------------------------------
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

    

   