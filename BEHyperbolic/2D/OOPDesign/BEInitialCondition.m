classdef (ConstructOnLoad) BEInitialCondition
   % Class help goes here
   properties (SetAccess = private, GetAccess = public)
      u_t0
      dudt_t0
      mu
      theta
   end 
 
   methods      
      function this = BEInitialCondition( u_t0 ,dudt_t0, mu, theta )
         % Method help here
            this.u_t0 = u_t0;
            this.dudt_t0 = dudt_t0;
            this.mu = mu;
            this.theta = theta;
      end
   end

end