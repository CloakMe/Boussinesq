classdef (ConstructOnLoad) BEDiscretizationParameters
   % Class help goes here
    properties (SetAccess = private, GetAccess = public)
        x
        y
        h
        tau
        tEnd
        estep
    end 
     
    properties (SetAccess = public, GetAccess = public)
        order
    end 
 
    methods

        function this = BEDiscretizationParameters( x, y ,h, order ,tau ,tEnd, estep )
            % Method help here
            this.x = x;
            this.y = y;
            this.h = h;
            this.order = order;
            this.tau = tau;
            this.tEnd = tEnd;
            this.estep = estep;
        end
    end

end