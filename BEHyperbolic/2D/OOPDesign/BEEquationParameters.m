classdef (ConstructOnLoad) BEEquationParameters
   % Class help goes here
   properties (SetAccess = private, GetAccess = public)
      alpha
      beta1
      beta2
      beta
      c
   end 
 
   methods
      
      function this = BEEquationParameters( alpha, beta1, beta2, c )
         % Method help here
            this.alpha = alpha;
            this.beta = beta1/beta2;
            this.c = c;
            this.beta1 = beta1;
            this.beta2 = beta2;
      end
   end

end