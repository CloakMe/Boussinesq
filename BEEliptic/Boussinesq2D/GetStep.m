function step = GetStep(h)

   step = 1;

   if(0.5>h && h>=0.25) 
       step = 2; 
   end
   if(0.25>h && h>=0.125) 
       step = 4; 
   end
   if(0.125>h && h>=0.0625) 
       step = 8; 
   end
   if(0.0625>h ) 
       step = 16; 
   end
end