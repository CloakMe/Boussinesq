function stopSwitch = Stop(figHandle,stopSwitch)
   pause(0.001);
   is_x = get(figHandle,'currentcharacter');
   if (is_x =='x')
       if(stopSwitch==1)
           stopSwitch=2;
       else
          stopSwitch=1;
       end
   end
end