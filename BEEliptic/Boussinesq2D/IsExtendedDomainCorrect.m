function correct = IsExtendedDomainCorrect(compBox)

    if(compBox.x_st>=compBox.x_st2 && compBox.x_end <= compBox.x_end2 &&...
            compBox.y_st>=compBox.y_st2 && compBox.y_end <= compBox.y_end2)
        
         if(compBox.x_st==compBox.x_st2 && compBox.x_end == compBox.x_end2 &&...
                 compBox.y_st==compBox.y_st2 && compBox.y_end == compBox.y_end2)
             
             correct = 0;
         else
            correct = 1;
         end
    else
        correct = 0;
    end
end