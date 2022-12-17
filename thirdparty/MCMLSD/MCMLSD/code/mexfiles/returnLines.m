function [begin_points,end_points,lslength,startpositions,endpositions]=returnLines(points_proj,final_labels)
count = 0;
on_seg = false;
begin_pointsT = zeros(20,2);
counter_begin_points=1;
end_pointsT = zeros(20,2);
counter_end_points=1;
lslengthT = zeros(20,1);
counter_lslength=1;
startpos = 0;
startpositionsT = zeros(20,1);
counter_startpositions=1;
endpositionsT = zeros(20,1);
counter_endpositions=1;

for i = 1:length(final_labels);
   if ~on_seg
       if final_labels(i) == 1
           begin_pointsT(counter_begin_points,:)=points_proj(i,:);
           counter_begin_points=counter_begin_points+1;
           startpositionsT(counter_startpositions) = i;
           counter_startpositions=counter_startpositions+1;
           startpos = i;
           on_seg = true;
           count = count + 1;
       end
   elseif i == length(final_labels)
       end_pointsT(counter_end_points,:) =points_proj(i,:);
       counter_end_points=counter_end_points+1;
       endpositionsT(counter_endpositions) = i;
       counter_endpositions=counter_endpositions+1;
       lslengthT(counter_lslength) = i-startpos;
       counter_lslength=counter_lslength+1;
       on_seg = false;
   elseif final_labels(i) == 0
       end_pointsT(counter_end_points,:) = points_proj(i-1,:);
       counter_end_points=counter_end_points+1;
       endpositionsT(counter_endpositions) = i-1;
       counter_endpositions=counter_endpositions+1;
       lslengthT(counter_lslength) = i-startpos;
       counter_lslength=counter_lslength+1;
       on_seg = false;
   end
end

begin_points = begin_pointsT(1:counter_begin_points-1,:);
end_points =end_pointsT(1:counter_end_points-1,:);
lslength =lslengthT(1:counter_lslength-1); 
startpositions =startpositionsT(1:counter_startpositions);
endpositions =endpositionsT(1:counter_endpositions-1);

end
