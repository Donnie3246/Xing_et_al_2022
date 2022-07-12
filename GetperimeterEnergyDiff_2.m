function out=GetperimeterEnergyDiff_2(cells,sigma, x_select,y_select,c, nb_c)
global LAM_PERIMETER TARGET_PERIMETER_SIZE
e_perimeter=0;                
if c>0
    Cell_area_c=sigma==c;
    perimeter_diff_c=perimeter_local_1(Cell_area_c,y_select,x_select);
    e_perimeter = e_perimeter + LAM_PERIMETER*(perimeter_diff_c)*(perimeter_diff_c+2* cells.perimeter(c)-2*cells.target_perimeter);
%     Cell_area_2=Cell_area_c; 
%     Cell_area_2(y_select,x_select)=double(~Cell_area_c(y_select,x_select));

%      if  bweuler(Cell_area_2)<1   %||A.NumObjects>1   %连通性计算，低于1才能
%          e_perimeter=e_perimeter+200000;
%      end    
end

if nb_c>0
    Cell_area_nbc=sigma==nb_c;
    perimeter_diff_nbc=perimeter_local_2(Cell_area_nbc,y_select,x_select);
    e_perimeter = e_perimeter + LAM_PERIMETER*(perimeter_diff_nbc)*(perimeter_diff_nbc+2*cells.perimeter(nb_c)-2*cells.target_perimeter);
%     Cell_area_2=Cell_area_nbc; 
%     Cell_area_2(y_select,x_select)=double(~Cell_area_nbc(y_select,x_select));

%     if  bweuler(Cell_area_2)<1%||A.NumObjects>1
%     e_perimeter=e_perimeter+200000;
%     end
end

out = e_perimeter;
end 

function perimeter_diff=perimeter_local_1(Cell_area,y_select,x_select)
perimeter_diff=sum(sum(2*Cell_area(y_select-1:y_select+1,x_select-1:x_select+1)-1))-1;
end
function perimeter_diff=perimeter_local_2(Cell_area,y_select,x_select)
perimeter_diff=-sum(sum(2*Cell_area(y_select-1:y_select+1,x_select-1:x_select+1)-1))-1;
end