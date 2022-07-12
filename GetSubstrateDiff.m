function out=GetSubstrateDiff( Substrate,nb_c,y_select,x_select)
e_area=0; 
if nb_c==0
    e_area = e_area-Substrate(y_select,x_select);
end
if nb_c>0
    e_area = e_area+Substrate(y_select,x_select);
end
out = e_area;
end