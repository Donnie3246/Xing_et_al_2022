function out=GetAct_PodosomeEnergyDiff_real2(Substrate,Cell_area, y_select,x_select,c, nb_c)
%更新:死皮会直接阻碍8联通区域的运动   '08-May-2021'
global LAM_ACT MAX_ACT  PEG_obstacle FN_obstacle 
e_act=0;
    lam_shrink=1;
    lam_stretch=1; 
    L_neighood=2;

if  Substrate(y_select,x_select)==0    
%     lam_shrink=1;
    energy_obstacle=PEG_obstacle;  
else
    energy_obstacle=FN_obstacle;
end

if c>0  %收缩
        A1=Cell_area{c,2}(y_select-1:y_select+1,x_select-1:x_select+1);
        A2=A1(A1>=0);GM_act_c=prod(A2)^(1/length(A2));
        if  Substrate(y_select,x_select)==0 

            lam_shrink=1;
            if max(max(Cell_area{c,3}(y_select-L_neighood:y_select+L_neighood,x_select-L_neighood:x_select+L_neighood)))>=1
                lam_shrink=1;
                energy_obstacle=0;
            end

            lam_stretch=lam_stretch*1;
        end
%         lam_shrink
        e_act = e_act +lam_shrink*LAM_ACT/MAX_ACT*GM_act_c-energy_obstacle;
end
if nb_c>0   %扩张
        A1=Cell_area{nb_c,2}(y_select-1:y_select+1,x_select-1:x_select+1);
        A2=A1(A1>=0);GM_act_nbc=prod(A2)^(1/length(A2));
   if  Substrate(y_select,x_select)==0    
	  
        lam_stretch=1;
        if max(max(Cell_area{nb_c,3}(y_select-L_neighood:y_select+L_neighood,x_select-L_neighood:x_select+L_neighood)))>=1
            lam_stretch=1;
            energy_obstacle=0;
        end
        lam_stretch=lam_stretch*1;
    
   end

        e_act = e_act - lam_stretch*LAM_ACT/MAX_ACT*GM_act_nbc+energy_obstacle;
end   

      out = e_act;
end