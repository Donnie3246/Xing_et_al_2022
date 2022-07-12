close all
clear;tic
%修改周长计算方法，与2015文献相同。方法类似黏附能计算
%将两种绘图方式组合至一起
%添加运动轨迹，速度矢量。%2021.05.28
%添加序参数。%2021.06.01
%添加边框 %2021.06.05
%修改初始种子及边界条件 %2021.06.09
%跨越式迁移，纯actin %2021.-0.9
%跨越式迁移新版，锚定效应 %2022.06.10
%边界间距为20pixel的条带
%引入伪足小体效应
%改为影响细胞Act值 %2022.06.13
%条带 %2022.06.14
%更新优化，添加颜色 %2022.06.20
%% Global parameters
global XMAX YMAX   STEP_MAX N_cell ; 
global J_LL J_LM LAM_AREA TARGET_CELL_SIZE LAM_PERIMETER    ...
       TARGET_PERIMETER_SIZE LAM_ACT MAX_ACT PEG_obstacle FN_obstacle flag_podosome;
%% Parameter Setting

% 绘图模式 % 
Pic_Type=2; %1：彩色分块图； 2：Act生长图
intervel_step=4; %绘图间隔(MCS)
N_step=1600;  %模拟总步长(MCS)
real_XMAX=300;  %图片宽度（pixel）
real_YMAX=160;  %图片高度（pixel）
  
%细胞基础参数 %   
LAM_AREA=2;                       %面积模量
LAM_PERIMETER=0.4;                %周长模量
TARGET_CELL_SIZE=320;             %细胞面积(/pixel)
TARGET_PERIMETER_SIZE= 250;       %细胞周长(8连通域下的周长)
J_LL=20;       %细胞间黏附能（无用）
J_LM=20;       %细胞-ECM黏附能
LAM_ACT=100;   %细胞骨架基础能量
MAX_ACT=25;    %FN上细胞骨架最大活性

%伪足小体参数 %  
PEG_obstacle=20.0;  %PEG能量阻碍（正值表示阻碍）
FN_obstacle=-30;    %FN能量吸附（负值表示吸附）
persist_time=5;     %伪足小体寿命
PEG_act=0;          %PEG处新生骨架的活性
PODO_act=45;        %PEG处有伪足小体激发的细胞骨架的活性，须大于MAX_ACT
P_Grow0=0.05*16;       %伪足小体生成概率
   threshold_podo=0.75; %成伪足小体的act活性阈值（threshold_podo*MAX_ACT）
   
P_Decay0=0.8;       %伪足小体降解（寿命减1）的概率
P_Grow=P_Grow0;

% 跨越辅助参数 %
%跨越前
FN_cell_area0=20;   %跨越时留在FN区域的细胞面积（锚定面积）
e_back0=200;        %减少锚定面积的耗能
%跨越后
FN_cell_area1=30;       %已经跨越到新FN的细胞面积
decrease_rate=1/4;      %跨越后负反馈导致伪足小体生成概率的降低比例
T_depolymerize0=100;    %跨越后负反馈导致伪足小体生成概率降低的持续时间(MCS)

% 边界条件限制 %
Distance=36;  %条纹间距（pixel）
width_bar=19;  %条纹宽度（pixel）
N_bar=3;      %条纹数量（pixel）
Distance=Distance+width_bar;

% 上色 % 
color_FN=[1,0.9,1];   %FN颜色
color_podo=[1,0,0];   %PEG颜色
cool0=cool(PODO_act*2);   
color_act=cool0(1:PODO_act,:);      %act颜色
my_map=[1 1 1;color_act;color_podo;color_FN];%flipud

%% 无需修改参数
merge=1;   %不可接触边界
XMAX=real_XMAX+merge*2; % domain size in x axis 
YMAX=real_YMAX+merge*2; % domain size in y axis %注意最后的圆直径减了3 
center_xc=(XMAX+1)/2; center_yc=(YMAX+1)/2; 
N_cell=1;            %细胞数,不能修改
STEP_MAX=real_XMAX*real_YMAX*N_step; % time step number
TEMPERATURE=20;
flag_podosome=1;
originL_FN=2;  %初始基底的编号

Name=['rotate',num2str(N_cell),'_',num2str(real_XMAX),'×',num2str(real_YMAX),'_LA',num2str(LAM_AREA),'_LP',num2str(LAM_PERIMETER),'_LAt',num2str(LAM_ACT),'_MA',num2str(MAX_ACT),...
        '_TA',num2str(TARGET_CELL_SIZE), '_TP',num2str(TARGET_PERIMETER_SIZE),'_JC',num2str(J_LL),'_JM',num2str(J_LM),'.tif'];
%% Boundary Condition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y]=meshgrid(1:XMAX,1:YMAX);
sigma=zeros(YMAX,XMAX);
sigma_order=find(((X-center_xc).^2+(Y-center_yc).^2)<(width_bar/2-1)^2);
order=randi(length(sigma_order),[N_cell,1]);
sigma(sigma_order(order))=1:N_cell;

Substrate_label=zeros(YMAX,XMAX);
sigma_orderS={};
for i=1:N_bar
         sigma_orderS{i}=find(((Y-center_yc-(i-2)*Distance).^2)<width_bar^2/4);
         Substrate_label(sigma_orderS{i})=i;
end

sigma_order_out=[];
for i=1:N_bar
        sigma_order_out=union(sigma_order_out,sigma_orderS{i});   
end
Substrate=zeros(YMAX,XMAX);
Substrate(sigma_order_out)=1;
% imshow(Substrate_label/7)
%% Initial Condition
%sigma = GetInitialCondition();   %初始化细胞及ECM分布,赋给sigma
Num_cell=max(sigma(:));
cells = SetInitialParameters(sigma);
Cell_area{Num_cell,2}=[];
im{N_step,1}=[];
XY{Num_cell,1}=[];
k_T=1;v_intervel=10;
V{Num_cell,1}=[];
order_parameter=zeros(N_step*intervel_step-v_intervel,1);
V_mean=zeros(N_step*intervel_step-v_intervel,1);
Omage=zeros(N_step*intervel_step-v_intervel,1);

%%
for i=1:Num_cell
Cell_area{i,1}=sigma==i;
Cell_area{i,2}=Cell_area{i,1}-1;
Cell_area{i,2}(Cell_area{i,2}<0)=nan;
Cell_area{i,3}=(1-Substrate).*Cell_area{i,1}*nan;
cells.perimeter(i)=0;
 for n_y=2:size(Cell_area{i,1},1)-1
    for n_x=2:size(Cell_area{i,1},2)-1
        if Cell_area{i,1}(n_y,n_x)==1
           cells.perimeter(i)=cells.perimeter(i)+8-sum(sum(Cell_area{i,1}(n_y-1:n_y+1,n_x-1:n_x+1 )))+1;
        end
    end
 end
   [yi,xi]=find(Cell_area{i,1}>0);
   XY{1,1}(i,1)=mean(xi);
   XY{1,1}(i,2)=mean(yi);
end

%% Setting for Model Variables and Parameters


%% Time Step
figure(1)
tic 
cells.target_perimeter=TARGET_PERIMETER_SIZE;
sigma_Act=zeros(size(sigma)); sigma_Act(sigma==0)=nan;
Pseudopodia =zeros(size(sigma))*nan;
FN_cell = Substrate.*sigma;
FN_cell_area=sum(FN_cell(:));
FN_cell_area_other=0;
otherL_FN=0;
flag_podosome1=1;
T_depolymerize=0;

%计算细胞连通域
[yy,xx]=find(sigma==1);
Cell_box=[min(xx),max(xx),min(yy),max(yy)];

sigma_C=sigma;  
[L,num] = bwlabel(sigma(Cell_box(3):Cell_box(4),Cell_box(1):Cell_box(2)),8);
%%
figure(1)

for step = 1:STEP_MAX
      x_select = randi([3+merge,real_XMAX-2+merge]);
      y_select = randi([3+merge,real_YMAX-2+merge]);
      c = sigma(y_select,x_select);
 
    nb_indx=randi(8);
    switch nb_indx    
        case 1
            nb_x=-1; nb_y=-1;
        case 2  
            nb_x=-1; nb_y=0;
        case 3
            nb_x=-1; nb_y=1;
        case 4
            nb_x=0; nb_y=-1;
        case 5
            nb_x=0; nb_y=1;
        case 6
            nb_x=1; nb_y=-1;
        case 7
            nb_x=1; nb_y=0;
        case 8
            nb_x=1; nb_y=1;
    end
    nb_c = sigma(y_select+nb_y,x_select+nb_x);
    
    if FN_cell_area<FN_cell_area0
        flag_podosome=0;
    else
        flag_podosome=1;
    end
    if FN_cell_area_other>=FN_cell_area1
        flag_podosome1=0;
        if T_depolymerize==0
            T_depolymerize=T_depolymerize0; % 单位MCS
        end   
    end
    if flag_podosome1==0
        
         if mod(step,real_YMAX*real_XMAX)==0 && T_depolymerize>0
             T_depolymerize=T_depolymerize-1;
         end
         if  T_depolymerize==0
            FN_cell_area_other=0;
            originL_FN=otherL_FN;
            otherL_FN=0;
            flag_podosome1=1;
         end       
        P_Grow=P_Grow0*decrease_rate;
    else
        P_Grow=P_Grow0;
    end
    
    if FN_cell_area<1
        cells.target_area=TARGET_CELL_SIZE/2;  cells.target_perimeter=6*(cells.target_area*pi)^0.5;
            P_Grow=0;
    else
        cells.target_area=TARGET_CELL_SIZE;   cells.target_perimeter= TARGET_PERIMETER_SIZE;
%             P_Grow=P_Grow0;
    end

    if c ~= nb_c
               e_back=0;
            if flag_podosome==0&&c>0
                if Substrate(y_select,x_select)==1
                    e_back=e_back0;
                end   
            end
        %连通性能量判定 
          if c > 0  %只有收缩（删除pixel）才会导致细胞连通性变化
              sigma_C=Cell_area{c,1};    sigma_C(y_select,x_select)=0;
              [L,num] = bwlabel(sigma_C(Cell_box(3):Cell_box(4),Cell_box(1):Cell_box(2)),8);
          end
          if  num>1  %细胞连通性>1时，附加断裂能
              e_connectivity=10^10;  %断裂能很大
          else
              e_connectivity=0;
          end
            
        e_all = e_connectivity+e_back+GetAdhesionEnergyDiff(cells, sigma, x_select, y_select, c, nb_c)...
            + GetSizeEnergyDiff(cells, c, nb_c)+ GetperimeterEnergyDiff_2(cells,sigma,x_select,y_select,c, nb_c)+...
            GetAct_PodosomeEnergyDiff_real2(Substrate,Cell_area, y_select,x_select,c, nb_c);

        if e_all >= 0
            prob = exp(-e_all/TEMPERATURE);
        elseif e_all < 0
            prob = 1.0;
        end
        % Update of the state
        if prob >= rand
            if c > 0 %收缩(删除pixel)
                cells.area(c) = cells.area(c)-1;
                perimeter_diff_c=perimeter_local_1(Cell_area{c,1},y_select,x_select);
                cells.perimeter(c)=cells.perimeter(c)+perimeter_diff_c;
                Cell_area{c,1}(y_select,x_select)=0;
                Cell_area{c,2}(y_select,x_select)=nan;
                
                % 连通性计算
                if x_select==Cell_box(1)
                    if sum(Cell_area{c,1}(:,x_select)>0)==0
                        Cell_box(1)=x_select+1;
                    end
                elseif x_select==Cell_box(2)
                    if sum(Cell_area{c,1}(:,x_select)>0)==0
                        Cell_box(2)=x_select-1;
                    end
                end               
                if y_select==Cell_box(3)
                    if sum(Cell_area{c,1}(y_select,:)>0)==0
                        Cell_box(3)=y_select+1;
                    end
                elseif y_select==Cell_box(4)
                    if sum(Cell_area{c,1}(y_select,:)>0)==0
                        Cell_box(4)=y_select-1;
                    end
                end              
                
                   if Substrate(y_select,x_select)==0
                       Cell_area{c,3}(y_select,x_select)=nan; 
                       Pseudopodia(y_select,x_select)=nan;
                   else
                       FN_cell(y_select,x_select)=0;
                       FN_cell_area=FN_cell_area-1;
                       %跨越反馈
                       if Substrate_label(y_select,x_select)==otherL_FN
                           FN_cell_area_other=FN_cell_area_other-1;
                       end
                   end           
            end

            if nb_c > 0    %铺展(复制pixel)    
              cells.area(nb_c) = cells.area(nb_c)+1;
              perimeter_diff_nbc=perimeter_local_2(Cell_area{nb_c,1},y_select,x_select);
              cells.perimeter(nb_c)=cells.perimeter(nb_c)+perimeter_diff_nbc;
              Cell_area{nb_c,1}(y_select,x_select)=1;
              % 连通性计算
              Cell_box=[min(Cell_box(1),x_select),max(Cell_box(2),x_select),min(Cell_box(3),y_select),max(Cell_box(4),y_select)];
                
              if max(max(Pseudopodia(y_select-2:y_select+2,x_select-2:x_select+2)))>0
                     copy_act=PODO_act;
              elseif max(max(FN_cell(y_select-1:y_select+1,x_select-1:x_select+1)))>0
                     copy_act=max(Cell_area{nb_c,2}(y_select+nb_y,x_select+nb_x),MAX_ACT);
              else
                     copy_act=PEG_act;
              end

              Cell_area{nb_c,2}(y_select,x_select)=copy_act;
              sigma_Act(y_select,x_select)=copy_act;
                                   
                   if Substrate(y_select,x_select)==0
                       Cell_area{nb_c,3}(y_select,x_select)=0;  
                       Pseudopodia(y_select,x_select)=0;                    
                   else

                       FN_cell(y_select,x_select)=Substrate_label(y_select,x_select);
                       FN_cell_area=FN_cell_area+1;     
                       %跨越反馈
                       if Substrate_label(y_select,x_select)~=originL_FN
                           otherL_FN=Substrate_label(y_select,x_select);
                           FN_cell_area_other=FN_cell_area_other+1;
                       end
                       
                   end
            else
                  sigma_Act(y_select,x_select)=nan;
                  
            end 
                  sigma(y_select,x_select) = nb_c;           
 
        end        
    end
    

    if mod(step,real_YMAX*real_XMAX)==0
        k_T= k_T+1;
        for i=1:Num_cell
             cellact=(Cell_area{i,2}>0);
             Cell_area{i,2}(cellact)=Cell_area{i,2}(cellact)-1;
             
             %%%%%
             p_grow=P_Grow;
             cell_PEG_grow=find(Cell_area{i,3}==0);
             R_grow=(rand(size(cell_PEG_grow))-p_grow<0);
             Act_grow=Cell_area{i,2}(cell_PEG_grow)/MAX_ACT;
             Act_grow=1./(1+exp(-1000*(Act_grow-threshold_podo)));
             Act_grow=Act_grow>rand(size(Act_grow));
             Cell_area{i,3}(cell_PEG_grow)=R_grow.*Act_grow*persist_time;  
             
             p_decay=P_Decay0;
             cell_PEG_decay=find(Cell_area{i,3}>0);
             R_decay=(rand(size(cell_PEG_decay))-p_decay<0);
             Cell_area{i,3}(cell_PEG_decay)=Cell_area{i,3}(cell_PEG_decay)-1*R_decay;  
             Cell_area{i,3}( Cell_area{i,3}<0)=0;
             Pseudopodia=Cell_area{i,3};
             
             
             [yi,xi]=find(Cell_area{i,1}>0);
             XY{k_T,1}(i,1)=mean(xi);  XY{k_T,1}(i,2)=mean(yi);
             if k_T>v_intervel
                 V{k_T-v_intervel,1}(i,1)=(XY{k_T,1}(i,1)-XY{k_T-v_intervel,1}(i,1))/v_intervel; 
                 V{k_T-v_intervel,1}(i,2)=(XY{k_T,1}(i,2)-XY{k_T-v_intervel,1}(i,2))/v_intervel; 
             end   
              
        end
        if k_T>v_intervel
               Rx=XY{k_T,1}(:,1)-center_xc; Ry=XY{k_T,1}(:,2)-center_yc;
               AM=Rx.*V{k_T-v_intervel,1}(:,2)-Ry.*V{k_T-v_intervel,1}(:,1);
               V_square=V{k_T-v_intervel,1}(:,2).^2+V{k_T-v_intervel,1}(:,1).^2;
               V_mean(k_T-v_intervel)=mean(V_square.^0.5);
               R_square=Rx.^2+Ry.^2;
               Omage(k_T-v_intervel)=mean(AM./(R_square+eps));
               order_parameter(k_T-v_intervel)=mean(AM./(R_square.^0.5+eps)./(V_square.^0.5+eps));
        end
           
           sigma_Act(sigma_Act>0)=sigma_Act(sigma_Act>0)-1;
    end
  
    %%% Output Figure
    if mod(step,real_YMAX*real_XMAX*intervel_step)==0
       disp(step/(real_YMAX*real_XMAX));
       sigma_Act1=sigma_Act;sigma_Act1(isnan(sigma_Act1))=-200;
       sigma1=sigma;sigma1(sigma1==0)=-200;
       sigma_Act2=sigma_Act1;
       sigma_Act2(Pseudopodia>0)=PODO_act+2;
       sigma_Act3=sigma_Act2;
       sigma_Act3(((Substrate==1).*isnan(sigma_Act))>0)=PODO_act+4;
       
           
   if Pic_Type==1
       OutputCellFigure(cells, sigma1)
       
          if step==real_YMAX*real_XMAX*intervel_step
           jet1=[1 1 1;hsv];
           colormap(jet1)% default
           colorbar
           caxis([-2 N_cell])
           axis equal
           axis off
           pause(0.1)                   

           end 
   elseif Pic_Type==2
      
        OutputCellActFigure(sigma, sigma_Act3)
       
        if step==real_YMAX*real_XMAX*intervel_step || step==real_YMAX*real_XMAX*intervel_step*2

           colormap(my_map)
           colorbar
           caxis([-2 PODO_act+5])
           axis equal
           axis off
           pause(0.1)
        end 
   end
   
   %       h1=plot(XY{k_T,1}(:,1),XY{k_T,1}(:,2),'k.');

      hold on
%       for i_d=1:3
%           tt=0:1:XMAX;
%           yy1=ones(1,XMAX+1).*(center_yc+(i_d-2)*Distance+R_cir);          
%           yy2=ones(1,XMAX+1).*(center_yc+(i_d-2)*Distance-R_cir);
%           plot(tt,yy1,'k--',tt,yy2,'k--')
%       end
        axis equal

  CurrFrame = getframe;   % 获取像素，否则无法显示动画     
  im{step/(real_YMAX*real_XMAX*intervel_step)} = frame2im(CurrFrame);  

   end
   
end
%%
hold on
axis equal
for i=1:Num_cell    
    for j=1:k_T    
       YX{i,1}(j,1)=XY{j,1}(i,1);
       YX{i,1}(j,2)=XY{j,1}(i,2);
    end
end

for i=1
  h1=plot(YX{i,1}(25:end,1),YX{i,1}(25:end,2),'Color',[i/Num_cell rand 1-i/Num_cell],'linewidth',1.5);   
end

h1=plot(XY{k_T,1}(40:end,1),XY{k_T,1}(40:end,2),'k.');     

%%
for i=1:N_step
    im_X(i)=size(im{i},1);
    im_Y(i)=size(im{i},2);
end
figure(200)
plot(1:N_step,im_X,'b',1:N_step,im_Y,'r')
%%
tic
StratF=21;
imwrite(im{StratF},Name);  
for i=StratF+1:N_step   
    if ~isempty(im{i})
    imwrite(im{i},Name,'WriteMode','append')
    else
        break;
    end
end
toc
% sum(t2)
%%
function perimeter_diff=perimeter_local_1(Cell_area,y_select,x_select)
perimeter_diff=sum(sum(2*Cell_area(y_select-1:y_select+1,x_select-1:x_select+1)-1))-1;
end
function perimeter_diff=perimeter_local_2(Cell_area,y_select,x_select)
perimeter_diff=-sum(sum(2*Cell_area(y_select-1:y_select+1,x_select-1:x_select+1)-1))-1;
end 