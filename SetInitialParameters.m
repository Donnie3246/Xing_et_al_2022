function cells=SetInitialParameters(sigma)
global TARGET_CELL_SIZE N_cell
cell_numb=N_cell;
for i=1:cell_numb
    cells.area(i,1)=sum(sum(sigma==i));
end
cells.target_area(1:cell_numb,1)=TARGET_CELL_SIZE;
cells.type=randi(1,cell_numb,1); % 1: light cell, 2: dark cell
end
