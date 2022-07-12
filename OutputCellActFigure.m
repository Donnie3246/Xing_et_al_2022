
function OutputCellActFigure(sigma, sigma_Act)

h1=pcolor(sigma_Act+1);
 set(h1, 'EdgeColor', 'none');
 %jet1=[1 1 1;cool];
 %colormap(jet1)% default
%   colorbar
%  caxis([-1 100])
% 
% hold on
% [y1, x1] = find(sigma ~= circshift(sigma,[1 0]));
% [y2, x2] = find(sigma ~= circshift(sigma,[0 1]));
% line_x=horzcat([x1 x1+1]',[x2 x2]');
% line_y=horzcat([y1 y1]',[y2 y2+1]');
% h2=plot(line_x,line_y,'-','linewidth',0.1);
% set(h2, 'Color', 'k'); 
% hold off

axis off

 %pause(0.0001)
end