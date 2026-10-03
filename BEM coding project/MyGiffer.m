function MyGiffer(xx,yy,M,V, GifName)
% versione modificata di MyGiffer fornito nella pagina del corso, ha anche
% il plot del poligono
figure('Position',[100 100 700 700]);   % dimensioni fisse
axes('Position',[0 0 1 1]);             % occupa tutta la figura
% p = polyshape(real(V), imag(V));
% 
% [inside,boundary] = isinterior(p, xx, yy);

pcolor(xx,yy,(real(M)));
hold on
colormax = max(max(abs(M)));

clim([-colormax,colormax]);
plot_polygon(V)

axis tight;  axis equal;  shading interp;
set(gca, 'visible', 'off')
f = getframe;
[im,map] = rgb2ind(f.cdata,256,'nodither');
for t_count = 1:25  % loop over 25 frames
    t = t_count*2*pi/25;
    
    pcolor(xx,yy,(real(exp(-t*1i)*M)))
    hold on
    plot_polygon(V)
    
    axis tight;    axis equal;   shading interp;
    f = getframe; %(gcf);
    im(:,:,1,t_count) = rgb2ind(f.cdata,map,'nodither');
end
imwrite(im,map,GifName,'gif','DelayTime',0,'LoopCount',inf);
hold off;
close;
end     
