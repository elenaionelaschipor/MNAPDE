function MyGiffer(xx,yy,M,V, GifName)
% Generate time-harmonic animation of complex field M on points xx/yy, 
% export as gif in file GifName --- not very reliable
% Do not touch mouse & keyboard while running!
figure('Position',[100 100 700 700]);   % dimensioni fisse
axes('Position',[0 0 1 1]);             % occupa tutta la figura
% p = polyshape(real(V), imag(V));
% 
% [inside,boundary] = isinterior(p, xx, yy);

pcolor(xx,yy,(real(M)));
hold on
colormax = max(max(abs(M)));

% plot(p, "FaceColor",[1,1,1], "EdgeColor", [1,1,1], "FaceAlpha",1)

clim([-colormax,colormax]);
plot_polygon(V)

axis tight;  axis equal;  shading interp;
set(gca, 'visible', 'off')
f = getframe; %(gcf); % gcf sometimes avoids errors but gives ugly gray box
[im,map] = rgb2ind(f.cdata,256,'nodither');
for t_count = 1:25  % loop over 25 frames
    t = t_count*2*pi/25;
    
    pcolor(xx,yy,(real(exp(-t*1i)*M)))
    hold on
    plot_polygon(V)

    % plot(p, "FaceColor",[1,1,1], "EdgeColor", [1,1,1], "FaceAlpha",1)
    
    axis tight;    axis equal;   shading interp;
    f = getframe; %(gcf);
    im(:,:,1,t_count) = rgb2ind(f.cdata,map,'nodither');
end
imwrite(im,map,GifName,'gif','DelayTime',0,'LoopCount',inf);
hold off;
close;
end     
