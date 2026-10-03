function plot_polygon(V_collection)
for V_ = V_collection
    V = V_{:};
    % disp(V)
    p = polyshape(real(V), imag(V));
    plot(p, "FaceColor",[1,1,1], "EdgeColor", [1,1,1], "FaceAlpha",1)
    
    hold on
    axis equal
end
end