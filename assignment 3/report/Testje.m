x = 1:10000;
y = cosd(x);
y2 = sind(x);

xi = x(1);
yi = y(1);
yi2 = y2(1);
figure(1); clf;
h = plot(xi, yi);
hold on;
h2 = plot(xi, yi2);

for k = 200:10000
    xi = x(1:k);
    yi = y(1:k);
    yi2 = y2(1:k);
    set(h,'Xdata',xi,'YData',yi);
    set(h2,'Xdata',xi,'YData',yi2);
    drawnow;
end;