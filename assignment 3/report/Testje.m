x = 1:100;
y = log(x);
DELAY = 0.05;
for i = 1:numel(x)
    clf()
    plot(x,y);
    hold on;
    plot(x(i),y(i),'r*');
    hold on
    pause(DELAY);
    drawnow
end