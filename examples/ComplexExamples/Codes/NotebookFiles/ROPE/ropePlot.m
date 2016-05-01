function ropePlot()

x = hdf5read('ropePlotx.h5','x');
y = hdf5read('ropePloty.h5','y');
t = hdf5read('ropePlott.h5','t');

plot3(x(1,:),t(1)*ones(41),y(1,:),'k');
set(gca,'Zdir','reverse')
set(gca,'Xdir','reverse')
hold on
for i=2:18:3350
    plot3(x(i,:),t(i)*ones(41),y(i,:),'k');
end
plot3(x(end,:),t(end)*ones(41),y(end,:),'k');
axis([-1.8 1.8 -1 4 -1 2])

plot3(x(:,end),t,y(:,end),'k','LineWidth',2);
plot3(x(:,1),t,y(:,1),'k');

hold off
end
