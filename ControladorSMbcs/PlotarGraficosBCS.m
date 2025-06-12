
 bcs=ESP();

figure;

bcs.envelope.fig.qH()
hold on
% plot(out.X(:,3)*3600,out.Y(:,2),'k.')
scatter(out.X(:,3)*3600,out.Y(:,2),30,out.time/3600);
plot(out.X(1,3)*3600,out.Y(1,2),'o','MarkerFaceColor',[0,1,0],'MarkerEdgeColor',[0,0,0])
plot(out.X(end,3)*3600,out.Y(end,2),'o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0])
text(out.X(1,3)*3600,out.Y(1,2),'t_0','HorizontalAlignment','left')
text(out.X(end,3)*3600,out.Y(end,2),'t_f','HorizontalAlignment','left')
ylabel('H/(m)');
xlabel('q/(m^3/h)');
h = colorbar;
ylabel(h, 't/(h)');
set(gca,'FontSize',34,'FontName','Time New Roman');



figure;
set(gca,'FontSize',28,'FontName','Time New Roman');
bcs.envelope.fig.fH()
hold on
% plot(out.X(:,4),out.Y(:,2),'k.')
scatter(out.X(:,4),out.Y(:,2),30,out.time/3600);
plot(out.X(1,4),out.Y(1,2),'o','MarkerFaceColor',[0,1,0],'MarkerEdgeColor',[0,0,0])
plot(out.X(end,4),out.Y(end,2),'o','MarkerFaceColor',[1,0,0],'MarkerEdgeColor',[0,0,0])
text(out.X(1,4),out.Y(1,2),'t_0','HorizontalAlignment','left')
text(out.X(end,4),out.Y(end,2),'t_f','HorizontalAlignment','left')
ylabel('H/(m)');
xlabel('f/(Hz)');
h = colorbar;
ylabel(h, 't/(h)');


figure;%Pin

plot(out.time/3600,out.URef*1e-5,'k','LineWidth',3);
hold on;
plot(out.time/3600,out.Y(:,1)*1e-5,'b--','LineWidth',3);
legend('Referência','Planta');
ylabel('p_{in}(t)/(bar)')
xlabel('t/(h)')
set(gca,'FontSize',32,'FontName','Time New Roman');
axis([0 (out.time(end)/3600) -inf inf])

figure;%U

plot(out.time/3600,out.U(:,1),'LineWidth',3);
 hold on;
% plot(out.time/3600,out.U(:,3),'LineWidth',3);
% plot(out.time/3600,out.U(:,2),'LineWidth',3);
% legend('u_{total}','u_{eq}','u_{chave}');
plot([0 out.time(end)]/3600,[35 35],'LineWidth',3);
plot([0 out.time(end)]/3600,[65 65],'LineWidth',3);
ylabel('f_{ref}(t)/(Hz)')
xlabel('t/(h)');
set(gca,'FontSize',34,'FontName','Time New Roman');
axis([0 (out.time(end)/3600) -inf inf])

figure;%Superficie
plot(out.time/3600,out.Sn./60^2,'r--','LineWidth',3);
% legend('Referência','Planta');
ylabel('Superfície/(Pa/min^2)')
xlabel('t/(h)')
axis([0 (out.time(end)/3600) -inf inf])
set(gca,'FontSize',34,'FontName','Time New Roman');

% figure
% scatter(out.edeie(:,1),out.edeie(:,2),200,out.time);
% colorbar
% hold on
% plot(out.edeie(:,1),out.edeie(:,2));
% 
% 
% figure
% scatter(out.Y(:,1),out.dpin(:,1),200,out.time)
% hold on
% plot(out.Y(:,1),out.dpin(:,1));
% colorbar

[rmse,vaf] = rmse_vaf(out.URef*.1e-5,out.Y(:,1)*.1e-5,length(out.Y(:,1)));

%   save('caso2.mat','out','bcs','rmse')