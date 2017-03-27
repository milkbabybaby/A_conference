load weight_circle

 figure (6);
  plot(weight_circle(:,1),'r-*','MarkerSize',6,'LineWidth', 2)
  hold on;
  plot(weight_circle(:,2),'b--x','MarkerSize',6,'LineWidth', 2)
  hold on;
  plot(weight_circle(:,3),'g-o','MarkerSize',6,'LineWidth', 2)
  hold on;
    plot(weight_stem(:,4),'m-.','MarkerSize',6,'LineWidth', 2)
    hold on;
     plot(weight_stem(:,5),'k:','MarkerSize',6,'LineWidth',2)


%set(gca,'TickLabelInterpreter','latex');
%  legend('A_{12}','A_{23}','A_{34}','A_{45}','A_{56}','A_{61}','Location','East');   
  set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
 xlabel('Iterations','FontName','Times New Roman','FontWeight','bold');
 ylabel('Edge weight ','FontName','Times New Roman','FontWeight','bold');
  export_fig weight_changes_1.eps -painters -transparent 