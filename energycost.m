

x=[1,2,3,4,5]


y=[6.9269*10^6,1.9722*10^4,900.7143,1.9722*10^4, 6.9269*10^6]


x2=[ 6,7,8]

y2=[6.927*10^6,1.9723*10^4,900.7140]

figure (1);
 bar(x,log10(y),'g','BarWidth',0.5)
  
  hold on;
   bar(x2,log10(y2),'y','BarWidth',0.5)
  set(gca, 'XTick', [1 2 3 4 5 6 7 8 ])
  set(gca,'XTickLabel',{'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'})
   legend('stem','circle')
   set(gca, 'LineWidth', 1.5);
      set(gca,'FontName','Times New Roman','FontWeight','bold')
 xlabel('Network','FontName','Times New Roman','FontWeight','bold');
 ylabel('lg E(A^*) ','FontName','Times New Roman','FontWeight','bold');

  export_fig energy.eps -painters -transparent 