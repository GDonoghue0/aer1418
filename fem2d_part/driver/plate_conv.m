
h = [0.12, 0.06, 0.03, 0.015];
N = length(h);
A1_ref = 0.00455889016913;
A2_ref = 0.00456105413299;
A3_ref = 0.00456145212005;

% A1_out = zeros(N,1);
% A2_out = zeros(N,1);
% A3_out = zeros(N,1);
% ndof = zeros(N,1);
% 
% for i = 1:N
%     [A1_out(i), ndof(i)] = plate(h(i),1);
%     A2_out(i) = plate(h(i),2, false);
%     A3_out(i) = plate(h(i),2, true);
% end

close all
figure(1)
semilogx(ndof,A1_out,'-o',ndof,A2_out,'-o',ndof,A3_out,'-o',ndof,A3_ref*ones(N,1),'-o','Linewidth',1.5)
title('Output Convergence')
xlabel('ndofs')
ylabel('$\ell^\circ(u_h)$','interpreter','latex')
legend('A1','A2','A3','Reference','location','best')

A1_err = abs(A1_ref-A1_out);
A2_err = abs(A2_ref-A2_out);
A3_err = abs(A3_ref-A3_out);

figure(2)
loglog(ndof,A1_ref-A1_out,'-o',ndof,A2_ref-A2_out,'-o',ndof,abs(A3_ref-A3_out),'-o','Linewidth',1.5)
title('Output Error Convergence')
xlabel('ndofs')
ylabel('$|\ell^\circ(u) - \ell^\circ(u_h)|$','interpreter','latex')
legend('A1','A2','A3','location','best')

A1_conv = log(A1_err(3,:)./A1_err(1,:)) ./ log(ndof(3)/ndof(1))
A2_conv = log(A2_err(3,:)./A2_err(1,:)) ./ log(ndof(3)/ndof(1))
A3_conv = log(A3_err(3,:)./A3_err(1,:)) ./ log(ndof(3)/ndof(1))

