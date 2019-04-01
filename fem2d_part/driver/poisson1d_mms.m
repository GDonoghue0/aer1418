

h = 1./2.^(1:5);
p = 1:2;

L2_err = zeros(5,2);
sH1_err = zeros(5,2);
out_err = zeros(5,2);
for i = 1:5
    for j = 1:2
        [L2_err(i,j), sH1_err(i,j), out_err(i,j)] = poisson1d(h(i),p(j));
    end
end

L2_conv = log(L2_err(end,:)./L2_err(end-1,:)) ./ log(h(end)/h(end-1));
sH1_conv = log(sH1_err(end,:)./sH1_err(end-1,:)) ./ log(h(end)/h(end-1));
out_conv = log(out_err(end,:)./out_err(end-1,:)) ./ log(h(end)/h(end-1));

conv = [L2_conv; sH1_conv; out_conv];
conv_expected = [2 3; 1 2; 2 4];


figure(2)
loglog(h, L2_err(:,1),'-o', h, L2_err(:,2),'-o','LineWidth',1.5)
title('L^2 Norm of Poisson Solution Error')
xlabel('$h$','interpreter','latex')
ylabel('$\|u - u_h\|_{L^2(\Omega)}$','interpreter','latex')
legend('p = 1','p = 2','location','northwest')

figure(3)
loglog(h, sH1_err(:,1),'-o', h, sH1_err(:,2),'-o','LineWidth',1.5)
title('H^1 Seminorm of Poisson Solution Error')
xlabel('$h$','interpreter','latex')
ylabel('$|u - u_h|_{H^1(\Omega)}$','interpreter','latex')
legend('p = 1','p = 2','location','northwest')

figure(4)
loglog(h, out_err(:,1),'-o', h, out_err(:,2),'-o','LineWidth',1.5)
title('Output Error')
xlabel('$h$','interpreter','latex')
ylabel('$|\ell^\circ(u) - \ell^\circ(u_h)|$','interpreter','latex')
legend('p = 1','p = 2','location','northwest')
