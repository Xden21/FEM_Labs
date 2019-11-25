%% mesh density

x = [2 1 0.5 0.33 0.25 .1];

bag_2d = [0.853 0.855 0.857 0.858 0.8585 0.8591];
bag_3d = [0.856 0.856 0.8597 0.8597 0.8597 0.8597];

l_2d = [0.35581 0.361066 0.364822 0.365948 0.366369  0.367247];
l_3d = [0.444866 0.444866 0.452022 0.452022 0.452022 0.452022];

f_2d = [20.9326 21.2091 21.5636 21.6741 21.7057  21.8137];
f_3d = [19.8196 20.2488 20.4949 20.4949 20.4949 20.4949];

subplot(1,3,1);
plot(x, bag_2d, 'r');hold on;plot(x, bag_3d, 'b');hold off;
title('Induction sensitivity to discretisation')
xlabel('characteristic size in the airgap [ag]')
ylabel('b_{ag} [ T ]');
legend('2D case', '3D case');

subplot(1,3,2)
plot(x, l_2d, 'r');hold on;plot(x, l_3d, 'b');hold off;
title('Inductance sensitivity to discretisation')
xlabel('characteristic size in the airgap [ag]')
ylabel('L [mH]');
legend('2D case', '3D case');

subplot(1,3,3);
plot(x, f_2d, 'r');hold on;plot(x, f_3d, 'b');hold off;
title('Attraction force sensitivity to discretisation')
xlabel('characteristic size in the airgap [ag]')
ylabel('F [N]');
legend('2D case', '3D case');

%% truncation

% Dirichlet
bag_r40_2d_dir = 0.853;
bag_r60_2d_dir = 0.857;
bag_inf_2d_dir = 0.858;

l_r40_2d_dir = 0.361066;
l_r60_2d_dir = 0.361782;
l_inf_2d_dir = 0.362262;

f_r40_2d_dir = 21.2091;
f_r60_2d_dir = 21.2163;
f_inf_2d_dir = 21.2226;

bag_r40_3d_dir = 0;
bag_r60_3d_dir = 0;
bag_inf_3d_dir = 0.0;

l_r40_3d_dir = 0.0;
l_r60_3d_dir = 0.0;
l_inf_3d_dir = 0.0;

f_r40_3d_dir = 0.0;
f_r60_3d_dir = 0.0;
f_inf_3d_dir = 0.0;

%Neumann
bag_r40_2d_neu = 0.862;
bag_r60_2d_neu = 0.860;
bag_inf_2d_neu = 0.857;

l_r40_2d_neu = 0.363938;
l_r60_2d_neu = 0.36303;
l_inf_2d_neu = 0.362271;

f_r40_2d_neu = 21.2402;
f_r60_2d_neu = 21.2299;
f_inf_2d_neu = 21.2215;

bag_r40_3d_neu = 0.0;
bag_r60_3d_neu = 0.0;
bag_inf_3d_neu = 0.0;

l_r40_3d_neu = 0.0;
l_r60_3d_neu = 0.0;
l_inf_3d_neu = 0.0;

f_r40_3d_neu = 0.0;
f_r60_3d_neu = 0;
f_inf_3d_neu = 0;

figure;
subplot(1,3,1);
plot([1 2 3], [bag_r40_2d_dir bag_r60_2d_dir bag_inf_2d_dir]);
hold on;
plot([1 2 3], [bag_r40_2d_neu bag_r60_2d_neu bag_inf_2d_neu]);
xticks([1 2 3])
xticklabels({'40', '60', 'inf'})
title('Induction sensitivity to truncation');
xlabel('Outer diameter [mm]')
ylabel('Inductione [ T ]')
legend('2d Dirichlet', '2d Neumann')

subplot(1,3,2);
plot([1 2 3], [l_r40_2d_dir l_r60_2d_dir l_inf_2d_dir]);
hold on;
plot([1 2 3], [l_r40_2d_neu l_r60_2d_neu l_inf_2d_neu]);
xticks([1 2 3])
xticklabels({'40', '60', 'inf'})
title('Inductance sensitivity to truncation');
xlabel('Outer diameter [mm]')
ylabel('Inductance [ mH ]')
legend('2d Dirichlet', '2d Neumann')

subplot(1,3,3);
plot([1 2 3], [f_r40_2d_dir f_r60_2d_dir f_inf_2d_dir]);
hold on;
plot([1 2 3], [f_r40_2d_neu f_r60_2d_neu f_inf_2d_neu]);
xticks([1 2 3])
xticklabels({'40', '60', 'inf'})
title('Attraction force sensitivity to truncation');
xlabel('Outer diameter [mm]')
ylabel('Attracion force [ N ]')
legend('2d Dirichlet', '2d Neumann')

%% Saturation
x = [1 2 4 6 8 10];

bag_2d_lin = [0.171 0.341 0.683 1.024 1.364 1.708];
bag_3d_lin = [0.171 0.342 0.684 1.027 1.369 1.711];
 
bag_2d_sat = [0.159 0.322 0.649 0.968 1.211 1.322];
bag_3d_sat = [0.166 0.337 0.678 0.985 1.137 1.214];

l_2d_lin = [0.361066 0.361066 0.361066 0.361066 0.361066 0.361066];
l_3d_lin = [0.444866 0.444866 0.444866 0.444866 0.444866 0.444866];
 
l_2d_sat = [0.337167 0.340843 0.343458 0.341533 0.319229 0.27695];
l_3d_sat = [0.432302 0.438406 0.440846 0.425683 0.364579 0.307602];

f_2d_lin = [0.848365 3.39346 13.5738 30.5412 54.2954 84.8365];
f_3d_lin = [0.792783 3.17113 12.6845 28.5402 50.7381 79.2783];

f_2d_sat = [0.737342 3.01478 12.2548 27.2833 42.5532 50.6036];
f_3d_sat = [0.744779 3.06881 12.4391 26.1809 34.689 39.3033];

figure;
subplot(1,3,1);
plot(x, bag_2d_lin, 'r');hold on;plot(x, bag_3d_lin, 'g');plot(x, bag_2d_sat, 'b');plot(x, bag_3d_sat, 'm');hold off;
title('Induction sensitivity to saturation')
xlabel('Current RMS [I]')
ylabel('b_{ag} [ T ]');
legend('2D case (lin)', '3D case(lin)','2D case (sat)', '3D case(sat)');

subplot(1,3,2)
plot(x, l_2d_lin, 'r');hold on;plot(x, l_3d_lin, 'g');plot(x, l_2d_sat, 'b');plot(x, l_3d_sat, 'm');hold off;
title('Inductance sensitivity to saturation')
xlabel('Current RMS [I]')
ylabel('L [mH]');
legend('2D case (lin)', '3D case (lin)','2D case (sat)', '3D case(sat)');

subplot(1,3,3);
plot(x, f_2d_lin, 'r');hold on;plot(x, f_3d_lin, 'g');plot(x, f_2d_sat, 'b');plot(x, f_3d_sat, 'm');hold off;
title('Attraction force sensitivity to saturation')
xlabel('Current RMS [I]')
ylabel('F [N]');
legend('2D case (lin)', '3D case (lin)','2D case (sat)', '3D case(sat)');
