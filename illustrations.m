

t = 0:0.01:6;
R = smoothstep(t-1) .* (1-smoothstep(t-4));
G = smoothstep(t-1.1) .* (1-smoothstep(t-4.1));
B = smoothstep(t-1.2) .* (1-smoothstep(t-4.2));

figure(1);
hold on
plot(t*10, R, 'r', 'linewidth', 2);
plot(t*10, G, 'g', 'linewidth', 2);
plot(t*10, B, 'b', 'linewidth', 2);
hold off
legend('R', 'G', 'B');
title('Channel intensities along an image line');
axis([0, 60, -.1, 1.1]);
xlabel('x');


R = smoothdelta(t-1) - smoothdelta(t-4);
G = smoothdelta(t-1.1) - smoothdelta(t-4.1);
B = smoothdelta(t-1.2) - smoothdelta(t-4.2);

figure(2);
hold on
plot(t*10, R, 'r', 'linewidth', 2);
plot(t*10, G, 'g', 'linewidth', 2);
plot(t*10, B, 'b', 'linewidth', 2);
hold off
legend('R_x', 'G_x', 'B_x');
title('x derivative of channel intensities along an image line');
axis([0, 60, -1.6, 1.6]);
xlabel('x');
