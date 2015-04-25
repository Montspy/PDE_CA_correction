
clear
close all

alphas = [50];
betas =  [0.5];
tfs =    [20];
filenames = {'wood'};

for i = 1:length(alphas)
	a = alphas(i);
	b = betas(i);
	tf = tfs(i)
	f = filenames{i};
	disp(fprintf('alpha=%f, beta=%f, tf=%f, file=%s', a, b, tf, f));
	[J, tt, E] = method2_fun(a, b, tf, [f, '.jpg']);
    
	figure(1);
	image(J);
	title('Corrected image');
	saveas(gcf, ['Results/', f, '/auto_method2_newCFL_', num2str(tf), 's_a', num2str(a), '_b', num2str(b), '.fig']);
	
	figure(2);
	plot(tt, E);
	saveas(gcf, ['Results/', f, '/auto_method2_newCFL_', num2str(tf), 's_a', num2str(a), '_b', num2str(b), '_Energy.fig']);
	
	fileID = fopen(['Results/', f, '/auto_method2_newCFL_', num2str(tf), 's_a', num2str(a), '_b', num2str(b), '_Score.txt'], 'w+');
	fprintf(fileID,'%f\n', score_image(J));
	fclose(fileID);
    close all
end

clear