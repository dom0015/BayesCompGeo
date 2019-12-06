for w=[1 2 3 12 23 13 123]
    %load(['res_okno' w '/SAMPLES_nic_2_4_1000_col.mat']);
    t=num2str(w);
    file=dir(['res_okno' t '/SAMPLES*.*']);
    disp(file.name)
    load(['res_okno' t '/' file.name])
    gamma=1;            % prior standard deviation
    mu0=[6 6]';       % prior mean
    u_real=[5 7]';
    visualization( mu0, gamma, u_real )
    plot(SAMPLES(:,1),SAMPLES(:,2),'.');
    title(t)
    xlim([3 9]); ylim([3 9])
end