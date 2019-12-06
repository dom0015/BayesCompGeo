SAMPLES = [];

no_confs=1;
no_win=2;

no_measurements=9*no_confs;%no_confs*(no_win+1);
N=5e5; %1e6;
period_saving=1e3;%1e4;%5e3;
C=N/period_saving;
s=10;

%label=['eight2' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['eight1234' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['four2' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['four1234' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['two1' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['two2' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['two3' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['two4' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['two13' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['two24' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['two1234' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];

%label=['grf1' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['grf3' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['grf13' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['grf_sikma1' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['grf2_sikma' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['grf3_sikma' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['grf4_sikma' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['grf_sikma1234' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];

%label=['boreholes1' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['boreholes2' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['boreholes3' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['boreholes4' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['boreholes1234' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];

%label=['boreholes_grf1' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['boreholes_grf2' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['boreholes_grf3' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['boreholes_grf4' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
%label=['boreholes_grf1234' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];
label=['boreholes_grf4_em' '_' num2str(s) '_' num2str(no_measurements) '_' num2str(N) '_col'];

for i=1:C
    load(['res/SAMPLES_' num2str(i) '_' label '.mat'],'SAMPLES_save');
    SAMPLES = [SAMPLES; SAMPLES_save];
	disp(i)
end


