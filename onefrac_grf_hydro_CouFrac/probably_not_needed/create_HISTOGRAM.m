


boundX=[-10 -1];
boundY=[-10 -1];
S=[SAMPLES(:,2),SAMPLES(:,1)];
SAMPLES_size=size(SAMPLES);

H=hist2d(S,boundY(1):0.1:boundY(2),boundX(1):0.1:boundX(2));
HX=hist(SAMPLES(:,1),boundX(1):0.1:boundX(2));
HY=hist(SAMPLES(:,2),boundY(1):0.1:boundY(2));

%figure; imagesc( boundX, boundY, H);

%label='PRIOR_twofrac'

save(['TWOFRAC_' label '.mat'],'H','boundX','boundY','SAMPLES_size','label','HX','HY')