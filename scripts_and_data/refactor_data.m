d = dir('.');
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..','allres'})) = [];

nfolds=length(nameFolds);

for iindex=1:nfolds
    disp(nameFolds{iindex})
cd(nameFolds{iindex});
cd('res');
conffileList = dir('CONFIGURATION0*');
disp(conffileList.name)

load(conffileList.name);

cd('..');

cd('..');
viz_corrplot_all(Nbatch,label,u_real,gamma,mu0,nameFolds{iindex})
end