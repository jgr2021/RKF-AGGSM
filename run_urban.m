function report=run_urban
% Replay the 14 manuscript methods on Odaiba-Trimble, scoring epochs 51:end.
here=fileparts(mfilename('fullpath'));setup_journal('urban');
t=load(fullfile(here,'data/urban_configuration.mat'));c=t.config;
b=load(fullfile(here,'data/Odaiba_trimble.mat'));e=load(fullfile(here,'data/urban_expected.mat'));
N=size(b.z,1);X=zeros(N,2,14);rmse=zeros(1,14);mae=rmse;seconds=rmse;delta=rmse;
for j=1:14
 [x,P,tm]=tokyo_call(j,b,c);
 assert(isreal(x)&&all(isfinite(x(:)))&&all(isfinite(P(:))));
 X(:,:,j)=x;error=sqrt(sum((x(51:end,:)-b.pos(51:end,:)).^2,2));
 rmse(j)=sqrt(mean(error.^2));mae(j)=mean(error);seconds(j)=tm.total_seconds;
 delta(j)=max(abs(x-e.expected_X(:,:,j)),[],'all');
 fprintf('%s: RMSE %.6f m, archived trajectory difference %.3g\n',e.names{j},rmse(j),delta(j));
end
report=struct('names',{e.names},'rmse',rmse,'mae',mae,'seconds',seconds,'archive_max_abs_difference',max(delta));
out=fullfile(here,'generated');if ~exist(out,'dir'),mkdir(out);end
save(fullfile(out,'urban_replay.mat'),'X','report','-v7');
f=fopen(fullfile(out,'urban_replay.json'),'w');fprintf(f,'%s',jsonencode(report,PrettyPrint=true));fclose(f);
end
