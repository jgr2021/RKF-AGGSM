function report=run_synthetic(nRuns,conditionIds)
% Replay frozen manuscript inputs. Metric: mean of per-trajectory RMSEs.
% run_synthetic(3) is a smoke/replay check; run_synthetic(1000) is full replay.
if nargin<1,nRuns=3;end
if nargin<2,conditionIds=1:32;end
assert(nRuns>=1 && nRuns<=1000 && nRuns==fix(nRuns));
assert(all(ismember(conditionIds,1:32)));
here=fileparts(mfilename('fullpath'));setup_journal('synthetic');
meta=jsondecode(fileread(fullfile(here,'parameters.json')));report=struct([]);
for id=conditionIds
 if iscell(meta.conditions),q=meta.conditions{id};else,q=meta.conditions(id);end
 tag=sprintf('%s_%s',q.family,q.paper_condition);
 b=load(fullfile(here,q.bank));expected=load(fullfile(here,'data',[tag '_expected.mat']));
 scores=nan(nRuns,14);seconds=nan(nRuns,14);
 for r=1:nRuns
  for j=1:14
   rng(b.seed_base+20*r+b.seed_offsets(j),'twister');
   [x,P,tm]=call_method(j,b.observations(:,:,r),b.truth(:,:,r),b);
   assert(isreal(x)&&all(isfinite(x(:)))&&all(isfinite(P(:))),'Nonfinite filter output');
   scores(r,j)=sqrt(mean(sum((x(:,1:2)-b.truth(:,1:2,r)).^2,2)));
   seconds(r,j)=tm.total_seconds;
  end
 end
 delta=max(abs(scores-expected.trajectory_rmse(1:nRuns,:)),[],'all');
 rec=struct('condition',tag,'runs',nRuns,'rmse',mean(scores,1),'archive_max_abs_difference',delta);
 if isempty(report),report=rec;else,report(end+1)=rec;end %#ok<AGROW>
 fprintf('%s: %d runs, max archive RMSE difference %.3g\n',tag,nRuns,delta);
 out=fullfile(here,'generated');if ~exist(out,'dir'),mkdir(out);end
 save(fullfile(out,[tag '_replay.mat']),'scores','seconds','rec','-v7');
end
f=fopen(fullfile(here,'generated/synthetic_replay.json'),'w');fprintf(f,'%s',jsonencode(report,PrettyPrint=true));fclose(f);
end
