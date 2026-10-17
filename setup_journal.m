function setup_journal(mode)
% Install only this release's paths; avoid mixing legacy filter versions.
if nargin<1,mode='synthetic';end
here=fileparts(mfilename('fullpath'));
dirs={fullfile(here,'shared'),fullfile(here,'synthetic'),fullfile(here,'synthetic/filters'),fullfile(here,'urban'),fullfile(here,'urban/filters')};
for k=1:numel(dirs)
 if any(strcmp(strsplit(path,pathsep),dirs{k})),rmpath(dirs{k});end
end
addpath(fullfile(here,'shared'),'-begin');
assert(any(strcmp(mode,{'synthetic','urban'})),'Unknown mode');
addpath(fullfile(here,mode),fullfile(here,mode,'filters'),'-begin');
maxNumCompThreads(1);
end
