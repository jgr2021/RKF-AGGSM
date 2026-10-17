function validate_release
here=fileparts(mfilename('fullpath'));
synthetic=run_synthetic(3);urban=run_urban;
assert(max([synthetic.archive_max_abs_difference])<1e-7,'Synthetic archive replay mismatch');
assert(urban.archive_max_abs_difference<1e-7,'Urban archive replay mismatch');
report=struct('matlab',version,'synthetic_conditions',32,'trajectories_per_condition',3,'methods',14,'synthetic_max_abs_difference',max([synthetic.archive_max_abs_difference]),'urban_max_abs_difference',urban.archive_max_abs_difference,'passed',true);
f=fopen(fullfile(here,'generated/validation.json'),'w');fprintf(f,'%s',jsonencode(report,PrettyPrint=true));fclose(f);
disp(report);
end
