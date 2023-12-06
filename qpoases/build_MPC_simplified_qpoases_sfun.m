clear mex;
USE_CPP = 1;

DYNOPT_PATH = getenv("DYNOPT_DIR");

solver = 'qpoases';

sfunc_file = ["MPC_simplified_qpoases_sfun.c"];
prob_files = ["MPC_simplified_qpoases_prep.c", "MPC_simplified_qpoases.c"];

req_files = ["qp_pcond.c", "qp_qpOASES.c", "qp_solver.c", "math_lib.c"];

req_files_path = [];
for k=1:length(req_files)
	req_files_path = [req_files_path, fullfile(DYNOPT_PATH, req_files(k))];
end

lib_path = fullfile(DYNOPT_PATH, solver);
req_files_path = mat2cell(req_files_path, 1, ones(1,numel(req_files_path)));
problem_files = mat2cell([sfunc_file, prob_files], 1, ones(1,numel([sfunc_file, prob_files])));
if USE_CPP
	mex('-v', ['-I' DYNOPT_PATH], ['-I', fullfile(DYNOPT_PATH, solver, 'include')], ['-L', lib_path], problem_files{:}, req_files_path{:}, fullfile(DYNOPT_PATH,"*.cpp"), ['-l', 'qpOASES'], '-outdir', './');
else
	mex('-v', ['-I' DYNOPT_PATH], ['-I', fullfile(DYNOPT_PATH, solver, 'include')], ['-L', lib_path], problem_files{:}, req_files_path{:}, ['-l', 'qpOASES'], '-outdir', './');
end

prob_files_with_path = ["MPC_simplified_qpoases_prep.c", "MPC_simplified_qpoases.c"];
req_files = arrayfun(@(x) "$(DYNOPT_DIR)/" + x, req_files);
tmp_files = arrayfun(@(x) x+" ", [prob_files_with_path, req_files]);
problem_info = struct(...
	'INCLUDES', ['-I$(DYNOPT_DIR) -I$(DYNOPT_DIR)/', solver ,'/include/'],...
	'LIB_PATH', ['-L$(DYNOPT_DIR)/', solver],...
	'LDFLAGS', ['-lqpOASES'],...
	'SOURCES', strjoin(tmp_files));
