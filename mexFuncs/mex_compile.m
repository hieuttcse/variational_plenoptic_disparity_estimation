mex -g -I'./libs/' ./libs/lf4d_matlab.c ./libs/mat_matlab.c ./lf_downscale.c -output lf_downscale

mex -g -I'./libs/' ./libs/lf4d_matlab.c ./libs/mat_matlab.c ./img_resample.c -output img_resample
mex -g -I'./libs/' ./libs/lf4d_matlab.c ./libs/mat_matlab.c ./lf_warp_with_z.c -output lf_warp_with_z
% %
% mex -g -I'./libs/' ./libs/mat_matlab.c ./warp_solver_gs.c -output warp_solver_gs
% %
% % mex -g -I'./libs/' ./libs/mat_matlab.c ./warp_solver_gs_robust_data.c -output warp_solver_gs_robust_data
mex -g -I'./libs/' ./libs/mat_matlab.c ./warp_solver_gs_robust_data_flow_driven.c -output warp_solver_gs_robust_data_flow_driven
%
mex -g -I'./libs/' ./libs/mat_matlab.c ./warp_solver_gs_color_robust_data_flow_driven.c -output warp_solver_gs_color_robust_data_flow_driven
%
%
% mex -g  -I'./libs/' ./libs/lf4d_matlab.c ./libs/mat_matlab.c ./libs/commons.c ./solver_gs.c -output solver_gs
%
% mex -g  -I'./libs/' ./libs/lf4d_matlab.c ./libs/mat_matlab.c ./libs/commons.c ./hole_filling.c -output hole_filling
%
% %
% mex -g -I'./libs/' ./libs/lf4d_matlab.c ./libs/mat_matlab.c ./img_warp_with_z.c -output img_warp_with_z



% mex -g -I'./libs/' ./libs/mat_matlab.c ./warp_solver_gs_robust_data_flow_driven_image_driven.c -output warp_solver_gs_robust_data_flow_driven_image_driven
%
% mex -g -I'./libs/' ./libs/mat_matlab.c ./warp_solver_gs_robust_data_flow_driven_image_driven_data_weight.c -output warp_solver_gs_robust_data_flow_driven_image_driven_data_weight
