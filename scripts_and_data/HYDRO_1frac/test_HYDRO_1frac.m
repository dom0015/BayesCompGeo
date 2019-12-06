load('uloha2.mat')
hydro_problem.mat_frac=1e-7;
[PRESSURE,u0,GRAD,Q,PRESSURE_diff] = tocouple_handle(D,hydro_problem);
    
u_real_approx=[
    0.2609
    0.1332
    0.1102
    0.0612
    0.1902
    0.0379
    0.1439
    0.0269
    0.0952
    0.0205];

[data_generator, ~] = set_fracture(10,[linspace(0,1,81)' zeros(81,1)]);
wrapper_hydro(-6,randn(10,1),hydro_problem,data_generator)
Q(2:end)'

