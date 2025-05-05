#ifndef SUBROUTINE_H
    #define SUBROUTINE_H
    
    /* `````` IN initilize.c File ``````` */
    void global_variable_initializer();
    void allocate_arrays();
    void load_initial_position();
    void initialise_kind_random();
    void load_initial_velocity();
    void load_final_config();
    void Initialize_Program(int mode);
    void deallocate_arrays();


    /* `````` IN utils.c File ``````` */
    void md_step_generator();
    double rand_01();
    double random_normal();
    

   /* `````` IN cell_list.c File ``````` */
    double particle_distance(int i, int j);
    int getcell_index(int i, int j, int k);
    void maps();
    void construct_cell_list();
    void construct_verlet_list();
    void force_calculation();
    void Initilize_force_using_cell();

  
    /* `````` IN time_evolve.c File ``````` */
    void Temperature();
    void Apply_PBC();
    void integrate_nve();
    void Initilize_NoseHoover();
    void U1_propagator();
    void U2_propagator();
    void U3_propagator();
    void U4_propagator(int start, int end);
    void integrate_nvt();   
    void compute_H_eff();


    /* `````` IN file_handle.c File ``````` */
    void create_path(const char *path);
    int directory_exists(const char *path);
    void clean_directory(const char *path);
    void save_data_infile(int step, const char *path);
    void save_final_config(const char *path);

#endif
