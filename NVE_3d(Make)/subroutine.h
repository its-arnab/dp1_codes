#ifndef SUBROUTINE_H
    #define SUBROUTINE_H
    
    /* `````` IN initilize.c File ``````` */
    void global_variable_initializer();
    void allocate_arrays();
    void load_initial_position();
    double rand_01();
    void load_initial_velocity();
    void load_final_config();
    void Initilize_Program(int mode);
    void deallocate_arrays();


   /* `````` IN cell_list.c File ``````` */
   double particle_distance(int i, int j);
   int getcell_index(int i, int j, int k);
   void maps();
   void construct_cell_list();
   void construct_verlet_list();

  
    /* `````` IN time_evolve.c File ``````` */
   void force_calculation();
   void Temperature();
   void Initilize_force_using_cell();
   void integrate();

#endif
