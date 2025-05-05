#include "global.h"


/*~~~~~~ Function to create new directory ~~~~~~~~~*/
void create_path(const char *path) {
    char tmp[512];
    snprintf(tmp, sizeof(tmp), "%s", path);
    size_t len = strlen(tmp);

    for (size_t i = 1; i < len; i++) {
        if (tmp[i] == '/') {
            tmp[i] = '\0';
            mkdir(tmp, 0777);
            tmp[i] = '/';
        }
    }
    mkdir(tmp, 0777);
}


/*~~~~~~ Function to check if directory exists ~~~~~~~~~*/
int directory_exists(const char *path) {
    struct stat st;
    return (stat(path, &st) == 0 && S_ISDIR(st.st_mode));
}


/*~~~~~~ Function to clean directory ~~~~~~~~~*/
void clean_directory(const char *path) {
    char cmd[512];
    snprintf(cmd, sizeof(cmd), "rm -rf %s/*", path);
    int status = system(cmd);
    if (status != 0) {
        fprintf(stderr, "Failed to clean directory: %s\n", path);
    }
}


/*~~~~~~ Function to save data in file ~~~~~~~~~*/
void save_data_infile(int step, const char *path)
{
    char filename[256];
    snprintf(filename, sizeof(filename), 
    "%s/step_%07d", path, step); // "step_0001.txt", "step_0002.txt", ...

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Error opening file");
        exit(1);
    }

    for (int i = 0; i < N; i++)
    fprintf(file, "%d %.10lf %.10lf %.10lf\n", pt_kind[i], xu[i], yu[i], zu[i]);


    fclose(file);
}


/*~~~~~~~~ Save final config ~~~~~~~~~*/
void save_final_config(const char *path)
{
    FILE *fid3;
    char final_config_path[256];
    snprintf(final_config_path, sizeof(final_config_path), "%s/final_config", path);
    fid3 = fopen(final_config_path, "w");
    if( fid3 == NULL  ) exit(1);

    for (int i = 0; i < M_chain; i++)
    {
        fprintf(fid3, "%.12lf \t %.12lf \t %.12lf \n", eta[i], p_eta[i], Q[i]);
    }

    for( int i = 0; i < N; i++ )
    {
        fprintf(fid3, "%.12lf\t %.12lf\t %.12lf\t %.12lf\t %.12lf\t %.12lf\t %d\n", 
                x[i], y[i], z[i], vx[i], vy[i], vz[i], pt_kind[i]);
    }
    fclose(fid3); 
}