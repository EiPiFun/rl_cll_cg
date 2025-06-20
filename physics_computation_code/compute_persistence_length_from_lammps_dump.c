#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float compute_vector_length(float vector_x, float vector_y, float vector_z){

    float vector_length;
    vector_length = sqrt(vector_x*vector_x+vector_y*vector_y+vector_z*vector_z);

    return(vector_length);

}

float compute_angle_cosine_between_two_vector(float vector_1_x, float vector_1_y, float vector_1_z, float vector_2_x, float vector_2_y, float vector_2_z){

    float dot_product;
    float vector_size_1;
    float vector_size_2;
    float angle_cosine;

    dot_product = vector_1_x*vector_2_x+vector_1_y*vector_2_y+vector_1_z*vector_2_z;
    vector_size_1 = compute_vector_length(vector_1_x,vector_1_y,vector_1_z);
    vector_size_2 = compute_vector_length(vector_2_x,vector_2_y,vector_2_z);
    angle_cosine = dot_product/vector_size_1/vector_size_2;

    return(angle_cosine);
}

void compute_one_frame_autocorrelation(int considered_backbone_atom_number, int backbone_atom_interval_number, int frame_count, float *backbone_atom_coordinate_array_pointer, float *per_frame_average_bond_length_array_pointer, float *per_frame_autocorrelation_array_pointer){

    float (*bond_vector)[3] = malloc(sizeof(float)*3*(considered_backbone_atom_number-1));

    float *autocorrelation_array = malloc(sizeof(float)*(considered_backbone_atom_number-1));

    float average_bond_length = 0.0;
    float fitter_persistence_bond_number;
    float fitted_persistence_length;

    float temp_vector_1[3];
    float temp_vector_2[3];
    float average_angle_cosine;

    for(int i=0;i<considered_backbone_atom_number-1;i++){
        for(int j=0;j<3;j++){
            temp_vector_1[j] = backbone_atom_coordinate_array_pointer[3*((1+backbone_atom_interval_number)*i+1)+j]-backbone_atom_coordinate_array_pointer[3*(1+backbone_atom_interval_number)*i+j];
            bond_vector[i][j] = temp_vector_1[j];
        }
        average_bond_length = compute_vector_length(temp_vector_1[0],temp_vector_1[1],temp_vector_1[2])+average_bond_length;
    }
    average_bond_length = average_bond_length/(considered_backbone_atom_number-1);

    per_frame_average_bond_length_array_pointer[frame_count] = average_bond_length;

    for(int i=0;i<considered_backbone_atom_number-1;i++){

        average_angle_cosine = 0.0;
        for(int j=0;j<considered_backbone_atom_number-1-i;j++){
            for(int k=0;k<3;k++){
                temp_vector_1[k] = bond_vector[j][k];
                temp_vector_2[k] = bond_vector[j+i][k];
            }
            average_angle_cosine = compute_angle_cosine_between_two_vector(temp_vector_1[0],temp_vector_1[1],temp_vector_1[2],temp_vector_2[0],temp_vector_2[1],temp_vector_2[2])+average_angle_cosine;   
        }
        average_angle_cosine = average_angle_cosine/(considered_backbone_atom_number-1-i);
        per_frame_autocorrelation_array_pointer[i+(considered_backbone_atom_number-1)*frame_count] = average_angle_cosine;
    }
    free(bond_vector);
    free(autocorrelation_array);
}

int main(int argc, char *argv[]){

    char *target_lammps_dump_file_location = argv[1];

    int backbone_atom_type;
    sscanf(argv[2],"%d",&backbone_atom_type);

    int backbone_atom_interval_number;
    sscanf(argv[3], "%d", &backbone_atom_interval_number);

    char *output_persistence_length_file_location = argv[4];

    printf("\nTarget lammps dump file: %s\nBackbone atom type: %d\nBackbone atom interval number: %d\nOutput persistence length file: %s\n", target_lammps_dump_file_location,backbone_atom_type,backbone_atom_interval_number,output_persistence_length_file_location);

    FILE *target_lammps_dump_file_pointer;

    FILE *output_persistence_length_file_pointer;

    char temp_row[99];

    int frame_count = 0;
    int total_frame_number;

    int temp_atom_id;
    int atom_count = 0;
    int one_frame_atom_number;
    int one_frame_backbone_atom_number;
    int one_frame_considered_backbone_atom_number;
    int temp_atom_type;

    int temp_residue_id;

    float box_limit[3][2];
    float box_size[3];

    target_lammps_dump_file_pointer = fopen(target_lammps_dump_file_location, "r");
    for(int i=0;i<4;i++){
        fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
    }
    sscanf(temp_row, "%d", &one_frame_atom_number);
    for(int i=0;i<5;i++){
        fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
    }
    while(!feof(target_lammps_dump_file_pointer)){
        if(sscanf(temp_row, "%f%f%f", &box_limit[0][1],&box_limit[1][1],&box_limit[2][1])==2){
            break;
        }else if(sscanf(temp_row, "%d%d%d\n", &temp_atom_id,&temp_residue_id,&temp_atom_type)==3){
            if(temp_atom_type==backbone_atom_type){
                atom_count = atom_count+1;
            }
        }
        fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
    }
    fclose(target_lammps_dump_file_pointer);

    one_frame_backbone_atom_number = atom_count;
    one_frame_considered_backbone_atom_number = one_frame_backbone_atom_number/(backbone_atom_interval_number+1);

    target_lammps_dump_file_pointer = fopen(target_lammps_dump_file_location, "r");
    for(int i=0;i<5;i++){
        fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
    }
    while(!feof(target_lammps_dump_file_pointer)){
        if(sscanf(temp_row, "%f%f%f", &box_limit[0][1],&box_limit[1][1],&box_limit[2][1])==2){
            frame_count = frame_count+1;
            for(int i=1;i<3;i++){
                fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
            }
        }
        fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
    }
    fclose(target_lammps_dump_file_pointer);

    total_frame_number = frame_count;

    printf("\nOne frame atom number: %d\nOne frame backbone atom number: %d\nOne frame considered atom number: %d\nTotal frame number: %d\n\n", one_frame_atom_number,one_frame_backbone_atom_number,one_frame_considered_backbone_atom_number,total_frame_number);

    float temp_atom_charge;
    float temp_atom_coordinate[3];

    float (*backbone_atom_coordinate_array)[3] = malloc(sizeof(float)*3*(one_frame_backbone_atom_number+1));
    float *per_frame_average_bond_length_array = malloc(sizeof(float)*total_frame_number);
    float (*per_frame_autocorrelation_array)[one_frame_considered_backbone_atom_number-1] = malloc(sizeof(float)*(one_frame_considered_backbone_atom_number-1)*total_frame_number);
    float *average_autocorrelation_array = malloc(sizeof(float)*(one_frame_considered_backbone_atom_number-1));

    target_lammps_dump_file_pointer = fopen(target_lammps_dump_file_location, "r");

    for(int i=0;i<9;i++){
        fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
    }

    atom_count = 0;
    frame_count = 0;

    while(!feof(target_lammps_dump_file_pointer)){

        if(sscanf(temp_row, "%f%f%f", &box_limit[0][1],&box_limit[1][1],&box_limit[2][1])==2){

            compute_one_frame_autocorrelation(one_frame_considered_backbone_atom_number,backbone_atom_interval_number,frame_count,&backbone_atom_coordinate_array[0][0],&per_frame_average_bond_length_array[0],&per_frame_autocorrelation_array[0][0]);
            frame_count = frame_count+1;
            atom_count = 0;

            for(int i=1;i<3;i++){
                fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
            }

        }else if(sscanf(temp_row, "%d%d%d%f%f%f%f\n", &temp_atom_id,&temp_residue_id,&temp_atom_type,&temp_atom_charge,&temp_atom_coordinate[0],&temp_atom_coordinate[1],&temp_atom_coordinate[2])==7){

            if(temp_atom_type==backbone_atom_type){
                /* Lammps units "real" distance unit is 1.0 Angstroms or 0.1 nm */
                for(int i=0;i<3;i++){
                    backbone_atom_coordinate_array[atom_count][i] = 0.1*temp_atom_coordinate[i];
                }
                atom_count = atom_count+1;
            }

        }
        fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
    }

    fclose(target_lammps_dump_file_pointer);

    compute_one_frame_autocorrelation(one_frame_considered_backbone_atom_number,backbone_atom_interval_number,frame_count,&backbone_atom_coordinate_array[0][0],&per_frame_average_bond_length_array[0],&per_frame_autocorrelation_array[0][0]);

    float average_bond_length = 0.0;
    for(int i=0;i<total_frame_number;i++){
        average_bond_length = per_frame_average_bond_length_array[i]+average_bond_length;
    }
    average_bond_length = average_bond_length/total_frame_number;

    float temp_autocorrelation;
    for(int i=0;i<one_frame_considered_backbone_atom_number-1;i++){
        temp_autocorrelation = 0.0;
        for(int j=0;j<total_frame_number;j++){
            temp_autocorrelation = per_frame_autocorrelation_array[j][i]+temp_autocorrelation;

        }
        temp_autocorrelation = temp_autocorrelation/(total_frame_number);
        average_autocorrelation_array[i] = temp_autocorrelation;
    }

    output_persistence_length_file_pointer = fopen(output_persistence_length_file_location, "w");
    fprintf(output_persistence_length_file_pointer, "%s\n", "# Bond interval number and autocorrelation, %6d%12.6f\n# Last line is average bond length and persistence length in nm, %.6f %.6f");
    for(int i=0;i<one_frame_considered_backbone_atom_number-1;i++){
        fprintf(output_persistence_length_file_pointer, "%6d%12.6f\n", i, average_autocorrelation_array[i]);
    }

    float temp_xx_dot = 0.0;
    float temp_xy_dot = 0.0;
    int fit_considered_range = one_frame_considered_backbone_atom_number-1;
    for(int i=0;i<fit_considered_range;i++){
        temp_xx_dot = i*i+temp_xx_dot;
        temp_autocorrelation = average_autocorrelation_array[i];
        if(temp_autocorrelation>0.0){
            temp_xy_dot = i*log(temp_autocorrelation)+temp_xy_dot;
        }else{
            temp_xy_dot = i*log(0.01)+temp_xy_dot;
        }
    }

    float fitted_k = temp_xy_dot/temp_xx_dot;
    float fitted_persistence_bond_number = -1.0/fitted_k;
    float fitted_persistence_length = -average_bond_length/fitted_k;
    fprintf(output_persistence_length_file_pointer, "%.6f %.6f\n", average_bond_length, fitted_persistence_length);

    fclose(output_persistence_length_file_pointer);

    free(backbone_atom_coordinate_array);
    free(per_frame_average_bond_length_array);
    free(per_frame_autocorrelation_array);
    free(average_autocorrelation_array);

    return(0);

}
