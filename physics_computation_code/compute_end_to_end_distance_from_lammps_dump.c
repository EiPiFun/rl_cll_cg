#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

float compute_vector_length(float vector_x, float vector_y, float vector_z){

    float vector_length;
    vector_length = sqrt(vector_x*vector_x+vector_y*vector_y+vector_z*vector_z);

    return(vector_length);

}

int main(int argc, char *argv[]){

    char *target_lammps_dump_file_location = argv[1];

    int backbone_atom_type;
    sscanf(argv[2],"%d",&backbone_atom_type);

    char *output_end_to_end_distance_file_location = argv[3];

    printf("\nTarget lammps dump file: %s\nBackbone atom type: %d\nOutput end to end distance file: %s\n", target_lammps_dump_file_location,backbone_atom_type,output_end_to_end_distance_file_location);

    FILE *target_lammps_dump_file_pointer;

    FILE *output_end_to_end_distance_file_pointer;

    char temp_row[99];

    int frame_count = 0;
    int total_frame_number;

    int temp_atom_id;
    int atom_count = 0;
    int one_frame_atom_number;
    int one_frame_backbone_atom_number;
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

    printf("\nOne frame atom number: %d\nOne frame backbone atom number: %d\nTotal frame number: %d\n\n", one_frame_atom_number,one_frame_backbone_atom_number,total_frame_number);

    float temp_atom_charge;
    float temp_atom_coordinate[3];
    float temp_vector[3];
    float temp_end_to_end_distance;

    float (*backbone_atom_coordinate_array)[3] = malloc(sizeof(float)*3*(one_frame_backbone_atom_number+1));
    float *per_frame_end_to_end_distance_array = malloc(sizeof(float)*total_frame_number);

    target_lammps_dump_file_pointer = fopen(target_lammps_dump_file_location, "r");
    output_end_to_end_distance_file_pointer = fopen(output_end_to_end_distance_file_location, "w");
    fprintf(output_end_to_end_distance_file_pointer, "%s\n", "# frame count and end to end distance in nm, the last line is the average, %d %.6f");

    for(int i=0;i<9;i++){
        fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
    }

    atom_count = 0;
    frame_count = 0;

    while(!feof(target_lammps_dump_file_pointer)){

        if(sscanf(temp_row, "%f%f%f", &box_limit[0][1],&box_limit[1][1],&box_limit[2][1])==2){

            for(int i=0;i<3;i++){
                temp_vector[i] = backbone_atom_coordinate_array[0][i]-backbone_atom_coordinate_array[one_frame_backbone_atom_number-1][i];
            }
            temp_end_to_end_distance = compute_vector_length(temp_vector[0],temp_vector[1],temp_vector[2]);
            per_frame_end_to_end_distance_array[frame_count] = temp_end_to_end_distance;

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

    for(int i=0;i<3;i++){
        temp_vector[i] = backbone_atom_coordinate_array[0][i]-backbone_atom_coordinate_array[one_frame_backbone_atom_number-1][i];
    }
    temp_end_to_end_distance = compute_vector_length(temp_vector[0],temp_vector[1],temp_vector[2]);
    per_frame_end_to_end_distance_array[frame_count] = temp_end_to_end_distance;

    float average_end_to_end_distance = 0.0;
    for(int i=0;i<total_frame_number;i++){
        average_end_to_end_distance = per_frame_end_to_end_distance_array[i]+average_end_to_end_distance;
    }
    average_end_to_end_distance = average_end_to_end_distance/total_frame_number;

    printf("Average end to end distance: %.6f\n",average_end_to_end_distance);

    for(int i=0;i<total_frame_number;i++){
        fprintf(output_end_to_end_distance_file_pointer,"%d %.6f\n",i+1,per_frame_end_to_end_distance_array[i]);
    }
    fprintf(output_end_to_end_distance_file_pointer,"%d %.6f\n",-1,average_end_to_end_distance);

    fclose(output_end_to_end_distance_file_pointer);

    free(backbone_atom_coordinate_array);
    free(per_frame_end_to_end_distance_array);

    return(0);

}
