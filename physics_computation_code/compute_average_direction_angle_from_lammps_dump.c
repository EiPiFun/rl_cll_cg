#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define pi acos(-1)

float compute_vector_length(float vector_x, float vector_y, float vector_z){

    float vector_length;
    vector_length = sqrt(vector_x*vector_x+vector_y*vector_y+vector_z*vector_z);

    return(vector_length);

}

void periodic_vector_processing(float *original_vector_pointer, float minimum_box_size, float *box_size_pointer){

    /* temp_distance[6] is the orignal distance, temp_distance[7] is the shortest distance */
    float temp_distance[8];
    float translated_vector[6][3];
    int shortest_index = 6;

    temp_distance[6] = compute_vector_length(original_vector_pointer[0],original_vector_pointer[1],original_vector_pointer[2]);

    while(temp_distance[6]>0.5*minimum_box_size){

        translated_vector[0][0] = original_vector_pointer[0]+box_size_pointer[0]; translated_vector[0][1] = original_vector_pointer[1]                    ; translated_vector[0][2] = original_vector_pointer[2];
        translated_vector[1][0] = original_vector_pointer[0]-box_size_pointer[0]; translated_vector[1][1] = original_vector_pointer[1]                    ; translated_vector[1][2] = original_vector_pointer[2];
        translated_vector[2][0] = original_vector_pointer[0]                    ; translated_vector[2][1] = original_vector_pointer[1]+box_size_pointer[1]; translated_vector[2][2] = original_vector_pointer[2];
        translated_vector[3][0] = original_vector_pointer[0]                    ; translated_vector[3][1] = original_vector_pointer[1]-box_size_pointer[1]; translated_vector[3][2] = original_vector_pointer[2];
        translated_vector[4][0] = original_vector_pointer[0]                    ; translated_vector[4][1] = original_vector_pointer[1]                    ; translated_vector[4][2] = original_vector_pointer[2]+box_size_pointer[2];
        translated_vector[5][0] = original_vector_pointer[0]                    ; translated_vector[5][1] = original_vector_pointer[1]                    ; translated_vector[5][2] = original_vector_pointer[2]-box_size_pointer[2];

        temp_distance[0] = compute_vector_length(translated_vector[0][0],translated_vector[0][1],translated_vector[0][2]);
        temp_distance[1] = compute_vector_length(translated_vector[1][0],translated_vector[1][1],translated_vector[1][2]);
        temp_distance[2] = compute_vector_length(translated_vector[2][0],translated_vector[2][1],translated_vector[2][2]);
        temp_distance[3] = compute_vector_length(translated_vector[3][0],translated_vector[3][1],translated_vector[3][2]);
        temp_distance[4] = compute_vector_length(translated_vector[4][0],translated_vector[4][1],translated_vector[4][2]);
        temp_distance[5] = compute_vector_length(translated_vector[5][0],translated_vector[5][1],translated_vector[5][2]);

        temp_distance[7] = temp_distance[6];
        for(int i=0;i<6;i++){
            if(temp_distance[i]<temp_distance[7]){
                temp_distance[7] = temp_distance[i];
                shortest_index = i;
            }
        }

        if(shortest_index!=6){
            for(int i=0;i<3;i++){
                original_vector_pointer[i] = translated_vector[shortest_index][i];
            }
        }
        temp_distance[6] = compute_vector_length(original_vector_pointer[0],original_vector_pointer[1],original_vector_pointer[2]);
    }
}

float compute_arccos_direction_angle(float *vector){

    /*float baseline_vector = {1.0,0.0,0.0};*/
    float angle_cosine = fabs(vector[0])/sqrt(vector[0]*vector[0]+vector[1]*vector[1]);

    if(angle_cosine>1.0){
        angle_cosine = 1.0;
    }

    float angle_size;
    if(vector[0]*vector[1]>0.0){
        angle_size = acos(angle_cosine);
    }else{
        angle_size = pi-acos(angle_cosine);
    }
    
    return(angle_size);
}

float compute_arcsin_direction_angle(float *vector){

    /*float baseline_vector = {1.0,0.0,0.0};*/
    float angle_sine = fabs(vector[1])/sqrt(vector[0]*vector[0]+vector[1]*vector[1]);

    if(angle_sine>1.0){
        angle_sine = 1.0;
    }

    float angle_size;
    if(vector[0]*vector[1]>0.0){
        angle_size = asin(angle_sine);
    }else{
        angle_size = -asin(angle_sine);
    }
    
    return(angle_size);
}

int main(int argc, char *argv[]){

    char *target_lammps_dump_file_location = argv[1];

    int target_first_atom_type;
    sscanf(argv[2], "%d", &target_first_atom_type);
    int target_second_atom_type;
    sscanf(argv[3], "%d", &target_second_atom_type);

    char *output_average_direction_angle_file_location = argv[4];

    printf("\nTarget lammps dump file: %s\nTarget first atom type: %d\nTarget second atom type: %d\nOutput average direction angle file: %s\n", target_lammps_dump_file_location,target_first_atom_type,target_second_atom_type,output_average_direction_angle_file_location);

    FILE *target_lammps_dump_file_pointer;
    FILE *output_average_direction_angle_file_pointer;

    char temp_row[99];

    int frame_count = 0;
    int total_frame_number;

    int temp_atom_id;
    int one_frame_atom_number;
    int temp_atom_type;

    int temp_residue_id;
    int residue_count = 0;
    int one_frame_residue_number;
    int last_residue_id = -1;

    float box_limit[3][2];
    float box_size[3];
    float minimum_box_size;

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
            if(temp_residue_id!=last_residue_id){
                residue_count = residue_count+1;
            }
            last_residue_id = temp_residue_id;
        }
        fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
    }
    fclose(target_lammps_dump_file_pointer);

    one_frame_residue_number = residue_count;

    printf("\nOne frame atom number: %d\nOne frame residue number: %d\n", one_frame_atom_number,one_frame_residue_number);

    last_residue_id = -1;
    residue_count = -1;

    float temp_atom_charge;

    int temp_first_atom_residue_id = -1;
    int temp_second_atom_residue_id = -2;

    float temp_atom_coordinate[3];

    float temp_first_atom_coordinate[3];
    float temp_second_atom_coordinate[3];

    float temp_vector[3];

    float temp_direction_angle = 0.0;
    float (*direction_angle_array)[2] = malloc(sizeof(float)*2*(one_frame_residue_number+1));

    float minimum_direction_angle[3];
    float maximum_direction_angle[3];
    float average_direction_angle[3];
    float direction_angle_range[2];
    float order_parameter[2];

    frame_count = 0;

    target_lammps_dump_file_pointer = fopen(target_lammps_dump_file_location, "r");

    for(int i=0;i<5;i++){
        fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
    }
    for(int i=0;i<3;i++){
        fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
        sscanf(temp_row, "%f %f", &box_limit[i][0],&box_limit[i][1]);
        box_size[i] = box_limit[i][1]-box_limit[i][0];
    }
    minimum_box_size = box_size[0];
    for(int i=1;i<3;i++){
        if(box_size[i]<minimum_box_size){
            minimum_box_size = box_size[i];
        }
    }
    fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);

    output_average_direction_angle_file_pointer = fopen(output_average_direction_angle_file_location, "w");
    fprintf(output_average_direction_angle_file_pointer, "%s\n", "# Frame count, minimum, maximum and average direction angle, %d %.6f %.6f %.6f");

    while(!feof(target_lammps_dump_file_pointer)){

        if(fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row)>0){

            if(sscanf(temp_row, "%f%f%f", &box_limit[0][1],&box_limit[1][1],&box_limit[2][1])==2){

                frame_count = frame_count+1;
                residue_count = -1;

                for(int i=0;i<2;i++){
                    minimum_direction_angle[i] = direction_angle_array[0][i];
                    maximum_direction_angle[i] = direction_angle_array[0][i];
                    average_direction_angle[i] = 0.0;
                    for(int j=0;j<one_frame_residue_number;j++){
                        temp_direction_angle = direction_angle_array[j][i];
                        average_direction_angle[i] = temp_direction_angle+average_direction_angle[i];
                        if(temp_direction_angle<minimum_direction_angle[i]){
                            minimum_direction_angle[i] = temp_direction_angle;
                        }
                        if(temp_direction_angle>maximum_direction_angle[i]){
                            maximum_direction_angle[i] = temp_direction_angle;
                        }
                    }
                    direction_angle_range[i] = maximum_direction_angle[i]-minimum_direction_angle[i];
                    average_direction_angle[i] = average_direction_angle[i]/one_frame_residue_number;
                }

                for(int i=0;i<2;i++){
                    order_parameter[i] = 0.0;
                    for(int j=0;j<one_frame_residue_number;j++){
                        temp_direction_angle = direction_angle_array[j][i]-average_direction_angle[i];
                        /* order_parameter[i] = 1.5*cos(temp_direction_angle)*cos(temp_direction_angle)-0.5+order_parameter[i]; */
                        order_parameter[i] = 1.5*fabs(cos(temp_direction_angle))-0.5+order_parameter[i];
                    }
                    order_parameter[i] = order_parameter[i]/one_frame_residue_number;
                }

                temp_direction_angle = fabs(average_direction_angle[0]-average_direction_angle[1]);

                minimum_direction_angle[2] = minimum_direction_angle[0];
                maximum_direction_angle[2] = maximum_direction_angle[0];
                average_direction_angle[2] = average_direction_angle[0];
                if(fabs(average_direction_angle[1])<0.16*pi){
                    if(order_parameter[1]>order_parameter[0]){
                        minimum_direction_angle[2] = minimum_direction_angle[1];
                        maximum_direction_angle[2] = maximum_direction_angle[1];
                        average_direction_angle[2] = average_direction_angle[1];
                    }
                }

                fprintf(output_average_direction_angle_file_pointer, "%d %.6f %.6f %.6f\n",frame_count,minimum_direction_angle[2],maximum_direction_angle[2],average_direction_angle[2]);

                box_limit[0][0] = box_limit[0][1];
                box_limit[0][1] = box_limit[1][1];
                box_size[0] = box_limit[0][1]-box_limit[0][0];
                for(int i=1;i<3;i++){
                    fscanf(target_lammps_dump_file_pointer, "%[^\n]%*c", temp_row);
                    sscanf(temp_row, "%f %f", &box_limit[i][0],&box_limit[i][1]);
                    box_size[i] = box_limit[i][1]-box_limit[i][0];
                }
                minimum_box_size = box_size[0];
                for(int i=1;i<3;i++){
                    if(box_size[i]<minimum_box_size){
                        minimum_box_size = box_size[i];
                    }
                }

            }else if(sscanf(temp_row, "%d%d%d%f%f%f%f\n", &temp_atom_id,&temp_residue_id,&temp_atom_type,&temp_atom_charge,&temp_atom_coordinate[0],&temp_atom_coordinate[1],&temp_atom_coordinate[2])==7){

                if(temp_residue_id!=last_residue_id){
                    residue_count = residue_count+1;
                }
                last_residue_id = temp_residue_id;

                if(temp_atom_type==target_first_atom_type){
                    temp_first_atom_residue_id = temp_residue_id;
                    for(int i=0;i<3;i++){
                        temp_first_atom_coordinate[i] = temp_atom_coordinate[i];
                    }
                    if(temp_first_atom_residue_id==temp_second_atom_residue_id){
                        for(int i=0;i<3;i++){
                            temp_vector[i] = temp_first_atom_coordinate[i]-temp_second_atom_coordinate[i];
                        }
                        periodic_vector_processing(&temp_vector[0],minimum_box_size,&box_size[0]);
                        temp_direction_angle = compute_arccos_direction_angle(&temp_vector[0]);
                        direction_angle_array[residue_count][0] = temp_direction_angle;
                        temp_direction_angle = compute_arcsin_direction_angle(&temp_vector[0]);
                        direction_angle_array[residue_count][1] = temp_direction_angle;
                    }
                }else if(temp_atom_type==target_second_atom_type){
                    temp_second_atom_residue_id = temp_residue_id;
                    for(int i=0;i<3;i++){
                        temp_second_atom_coordinate[i] = temp_atom_coordinate[i];
                    }
                    if(temp_first_atom_residue_id==temp_second_atom_residue_id){
                        for(int i=0;i<3;i++){
                            temp_vector[i] = temp_first_atom_coordinate[i]-temp_second_atom_coordinate[i];
                        }
                        periodic_vector_processing(&temp_vector[0],minimum_box_size,&box_size[0]);
                        temp_direction_angle = compute_arccos_direction_angle(&temp_vector[0]);
                        direction_angle_array[residue_count][0] = temp_direction_angle;
                        temp_direction_angle = compute_arcsin_direction_angle(&temp_vector[0]);
                        direction_angle_array[residue_count][1] = temp_direction_angle;
                    }
                }
            }
        }else{
            fscanf(target_lammps_dump_file_pointer, "%c", temp_row);
        }
    }

    fclose(target_lammps_dump_file_pointer);

    frame_count = frame_count+1;

    for(int i=0;i<2;i++){
        minimum_direction_angle[i] = direction_angle_array[0][i];
        maximum_direction_angle[i] = direction_angle_array[0][i];
        average_direction_angle[i] = 0.0;
        for(int j=0;j<one_frame_residue_number;j++){
            temp_direction_angle = direction_angle_array[j][i];
            average_direction_angle[i] = temp_direction_angle+average_direction_angle[i];
            if(temp_direction_angle<minimum_direction_angle[i]){
                minimum_direction_angle[i] = temp_direction_angle;
            }
            if(temp_direction_angle>maximum_direction_angle[i]){
                maximum_direction_angle[i] = temp_direction_angle;
            }
        }
        direction_angle_range[i] = maximum_direction_angle[i]-minimum_direction_angle[i];
        average_direction_angle[i] = average_direction_angle[i]/one_frame_residue_number;
    }

    for(int i=0;i<2;i++){
        order_parameter[i] = 0.0;
        for(int j=0;j<one_frame_residue_number;j++){
            temp_direction_angle = direction_angle_array[j][i]-average_direction_angle[i];
            /* order_parameter[i] = 1.5*cos(temp_direction_angle)*cos(temp_direction_angle)-0.5+order_parameter[i]; */
            order_parameter[i] = 1.5*fabs(cos(temp_direction_angle))-0.5+order_parameter[i];
        }
        order_parameter[i] = order_parameter[i]/one_frame_residue_number;
    }

    temp_direction_angle = fabs(average_direction_angle[0]-average_direction_angle[1]);

    minimum_direction_angle[2] = minimum_direction_angle[0];
    maximum_direction_angle[2] = maximum_direction_angle[0];
    average_direction_angle[2] = average_direction_angle[0];
    if(fabs(average_direction_angle[1])<0.16*pi){
        if(order_parameter[1]>order_parameter[0]){
            minimum_direction_angle[2] = minimum_direction_angle[1];
            maximum_direction_angle[2] = maximum_direction_angle[1];
            average_direction_angle[2] = average_direction_angle[1];
        }
    }

    fprintf(output_average_direction_angle_file_pointer, "%d %.6f %.6f %.6f\n",frame_count,minimum_direction_angle[2],maximum_direction_angle[2],average_direction_angle[2]);

    fclose(output_average_direction_angle_file_pointer);

    total_frame_number = frame_count;

    printf("Total frame number: %d\n", total_frame_number);

    free(direction_angle_array);

    return(0);
}