import numpy
import os
import subprocess
import sys

from coarse_grained_cellulose_post_processing import CoarseGrainedCellulosePostProcessing

lammps_executive_name = 'lmp'
lammps_backend_switch = 'omp'
coarse_grained_cellulose_computation_timeout = 600.0

relative_hb_r_in = 1.00
relative_hb_r_cut = 1.30
relative_lj_r_in = 1.60
relative_lj_r_cut = 1.70

from coarse_grained_cellulose_post_processing import default_fitted_elastic_modulus

from coarse_grained_cellulose_post_processing import default_elastic_modulus_match_degree

from coarse_grained_cellulose_post_processing import default_average_bond_length
from coarse_grained_cellulose_post_processing import default_average_persistence_length
from coarse_grained_cellulose_post_processing import default_average_end_to_end_distance

from coarse_grained_cellulose_post_processing import default_persistence_length_match_degree
from coarse_grained_cellulose_post_processing import default_end_to_end_distance_match_degree

from coarse_grained_cellulose_post_processing import default_system_stable_enough_switch

from coarse_grained_cellulose_post_processing import default_fracture_average_direction_angle_match_degree
from coarse_grained_cellulose_post_processing import default_fracture_rotation_angle_match_degree

from coarse_grained_cellulose_post_processing import default_max_stress_match_degree
from coarse_grained_cellulose_post_processing import default_max_strain_energy_match_degree
from coarse_grained_cellulose_post_processing import default_equivalent_strain_match_degree

class CoarseGrainedCelluloseComputation():

    def __init__(self,working_directory):

        self.working_directory = working_directory

        self.physics_computation_1_inputs_directory = working_directory+'physics_computation_1/inputs/'
        self.physics_computation_1_results_directory = working_directory+'physics_computation_1/results/'
        self.physics_computation_2_inputs_directory = working_directory+'physics_computation_2/inputs/'
        self.physics_computation_2_results_directory = working_directory+'physics_computation_2/results/'
        self.physics_computation_3_inputs_directory = working_directory+'physics_computation_3/inputs/'
        self.physics_computation_3_results_directory = working_directory+'physics_computation_3/results/'
        self.physics_computation_4_inputs_directory = working_directory+'physics_computation_4/inputs/'
        self.physics_computation_4_results_directory = working_directory+'physics_computation_4/results/'
        self.physics_computation_5_inputs_directory = working_directory+'physics_computation_5/inputs/'
        self.physics_computation_5_results_directory = working_directory+'physics_computation_5/results/'
        self.physics_computation_6_inputs_directory = working_directory+'physics_computation_6/inputs/'
        self.physics_computation_6_results_directory = working_directory+'physics_computation_6/results/'
        self.physics_computation_7_inputs_directory = working_directory+'physics_computation_7/inputs/'
        self.physics_computation_7_results_directory = working_directory+'physics_computation_7/results/'

        self.physics_computation_1_in_file_location = self.physics_computation_1_inputs_directory+'coarse_grained_cellulose.in'
        self.physics_computation_2_in_file_location = self.physics_computation_2_inputs_directory+'coarse_grained_cellulose.in'
        self.physics_computation_3_in_file_location = self.physics_computation_3_inputs_directory+'coarse_grained_cellulose.in'
        self.physics_computation_4_in_file_location = self.physics_computation_4_inputs_directory+'coarse_grained_cellulose.in'
        self.physics_computation_5_in_file_location = self.physics_computation_5_inputs_directory+'coarse_grained_cellulose.in'
        self.physics_computation_6_in_file_location = self.physics_computation_6_inputs_directory+'coarse_grained_cellulose.in'
        self.physics_computation_7_in_file_location = self.physics_computation_7_inputs_directory+'coarse_grained_cellulose.in'

        self.lammps_physics_computation_1_command = (lammps_executive_name,'-in',self.physics_computation_1_in_file_location)
        self.lammps_physics_computation_2_command = (lammps_executive_name,'-in',self.physics_computation_2_in_file_location)
        self.lammps_physics_computation_3_command = (lammps_executive_name,'-in',self.physics_computation_3_in_file_location)
        self.lammps_physics_computation_4_command = (lammps_executive_name,'-in',self.physics_computation_4_in_file_location)
        self.lammps_physics_computation_5_command = (lammps_executive_name,'-in',self.physics_computation_5_in_file_location)
        self.lammps_physics_computation_6_command = (lammps_executive_name,'-in',self.physics_computation_6_in_file_location)
        self.lammps_physics_computation_7_command = (lammps_executive_name,'-in',self.physics_computation_7_in_file_location)

        self.post_processing_results_location = working_directory+'post_processing_results.txt'

    def generate_in_file(self,bc_1,bc_2,bc_3,ac_1,ac_2,ac_3,ac_4,ac_5,ac_6,ic_1,ic_2,ic_3,pc_1_2_1_hb,pc_1_2_2_hb,pc_1_3_1_hb,pc_1_3_2_hb,pc_1_1_1,pc_1_1_2,pc_1_2_1,pc_1_2_2,pc_1_3_1,pc_1_3_2,pc_2_2_1,pc_2_2_2,pc_2_3_1,pc_2_3_2,pc_3_3_1,pc_3_3_2):

        inputs_directories = (self.physics_computation_1_inputs_directory,self.physics_computation_2_inputs_directory,self.physics_computation_3_inputs_directory,self.physics_computation_4_inputs_directory,self.physics_computation_5_inputs_directory,self.physics_computation_6_inputs_directory,self.physics_computation_7_inputs_directory)
        in_file_locations = (self.physics_computation_1_in_file_location,self.physics_computation_2_in_file_location,self.physics_computation_3_in_file_location,self.physics_computation_4_in_file_location,self.physics_computation_5_in_file_location,self.physics_computation_6_in_file_location,self.physics_computation_7_in_file_location)

        for i in (0,1,2,3,4,5,6):

            constant_head = open(inputs_directories[i]+'in_file_constant_head.txt').read()
            constant_tail = open(inputs_directories[i]+'in_file_constant_tail.txt').read()

            coeff_information = '\n\
bond_coeff 1 {0} 5.242 # CL1-CL1\n\
bond_coeff 2 {1} 1.955 # CL2-CL1\n\
bond_coeff 3 {2} 2.250 # CL3-CL1\n\
\n\
angle_coeff 1 {3} 166.2 # CL1-CL1-CL1\n\
angle_coeff 2 {4}  77.0 # CL2-CL1-CL1_1\n\
angle_coeff 3 {5}  71.4 # CL3-CL1-CL1_1\n\
angle_coeff 4 {6} 112.5 # CL2-CL1-CL1_2\n\
angle_coeff 5 {7}  97.3 # CL3-CL1-CL1_2\n\
angle_coeff 6 {8} 147.5 # CL2-CL1-CL3\n\
\n\
improper_coeff 1 {9} 176.7 # CL1-CL1-CL1-CL1\n\
improper_coeff 2 {10}   7.8 # CL2-CL1-CL1-CL3\n\
improper_coeff 3 {11}   9.1 # CL3-CL1-CL1-CL2\n\
\n\
pair_coeff 1 2 hbond/dreiding/lj 3 i {12} {13} 4 {14} {15} 120.0 # CL1-CL2\n\
pair_coeff 1 3 hbond/dreiding/lj 2 i {16} {17} 4 {18} {19} 120.0 # CL1-CL3\n\
pair_coeff 1 1 lj/smooth {20} {21} {22} {23}                     # CL1-CL1\n\
pair_coeff 1 2 lj/smooth {24} {25} {26} {27}                     # CL1-CL2\n\
pair_coeff 1 3 lj/smooth {28} {29} {30} {31}                     # CL1-CL3\n\
pair_coeff 2 2 lj/smooth {32} {33} {34} {35}                     # CL2-CL2\n\
pair_coeff 2 3 lj/smooth {36} {37} {38} {39}                     # CL2-CL3\n\
pair_coeff 3 3 lj/smooth {40} {41} {42} {43}                     # CL3-CL3\n\
\n'.format(bc_1,bc_2,bc_3,ac_1,ac_2,ac_3,ac_4,ac_5,ac_6,ic_1,ic_2,ic_3,pc_1_2_1_hb,pc_1_2_2_hb,relative_hb_r_in*pc_1_2_2_hb,relative_hb_r_cut*pc_1_2_2_hb,pc_1_3_1_hb,pc_1_3_2_hb,relative_hb_r_in*pc_1_3_2_hb,relative_hb_r_cut*pc_1_3_2_hb,pc_1_1_1,pc_1_1_2,relative_lj_r_in*pc_1_1_2,relative_lj_r_cut*pc_1_1_2,pc_1_2_1,pc_1_2_2,relative_lj_r_in*pc_1_2_2,relative_lj_r_cut*pc_1_2_2,pc_1_3_1,pc_1_3_2,relative_lj_r_in*pc_1_3_2,relative_lj_r_cut*pc_1_3_2,pc_2_2_1,pc_2_2_2,relative_lj_r_in*pc_2_2_2,relative_lj_r_cut*pc_2_2_2,pc_2_3_1,pc_2_3_2,relative_lj_r_in*pc_2_3_2,relative_lj_r_cut*pc_2_3_2,pc_3_3_1,pc_3_3_2,relative_lj_r_in*pc_3_3_2,relative_lj_r_cut*pc_3_3_2)

            in_file_content = constant_head+coeff_information+constant_tail
            in_file = open(in_file_locations[i],'w')
            in_file.write(in_file_content)
            in_file.close()

    def physics_computation_1(self,queue):

        try:
            subprocess.run(self.lammps_physics_computation_1_command,cwd=self.physics_computation_1_results_directory,stdout=subprocess.DEVNULL,check=True)
            queue.get()
            queue.put(True)
            return True
        except:
            queue.get()
            queue.put(False)
            return False

    def physics_computation_2(self,queue):

        try:
            subprocess.run(self.lammps_physics_computation_2_command,cwd=self.physics_computation_2_results_directory,stdout=subprocess.DEVNULL,check=True)
            queue.get()
            queue.put(True)
            return True
        except:
            queue.get()
            queue.put(False)
            return False

    def physics_computation_3(self,queue):

        try:
            subprocess.run(self.lammps_physics_computation_3_command,cwd=self.physics_computation_3_results_directory,stdout=subprocess.DEVNULL,check=True)
            queue.get()
            queue.put(True)
            return True
        except:
            queue.get()
            queue.put(False)
            return False

    def physics_computation_4(self,queue):

        try:
            subprocess.run(self.lammps_physics_computation_4_command,cwd=self.physics_computation_4_results_directory,stdout=subprocess.DEVNULL,check=True)
            queue.get()
            queue.put(True)
            return True
        except:
            queue.get()
            queue.put(False)
            return False

    def physics_computation_5(self,queue):

        try:
            subprocess.run(self.lammps_physics_computation_5_command,cwd=self.physics_computation_5_results_directory,stdout=subprocess.DEVNULL,check=True)
            queue.get()
            queue.put(True)
            return True
        except:
            queue.get()
            queue.put(False)
            return False

    def physics_computation_6(self,queue):

        try:
            subprocess.run(self.lammps_physics_computation_6_command,cwd=self.physics_computation_6_results_directory,stdout=subprocess.DEVNULL,check=True)
            queue.get()
            queue.put(True)
            return True
        except:
            queue.get()
            queue.put(False)
            return False

    def physics_computation_7(self,queue):

        try:
            subprocess.run(self.lammps_physics_computation_7_command,cwd=self.physics_computation_7_results_directory,stdout=subprocess.DEVNULL,check=True)
            queue.get()
            queue.put(True)
            return True
        except:
            queue.get()
            queue.put(False)
            return False

    def post_processing(self):

        post_processing_results_array = CoarseGrainedCellulosePostProcessing(self.working_directory).post_processing_bundle_v3()

        numpy.savetxt(self.post_processing_results_location,post_processing_results_array,fmt='%.6e',delimiter=',')

        return post_processing_results_array

if __name__=='__main__':

    import multiprocessing

    import pathlib
    import time
    import shutil

    current_user_main_directory = str(pathlib.Path.home())+'/'
    try:
        coarse_grained_cellulose_computation_home_directory = os.getenv('CLL_CG_HOME')+'/'
    except:
        coarse_grained_cellulose_computation_home_directory = current_user_main_directory

    lammps_results_overwrite_switch = False

    code_directory = sys.path[0]+'/'
    physics_computation_1_code_inputs_directory = code_directory+'physics_computation_1/inputs/'
    physics_computation_2_code_inputs_directory = code_directory+'physics_computation_2/inputs/'
    physics_computation_3_code_inputs_directory = code_directory+'physics_computation_3/inputs/'
    physics_computation_4_code_inputs_directory = code_directory+'physics_computation_4/inputs/'
    physics_computation_5_code_inputs_directory = code_directory+'physics_computation_5/inputs/'
    physics_computation_6_code_inputs_directory = code_directory+'physics_computation_6/inputs/'
    physics_computation_7_code_inputs_directory = code_directory+'physics_computation_7/inputs/'

    run_mode = 'normal'

    if len(sys.argv) > 1:
        run_mode = sys.argv[1]

    if run_mode == 'debug':
        working_directory = sys.path[0]+'/debug/'
    else:
        working_directory_base = coarse_grained_cellulose_computation_home_directory+'coarse_grained_cellulose_parameters_data/physics_computations/manual/coarse_grained_cellulose_computations/coarse_grained_cellulose_computation'
        if lammps_results_overwrite_switch == True:
            working_directory = working_directory_base+'/'
        else:
            current_time = time.strftime('%Y-%m-%d-%H:%M:%S')
            working_directory = working_directory_base+'-'+current_time+'/'

    try:
        pathlib.Path(working_directory).mkdir(parents=True,exist_ok=True)
        shutil.rmtree(working_directory)
    except:
        current_time = time.strftime('%Y-%m-%d-%H:%M:%S')
        working_directory = working_directory+current_time+'/'

    physics_computation_1_inputs_directory = working_directory+'physics_computation_1/inputs/'
    physics_computation_1_results_directory = working_directory+'physics_computation_1/results/'
    physics_computation_2_inputs_directory = working_directory+'physics_computation_2/inputs/'
    physics_computation_2_results_directory = working_directory+'physics_computation_2/results/'
    physics_computation_3_inputs_directory = working_directory+'physics_computation_3/inputs/'
    physics_computation_3_results_directory = working_directory+'physics_computation_3/results/'
    physics_computation_4_inputs_directory = working_directory+'physics_computation_4/inputs/'
    physics_computation_4_results_directory = working_directory+'physics_computation_4/results/'
    physics_computation_5_inputs_directory = working_directory+'physics_computation_5/inputs/'
    physics_computation_5_results_directory = working_directory+'physics_computation_5/results/'
    physics_computation_6_inputs_directory = working_directory+'physics_computation_6/inputs/'
    physics_computation_6_results_directory = working_directory+'physics_computation_6/results/'
    physics_computation_7_inputs_directory = working_directory+'physics_computation_7/inputs/'
    physics_computation_7_results_directory = working_directory+'physics_computation_7/results/'

    for temp_directory in (working_directory,physics_computation_1_inputs_directory,physics_computation_1_results_directory,physics_computation_2_inputs_directory,physics_computation_2_results_directory,physics_computation_3_inputs_directory,physics_computation_3_results_directory,physics_computation_4_inputs_directory,physics_computation_4_results_directory,physics_computation_5_inputs_directory,physics_computation_5_results_directory,physics_computation_6_inputs_directory,physics_computation_6_results_directory,physics_computation_7_inputs_directory,physics_computation_7_results_directory):
        pathlib.Path(temp_directory).mkdir(parents=True,exist_ok=True)

    for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data'):
        shutil.copyfile(physics_computation_1_code_inputs_directory+file_name,physics_computation_1_inputs_directory+file_name)

    for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data','compute_persistence_length_from_lammps_dump','compute_end_to_end_distance_from_lammps_dump'):
        shutil.copyfile(physics_computation_2_code_inputs_directory+file_name,physics_computation_2_inputs_directory+file_name)

    for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data','compute_persistence_length_from_lammps_dump','compute_end_to_end_distance_from_lammps_dump'):
        shutil.copyfile(physics_computation_3_code_inputs_directory+file_name,physics_computation_3_inputs_directory+file_name)

    for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data','compute_persistence_length_from_lammps_dump','compute_end_to_end_distance_from_lammps_dump'):
        shutil.copyfile(physics_computation_4_code_inputs_directory+file_name,physics_computation_4_inputs_directory+file_name)

    for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data','compute_average_direction_angle_from_lammps_dump'):
        shutil.copyfile(physics_computation_5_code_inputs_directory+file_name,physics_computation_5_inputs_directory+file_name)

    for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data','compute_average_direction_angle_from_lammps_dump'):
        shutil.copyfile(physics_computation_6_code_inputs_directory+file_name,physics_computation_6_inputs_directory+file_name)

    for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data','compute_average_direction_angle_from_lammps_dump'):
        shutil.copyfile(physics_computation_7_code_inputs_directory+file_name,physics_computation_7_inputs_directory+file_name)

    coarse_grained_cellulose_computation_class = CoarseGrainedCelluloseComputation(working_directory)
    coarse_grained_cellulose_computation_class.generate_in_file(60.00,45.01,97.07,25.00,19.44,28.38,10.90,12.96,12.06,0.16,0.21,0.18,1.4,6.4,1.4,6.0,0.8,5.1,0.8,4.4,0.8,3.9,0.8,5.1,0.8,3.8,0.8,5.1)

    coarse_grained_cellulose_computation_1_queue = multiprocessing.Queue()
    coarse_grained_cellulose_computation_2_queue = multiprocessing.Queue()
    coarse_grained_cellulose_computation_3_queue = multiprocessing.Queue()
    coarse_grained_cellulose_computation_4_queue = multiprocessing.Queue()
    coarse_grained_cellulose_computation_5_queue = multiprocessing.Queue()
    coarse_grained_cellulose_computation_6_queue = multiprocessing.Queue()
    coarse_grained_cellulose_computation_7_queue = multiprocessing.Queue()
    coarse_grained_cellulose_computation_queue_array = (coarse_grained_cellulose_computation_1_queue,coarse_grained_cellulose_computation_2_queue,coarse_grained_cellulose_computation_3_queue,coarse_grained_cellulose_computation_4_queue,coarse_grained_cellulose_computation_5_queue,coarse_grained_cellulose_computation_6_queue,coarse_grained_cellulose_computation_7_queue)
    for queue in coarse_grained_cellulose_computation_queue_array:
        queue.put(False)

    coarse_grained_cellulose_computation_1_process = multiprocessing.Process(target=coarse_grained_cellulose_computation_class.physics_computation_1,daemon=True,args=(coarse_grained_cellulose_computation_1_queue,))
    coarse_grained_cellulose_computation_2_process = multiprocessing.Process(target=coarse_grained_cellulose_computation_class.physics_computation_2,daemon=True,args=(coarse_grained_cellulose_computation_2_queue,))
    coarse_grained_cellulose_computation_3_process = multiprocessing.Process(target=coarse_grained_cellulose_computation_class.physics_computation_3,daemon=True,args=(coarse_grained_cellulose_computation_3_queue,))
    coarse_grained_cellulose_computation_4_process = multiprocessing.Process(target=coarse_grained_cellulose_computation_class.physics_computation_4,daemon=True,args=(coarse_grained_cellulose_computation_4_queue,))
    coarse_grained_cellulose_computation_5_process = multiprocessing.Process(target=coarse_grained_cellulose_computation_class.physics_computation_5,daemon=True,args=(coarse_grained_cellulose_computation_5_queue,))
    coarse_grained_cellulose_computation_6_process = multiprocessing.Process(target=coarse_grained_cellulose_computation_class.physics_computation_6,daemon=True,args=(coarse_grained_cellulose_computation_6_queue,))
    coarse_grained_cellulose_computation_7_process = multiprocessing.Process(target=coarse_grained_cellulose_computation_class.physics_computation_7,daemon=True,args=(coarse_grained_cellulose_computation_7_queue,))

    coarse_grained_cellulose_computation_7_process.start()

    coarse_grained_cellulose_computation_process_array = (coarse_grained_cellulose_computation_1_process,coarse_grained_cellulose_computation_2_process,coarse_grained_cellulose_computation_3_process,coarse_grained_cellulose_computation_4_process)
    for process in coarse_grained_cellulose_computation_process_array:
        process.start()
    for process in coarse_grained_cellulose_computation_process_array:
        process.join(timeout=coarse_grained_cellulose_computation_timeout)

    coarse_grained_cellulose_computation_process_array = (coarse_grained_cellulose_computation_5_process,coarse_grained_cellulose_computation_6_process)
    for process in coarse_grained_cellulose_computation_process_array:
        process.start()
    for process in coarse_grained_cellulose_computation_process_array:
        process.join(timeout=coarse_grained_cellulose_computation_timeout)
    
    coarse_grained_cellulose_computation_7_process.join(timeout=coarse_grained_cellulose_computation_timeout)

    coarse_grained_cellulose_computation_return_array = [False,False,False,False,False,False,False]
    for i in (0,1,2,3,4,5,6):
        coarse_grained_cellulose_computation_return_array[i] = coarse_grained_cellulose_computation_queue_array[i].get()
    for queue in coarse_grained_cellulose_computation_queue_array:
        queue.close()
    for process in coarse_grained_cellulose_computation_process_array:
        if process.is_alive() == True:
            process.kill()
        process.join()
        process.close()

    coarse_grained_cellulose_computation_return = True
    for i in (0,1,2,3,4,5,6):
        if coarse_grained_cellulose_computation_return_array[i] == False:
            coarse_grained_cellulose_computation_return = False

    if coarse_grained_cellulose_computation_return == True:
        fitted_elastic_modulus,elastic_modulus_match_degree,average_persistence_length,persistence_length_match_degree,average_end_to_end_distance,end_to_end_distance_match_degree,system_stable_enough_switch,fracture_average_direction_angle_match_degree,fracture_rotation_angle_match_degree,max_stress_match_degree,max_strain_energy_match_degree,equivalent_strain_match_degree = coarse_grained_cellulose_computation_class.post_processing()
    else:
        fitted_elastic_modulus = default_fitted_elastic_modulus
        elastic_modulus_match_degree = default_elastic_modulus_match_degree
        average_persistence_length = default_average_persistence_length
        persistence_length_match_degree = default_persistence_length_match_degree
        average_end_to_end_distance = default_average_end_to_end_distance
        end_to_end_distance_match_degree = default_end_to_end_distance_match_degree
        system_stable_enough_switch = default_system_stable_enough_switch
        fracture_average_direction_angle_match_degree = default_fracture_average_direction_angle_match_degree
        fracture_rotation_angle_match_degree = default_fracture_rotation_angle_match_degree
        max_stress_match_degree = default_max_stress_match_degree
        max_strain_energy_match_degree = default_max_strain_energy_match_degree
        equivalent_strain_match_degree = default_equivalent_strain_match_degree

    print(fitted_elastic_modulus)
    print(elastic_modulus_match_degree)
    print(average_persistence_length)
    print(persistence_length_match_degree)
    print(average_end_to_end_distance)
    print(end_to_end_distance_match_degree)
    print(system_stable_enough_switch)
    print(fracture_average_direction_angle_match_degree)
    print(fracture_rotation_angle_match_degree)
    print(max_stress_match_degree)
    print(max_strain_energy_match_degree)
    print(equivalent_strain_match_degree)

