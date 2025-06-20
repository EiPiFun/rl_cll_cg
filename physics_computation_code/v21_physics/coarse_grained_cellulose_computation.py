import numpy
import os
import subprocess
import sys

from coarse_grained_cellulose_post_processing import CoarseGrainedCellulosePostProcessing

lammps_executive_name = 'lmp'
lammps_backend_switch = 'omp'
coarse_grained_cellulose_computation_timeout = 600.0

from coarse_grained_cellulose_post_processing import default_fitted_elastic_modulus

from coarse_grained_cellulose_post_processing import default_elastic_modulus_match_degree

from coarse_grained_cellulose_post_processing import default_average_bond_length
from coarse_grained_cellulose_post_processing import default_average_persistence_length
from coarse_grained_cellulose_post_processing import default_average_end_to_end_distance

from coarse_grained_cellulose_post_processing import default_persistence_length_match_degree
from coarse_grained_cellulose_post_processing import default_end_to_end_distance_match_degree

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

        self.physics_computation_1_in_file_location = self.physics_computation_1_inputs_directory+'coarse_grained_cellulose.in'
        self.physics_computation_2_in_file_location = self.physics_computation_2_inputs_directory+'coarse_grained_cellulose.in'
        self.physics_computation_3_in_file_location = self.physics_computation_3_inputs_directory+'coarse_grained_cellulose.in'
        self.physics_computation_4_in_file_location = self.physics_computation_4_inputs_directory+'coarse_grained_cellulose.in'

        self.lammps_physics_computation_1_command = (lammps_executive_name,'-in',self.physics_computation_1_in_file_location)
        self.lammps_physics_computation_2_command = (lammps_executive_name,'-in',self.physics_computation_2_in_file_location)
        self.lammps_physics_computation_3_command = (lammps_executive_name,'-in',self.physics_computation_3_in_file_location)
        self.lammps_physics_computation_4_command = (lammps_executive_name,'-in',self.physics_computation_4_in_file_location)

        self.post_processing_results_location = working_directory+'post_processing_results.txt'

    def generate_in_file(self,bc_1,bc_2,bc_3,ac_1,ac_2,ac_3,ac_4,ac_5,ac_6,ic_1,ic_2,ic_3):

        inputs_directories = (self.physics_computation_1_inputs_directory,self.physics_computation_2_inputs_directory,self.physics_computation_3_inputs_directory,self.physics_computation_4_inputs_directory)
        in_file_locations = (self.physics_computation_1_in_file_location,self.physics_computation_2_in_file_location,self.physics_computation_3_in_file_location,self.physics_computation_4_in_file_location)

        for i in (0,1,2,3):

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
\n'.format(bc_1,bc_2,bc_3,ac_1,ac_2,ac_3,ac_4,ac_5,ac_6,ic_1,ic_2,ic_3)

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

    def post_processing(self):

        post_processing_results_array = CoarseGrainedCellulosePostProcessing(self.working_directory).post_processing_bundle_v1()

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

    for temp_directory in (working_directory,physics_computation_1_inputs_directory,physics_computation_1_results_directory,physics_computation_2_inputs_directory,physics_computation_2_results_directory,physics_computation_3_inputs_directory,physics_computation_3_results_directory,physics_computation_4_inputs_directory,physics_computation_4_results_directory):
        pathlib.Path(temp_directory).mkdir(parents=True,exist_ok=True)

    for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data'):
        shutil.copyfile(physics_computation_1_code_inputs_directory+file_name,physics_computation_1_inputs_directory+file_name)

    for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data','compute_persistence_length_from_lammps_dump','compute_end_to_end_distance_from_lammps_dump'):
        shutil.copyfile(physics_computation_2_code_inputs_directory+file_name,physics_computation_2_inputs_directory+file_name)

    for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data','compute_persistence_length_from_lammps_dump','compute_end_to_end_distance_from_lammps_dump'):
        shutil.copyfile(physics_computation_3_code_inputs_directory+file_name,physics_computation_3_inputs_directory+file_name)

    for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data','compute_persistence_length_from_lammps_dump','compute_end_to_end_distance_from_lammps_dump'):
        shutil.copyfile(physics_computation_4_code_inputs_directory+file_name,physics_computation_4_inputs_directory+file_name)

    coarse_grained_cellulose_computation_class = CoarseGrainedCelluloseComputation(working_directory)
    coarse_grained_cellulose_computation_class.generate_in_file(60.00,45.01,97.07,25.00,19.44,28.38,10.90,12.96,12.06,0.16,0.21,0.18)

    coarse_grained_cellulose_computation_1_queue = multiprocessing.Queue()
    coarse_grained_cellulose_computation_2_queue = multiprocessing.Queue()
    coarse_grained_cellulose_computation_3_queue = multiprocessing.Queue()
    coarse_grained_cellulose_computation_4_queue = multiprocessing.Queue()
    coarse_grained_cellulose_computation_queue_array = (coarse_grained_cellulose_computation_1_queue,coarse_grained_cellulose_computation_2_queue,coarse_grained_cellulose_computation_3_queue,coarse_grained_cellulose_computation_4_queue)
    for queue in coarse_grained_cellulose_computation_queue_array:
        queue.put(False)

    coarse_grained_cellulose_computation_1_process = multiprocessing.Process(target=coarse_grained_cellulose_computation_class.physics_computation_1,daemon=True,args=(coarse_grained_cellulose_computation_1_queue,))
    coarse_grained_cellulose_computation_2_process = multiprocessing.Process(target=coarse_grained_cellulose_computation_class.physics_computation_2,daemon=True,args=(coarse_grained_cellulose_computation_2_queue,))
    coarse_grained_cellulose_computation_3_process = multiprocessing.Process(target=coarse_grained_cellulose_computation_class.physics_computation_3,daemon=True,args=(coarse_grained_cellulose_computation_3_queue,))
    coarse_grained_cellulose_computation_4_process = multiprocessing.Process(target=coarse_grained_cellulose_computation_class.physics_computation_4,daemon=True,args=(coarse_grained_cellulose_computation_4_queue,))

    coarse_grained_cellulose_computation_4_process.start()
    
    coarse_grained_cellulose_computation_process_array = (coarse_grained_cellulose_computation_1_process,coarse_grained_cellulose_computation_2_process,coarse_grained_cellulose_computation_3_process)
    for process in coarse_grained_cellulose_computation_process_array:
        process.start()
        process.join(timeout=coarse_grained_cellulose_computation_timeout)

    coarse_grained_cellulose_computation_4_process.join()

    coarse_grained_cellulose_computation_return_array = [False,False,False,False]
    for i in (0,1,2,3):
        coarse_grained_cellulose_computation_return_array[i] = coarse_grained_cellulose_computation_queue_array[i].get()
    for queue in coarse_grained_cellulose_computation_queue_array:
        queue.close()
    for process in coarse_grained_cellulose_computation_process_array:
        if process.is_alive() == True:
            process.kill()
        process.join()
        process.close()

    coarse_grained_cellulose_computation_return = True
    for i in (0,1,2,3):
        if coarse_grained_cellulose_computation_return_array[i] == False:
            coarse_grained_cellulose_computation_return = False

    if coarse_grained_cellulose_computation_return == True:
        fitted_elastic_modulus,elastic_modulus_match_degree,average_bond_length,average_persistence_length,persistence_length_match_degree,average_end_to_end_distance,end_to_end_distance_match_degree = coarse_grained_cellulose_computation_class.post_processing()
    else:
        fitted_elastic_modulus = default_fitted_elastic_modulus
        elastic_modulus_match_degree = default_elastic_modulus_match_degree
        average_bond_length = default_average_bond_length
        average_persistence_length = default_average_persistence_length
        persistence_length_match_degree = default_persistence_length_match_degree
        average_end_to_end_distance = default_average_end_to_end_distance
        end_to_end_distance_match_degree = default_end_to_end_distance_match_degree

    print(fitted_elastic_modulus)
    print(elastic_modulus_match_degree)
    print(average_bond_length)
    print(average_persistence_length)
    print(persistence_length_match_degree)
    print(average_end_to_end_distance)
    print(end_to_end_distance_match_degree)

