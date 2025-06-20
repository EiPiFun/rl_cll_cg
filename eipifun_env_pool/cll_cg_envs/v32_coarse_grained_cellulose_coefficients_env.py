import gymnasium as gym
import multiprocessing
import numpy
import os
import pathlib
import shutil
import sys
import time

code_version = 'v31'
physics_version = 'v32'

coefficient_precision_format = '%.4g'

current_user_main_directory = str(pathlib.Path.home())+'/'
try:
    coarse_grained_cellulose_computation_home_directory = os.getenv('CLL_CG_HOME')+'/'
except:
    coarse_grained_cellulose_computation_home_directory = current_user_main_directory

coarse_grained_cellulose_computation_code_directory = coarse_grained_cellulose_computation_home_directory+'coarse_grained_cellulose_coefficients_code/physics_computation_code/'+code_version+'_physics/'
coarse_grained_cellulose_physics_computation_1_code_inputs_directory = coarse_grained_cellulose_computation_code_directory+'physics_computation_1/inputs/'
coarse_grained_cellulose_physics_computation_2_code_inputs_directory = coarse_grained_cellulose_computation_code_directory+'physics_computation_2/inputs/'
coarse_grained_cellulose_physics_computation_3_code_inputs_directory = coarse_grained_cellulose_computation_code_directory+'physics_computation_3/inputs/'
coarse_grained_cellulose_physics_computation_4_code_inputs_directory = coarse_grained_cellulose_computation_code_directory+'physics_computation_4/inputs/'

sys.path.append(coarse_grained_cellulose_computation_code_directory)

from coarse_grained_cellulose_computation import CoarseGrainedCelluloseComputation
from coarse_grained_cellulose_computation import lammps_executive_name
from coarse_grained_cellulose_computation import lammps_backend_switch
from coarse_grained_cellulose_computation import coarse_grained_cellulose_computation_timeout

from coarse_grained_cellulose_post_processing import default_fitted_elastic_modulus
from coarse_grained_cellulose_post_processing import default_elastic_modulus_match_degree
from coarse_grained_cellulose_post_processing import default_average_bond_length
from coarse_grained_cellulose_post_processing import default_average_persistence_length
from coarse_grained_cellulose_post_processing import default_persistence_length_match_degree
from coarse_grained_cellulose_post_processing import default_average_end_to_end_distance
from coarse_grained_cellulose_post_processing import default_end_to_end_distance_match_degree

class CoarseGrainedCelluloseCoefficientsEnv(gym.Env):

    def generate_coarse_grained_cellulose_computation_working_directory(self,lammps_results_overwrite_switch):
        coarse_grained_cellulose_computation_working_directory_base = coarse_grained_cellulose_computation_home_directory+'coarse_grained_cellulose_coefficients_data/physics_computations/train/'+physics_version+'_physics_computations/'+physics_version+'_physics_computation'
        if lammps_results_overwrite_switch == True:
            coarse_grained_cellulose_computation_working_directory = coarse_grained_cellulose_computation_working_directory_base+'/'
        else:
            current_time = time.strftime('%Y-%m-%d-%H:%M:%S')
            coarse_grained_cellulose_computation_working_directory = coarse_grained_cellulose_computation_working_directory_base+'-'+current_time+'/'
        return coarse_grained_cellulose_computation_working_directory

    def __init__(self):

        super().__init__()

        self.max_cpu_logical_processor_no = os.cpu_count()-1
        self.lammps_results_overwrite_switch = True

        self.pcs_pool = []
        self.coarse_grained_cellulose_data_pool = []
        self.action = None
        self.reward = 0.0
        self.observation = 1

        self.terminated = False
        self.truncated = False

        self.step_count = 0

        self.action_space = gym.spaces.Box(
            (-1.0*numpy.ones((3),dtype=numpy.float32)),
            (+1.0*numpy.ones((3),dtype=numpy.float32)),
        )

        self.observation_space = gym.spaces.Discrete(2)

        train_data_directory = coarse_grained_cellulose_computation_home_directory+'/coarse_grained_cellulose_coefficients_data/train_data/'+physics_version+'_train_data/'
        pathlib.Path(train_data_directory).mkdir(parents=True,exist_ok=True)

        self.coefficients_temp_pool_file_location = train_data_directory+'coefficients_temp_pool.txt'
        self.coarse_grained_cellulose_data_temp_pool_file_location = train_data_directory+'coarse_grained_cellulose_data_temp_pool.txt'

        self.coefficients_pool_file_location_base = train_data_directory+'coefficients_pool-'
        self.coarse_grained_cellulose_data_pool_file_location_base = train_data_directory+'coarse_grained_cellulose_data_pool-'

        self.coarse_grained_cellulose_coefficients_pool_header = 'step_count,bc_1,bc_2,bc_3,ac_1,ac_2,ac_3,ac_4,ac_5,ac_6,ic_1,ic_2,ic_3'
        self.coarse_grained_cellulose_data_pool_header = 'step_count,fitted_elastic_modulus,elastic_modulus_match_degree,average_bond_length,average_persistence_length,persistence_length_match_degree,average_end_to_end_distance,end_to_end_distance_match_degree,reward_value'

    def step(self, action):

        self.step_count = self.step_count+1
        self.observation = 1

        if action is not None:

            temp_bc_1 = float(coefficient_precision_format % (10.0*action[0]+70.0))
            temp_bc_2 = float(coefficient_precision_format % (1.2866*temp_bc_1))
            temp_bc_3 = float(coefficient_precision_format % (1.7531*temp_bc_1))
            temp_ac_1 = float(coefficient_precision_format % (5.0*action[1]+30.0))
            temp_ac_2 = float(coefficient_precision_format % (1.2524*temp_ac_1))
            temp_ac_3 = float(coefficient_precision_format % (0.9854*temp_ac_1))
            temp_ac_4 = float(coefficient_precision_format % (0.8809*temp_ac_1))
            temp_ac_5 = float(coefficient_precision_format % (0.8386*temp_ac_1))
            temp_ac_6 = float(coefficient_precision_format % (0.7240*temp_ac_1))
            temp_ic_1 = float(coefficient_precision_format % (0.35*action[2]+0.50))
            temp_ic_2 = float(coefficient_precision_format % (0.9139*temp_ic_1))
            temp_ic_3 = float(coefficient_precision_format % (0.9979*temp_ic_1))

        temp_coefficients = (self.step_count,temp_bc_1,temp_bc_2,temp_bc_3,temp_ac_1,temp_ac_2,temp_ac_3,temp_ac_4,temp_ac_5,temp_ac_6,temp_ic_1,temp_ic_2,temp_ic_3)
        self.pcs_pool.append(temp_coefficients)
        numpy.savetxt(self.coefficients_temp_pool_file_location,self.pcs_pool,header=self.coarse_grained_cellulose_coefficients_pool_header,fmt='%.6e',delimiter=',')

        temp_coarse_grained_cellulose_computation_working_directory = self.generate_coarse_grained_cellulose_computation_working_directory(self.lammps_results_overwrite_switch)

        try:
            pathlib.Path(temp_coarse_grained_cellulose_computation_working_directory).mkdir(parents=True,exist_ok=True)
            shutil.rmtree(temp_coarse_grained_cellulose_computation_working_directory)
        except:
            current_time = time.strftime('%Y-%m-%d-%H:%M:%S')
            temp_coarse_grained_cellulose_computation_working_directory = temp_coarse_grained_cellulose_computation_working_directory+current_time+'/'

        temp_coarse_grained_cellulose_physics_computation_1_inputs_directory = temp_coarse_grained_cellulose_computation_working_directory+'physics_computation_1/inputs/'
        temp_coarse_grained_cellulose_physics_computation_1_results_directory = temp_coarse_grained_cellulose_computation_working_directory+'physics_computation_1/results/'
        temp_coarse_grained_cellulose_physics_computation_2_inputs_directory = temp_coarse_grained_cellulose_computation_working_directory+'physics_computation_2/inputs/'
        temp_coarse_grained_cellulose_physics_computation_2_results_directory = temp_coarse_grained_cellulose_computation_working_directory+'physics_computation_2/results/'
        temp_coarse_grained_cellulose_physics_computation_3_inputs_directory = temp_coarse_grained_cellulose_computation_working_directory+'physics_computation_3/inputs/'
        temp_coarse_grained_cellulose_physics_computation_3_results_directory = temp_coarse_grained_cellulose_computation_working_directory+'physics_computation_3/results/'
        temp_coarse_grained_cellulose_physics_computation_4_inputs_directory = temp_coarse_grained_cellulose_computation_working_directory+'physics_computation_4/inputs/'
        temp_coarse_grained_cellulose_physics_computation_4_results_directory = temp_coarse_grained_cellulose_computation_working_directory+'physics_computation_4/results/'

        for temp_directory in (temp_coarse_grained_cellulose_computation_working_directory,temp_coarse_grained_cellulose_physics_computation_1_inputs_directory,temp_coarse_grained_cellulose_physics_computation_1_results_directory,temp_coarse_grained_cellulose_physics_computation_2_inputs_directory,temp_coarse_grained_cellulose_physics_computation_2_results_directory,temp_coarse_grained_cellulose_physics_computation_3_inputs_directory,temp_coarse_grained_cellulose_physics_computation_3_results_directory,temp_coarse_grained_cellulose_physics_computation_4_inputs_directory,temp_coarse_grained_cellulose_physics_computation_4_results_directory):
            pathlib.Path(temp_directory).mkdir(parents=True,exist_ok=True)

        for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data'):
            shutil.copyfile(coarse_grained_cellulose_physics_computation_1_code_inputs_directory+file_name,temp_coarse_grained_cellulose_physics_computation_1_inputs_directory+file_name)

        for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data','compute_persistence_length_from_lammps_dump','compute_end_to_end_distance_from_lammps_dump'):
            shutil.copyfile(coarse_grained_cellulose_physics_computation_2_code_inputs_directory+file_name,temp_coarse_grained_cellulose_physics_computation_2_inputs_directory+file_name)

        for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data','compute_persistence_length_from_lammps_dump','compute_end_to_end_distance_from_lammps_dump'):
            shutil.copyfile(coarse_grained_cellulose_physics_computation_3_code_inputs_directory+file_name,temp_coarse_grained_cellulose_physics_computation_3_inputs_directory+file_name)

        for file_name in ('in_file_constant_head.txt','in_file_constant_tail.txt','coarse_grained_cellulose.data','compute_persistence_length_from_lammps_dump','compute_end_to_end_distance_from_lammps_dump'):
            shutil.copyfile(coarse_grained_cellulose_physics_computation_4_code_inputs_directory+file_name,temp_coarse_grained_cellulose_physics_computation_4_inputs_directory+file_name)

        temp_coarse_grained_cellulose_computation_class = CoarseGrainedCelluloseComputation(temp_coarse_grained_cellulose_computation_working_directory)
        temp_coarse_grained_cellulose_computation_class.generate_in_file(temp_bc_1,temp_bc_2,temp_bc_3,temp_ac_1,temp_ac_2,temp_ac_3,temp_ac_4,temp_ac_5,temp_ac_6,temp_ic_1,temp_ic_2,temp_ic_3)

        temp_coarse_grained_cellulose_computation_1_queue = multiprocessing.Queue()
        temp_coarse_grained_cellulose_computation_2_queue = multiprocessing.Queue()
        temp_coarse_grained_cellulose_computation_3_queue = multiprocessing.Queue()
        temp_coarse_grained_cellulose_computation_4_queue = multiprocessing.Queue()
        temp_coarse_grained_cellulose_computation_queue_array = (temp_coarse_grained_cellulose_computation_1_queue,temp_coarse_grained_cellulose_computation_2_queue,temp_coarse_grained_cellulose_computation_3_queue,temp_coarse_grained_cellulose_computation_4_queue)
        for queue in temp_coarse_grained_cellulose_computation_queue_array:
            queue.put(False)

        temp_coarse_grained_cellulose_computation_1_process = multiprocessing.Process(target=temp_coarse_grained_cellulose_computation_class.physics_computation_1,daemon=True,args=(temp_coarse_grained_cellulose_computation_1_queue,))
        temp_coarse_grained_cellulose_computation_2_process = multiprocessing.Process(target=temp_coarse_grained_cellulose_computation_class.physics_computation_2,daemon=True,args=(temp_coarse_grained_cellulose_computation_2_queue,))
        temp_coarse_grained_cellulose_computation_3_process = multiprocessing.Process(target=temp_coarse_grained_cellulose_computation_class.physics_computation_3,daemon=True,args=(temp_coarse_grained_cellulose_computation_3_queue,))
        temp_coarse_grained_cellulose_computation_4_process = multiprocessing.Process(target=temp_coarse_grained_cellulose_computation_class.physics_computation_4,daemon=True,args=(temp_coarse_grained_cellulose_computation_4_queue,))

        temp_coarse_grained_cellulose_computation_4_process.start()

        temp_coarse_grained_cellulose_computation_process_array = (temp_coarse_grained_cellulose_computation_1_process,temp_coarse_grained_cellulose_computation_2_process,temp_coarse_grained_cellulose_computation_3_process)

        for process in temp_coarse_grained_cellulose_computation_process_array:
            process.start()
            process.join(timeout=coarse_grained_cellulose_computation_timeout)
        
        temp_coarse_grained_cellulose_computation_4_process.join(timeout=coarse_grained_cellulose_computation_timeout)

        temp_coarse_grained_cellulose_computation_return_array = [False,False,False,False]
        for i in (0,1,2,3):
            temp_coarse_grained_cellulose_computation_return_array[i] = temp_coarse_grained_cellulose_computation_queue_array[i].get()
        for queue in temp_coarse_grained_cellulose_computation_queue_array:
            queue.close()
        for process in temp_coarse_grained_cellulose_computation_process_array:
            if process.is_alive() == True:
                process.kill()
            process.join()
            process.close()

        temp_coarse_grained_cellulose_computation_return = True
        for i in (0,1,2,3):
            if temp_coarse_grained_cellulose_computation_return_array[i] == False:
                temp_coarse_grained_cellulose_computation_return = False

        if temp_coarse_grained_cellulose_computation_return == True:
            temp_fitted_elastic_modulus,temp_elastic_modulus_match_degree,temp_average_bond_length,temp_average_persistence_length,temp_persistence_length_match_degree,temp_average_end_to_end_distance,temp_end_to_end_distance_match_degree = temp_coarse_grained_cellulose_computation_class.post_processing()
        else:
            temp_fitted_elastic_modulus = default_fitted_elastic_modulus
            temp_elastic_modulus_match_degree = default_elastic_modulus_match_degree
            temp_average_bond_length = default_average_bond_length
            temp_average_persistence_length = default_average_persistence_length
            temp_persistence_length_match_degree = default_persistence_length_match_degree
            temp_average_end_to_end_distance = default_average_end_to_end_distance
            temp_end_to_end_distance_match_degree = default_end_to_end_distance_match_degree

        temp_reward_value = 5.0*(temp_elastic_modulus_match_degree+temp_persistence_length_match_degree+temp_end_to_end_distance_match_degree)

        temp_coarse_grained_cellulose_data = (self.step_count,temp_fitted_elastic_modulus,temp_elastic_modulus_match_degree,temp_average_bond_length,temp_average_persistence_length,temp_persistence_length_match_degree,temp_average_end_to_end_distance,temp_end_to_end_distance_match_degree,temp_reward_value)
        self.coarse_grained_cellulose_data_pool.append(temp_coarse_grained_cellulose_data)
        numpy.savetxt(self.coarse_grained_cellulose_data_temp_pool_file_location,self.coarse_grained_cellulose_data_pool,header=self.coarse_grained_cellulose_data_pool_header,fmt='%.6e',delimiter=',')

        self.reward = temp_reward_value

        self.terminated = True

        return self.observation,self.reward,self.terminated,self.truncated,{}

    def reset(self, seed=None, options=None):
        self.observation = 1
        return self.observation, {}

    def render(self):
        print("No graphical output")

    def close(self):
        current_time = time.strftime('%Y-%m-%d-%H:%M:%S')
        numpy.savetxt(self.coefficients_pool_file_location_base+current_time+'.txt',self.pcs_pool,header=self.coarse_grained_cellulose_coefficients_pool_header,fmt='%.6e',delimiter=',')
        numpy.savetxt(self.coarse_grained_cellulose_data_pool_file_location_base+current_time+'.txt',self.coarse_grained_cellulose_data_pool,header=self.coarse_grained_cellulose_data_pool_header,fmt='%.6e',delimiter=',')
        print('None of unexpectations met.')

