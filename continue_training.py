import gymnasium as gym
import os
import sys
import time
from stable_baselines3 import SAC

target_version = sys.argv[1]

one_train_loop_timestep = int(sys.argv[2])

total_train_loop_number = int(sys.argv[3])

current_user_main_directory = os.path.expanduser('~')+'/'
try:
    coarse_grained_cellulose_computation_home_directory = os.getenv('CLL_CG_HOME')+'/'
except:
    coarse_grained_cellulose_computation_home_directory = current_user_main_directory

trained_agent_directory = coarse_grained_cellulose_computation_home_directory+'coarse_grained_cellulose_coefficients_data/train_data/{0}_train_data/{0}_trained_agents/'.format(target_version)
trained_agent_file_name_base = '{}_coarse_grained_cellulose_sac_agent-'.format(target_version)

env_name = '{}CoarseGrainedCelluloseCoefficientsEnv-v0'.format(target_version)

for i in range(1,total_train_loop_number+1):

    print('\33c',end='\r')
    print('Target version: {}'.format(target_version))
    print('Reinforcement algorithm: SAC')
    print('Loop Count: '+str(i))

    env = gym.make(env_name)

    temp_file_name_list = sorted(os.listdir(trained_agent_directory))
    temp_trained_agent_file_name_list = []

    for file_name in temp_file_name_list:
        if file_name.__contains__('_sac_agent-'):
            temp_trained_agent_file_name_list.append(file_name)

    temp_latest_trained_agent_file_name = temp_trained_agent_file_name_list[-1]

    model = SAC('MlpPolicy',env).load(trained_agent_directory+temp_latest_trained_agent_file_name,env)
    model.learn(total_timesteps=one_train_loop_timestep)

    current_time = time.strftime('%Y-%m-%d-%H:%M:%S')

    model.save(trained_agent_directory+trained_agent_file_name_base+current_time+'.zip')

    env.close()
