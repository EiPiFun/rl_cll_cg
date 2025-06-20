import gymnasium as gym
import os
import pathlib
import sys
import time
from stable_baselines3 import SAC

target_version = sys.argv[1]

current_user_main_directory = str(pathlib.Path.home())+'/'
try:
    coarse_grained_cellulose_computation_home_directory = os.getenv('CLL_CG_HOME')+'/'
except:
    coarse_grained_cellulose_computation_home_directory = current_user_main_directory

trained_agent_directory = coarse_grained_cellulose_computation_home_directory+'coarse_grained_cellulose_coefficients_data/train_data/{0}_train_data/{0}_trained_agents/'.format(target_version)
trained_agent_file_name_base = '{}_coarse_grained_cellulose_sac_agent-'.format(target_version)

pathlib.Path(trained_agent_directory).mkdir(parents=True,exist_ok=True)

env_name = '{}CoarseGrainedCelluloseCoefficientsEnv-v0'.format(target_version)

print('\33c',end='\r')
print('Target version: {}'.format(target_version))
print('Reinforcement algorithm: SAC')

env = gym.make(env_name)

model = SAC('MlpPolicy',env).learn(total_timesteps=32)

current_time = time.strftime('%Y-%m-%d-%H:%M:%S')
model.save(trained_agent_directory+trained_agent_file_name_base+current_time+'.zip')

env.close()
