import gymnasium as gym
import sys
from stable_baselines3.common.env_checker import check_env

target_version = sys.argv[1]

env_name = '{}CoarseGrainedCelluloseCoefficientsEnv-v0'.format(target_version)

print('\33c',end='\r')
print('Target version: {}'.format(target_version))
print('Reinforcement algorithm: SAC')

env = gym.make(env_name)

check_env(env)

env.close()
