import numpy
import os
import pathlib
import sys

target_version = sys.argv[1]

current_user_main_directory = str(pathlib.Path.home())+'/'
try:
    coarse_grained_cellulose_computation_home_directory = os.getenv('CLL_CG_HOME')+'/'
except:
    coarse_grained_cellulose_computation_home_directory = current_user_main_directory

train_data_directory = coarse_grained_cellulose_computation_home_directory+'coarse_grained_cellulose_coefficients_data/train_data/{0}_train_data/'.format(target_version)

train_loop_statistic_info_pools_directory = coarse_grained_cellulose_computation_home_directory+'coarse_grained_cellulose_coefficients_data/train_loop_statistic_info_pools/'
pathlib.Path(train_loop_statistic_info_pools_directory).mkdir(parents=True,exist_ok=True)

print('\33c',end='\r')
print('Target version: {}'.format(target_version))

train_data_file_name_list = sorted(os.listdir(train_data_directory))
cellulose_data_pool_file_name_list = []
coefficients_pool_file_name_list = []

for train_data_file_name in train_data_file_name_list:
    if train_data_file_name.__contains__('cellulose_data_pool-'):
        cellulose_data_pool_file_name_list.append(train_data_file_name)
    elif train_data_file_name.__contains__('coefficients_pool-'):
        coefficients_pool_file_name_list.append(train_data_file_name)

total_cellulose_data_pool_file_number = len(cellulose_data_pool_file_name_list)
total_coefficients_pool_file_number = len(coefficients_pool_file_name_list)

if total_cellulose_data_pool_file_number > 0 and total_cellulose_data_pool_file_number == total_coefficients_pool_file_number:

    reward_data_column_number = -1
    train_data_statistic_info_pool = numpy.zeros((total_cellulose_data_pool_file_number,3))

    for i in range(0,total_cellulose_data_pool_file_number):

        temp_cellulose_data_pool_file_name = cellulose_data_pool_file_name_list[i]
        try:
            temp_cellulose_data_pool = numpy.loadtxt(train_data_directory+temp_cellulose_data_pool_file_name)
        except:
            temp_cellulose_data_pool = numpy.loadtxt(train_data_directory+temp_cellulose_data_pool_file_name,delimiter=',')

        if len(temp_cellulose_data_pool) > 0:
            if temp_cellulose_data_pool.ndim > 1:
                temp_reward_data_column = temp_cellulose_data_pool[:,reward_data_column_number]
                temp_total_step_number = len(temp_reward_data_column)
                temp_mean_reward_value = numpy.mean(temp_reward_data_column)
                temp_max_reward_value = numpy.max(temp_reward_data_column)
            else:
                temp_total_step_number = 1
                temp_mean_reward_value = temp_cellulose_data_pool[reward_data_column_number]
                temp_max_reward_value = temp_cellulose_data_pool[reward_data_column_number]
        else:
            temp_total_step_number = 0
            temp_mean_reward_value = 0.0
            temp_max_reward_value = 0.0

        train_data_statistic_info_pool[i,:] = (temp_total_step_number,temp_mean_reward_value,temp_max_reward_value)

    print('\nTotal train loop number: {}'.format(total_cellulose_data_pool_file_number))
    print('\nTotal step number|Mean reward|Max reward\nof each train loop\n')
    print(train_data_statistic_info_pool)
    numpy.savetxt(train_loop_statistic_info_pools_directory+target_version+'_train_loop_statistic_info_pool.txt',train_data_statistic_info_pool)

else:
    print('Improper train data files')
