import numpy
import os
import sys

target_version = sys.argv[1]

current_user_main_directory = os.path.expanduser('~')+'/'
try:
    coarse_grained_cellulose_computation_home_directory = os.getenv('CLL_CG_HOME')+'/'
except:
    coarse_grained_cellulose_computation_home_directory = current_user_main_directory

train_data_directory = coarse_grained_cellulose_computation_home_directory+'coarse_grained_cellulose_coefficients_data/train_data/{0}_train_data/'.format(target_version)

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

    step_count_column_number = 0
    reward_data_column_number = -1

    latest_cellulose_data_pool_file_name = cellulose_data_pool_file_name_list[-1]
    latest_coefficients_pool_file_name = coefficients_pool_file_name_list[-1]

    try:
        latest_cellulose_data_pool = numpy.loadtxt(train_data_directory+latest_cellulose_data_pool_file_name)
        latest_coefficients_pool = numpy.loadtxt(train_data_directory+latest_coefficients_pool_file_name)
    except:
        latest_cellulose_data_pool = numpy.loadtxt(train_data_directory+latest_cellulose_data_pool_file_name,delimiter=',')
        latest_coefficients_pool = numpy.loadtxt(train_data_directory+latest_coefficients_pool_file_name,delimiter=',')

    if len(latest_cellulose_data_pool) > 0:

        if latest_cellulose_data_pool.ndim > 1:

            step_count_column = latest_cellulose_data_pool[:,step_count_column_number]
            reward_data_column = latest_cellulose_data_pool[:,reward_data_column_number]

            total_record_number = len(reward_data_column)
            print('\nTotal step record number: {}'.format(total_record_number))

            min_step_count = int(numpy.min(step_count_column))
            max_step_count = int(numpy.max(step_count_column))
            print('Step count range: {0}-{1}'.format(min_step_count,max_step_count))

            mean_reward_value = numpy.mean(reward_data_column)
            print('Mean reward value: {}'.format(mean_reward_value))

            max_reward_value = numpy.max(reward_data_column)
            print('\nMax reward value: {}'.format(max_reward_value))

            max_reward_value_corresponding_cellulose_data_row_number = numpy.argmax(reward_data_column)
            max_reward_value_corresponding_cellulose_data_row_data = latest_cellulose_data_pool[max_reward_value_corresponding_cellulose_data_row_number,:]
            max_reward_value_corresponding_step_count = int(max_reward_value_corresponding_cellulose_data_row_data[step_count_column_number])
            max_reward_value_corresponding_cellulose_data = max_reward_value_corresponding_cellulose_data_row_data[step_count_column_number+1:]

            try:
                max_reward_value_corresponding_coefficients_row_data = latest_coefficients_pool[numpy.where(latest_coefficients_pool[:,step_count_column_number]==max_reward_value_corresponding_step_count)[0][0],:]
                max_reward_value_corresponding_coefficients = max_reward_value_corresponding_coefficients_row_data[step_count_column_number+1:]
                print('Max reward corresponding step count: {}'.format(max_reward_value_corresponding_step_count))
                print('Max reward corresponding cellulose data: ',*max_reward_value_corresponding_cellulose_data,sep=',')
                print('Max reward corresponding coefficients: ',*max_reward_value_corresponding_coefficients,sep=',')
                best_cellulose_data_number = min(10,len(reward_data_column))
                sorted_cellulose_data_pool = latest_cellulose_data_pool[(-reward_data_column).argsort()]
                max_several_reward_values_corresponding_cellulose_data_row_data_pool = sorted_cellulose_data_pool[0:best_cellulose_data_number,:]
                max_several_reward_values_corresponding_step_count_pool = max_several_reward_values_corresponding_cellulose_data_row_data_pool[:,step_count_column_number]
                max_several_reward_values_corresponding_cellulose_data_pool = max_several_reward_values_corresponding_cellulose_data_row_data_pool[:,step_count_column_number+1:]
                print('\nBest {} cellulose data:'.format(best_cellulose_data_number))
                for i in range(0,best_cellulose_data_number):
                    temp_step_count = max_several_reward_values_corresponding_step_count_pool[i]
                    temp_cellulose_data = max_several_reward_values_corresponding_cellulose_data_pool[i]
                    print(int(temp_step_count),':',*temp_cellulose_data,sep=',')
                print('\nBest {} cellulose data corresponding coefficients:'.format(best_cellulose_data_number))
                for i in range(0,best_cellulose_data_number):
                    temp_step_count = max_several_reward_values_corresponding_step_count_pool[i]
                    temp_max_reward_value_corresponding_coefficients_row_data = latest_coefficients_pool[numpy.where(latest_coefficients_pool[:,step_count_column_number]==temp_step_count)[0][0],:]
                    temp_max_reward_value_corresponding_coefficients = temp_max_reward_value_corresponding_coefficients_row_data[step_count_column_number+1:]
                    print(int(temp_step_count),':',*temp_max_reward_value_corresponding_coefficients,sep=',')
            except:
                print('cellulose data pool and coefficients pool are not matched.')

        else:

            print('\nOnly one step data')

            max_reward_value = latest_cellulose_data_pool[reward_data_column_number]
            print('\nMax reward value: {}'.format(max_reward_value))

            max_reward_value_corresponding_step_count = int(latest_cellulose_data_pool[step_count_column_number])
            max_reward_value_corresponding_cellulose_data = latest_cellulose_data_pool[step_count_column_number+1:]

            if max_reward_value_corresponding_step_count == latest_coefficients_pool[step_count_column_number]:

                max_reward_value_corresponding_coefficients = latest_coefficients_pool[step_count_column_number+1:]
                print('Step Count: {}'.format(max_reward_value_corresponding_step_count))
                print('Corresponding cellulose_data: ',*max_reward_value_corresponding_cellulose_data,sep=',')
                print('Corresponding coefficients: ',*max_reward_value_corresponding_coefficients,sep=',')

            else:
                print('cellulose data pool and coefficients pool are not matched.')

    else:
        print('Empty cellulose data file.')

else:
    print('Improper train data files')
