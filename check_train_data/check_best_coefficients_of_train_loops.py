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

best_cellulose_model_corresponding_coeffcients_pools_directory = coarse_grained_cellulose_computation_home_directory+'coarse_grained_cellulose_coefficients_data/best_cellulose_model_corresponding_coeffcients_pools/'
pathlib.Path(best_cellulose_model_corresponding_coeffcients_pools_directory).mkdir(parents=True,exist_ok=True)

print('\33c',end='\r')
print('Target version: {}'.format(target_version))

train_data_file_name_list = sorted(os.listdir(train_data_directory))
cellulose_data_pool_file_name_list = []
coeffcients_pool_file_name_list = []

for train_data_file_name in train_data_file_name_list:
    if train_data_file_name.__contains__('_cellulose_data_pool-'):
        cellulose_data_pool_file_name_list.append(train_data_file_name)
    elif train_data_file_name.__contains__('coefficients_pool-'):
        coeffcients_pool_file_name_list.append(train_data_file_name)

total_cellulose_data_pool_file_number = len(cellulose_data_pool_file_name_list)
total_coeffcients_pool_file_number = len(coeffcients_pool_file_name_list)

if total_cellulose_data_pool_file_number > 0 and total_cellulose_data_pool_file_number == total_coeffcients_pool_file_number:

    step_count_column_number = 0
    reward_data_column_number = -1

    best_cellulose_model_reward_value_pool = []
    best_cellulose_model_corresponding_step_count_pool = []
    best_cellulose_model_corresponding_cellulose_data_pool = []
    best_cellulose_model_corresponding_coeffcients_pool = []

    total_empty_cellulose_data_pool_file_number = 0
    total_unmatched_cellulose_data_pool_file_number = 0

    for i in range(0,total_cellulose_data_pool_file_number):

        temp_cellulose_data_pool_file_name = cellulose_data_pool_file_name_list[i]
        temp_coeffcients_pool_file_name = coeffcients_pool_file_name_list[i]
        try:
            temp_cellulose_data_pool = numpy.loadtxt(train_data_directory+temp_cellulose_data_pool_file_name)
            temp_coeffcients_pool = numpy.loadtxt(train_data_directory+temp_coeffcients_pool_file_name)
        except:
            temp_cellulose_data_pool = numpy.loadtxt(train_data_directory+temp_cellulose_data_pool_file_name,delimiter=',')
            temp_coeffcients_pool = numpy.loadtxt(train_data_directory+temp_coeffcients_pool_file_name,delimiter=',')

        if len(temp_cellulose_data_pool) > 0:

            if temp_cellulose_data_pool.ndim > 1:

                temp_reward_data_column = temp_cellulose_data_pool[:,reward_data_column_number]
                temp_max_reward_value = numpy.max(temp_reward_data_column)
                temp_max_reward_value_corresponding_cellulose_data_row_number = numpy.argmax(temp_reward_data_column)
                temp_max_reward_value_corresponding_cellulose_data_row_data = temp_cellulose_data_pool[temp_max_reward_value_corresponding_cellulose_data_row_number,:]
                temp_max_reward_value_corresponding_step_count = int(temp_max_reward_value_corresponding_cellulose_data_row_data[step_count_column_number])
                temp_max_reward_value_corresponding_cellulose_data = temp_max_reward_value_corresponding_cellulose_data_row_data[step_count_column_number+1:]

                try:
                    temp_max_reward_value_corresponding_coeffcients_row_data = temp_coeffcients_pool[numpy.where(temp_coeffcients_pool[:,step_count_column_number]==temp_max_reward_value_corresponding_step_count)[0][0],:]
                    temp_max_reward_value_corresponding_coeffcients = temp_max_reward_value_corresponding_coeffcients_row_data[step_count_column_number+1:]

                    best_cellulose_model_reward_value_pool.append(temp_max_reward_value)
                    best_cellulose_model_corresponding_step_count_pool.append(temp_max_reward_value_corresponding_step_count)

                    if len(best_cellulose_model_corresponding_cellulose_data_pool) == 0:
                        best_cellulose_model_corresponding_cellulose_data_pool = temp_max_reward_value_corresponding_cellulose_data
                    else:
                        best_cellulose_model_corresponding_cellulose_data_pool = numpy.vstack((best_cellulose_model_corresponding_cellulose_data_pool,temp_max_reward_value_corresponding_cellulose_data))

                    if len(best_cellulose_model_corresponding_coeffcients_pool) == 0:
                        best_cellulose_model_corresponding_coeffcients_pool = temp_max_reward_value_corresponding_coeffcients
                    else:
                        best_cellulose_model_corresponding_coeffcients_pool = numpy.vstack((best_cellulose_model_corresponding_coeffcients_pool,temp_max_reward_value_corresponding_coeffcients))
                except:
                    total_unmatched_cellulose_data_pool_file_number = total_unmatched_cellulose_data_pool_file_number+1

            else:

                temp_max_reward_value = temp_cellulose_data_pool[reward_data_column_number]
                temp_max_reward_value_corresponding_step_count = temp_cellulose_data_pool[step_count_column_number]
                temp_max_reward_value_corresponding_cellulose_data = temp_cellulose_data_pool[step_count_column_number+1:]

                if temp_max_reward_value_corresponding_step_count == temp_coeffcients_pool[step_count_column_number]:

                    temp_max_reward_value_corresponding_coeffcients = temp_coeffcients_pool[step_count_column_number+1:]

                    best_cellulose_model_reward_value_pool.append(temp_max_reward_value)
                    best_cellulose_model_corresponding_step_count_pool.append(temp_max_reward_value_corresponding_step_count)

                    if len(best_cellulose_model_corresponding_cellulose_data_pool) == 0:
                        best_cellulose_model_corresponding_cellulose_data_pool = temp_max_reward_value_corresponding_cellulose_data
                    else:
                        best_cellulose_model_corresponding_cellulose_data_pool = numpy.vstack((best_cellulose_model_corresponding_cellulose_data_pool,temp_max_reward_value_corresponding_cellulose_data))

                    if len(best_cellulose_model_corresponding_coeffcients_pool) == 0:
                        best_cellulose_model_corresponding_coeffcients_pool = temp_max_reward_value_corresponding_coeffcients
                    else:
                        best_cellulose_model_corresponding_coeffcients_pool = numpy.vstack((best_cellulose_model_corresponding_coeffcients_pool,temp_max_reward_value_corresponding_coeffcients))

                else:
                    total_unmatched_cellulose_data_pool_file_number = total_unmatched_cellulose_data_pool_file_number+1
        else:
            total_empty_cellulose_data_pool_file_number = total_empty_cellulose_data_pool_file_number+1

    print('\nTotal train loop number: {}'.format(total_cellulose_data_pool_file_number))
    print('\nTotal empty cellulose data pool file number: {}'.format(total_empty_cellulose_data_pool_file_number))
    print('\nTotal unmatched cellulose data pool file number: {}'.format(total_unmatched_cellulose_data_pool_file_number))
    print('\nBest cellulose model reward value pool:')
    print(numpy.array(best_cellulose_model_reward_value_pool))
    print('\nBest cellulose model corresponding step count pool:')
    print(numpy.array(best_cellulose_model_corresponding_step_count_pool))
    print('\nBest cellulose model corresponding cellulose data pool:')
    if best_cellulose_model_corresponding_cellulose_data_pool.ndim > 1:
        for i in range(0,len(best_cellulose_model_corresponding_cellulose_data_pool)):
           print(i+1,':',*best_cellulose_model_corresponding_cellulose_data_pool[i],sep=',')
    else:
        print(1,':',*best_cellulose_model_corresponding_cellulose_data_pool,sep=',')
    print('\nBest cellulose model corresponding coeffcients pool:')
    if best_cellulose_model_corresponding_coeffcients_pool.ndim > 1:
        for i in range(0,len(best_cellulose_model_corresponding_coeffcients_pool)):
            print(i+1,':',*best_cellulose_model_corresponding_coeffcients_pool[i],sep=',')
    else:
        print(1,':',*best_cellulose_model_corresponding_coeffcients_pool,sep=',')
    numpy.savetxt(best_cellulose_model_corresponding_coeffcients_pools_directory+target_version+'_best_cellulose_model_corresponding_coeffcients_pool.txt',best_cellulose_model_corresponding_coeffcients_pool)

else:
    print('Improper train data files')
