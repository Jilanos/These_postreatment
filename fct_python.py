import numpy as np
#VITRO
general_path = 'F:\\data_vitro\\iter_21\\PULSE_0_75\\'


general_path = 'C:\\Users\\PM263553\\Desktop\\These\\big_projects\\monkey_mircen\\Manip_26__21_06_2023\\tof_all_channels\\'

data = np.load(general_path+'data.npy')


print(data.shape)

def alerte_absence_de_signal(array, seuil_still = 3000, seuil_noise = 3, n_noise = 10):
    resultat = 1
    for i in range(array.shape[0]):
        mean, std, count = np.mean(array[i,:seuil_still]), np.std(array[i,:seuil_still]), 0
        for j in range(seuil_still,array.shape[1]):
            if array[i,j] > mean + seuil_noise*std:
                count += 1
        if count > n_noise:
            resultat *=1
        else:
            resultat *=0
            print('No signal on channel '+str(i+1))
    if resultat == 1:
        print('All channels are OK!!') 

alerte_absence_de_signal(data, seuil_still = 3000, seuil_noise = 10, n_noise = 10)




    