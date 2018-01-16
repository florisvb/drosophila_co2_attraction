import orchard.annotated_time_on_pad_analysis as ann
import pickle
import numpy as np

if __name__ == '__main__':

    data_dict_lengths = ann.load_lengths()
    labels = ['co2_15sccm', 'co2_60sccm', 'eth_60sccm', 'eth_200sccm']
    data = {label: np.array(data_dict_lengths[label])/30. for label in labels}

    f = open('co2_eth_time_on_pad.pickle', 'w')
    pickle.dump(data, f)
    f.close()