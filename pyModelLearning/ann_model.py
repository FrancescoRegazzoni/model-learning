import scipy.io as sio
import numpy as np
import configparser
import os

class ANNmodel:
    def __init__(self, path, relative = True):

        if relative:
            config = configparser.ConfigParser()
            script_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))

            print(script_path)

            config.read(script_path + '/options.ini')
            datapath = config['paths']['datapath']
            path = datapath + '/' + path

            print(path)

        data = sio.loadmat(path)

        self.num_states = data['N'][0,0]
        self.num_inputs = data['nU'][0,0]
        self.num_outputs = data['nY'][0,0]
        self.use_G = data['useG'][0,0] > 0
        self.initial_state = data['x0'][:, 0]

        if len(self.initial_state) != self.num_states:
            raise Exception('x0 has the wrong size')

        self.f_weights = data['W'][0]
        self.f_biases = data['T'][0]
        self.f_num_hidden_layers = len(self.f_weights) - 1

        if self.use_G:
            self.g_weights = data['W_G'][0]
            self.g_biases = data['T_G'][0]
            self.g_num_hidden_layers = len(self.g_weights) - 1

        self.rhs = lambda x, u: self.ANN(np.concatenate([u,x]), self.f_weights, self.f_biases)

        if self.use_G:
            self.obs = lambda x: self.ANN(x, self.g_weights, self.g_biases)
        else:
            self.obs = lambda x: x[:self.num_outputs]


    def ANN(self, input, weights, biases):
        y = input
        for i in range(len(weights)):
            y = np.matmul(weights[i], y) - biases[i][:,0]
            if i < len(weights) - 1:
                y = np.tanh(y)
        return y