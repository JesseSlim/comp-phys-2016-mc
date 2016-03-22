import numpy as np
import random
import math
from numba import jit

class Polymer:
    def __init__(self, bead_position, weight, size, theta, epsilon, sigma, T):
        self.bead_position = bead_position
        self.weight = weight
        self.size = size
        self.theta = theta
        self.epsilon = epsilon
        self.sigma = sigma
        self.T = T
        self.alive = True
        self.weight3 = 1
    def addBead(self):
        theta_weight_sum = 0                       # Sum of the weights
        while (theta_weight_sum == 0):
            theta_offset = random.uniform(0,1)
            theta_weight = np.zeros([self.theta,1])    # Weight of angles
            for i in range(0,self.theta):
                E = 0
                new_bead_x = self.bead_position[self.size-1, 0] + math.cos(2 * math.pi * (i+theta_offset) / self.theta);
                new_bead_y = self.bead_position[self.size-1, 1] + math.sin(2 * math.pi * (i+theta_offset) / self.theta);
                for j in range(0,self.size):
                    dx = self.bead_position[j, 0] - new_bead_x
                    dy = self.bead_position[j, 1] - new_bead_y
                    d = math.sqrt(dx**2+dy**2)
                    E += 4 * self.epsilon *((self.sigma / d)**12 - (self.sigma/d)**6)
                theta_weight[i] = math.exp(-E/self.T)
                theta_weight_sum += theta_weight[i]
        
        rand = random.uniform(0,1)
        self.weight = self.weight * theta_weight_sum
        theta_weight = np.divide(theta_weight, theta_weight_sum)
        for i in range(0,self.theta):
            rand -= theta_weight[i]
            if (rand < 0):
                direction = i
                new_bead = np.zeros([1,2])
                new_bead[0,0] = self.bead_position[self.size-1,0] + math.cos(2 * math.pi * (i+theta_offset) / self.theta)
                new_bead[0,1] = self.bead_position[self.size-1,1] + math.sin(2 * math.pi * (i+theta_offset) / self.theta)
                self.bead_position = np.append(self.bead_position, new_bead,axis=0)
                self.size += 1
                break
        if (self.size == 3):
            self.weight3 = self.weight;
        return self.weight;
    def getLength(self, l):
        if l < 2:
            l = self.size
        dx = (self.bead_position[0,0] - self.bead_position[l-1,0])**2
        dy = (self.bead_position[0,1] - self.bead_position[l-1,1])**2
        return math.sqrt(dx+dy)
    def getPosition (self,i):
        if (self.size > i):
            return self.bead_position[i,:]
        else:
            return np.zeros([1,2])