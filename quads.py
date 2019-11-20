#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np

class Quadrature:
    def __init__(self, mus, weights):
        self.mus = mus
        self.weights = weights

        assert len(self.mus) == len(self.weights)

    def plot(self):

