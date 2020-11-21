
import itertools
import numpy as np
import matplotlib.pyplot as plt


class Field:
    """ Class that contains all necessary info about the field and how is discretized in the simulation.
    """

    field_cond_fig = None

    veg_type = {'agricultural': -0.3,
                'thickets': 0,
                'Hallepo-pine': 0.4}
    #Just example, should be mofified

    veg_density = {'sparse': -0.4,
                   'normal': 0,
                   'dense': 0.3}
    #Just example, should be mofified

    def __init__(self, dimension, wind_velocity = 0, wind_direction = [0, 0], cell_states = None, cell_height = None,
                 cell_veg_type = None, cell_veg_density = None, cell_size = 10):
        self.dimension = dimension
        self.wind_velocity = wind_velocity
        self.wind_direction = wind_direction

        if cell_states is None:
            cell_states = np.full(dimension, 1) #maybe should change 1 for a constant like FUEL
        self.cell_states = cell_states

        if cell_height is None:
            cell_height = np.zeros(dimension)
        self.cell_heigh = cell_height

        if cell_veg_type is None:
            cell_veg_type = np.full(dimension, self.veg_type['thickets'])
        self.cell_veg_type = np.full(dimension, cell_veg_type)

        if cell_veg_density is None:
            cell_veg_density = np.full(dimension, self.veg_density['normal'])
        self.cell_veg_density = np.full(dimension, cell_veg_density)

        self.cell_size = cell_size
        self.original_state = self.cell_states

        self.field_cond_fig, ax1 = plt.subplots(2, 2, figsize=(8, 6))
        ax1[0, 0].set_title('Cell Heights')
        h_ax = ax1[0, 0].imshow(self.cell_heigh)
        self.field_cond_fig.colorbar(h_ax, ax=ax1[0, 0])
        ax1[0, 1].set_title('Wind direction and velocity')
        ax1[0, 1].text(0.95, 0.01, f'Wind velocity: {self.wind_velocity} m/s',
                       verticalalignment='bottom', horizontalalignment='right',
                       transform=ax1[0, 1].transAxes,
                       color='green', fontsize=15)
        ax1[0, 1].quiver([0], [0], [self.wind_direction[1]], [-self.wind_direction[0]], angles='xy',
                         scale_units='xy',
                         scale=1.6)
        ax1[0, 1].set_axis_off()
        ax1[0, 1].set_xlim(-1, 1)
        ax1[0, 1].set_ylim(-1, 1)
        ax1[1, 0].set_title('Vegetation type')
        ax1[1, 0].imshow(self.cell_veg_type)
        ax1[1, 1].set_title('Vegetation density')
        ax1[1, 1].imshow(self.cell_veg_density)


    def set_states(self, state_mat):
        '''
        TO DO add funtionality to make sure is adding correct type of data and size
        :param state_mat:
        :return:
        '''
        self.cell_states = state_mat

    def set_heights(self, heights_mat):
        self.cell_heigh = heights_mat

    def reset_state(self):
        self.cell_states = self.original_state

    def plot(self):
        plt.draw()
        plt.pause(0.1)

