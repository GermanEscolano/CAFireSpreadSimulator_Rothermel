import itertools
import random
import numpy as np
import matplotlib.pyplot as plt

from ca_classes import neighborhood, field_class


class Fire_simulation:
    neighborhood_obj = neighborhood.MooreNeighborhood(neighborhood.EdgeRule.IGNORE_MISSING_NEIGHBORS_OF_EDGE_CELLS)
    p_h = 0.36
    C1 = 0.045
    C2 = 0.131
    C3 = 0.3
    verbose = False
    period_count = 0
    first_sim = True


    def __init__(self, field, fire_origin, max_period_num = 100, plot=False):
        self.field = field
        self.fire_origin = fire_origin
        self.max_period_num = max_period_num
        self.plot = plot

    def run(self):
        still_fire = True
        self.period_count = 0
        if self.plot:
            fig2, ax2 = plt.subplots(figsize=(8, 6))
        while self.period_count < self.max_period_num and still_fire:
            if self.max_period_num > 29:
                if self.period_count % 10 == 0:
                    print(f'Period {self.period_count} / {self.max_period_num}') ################
            if self.plot:
                ax2.clear()
                ax2.set_title('Ongoing simulation')
                ax2.text(105, 10, f'Period: {self.period_count}')
                ax2.imshow(self.field.cell_states)
                plt.draw()
                plt.pause(0.3)
                if self.first_sim == True:
                    input()
                    self.first_sim = False

            self.evolve()
            self.period_count += 1
            if np.any(self.field.cell_states == 2):
                still_fire = True
            else:
                still_fire = False
        if self.plot:
            plt.close(fig2)
        return self.field.cell_states

    def evolve(self):
        new_cell_states = np.copy(self.field.cell_states)
        for coord in itertools.product(*[range(dim) for dim in self.field.dimension]):
            if Fire_simulation.verbose: print(f'Evaluated cell: {coord}')
            if self.field.cell_states[coord] == 1:
                prob_no_set_fire = self.get_cell_prob_no_burn(coord)
                if random.random() > prob_no_set_fire:
                    new_cell_states[coord] = 2
            elif self.field.cell_states[coord] == 2:
                new_cell_states[coord] = 3
        self.field.cell_states = new_cell_states
        return new_cell_states

    def get_cell_prob_no_burn(self, coord):
        coord_neigs = self.neighborhood_obj.calculate_cell_neighbor_coordinates(coord, self.field.dimension)

        if Fire_simulation.verbose: print(f'Neig cells: {coord_neigs}')

        list_prob_propagate_from_neig = []

        for coord_neig in coord_neigs:
            coord_neig = tuple(coord_neig)
            if self.field.cell_states[coord_neig] == 2:
                prob_propagate_cell_to_cell = self.get_prob_propagate_from_neig(coord, coord_neig)
                list_prob_propagate_from_neig.append(prob_propagate_cell_to_cell)
                if Fire_simulation.verbose: print(
                    f'Burning_neig: {coord_neig} ---> prob fire propagates = {prob_propagate_cell_to_cell}')

        cell_prob_no_burn = np.prod([1 - elem for elem in list_prob_propagate_from_neig])

        if Fire_simulation.verbose: print(f'Prob evaluated cell does not set on fire: {cell_prob_no_burn}')
        if Fire_simulation.verbose: print(f'------------')

        return cell_prob_no_burn

    def get_prob_propagate_from_neig(self, coord_orig, coord_neig):

        p_veg = self.field.cell_veg_type[coord_neig]

        p_den = self.field.cell_veg_density[coord_neig]

        V = self.field.wind_velocity
        wind_direction = self.field.wind_direction
        vector_from_neig_to_orig = np.array(coord_orig) - np.array(coord_neig)

        propagation_wind_angle = self.angle_between_vectors(vector_from_neig_to_orig, wind_direction)
        f_t = np.exp(self.C2 * V * (np.cos(propagation_wind_angle)))
        p_w = f_t * np.exp(self.C1 * V)

        height_neig = self.field.cell_heigh[coord_neig]
        height_orig = self.field.cell_heigh[coord_orig]

        p_heigh = np.exp(self.C3 * (height_orig - height_neig))

        return self.p_h * (1 + p_veg) * (1 + p_den) * p_w * p_heigh

    def start_fire(self):
        for coord in self.fire_origin:
            self.field.cell_states[coord] = 2

    def set_fire_parameters(self, p_h, C1, C2, C3):
        self.p_h = p_h
        self.C1 = C1
        self.C2 = C2
        self.C3 = C3

    @staticmethod
    def angle_between_vectors(v1, v2):
        if type(v1) is tuple or type(v1) is tuple:
            v1 = np.array(v1)
            v2 = np.array(v2)
        """ Returns the angle in radians between vectors 'v1' and 'v2'::"""
        return np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))


if __name__ == '__main__':
    field_obj = field_class.Field([11, 11])

    fire_simul_obj = Fire_simulation(field_obj, [(5, 5)])

    fire_simul_obj.set_fire_parameters(1, 0.3, 0.2, 0.1)

    fire_simul_obj.start_fire()

    print(fire_simul_obj.field.cell_states)

    fig, ax = plt.subplots(figsize=(8,6))

    ax.imshow(fire_simul_obj.field.cell_states)
    plt.draw()
    plt.pause(1)

    for _ in range(4):
        current_state = fire_simul_obj.evolve()
        ax.clear()
        ax.imshow(current_state)
        print(current_state)
        plt.draw()
        plt.pause(1)

    print('-----------------------------------')

    field_obj2 = field_class.Field([11, 11])

    fire_simul_obj2 = Fire_simulation(field_obj2, [(5, 5)], 4)

    fire_simul_obj2.set_fire_parameters(1, 0.3, 0.2, 0.1)

    fire_simul_obj2.start_fire()

    print(fire_simul_obj2.field.cell_states)

    final_state = fire_simul_obj2.run()

    fig2, ax2 = plt.subplots(figsize=(8, 6))

    ax2.imshow(final_state)
    plt.draw()
    plt.pause(3)


