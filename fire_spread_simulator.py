import numpy as np
import itertools

from ca_classes.field_class import Field
from ca_classes.fire_simulation_class import Fire_simulation
from ca_classes.MCE_class import MCE

########### External data ############

dimension = [40, 40]

test_cell_states = np.full(dimension, 1)
for coord in itertools.product(*[range(dim) for dim in dimension]):
    if coord[0] in list(range(40, 60)) and coord[1] in list(range(70, 80)):
        test_cell_states[coord] = 0

test_cell_height = np.zeros(dimension)
for i in range(dimension[0]):
    for j in range(dimension[1]):
        test_cell_height[i, j] = 1.6 * np.sqrt(i ** 2 + j ** 2)

#######################################

obj_field = Field(dimension=dimension,
                  wind_velocity=0,
                  wind_direction=[0, 0],
                  cell_states=None,
                  cell_height=None,
                  cell_veg_type=0,
                  cell_veg_density=0,
                  cell_size=15)

# obj_field.plot()

obj_fire_simul = Fire_simulation(field=obj_field,
                                 fire_origin=[(20, 20)],
                                 max_period_num=15,
                                 plot=True)

obj_MCE = MCE(ca_fire_simul=obj_fire_simul,
              rep_number=100)

obj_MCE.run()

obj_MCE.generate_report("sim_p_0.8_wind_height")

# obj_MCE.plot()


# obj_field = Field(dimension=[100, 100],
#                   wind_velocity=4,
#                   wind_direction=[0, 1],
#                   cell_states=None,
#                   cell_height=test_cell_height,
#                   cell_veg_type=0,
#                   cell_veg_density=0,
#                   cell_size=15)
#
# obj_fire_simul = Fire_simulation(field=obj_field,
#                                  fire_origin=[(20, 20)],
#                                  max_period_num=15,
#                                  plot=False)
#
# obj_MCE = MCE(ca_fire_simul=obj_fire_simul,
#               rep_number=100)
#
# obj_MCE.run()
