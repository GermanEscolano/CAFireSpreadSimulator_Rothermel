


import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import shutil
from decimal import Decimal

class MCE:

    rep_count = 0
    final_res_fig = None

    def __init__(self, ca_fire_simul, rep_number):
        self.ca_fire_simul = ca_fire_simul
        dimension = self.ca_fire_simul.field.dimension
        self.rep_number = rep_number
        self.running_avg = np.zeros(dimension)
        self.running_var = np.zeros(dimension)

    def run(self, running_avg_step = 0, running_var_step = 0, verbose = False, plot_results = False):
        self.rep_count = 0
        s_k = 0
        while self.rep_count < self.rep_number:
            print(f'Replication {self.rep_count} / {self.rep_number}')  ################
            self.ca_fire_simul.field.reset_state()
            self.ca_fire_simul.start_fire()
            current_cell_state = (self.ca_fire_simul.run()>1).astype(float)
            last_running_avg = self.running_avg

            self.running_avg = self.running_avg + (
                current_cell_state - self.running_avg) / (self.rep_count + 1)

            s_k = s_k + (current_cell_state - self.running_avg) * (
                current_cell_state - last_running_avg)
            self.running_var = s_k / (self.rep_count + 1)

            self.rep_count += 1

        self.final_res_fig, ax = plt.subplots(figsize=(8, 6))
        ax.set_title("Probability of burned cell")
        avg = ax.imshow(self.running_avg)
        self.final_res_fig.colorbar(avg, ax=ax)


        return self.running_avg

    def plot(self):
        plt.show()

    def generate_report(self, report_name = None):

        if report_name == None:
            report_name = "simulation_report" + '.csv'
        else:
            report_name = report_name + '.csv'

        # Path to be created
        dir_path = "/Users/gescolano/PycharmProjects/Base_Line_CA_Fire_Simulator/sim_resuls/" + report_name

        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
        os.makedirs(dir_path)

        file_name = "sim_parameters"

        file_path = dir_path + "/" + file_name + '.csv'
        with open(file_path, 'w+') as f:
            writer = csv.writer(f)
            writer.writerow(['Param', 'Value'])
            writer.writerow(("dimension", self.ca_fire_simul.field.dimension))
            writer.writerow(("cell_size", self.ca_fire_simul.field.cell_size))
            writer.writerow(("fire_origin", self.ca_fire_simul.fire_origin))
            writer.writerow(("max_period_num", self.ca_fire_simul.max_period_num))
            writer.writerow(("rep_number", self.rep_number))
            writer.writerow(("p_h", self.ca_fire_simul.p_h))
            # writer.writerow(("C1", self.ca_fire_simul.C1))
            # writer.writerow(("C2", self.ca_fire_simul.C2))
            # writer.writerow(("test", '%.3E' % Decimal('0.139')))
            # writer.writerow(("C3", self.ca_fire_simul.C3))

        self.ca_fire_simul.field.field_cond_fig.savefig(dir_path + "/" + "field_condition.png")

        self.final_res_fig.savefig(dir_path + "/" + "result.png")





if __name__ == '__main__':
    from ca_classes import field_class, fire_simulation_class

    field_obj = field_class.Field([300, 300])

    fire_simul_obj = fire_simulation_class.Fire_simulation(field_obj, [(5, 5)], 1000)

    fire_simul_obj.start_fire()

    MCE_obj = MCE(fire_simul_obj, 10)

    final_prob_burn = MCE_obj.run()

    print(final_prob_burn)

    final_state = fire_simul_obj.run()

    fig1, ax1 = plt.subplots(figsize=(8, 6))

    ax1.imshow(final_prob_burn)

    plt.show()
