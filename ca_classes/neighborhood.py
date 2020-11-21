"""
Copyright 2019 Richard Feistenauer

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import enum
import operator
import itertools
import math


class EdgeRule(enum.Enum):
    IGNORE_EDGE_CELLS = 0
    IGNORE_MISSING_NEIGHBORS_OF_EDGE_CELLS = 1
    FIRST_AND_LAST_CELL_OF_DIMENSION_ARE_NEIGHBORS = 2


class Neighborhood:
    def __init__(self, neighbors_relative, edge_rule: EdgeRule):
        """ Defines a neighborhood of a cell.
        :param neighbors_relative: List of relative coordinates for cell neighbors.
        :param edge_rule: EdgeRule to define, how cells on the edge of the grid will be handled.
        """
        self._rel_neighbors = neighbors_relative
        self.__edge_rule = edge_rule
        self.__grid_dimensions = []

    def calculate_cell_neighbor_coordinates(self, cell_coordinate, grid_dimensions):
        """ Get a list of absolute coordinates for the cell neighbors.
            The EdgeRule can reduce the returned neighbor count.
        :param cell_coordinate:  The coordinate of the cell.
        :param grid_dimensions:  The dimensions of the grid, to apply the edge the rule.
        :return: list of absolute coordinates for the cells neighbors.
        """
        self.__grid_dimensions = grid_dimensions
        return list(self._neighbors_generator(cell_coordinate))

    def get_id_of_neighbor_from_relative_coordinate(self, rel_coordinate):
        return self._rel_neighbors.index(rel_coordinate)

    def _neighbors_generator(self, cell_coordinate):
        if not self._does_ignore_edge_cell_rule_apply(cell_coordinate):
            for rel_n in self._rel_neighbors:
                yield from self._calculate_abs_neighbor_and_decide_validity(cell_coordinate, rel_n)

    def _calculate_abs_neighbor_and_decide_validity(self, cell_coordinate, rel_n):
        n = list(map(operator.add, rel_n, cell_coordinate))
        n_folded = self.__apply_edge_overflow(n)
        if n == n_folded or self.__edge_rule == EdgeRule.FIRST_AND_LAST_CELL_OF_DIMENSION_ARE_NEIGHBORS:
            yield n_folded

    def _does_ignore_edge_cell_rule_apply(self, coordinate):
        return self.__edge_rule == EdgeRule.IGNORE_EDGE_CELLS and self.__is_coordinate_on_an_edge(coordinate)

    def __is_coordinate_on_an_edge(self, coordinate):
        return any(0 == ci or ci == di-1 for ci, di in zip(coordinate, self.__grid_dimensions))

    def __apply_edge_overflow(self, n):
        return list(map(lambda ni, di: (ni + di) % di, n, self.__grid_dimensions))


class MooreNeighborhood(Neighborhood):
    """ Moore defined a neighborhood with a radius applied on a the non euclidean distance to other cells in the grid.
        Example:
            2 dimensions
            C = cell of interest
            N = neighbor of cell
            X = no neighbor of cell

                  Radius 1                     Radius 2
               X  X  X  X  X                N  N  N  N  N
               X  N  N  N  X                N  N  N  N  N
               X  N  C  N  X                N  N  C  N  N
               X  N  N  N  X                N  N  N  N  N
               X  X  X  X  X                N  N  N  N  N
    """

    def __init__(self, edge_rule: EdgeRule = EdgeRule.IGNORE_EDGE_CELLS, radius=1, dimension=2):
        super().__init__(tuple(_rel_neighbor_generator(dimension, radius, lambda rel_n: True)),
                         edge_rule)


class VonNeumannNeighborhood(Neighborhood):
    """ Von Neumann defined a neighborhood with a radius applied to Manhatten distance
        (steps between cells without diagonal movement).
        Example:
            2 dimensions
            C = cell of interest
            N = neighbor of cell
            X = no neighbor of cell

                  Radius 1                     Radius 2
               X  X  X  X  X                X  X  N  X  X
               X  X  N  X  X                X  N  N  N  X
               X  N  C  N  X                N  N  C  N  N
               X  X  N  X  X                X  N  N  N  X
               X  X  X  X  X                X  X  N  X  X
    """

    def __init__(self, edge_rule: EdgeRule = EdgeRule.IGNORE_EDGE_CELLS, radius=1, dimension=2):
        self.radius = radius
        super().__init__(tuple(_rel_neighbor_generator(dimension, radius, self.neighbor_rule)),
                         edge_rule)

    def neighbor_rule(self, rel_n):
        cross_sum = 0
        for ci in rel_n:
            cross_sum += abs(ci)
        return cross_sum <= self.radius


class RadialNeighborhood(Neighborhood):
    """ Neighborhood with a radius applied to euclidean distance + delta

        Example:
            2 dimensions
            C = cell of interest
            N = neighbor of cell
            X = no neighbor of cell

                  Radius 2                     Radius 3
            X  X  X  X  X  X  X          X  X  N  N  N  X  X
            X  X  N  N  N  X  X          X  N  N  N  N  N  X
            X  N  N  N  N  N  X          N  N  N  N  N  N  N
            X  N  N  C  N  N  X          N  N  N  C  N  N  N
            X  N  N  N  N  N  X          N  N  N  N  N  N  N
            X  X  N  N  N  X  X          X  N  N  N  N  N  X
            X  X  X  X  X  X  X          X  X  N  N  N  X  X
    """

    def __init__(self, edge_rule: EdgeRule = EdgeRule.IGNORE_EDGE_CELLS, radius=1, delta_=.25, dimension=2):
        self.radius = radius
        self.delta = delta_
        super().__init__(tuple(_rel_neighbor_generator(dimension, radius, self.neighbor_rule)),
                         edge_rule)

    def neighbor_rule(self, rel_n):
        cross_sum = 0
        for ci in rel_n:
            cross_sum += pow(ci, 2)

        return math.sqrt(cross_sum) <= self.radius + self.delta


class HexagonalNeighborhood(Neighborhood):
    """ Defines a Hexagonal neighborhood in a rectangular two dimensional grid:

        Example:
            Von Nexagonal neighborhood in 2 dimensions with radius 1 and 2
            C = cell of interest
            N = neighbor of cell
            X = no neighbor of cell

                  Radius 1                     Radius 2
               X   X   X   X   X           X   N   N   N   X
                 X   N   N   X               N   N   N   N
               X   N   C   N   X           N   N   C   N   N
                 X   N   N   X               N   N   N   N
               X   X   X   X   X           X   N   N   N   X


        Rectangular representation: Radius 1

          Row % 2 == 0            Row % 2 == 1
            N  N  X                 X  N  N
            N  C  N                 N  C  N
            N  N  X                 X  N  N

        Rectangular representation: Radius 2
          Row % 2 == 0            Row % 2 == 1
          X  N  N  N  X           X  N  N  N  X
          N  N  N  N  X           X  N  N  N  N
          N  N  C  N  N           N  N  C  N  N
          N  N  N  N  X           X  N  N  N  N
          X  N  N  N  X           X  N  N  N  X
    """

    def __init__(self, edge_rule: EdgeRule = EdgeRule.IGNORE_EDGE_CELLS, radius=1):
        neighbor_lists = [[(0, 0)],
                          [(0, 0)]]

        self.__calculate_hexagonal_neighborhood(neighbor_lists, radius)

        super().__init__(neighbor_lists, edge_rule)

    def __calculate_hexagonal_neighborhood(self, neighbor_lists, radius):
        for r in range(1, radius + 1):
            for i, n in enumerate(neighbor_lists):
                n = _grow_neighbours(n)
                n = self.__add_rectangular_neighbours(n, r, i % 2 == 1)
                n = sorted(n, key=(lambda ne: [ne[1], ne[0]]))
                n.remove((0, 0))
                neighbor_lists[i] = n

    def get_id_of_neighbor_from_relative_coordinate(self, rel_coordinate):
        raise NotImplementedError

    def _neighbors_generator(self, cell_coordinate):
        if not self._does_ignore_edge_cell_rule_apply(cell_coordinate):
            for rel_n in self._rel_neighbors[cell_coordinate[1] % 2]:
                yield from self._calculate_abs_neighbor_and_decide_validity(cell_coordinate, rel_n)

    @staticmethod
    def __add_rectangular_neighbours(neighbours, radius, is_odd):
        new_neighbours = []
        for x in range(0, radius + 1):
            if is_odd:
                x -= int(radius/2)
            else:
                x -= int((radius + 1) / 2)

            for y in range(-radius, radius + 1):
                new_neighbours.append((x, y))
        new_neighbours.extend(neighbours)
        return list(set(new_neighbours))


def _rel_neighbor_generator(dimension, range_, rule):
    for c in itertools.product(range(-range_, range_ + 1), repeat=dimension):
        if rule(c) and c != (0, ) * dimension:
            yield tuple(reversed(c))


def _grow_neighbours(neighbours):
    new_neighbours = neighbours[:]
    for n in neighbours:
        new_neighbours.append((n[0], n[1] - 1))
        new_neighbours.append((n[0] - 1, n[1]))
        new_neighbours.append((n[0] + 1, n[1]))
        new_neighbours.append((n[0], n[1] + 1))
    return list(set(new_neighbours))
