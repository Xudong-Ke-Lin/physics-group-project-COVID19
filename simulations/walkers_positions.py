import networkx as nx
from matplotlib import pyplot as plt
import random
from matplotlib.animation import FuncAnimation
import numpy as np
from matplotlib.patches import Rectangle



###################################################
# System dimensions
###################################################

Np = 20  # number of people

Nx = 20  # number of nodes in x-direction
Ny = 22  # number of nodes in y-direction

Lx = 20  # width of room (m)
Ly = 22  # height of room (m)
# reference for finite difference diffusion simulation
dx = Lx / Nx
dy = Ly / Ny

simuTime = 2000  # seconds
dt = 1  # seconds
simuSteps = round(simuTime / dt) # number of simulations steps


def dist(a, b):
    '''
     A function to evaluate the estimate of the distance from the a node to the target.
    The function takes two nodes arguments and must return a number.
    :param a: the posiiton of the start node
    :param b: the position of the target node
    :return: the distance between the nodes.
    '''
    (x1, y1) = a
    (x2, y2) = b
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2) ** 0.5


G = nx.grid_2d_graph(Ny, Nx)  # create 2D graph with  dimensions of (Ny, Nx)
disallowed_nodes = []  # store nodes where a walker is not allowed to go
for edge in G.edges:  # for each possible edge between nodes
    G[edge[0][0], edge[0][1]][edge[1][0], edge[1][1]]['weight'] = 0  # make each edge have zero weight


for l in range(4, 20, 4):  # spaces of 4 between benches
    for i in range(l, l + 2):  # height of bench is 2
        for k in range(0, 8):  # length of bench is 8
            for edge in G.edges((i, k)):

                # give weight of 10000000 for disallowed routes.
                G[edge[0][0], edge[0][1]][edge[1][0], edge[1][1]]['weight'] = 10000000
                # [0][0] is to select the first element of the first node that makes up the edge
                if not ((i, k) in disallowed_nodes):
                    disallowed_nodes.append((i, k))
                else:
                    continue
        # repeat process but for RHS benches
        for k in range(12, 20):
            for edge in G.edges((i, k)):
                G[edge[0][0], edge[0][1]][edge[1][0], edge[1][1]]['weight'] = 10000000
                if not ((i, k) in disallowed_nodes):
                    disallowed_nodes.append((i, k))
                else:
                    continue

# need to store nodes that you're allowed to traverse
allowed_nodes = list(G.nodes)  # use list to make copy and not a memory reference

# remove disallowed nodes from allowed node list
for node in disallowed_nodes:
    if node in allowed_nodes:
        allowed_nodes.remove(node)

start_nodes = [] # list to store the start node/ base position for a walker
target_nodes = [] # list to store the target / final position for a walker

for i in range(Np): # for each walker
    start_node = random.choice(allowed_nodes)  # randomise the initial position
    final_node = random.choice(allowed_nodes) # randomise the final position
    for edge in G.edges(start_node):
        G[edge[0][0], edge[0][1]][edge[1][0], edge[1][1]]['weight'] = 10000000

    start_nodes.append(start_node)
    target_nodes.append(final_node)

walkers_positions = np.ones((Np, simuSteps, 2)) # 1st index is person, 2nd index is time, 3rd is x coord, 4th is y coord
current_nodes = list(start_nodes) # current node for each walker - need to use a list or it copies a reference


random_starts = []  # create random start times for walkers
for i in range(0, Np):
    n = random.randint(1, int(simuSteps/4))
    random_starts.append(n)


counters = [0] * Np # store the amount of time steps elapsed since a walker reacher their taregt node
forward = [True] * Np # list to store whether a walker so be going to/from their target node


def gen_path(start, finish):
    '''
    Function to find the next node a walker should take to avoid people and get to their target node in the most optimal
    way.
    :param start: the start node of the walker
    :param finish:  the target node of the walker
    :return: the next node that the walker should take
    '''

    # other_nodes is all the nodes besides the current node
    other_nodes = list(current_nodes)
    other_nodes.pop(s)

    # make edges for those positions high weighting
    for node in other_nodes:  # for all nodes besides the walkers current node
        for edge in G.edges(node):

            # make the edges of these nodes untraversable so the walker cannot collide into another walker
            G[edge[0][0], edge[0][1]][edge[1][0], edge[1][1]]['weight'] = 10000000000 # used very large number

    # calculate optimal path
    optimal_path = nx.astar_path(G, start, finish, heuristic=dist,
                                 weight='weight')  # use a-star to find the opt. path
    new_node = optimal_path[1] # store the new node / the next node the walker should take
    return new_node # return this node


for k in range(0, simuSteps):  # for each step in the simulation
    for s in range(0, Np):  # for each walker

        if k >= random_starts[s]: # if some random time has passed for a walker then the walker should walk around

            if forward[s]: # if the walker should be walking to their final point

                if current_nodes[s] == target_nodes[s]: # and if their current position isn't the final point

                    if counters[s] > 40: # need to make random for each person and upon each time target node is reached
                        counters[s] = 0
                        forward[s] = False
                        new_node = gen_path(current_nodes[s], start_nodes[s])  # find the next node they should take
                    else:
                        counters[s] = counters[s] + 1
                        forward[s] = True
                        new_node = current_nodes[s]
                else:
                    new_node = gen_path(current_nodes[s], target_nodes[s])  # find the next node they should take
                    # look at all current positions besides our position

            else: # if the walker is going back to their seat / start position

                # compute reverse path, current_node -> people_start
                # move node to next point in reverse path using function
                if current_nodes[s] == start_nodes[s]:
                    forward[s] = True
                    new_node = gen_path(current_nodes[s], target_nodes[s])
                else:
                    forward[s] = False
                    new_node = gen_path(current_nodes[s], start_nodes[s])

            current_nodes[s] = new_node

            for node in current_nodes:
                # if nodes were originally in allowed nodes remove high weight (weight = 0)
                if node in allowed_nodes:
                    for edge in G.edges(node):
                        G[edge[0][0], edge[0][1]][edge[1][0], edge[1][1]]['weight'] = 0

        # add new walker position to numpy array
        walkers_positions[s, k, 0] = int(current_nodes[s][1]) # x coord
        walkers_positions[s, k, 1] = int(current_nodes[s][0]) # y coord


# fig, ax = plt.subplots()
#
# def update(num, walkers_positions, scatter):
#     '''
#
#     :param num: the current animation iteration number.
#     :param walkers_positions: the position tensor of the walkers.
#     :param scatter: the scatter plot object
#     :return: the updated scatter plot object
#     '''
#     scatter.set_offsets(np.c_[walkers_positions[:,num,0], walkers_positions[:,num,1]])
#     return scatter,
#
# colors = 255*np.random.rand(Np).tolist()  # list of colors use for scatter plot of walkers
#
# ax.set_title('Walkers in a public place')
# ax.set_xlabel('x [m]')
# ax.set_ylabel('y [m]')
# ax.set_xlim(0, Lx)
# ax.set_ylim(0, Ly)
# colors = np.random.rand(Np).tolist()  # list of colors use for scatter plot of walkers
# scatter = ax.scatter(walkers_positions[:,0,0], walkers_positions[:,0,1], c=colors)  # scatter seems to run slow compared to plot
# ax.add_patch(Rectangle((4, 0), 2, 8, facecolor='black'))
# ax.add_patch(Rectangle((8, 0), 2, 8, facecolor='black'))
# ax.add_patch(Rectangle((12, 0), 2, 8, facecolor='black'))
# ax.add_patch(Rectangle((16, 0), 2, 8, facecolor='black'))
# ax.add_patch(Rectangle((4,14), 2, 8, facecolor='black'))
# ax.add_patch(Rectangle((8,14), 2, 8, facecolor='black'))
# ax.add_patch(Rectangle((12,14), 2, 8, facecolor='black'))
# ax.add_patch(Rectangle((16,14), 2, 8, facecolor='black'))
#
# ani1 = FuncAnimation(fig, update, simuSteps, fargs=[walkers_positions, scatter],
#                               interval=1000, blit=True)
#
# aaa = bbb = ccc = ddd = eee = fff = ggg = 1
# ani1.save(f'walkers_{aaa}_{bbb}_{ccc}_{ddd}_{eee}_{fff}_{ggg}.mp4')
# print("Animation saved!")
