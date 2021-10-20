import networkx as nx
import random
import numpy as np
#import scipy as sc
import matplotlib.pyplot as plt

#import scipy.optimize
from sklearn import metrics
from scipy.optimize import curve_fit
from tqdm import tqdm

# create random graph
G = nx.bipartite.gnmk_random_graph(15, 10, 50, seed=123)
# get layout
top = nx.bipartite.sets(G)[0]
pos = nx.bipartite_layout(G, top)
# get adjacency matrix
A = nx.adjacency_matrix(G)
A = A.toarray()
# plot adjacency matrix
plt.title("Adjacency Matrix")
plt.imshow(A, cmap="Greys")
plt.show()
# plot graph visualisation
nx.draw(G, pos, with_labels=True)

# Remove 30% of the edges
proportion_edges = 0.3
# this is our test set
edge_subset = random.sample(G.edges(), int(proportion_edges * G.number_of_edges()))
# Create a copy of the graph and remove the edges
G_train = G.copy()
G_train.remove_edges_from(edge_subset)
# adjacency matrix
A_train = nx.adjacency_matrix(G_train)
A_train = A_train.toarray()
# visualise the train graph
plt.title("Train Graph")
nx.draw(G_train, pos, with_labels=True)
plt.show()


G_test = nx.Graph()
G_test.add_edges_from(edge_subset)
# visualise the test graph
plt.title("Test Graph")
nx.draw(G_test, pos, with_labels=True)
plt.show()

"""
Mean Average Precision
"""


def MAP(G_test, G_pred, thres=0):
    # calculate avePrecision for each node and its neighbors
    avePs = []

    # loop through every node
    for node in tqdm(G_test.nodes()):
        # get predicted edges sorted in ranking order
        rankedPredWeights = sorted(G_pred[node].items(), key=lambda x: -x[1]['weight'])
        # only include edges that exist i.e. predicted rank / weight > threshold
        rankedPred = filter(lambda x: x[1]['weight'] > thres, rankedPredWeights)
        # get the rank
        pred = [x[0] for x in rankedPred]
        # calculate rel (existence of predicted edge in the groundtruth/actual set of edges)
        # get groundtruth neighbors
        gt = set(G_test[node])
        rel = np.array([x in gt for x in pred])
        # calculate P accumulative average of precision
        predLength = len(pred)
        P = np.array([
            sum(rel[:i + 1]) / len(rel[:i + 1]) for i in range(predLength)
        ])
        # calculate aveP
        aveP = (rel @ P) / len(gt)
        # keep track of results
        avePs.append(aveP)
    MAPvalue = sum(avePs) / len(avePs)
    print("MAP: {}".format(MAPvalue))
    return MAPvalue


"""
Visualise Receiver Operating Charateristics Curve and Precision-Recall Curve
"""


def ROC_PRC(pred, G):
    y_score = [p[2] for p in pred]
    y_true = [G.has_edge(p[0], p[1]) for p in pred]

    fig, (ax1, ax2) = plt.subplots(1, 2)

    fpr, tpr, thresholds = metrics.precision_recall_curve(y_true, y_score)
    ax1.plot(fpr, tpr)
    ax1.set_title("Precision-Recall Curve")
    # ax1.set_xlabel("fpr")
    # ax1.set_ylabel("tpr")

    fpr, tpr, thresholds = metrics.roc_curve(y_true, y_score)
    ax2.plot(fpr, tpr)
    ax2.set_title("ROC Curve, AUC = {:.2f}".format(metrics.roc_auc_score(y_true, y_score)))
    # ax2.set_xlabel("fpr")
    # ax2.set_ylabel("tpr")

    plt.show()


"""
Kernel as Curve Fitting Problem
"""
# eigenvalue decomposition
V_train, U_train = np.linalg.eig(A_train)
# U.T * Atest * U
target_V = U_train.T @ A @ U_train
# take only the diagonals
target_V = np.diag(target_V)
# plot the pattern
plt.figure(figsize=(5, 3))
plt.xlabel("V_train")
plt.ylabel("V_test")
plt.scatter(V_train, target_V, c='b')
plt.show()


# odd path counting kernel function
class OddPathCountingKernel:
    def __init__(self):
        self.a1 = 0
        self.a3 = 0
        self.a5 = 0
        self.a7 = 0

    def func(self, V, a1, a3, a5, a7):
        return V * a1 + V ** 3 * a3 + V ** 5 * a5 + V ** 7 * a7

    def fit(self, V_train, target_V):
        # do curve fitting
        popt, pcov = curve_fit(self.func, V_train, target_V)
        self.a1, self.a3, self.a5, self.a7 = popt

    def pred(self, V_train):
        return self.func(V_train, self.a1, self.a3, self.a5, self.a7)


# sinh pseudokernel function
class SinhPseudokernel:
    def __init__(self):
        self.alpha = 0

        def func(self, V, alpha):
            return np.array([
                alpha * (np.exp(lamb) - np.exp(-lamb)) for lamb in V
            ])

    def fit(self, V_train, target_V):
        # do curve fitting
        popt, pcov = curve_fit(self.func, V_train, target_V)
        self.alpha, = popt

    def pred(self, V_train):
        return self.func(V_train, self.alpha)


# odd neumann pseudokernel function
class OddNeumannPseudokernel:
    def __init__(self):
        self.alpha = 0

    def func(self, V, alpha):
        return np.array([
            alpha * (1 / (1 - lamb) - 1 / (1 + lamb)) for lamb in V
        ])

    def fit(self, V_train, target_V):
        # do curve fitting
        popt, pcov = curve_fit(self.func, V_train, target_V)
        self.alpha, = popt

    def pred(self, V_train):
        return self.func(V_train, self.alpha)


# fit kernel function
for kernel in [OddPathCountingKernel(), SinhPseudokernel(), OddNeumannPseudokernel()]:

    print(kernel)

    # kernel = OddPathCountingKernel()
    kernel.fit(V_train, target_V)
    V_pred = kernel.pred(V_train)

    # assume our function is exponential V_train with alpha = 0.6
    plt.figure(figsize=(5, 3))
    plt.xlabel("V_train")
    plt.ylabel("V_test")
    plt.scatter(V_train, target_V, c='b', label="target_V")
    plt.scatter(V_train, V_pred, c='r', label="predicted_V")
    plt.legend()
    plt.show()

    # transformation
    Apred = U_train @ np.diag(V_pred) @ U_train.T
    Apred = Apred.real

    # make edges of prediction
    pred = [(i, j, Apred[i, j]) for i in range(Apred.shape[0]) for j in range(Apred.shape[1])]

    # create graph
    G_pred = nx.Graph()
    G_pred.add_weighted_edges_from(pred)

    # evaluate MAP and precision
    ROC_PRC(pred, G)
    MAP(G_test, G_pred)
