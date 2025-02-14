%matplotlib inline
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import deeptime
import pyemma
from pyemma.util.contexts import settings

#Loading data
pdb = 'chainP.pdb'
files = 'kras.dcd'
#files = []

#Feature selection
torsions_feat = pyemma.coordinates.featurizer(pdb)
torsions_feat.add_backbone_torsions(cossin=True, periodic=False)
torsions_data = pyemma.coordinates.load(files, features=torsions_feat)
labels = ['backbone\ntorsions']
#cossin=True means each angle will be returned as a pair of (sin(x), cos(x)). This is useful, if you calculate the mean (e.g TICA/PCA, clustering) in that space.
ED_feat = pyemma.coordinates.featurizer(pdb)
ED_feat.add_backbone_torsions(selstr="resid 1 to 86", cossin=True, periodic=False)
ED_data = pyemma.coordinates.load(files, features=ED_feat)
labels += ['ED\ntorsions']
sidechain_feat = pyemma.coordinates.featurizer(pdb)
sidechain_feat.add_sidechain_torsions(selstr="resid 1 to 86", cossin=True, periodic=False)
sidechain_data = pyemma.coordinates.load(files, features=sidechain_feat)
labels += ['sidechain\ntorsions']

def score_cv(data, dim, lag, number_of_splits=10, validation_fraction=0.5):
    """Compute a cross-validated VAMP2 score.

    We randomly split the list of independent trajectories into
    a training and a validation set, compute the VAMP2 score,
    and repeat this process several times.

    Parameters
    ----------
    data : list of numpy.ndarrays
        The input data.
    dim : int
        Number of processes to score; equivalent to the dimension
        after projecting the data with VAMP2.
    lag : int
        Lag time for the VAMP2 scoring.
    number_of_splits : int, optional, default=10
        How often do we repeat the splitting and score calculation.
    validation_fraction : int, optional, default=0.5
        Fraction of trajectories which should go into the validation
        set during a split.
    """
    # we temporarily suppress very short-lived progress bars
    with pyemma.util.contexts.settings(show_progress_bars=False):
        nval = int(len(data) * validation_fraction)
        scores = np.zeros(number_of_splits)
        for n in range(number_of_splits):
            ival = np.random.choice(len(data), size=nval, replace=False)
            vamp = pyemma.coordinates.vamp(
                [d for i, d in enumerate(data) if i not in ival], lag=lag, dim=dim)
            scores[n] = vamp.score([d for i, d in enumerate(data) if i in ival])
    return scores


dim = 10

fig, axes = plt.subplots(1, 3, figsize=(12, 3), sharey=True)
for ax, lag in zip(axes.flat, [5, 10, 20]):
    torsions_scores = score_cv(torsions_data, lag=lag, dim=dim)
    scores = [torsions_scores.mean()]
    errors = [torsions_scores.std()]
    ED_scores = score_cv(ED_data, lag=lag, dim=dim)
    scores += [ED_scores.mean()]
    errors += [ED_scores.std()]
    sidechain_scores = score_cv(sidechain_data, lag=lag, dim=dim)
    scores += [sidechain_scores.mean()]
    errors += [sidechain_scores.std()]
    ax.bar(labels, scores, yerr=errors, color=['C0', 'C1', 'C2'])
    ax.set_title(r'lag time $\tau$={:.1f}ns'.format(lag * 0.1))
    if lag == 80:
        # save for later
        vamp_bars_plot = dict(
            labels=labels, scores=scores, errors=errors, dim=dim, lag=lag)
axes[0].set_ylabel('VAMP2 score')
plt.savefig('feature_selection.png', dpi=300)
plt.close()

#detailed VAMP2 score analysis for selected feature
##sometimes the feature with lower score can have better performance
lags = [1, 5, 20, 80, 100, 150, 200]
dims = [i + 1 for i in range(10)]

fig, ax = plt.subplots()
for i, lag in enumerate(lags):
    scores_ = np.array([score_cv(sidechain_data, dim, lag)
                        for dim in dims])
    scores = np.mean(scores_, axis=1)
    errors = np.std(scores_, axis=1, ddof=1)
    color = 'C{}'.format(i)
    ax.fill_between(dims, scores - errors, scores + errors, alpha=0.3, facecolor=color)
    ax.plot(dims, scores, '--o', color=color, label='lag={:.1f}ns'.format(lag * 0.1))
ax.legend()
ax.set_xlabel('number of dimensions')
ax.set_ylabel('VAMP2 score')
plt.savefig('dimension_selection.png', dpi=300)
plt.close()


#Coordinate transform and discretization using TICA
tica = pyemma.coordinates.tica(sidechain_data, dim=8, lg=160)
tica_output = tica.get_output()
tica_concatenated = np.concatenate(tica_output)

#visualize the marginal and joint distributions of our TICA components by simple histograming
fig, axes = plt.subplots(1, 6, figsize=(30, 4))
pyemma.plots.plot_feature_histograms(
    tica_concatenated,
    ax=axes[0],
    feature_labels=['IC1', 'IC2', 'IC3', 'IC4', 'IC5', 'IC6', 'IC7','IC8'],
    ylog=True)
pyemma.plots.plot_density(*tica_concatenated[:, :2].T, ax=axes[1], logscale=True)
axes[1].set_xlabel('IC 1')
axes[1].set_ylabel('IC 2')
pyemma.plots.plot_density(*tica_concatenated[:, [2,3]].T, ax=axes[3], logscale=True)
axes[3].set_xlabel('IC 3')
axes[3].set_ylabel('IC 4')
pyemma.plots.plot_density(*tica_concatenated[:, [0,2]].T, ax=axes[2], logscale=True)
axes[2].set_xlabel('IC 1')
axes[2].set_ylabel('IC 3')
pyemma.plots.plot_density(*tica_concatenated[:, [1,2]].T, ax=axes[4], logscale=True)
axes[4].set_xlabel('IC 2')
axes[4].set_ylabel('IC 3')
pyemma.plots.plot_density(*tica_concatenated[:, [4,5]].T, ax=axes[5], logscale=True)
axes[5].set_xlabel('IC 5')
axes[5].set_ylabel('IC 6')
plt.savefig('dimension_selection-for-constructing_msm.png', dpi=300)
plt.close()

#Discretization
n_clustercenters = [5, 10, 30, 75, 200, 450]

scores = np.zeros((len(n_clustercenters), 5))
for n, k in enumerate(n_clustercenters):
    for m in range(5):
        with pyemma.util.contexts.settings(show_progress_bars=False):
            _cl = pyemma.coordinates.cluster_kmeans(
                tica_output, k=k, max_iter=50, stride=50)
            _msm = pyemma.msm.estimate_markov_model(_cl.dtrajs, 5)
            scores[n, m] = _msm.score_cv(
                _cl.dtrajs, n=1, score_method='VAMP2', score_k=min(10, k))

fig, ax = plt.subplots()
lower, upper = pyemma.util.statistics.confidence_interval(scores.T.tolist(), conf=0.9)
ax.fill_between(n_clustercenters, lower, upper, alpha=0.3)
ax.plot(n_clustercenters, np.mean(scores, axis=1), '-o')
ax.semilogx()
ax.set_xlabel('number of cluster centers')
ax.set_ylabel('VAMP-2 score')
plt.savefig('discretization.png', dpi=300)
plt.close()

#number of states k set to be 75
cluster = pyemma.coordinates.cluster_kmeans(
    tica_output, k=75, max_iter=50, stride=10, fixed_seed=1)
dtrajs_concatenated = np.concatenate(cluster.dtrajs)
fig, ax = plt.subplots(figsize=(4, 4))
pyemma.plots.plot_density(
    *tica_concatenated[:, :2].T, ax=ax, cbar=False, alpha=0.3)
ax.scatter(*cluster.clustercenters[:, :2].T, s=5, c='C1')
ax.set_xlabel('IC 1')
ax.set_ylabel('IC 2')
plt.savefig('number-of_states.png', dpi=300)
plt.close()

#MSM estimation and validation
its = pyemma.msm.its(cluster.dtrajs, lags=50, nits=10, errors='bayes')
pyemma.plots.plot_implied_timescales(its, units='ns', dt=0.1);
plt.savefig("implied_timescales.png", dpi=300)
plt.close()

#selected the lagtime=40*0.1 ns, construct the msm
msm = pyemma.msm.bayesian_markov_model(cluster.dtrajs, lag=40, dt_traj='0.1 ns')
print('fraction of states used = {:.2f}'.format(msm.active_state_fraction))
print('fraction of counts used = {:.2f}'.format(msm.active_count_fraction))

#ck-test
nstates = 5
cktest = msm.cktest(nstates, mlags=6)
pyemma.plots.plot_cktest(cktest, dt=0.1, units='ns');
plt.savefig("ck-test.png", dpi=300)
plt.close()

#MSM spectral analysis
def its_separation_err(ts, ts_err):
    """
    Error propagation from ITS standard deviation to timescale separation.
    """
    return ts[:-1] / ts[1:] * np.sqrt(
        (ts_err[:-1] / ts[:-1])**2 + (ts_err[1:] / ts[1:])**2)


nits = 15

timescales_mean = msm.sample_mean('timescales', k=nits)
timescales_std = msm.sample_std('timescales', k=nits)

fig, axes = plt.subplots(1, 2, figsize=(10, 4))


axes[0].errorbar(
    range(1, nits + 1),
    timescales_mean,
    yerr=timescales_std,
    fmt='.', markersize=10)
axes[1].errorbar(
    range(1, nits),
    timescales_mean[:-1] / timescales_mean[1:],
    yerr=its_separation_err(
        timescales_mean,
        timescales_std),
    fmt='.',
    markersize=10,
    color='C0')

for i, ax in enumerate(axes):
    ax.set_xticks(range(1, nits + 1))
    ax.grid(True, axis='x', linestyle=':')

axes[0].axhline(msm.lag * 0.1, lw=1.5, color='k')
axes[0].axhspan(0, msm.lag * 0.1, alpha=0.3, color='k')
axes[0].set_xlabel('implied timescale index')
axes[0].set_ylabel('implied timescales / ns')
axes[1].set_xticks(range(1, nits))
axes[1].set_xticklabels(
    ["{:d}/{:d}".format(k, k + 1) for k in range(1, nits + 2)],
    rotation=45)
axes[1].set_xlabel('implied timescale indices')
axes[1].set_ylabel('timescale separation')
fig.tight_layout()

#analyzing the stationary distribution and the free energy computed over the first two TICA coordinates.
fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, sharey=True)
pyemma.plots.plot_contour(
    *tica_concatenated[:, [1,2]].T,
    msm.pi[dtrajs_concatenated],
    ax=axes[0],
    mask=True,
    cbar_label='stationary distribution')
pyemma.plots.plot_free_energy(
    *tica_concatenated[:, [1,2]].T,
    weights=np.concatenate(msm.trajectory_weights()),
    ax=axes[1],
    legacy=False)
for ax in axes.flat:
    ax.set_xlabel('IC 2')
axes[0].set_ylabel('IC 3')
axes[0].set_title('Stationary distribution', fontweight='bold')
axes[1].set_title('Reweighted free energy surface', fontweight='bold')
plt.savefig("distribution_reweighted_fes.png", dpi=300)
plt.close()

#Perron cluster cluster analysis
msm.pcca(nstates)
fig, axes = plt.subplots(1, 5, figsize=(15, 3), sharex=True, sharey=True)
for i, ax in enumerate(axes.flat):
    pyemma.plots.plot_contour(
        *tica_concatenated[:, [1,2]].T,
        msm.metastable_distributions[i][dtrajs_concatenated],
        ax=ax,
        cmap='afmhot_r',
        mask=True,
        cbar_label='metastable distribution {}'.format(i + 1))
    ax.set_xlabel('IC 2')
axes[0].set_ylabel('IC 3')
plt.savefig("distribution_reweighted_fes.png", dpi=300)
plt.close()

#metastable assignments
metastable_traj = msm.metastable_assignments[dtrajs_concatenated]

fig, ax = plt.subplots(figsize=(5, 4))
_, _, misc = pyemma.plots.plot_state_map(
    *tica_concatenated[:, [1,2]].T, metastable_traj, ax=ax)
ax.set_xlabel('IC 2')
ax.set_ylabel('IC 3')
misc['cbar'].set_ticklabels([r'$\mathcal{S}_%d$' % (i + 1)
                             for i in range(nstates)])
plt.savefig("metastable_assignments.png", dpi=300)
plt.close()

#generate a number of representative sample structures for each macrostate
pcca_samples = msm.sample_by_distributions(msm.metastable_distributions, 10)
torsions_source = pyemma.coordinates.source(files, features=torsions_feat)
pyemma.coordinates.save_trajs(
    torsions_source,
    pcca_samples,
    outfiles=['./new-try/pcca{}_10samples.pdb'.format(n + 1)
              for n in range(msm.n_metastable)])

# calculate stationary distribution which encodes the free energy of the states.
print('state\tπ\t\tG/kT')
for i, s in enumerate(msm.metastable_sets):
    print('x_{} = {:f}'.format(i + 1, msm.pi[s].sum()))
    w.writelines('x_{} = {:f}'.format(i + 1, msm.pi[s].sum()))
    w.writelines('\n')
    prop_list.append('{:f}'.format(msm.pi[s].sum()))
w.close()

#extract mean first passage times (MFPTs) between PCCA++ metastable states
from itertools import product
from pandas import DataFrame
# Compute MFPT matrix
mfpt = np.zeros((nstates, nstates))
for i, j in product(range(nstates), repeat=2):
    mfpt[i, j] = msm.mfpt(
        msm.metastable_sets[i],
        msm.metastable_sets[j])

# Save as .dat file (space-separated values)
np.savetxt("mfpt.dat", mfpt, fmt="%.2f")

#Transition between states
A = msm.metastable_sets[0]
B = np.concatenate(msm.metastable_sets[1:])
print('MFPT 1 -> other: ({:6.1f} ± {:5.1f}) ns'.format(
    msm.sample_mean('mfpt', A, B), msm.sample_std('mfpt', A, B)))
print('MFPT other -> 1: ({:.1f} ± {:5.1f}) ns'.format(
    msm.sample_mean('mfpt', B, A), msm.sample_std('mfpt', B, A)))

a = msm.metastable_sets[4]
b = np.concatenate([msm.metastable_sets[i] for i in [0, 1, 2, 3]])
print('MFPT 5 -> other: ({:6.1f} ± {:5.1f}) ns'.format(
    msm.sample_mean('mfpt', a, b), msm.sample_std('mfpt', a, b)))
print('MFPT other -> 5: ({:.1f} ± {:5.1f}) ns'.format(
    msm.sample_mean('mfpt', b, a), msm.sample_std('mfpt', b, a)))
