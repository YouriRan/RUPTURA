"""The brain of the self-driving laboratory: batch Bayesian optimization of two objectives under a constraint.

Each objective and the constraint have a Gaussian-process surrogate on the unit cube. A batch is filled one setup at a
time. ParEGO (Knowles, IEEE Trans. Evol. Comput. 10 (2006) 50) draws a random weighting of the two objectives and
scalarizes them; the next setup maximizes the expected improvement of that scalarization (Jones et al., J. Global Optim.
13 (1998) 455) times the probability that the constraint holds (Gardner et al., ICML 2014), or, before any setup meets
the constraint, that probability alone. The surrogates then take their own prediction at that setup as if it had been
measured (Kriging believer; Ginsbourger et al., 2010), so that the next setup of the batch goes elsewhere.
"""

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize
from scipy.special import ndtr
from scipy.stats import qmc


class GaussianProcess:
    """Zero-mean GP on standardized targets with an ARD squared-exponential kernel.

    The hyperparameters maximize the marginal likelihood times a log-normal prior on the length scales around half the
    unit cube, as BoTorch does with its Gamma prior: with a few tens of points in seven or eight dimensions, the
    likelihood alone can explain the data with one very short length scale and call the other variables irrelevant.
    """

    def __init__(self, inputs, targets):
        self.inputs = np.asarray(inputs, dtype=float)
        targets = np.asarray(targets, dtype=float)
        self.mean = targets.mean()
        self.scale = targets.std() if targets.std() > 0.0 else 1.0
        self.targets = (targets - self.mean) / self.scale
        dimensions = self.inputs.shape[1]
        # log length scales, log signal variance, log noise variance
        start = np.concatenate([np.full(dimensions, np.log(0.5)), [0.0, np.log(1.0e-3)]])
        bounds = [(np.log(0.03), np.log(5.0))] * dimensions + [(np.log(0.05), np.log(20.0)), (np.log(1.0e-6), np.log(0.1))]
        fits = [minimize(self._negative_log_likelihood, guess, method="L-BFGS-B", bounds=bounds)
                for guess in (start, start + np.r_[np.full(dimensions, 1.0), 0.0, 0.0])]
        self._set(min(fits, key=lambda fit: fit.fun).x)

    def _set(self, parameters):
        dimensions = self.inputs.shape[1]
        self.lengths = np.exp(parameters[:dimensions])
        self.variance = np.exp(parameters[dimensions])
        self.noise = np.exp(parameters[dimensions + 1])
        matrix = self._kernel(self.inputs)
        matrix[np.diag_indices_from(matrix)] += self.noise + 1.0e-10
        self.factor = cho_factor(matrix, lower=True)
        self.weights = cho_solve(self.factor, self.targets)

    def _kernel(self, points):
        difference = (points[:, None, :] - self.inputs[None, :, :]) / self.lengths
        return self.variance * np.exp(-0.5 * np.sum(difference * difference, axis=-1))

    def _negative_log_likelihood(self, parameters):
        try:
            self._set(parameters)
        except np.linalg.LinAlgError:
            return 1.0e10
        log_lengths = parameters[: self.inputs.shape[1]]
        return (0.5 * self.targets @ self.weights + np.sum(np.log(np.diag(self.factor[0])))
                + 0.5 * np.sum((log_lengths - np.log(0.5)) ** 2))

    def predict(self, points):
        """Mean and standard deviation at the points, in the units of the targets."""
        cross = self._kernel(np.atleast_2d(points))
        variance = self.variance - np.sum(cross * cho_solve(self.factor, cross.T).T, axis=1)
        return self.mean + self.scale * (cross @ self.weights), self.scale * np.sqrt(np.maximum(variance, 1.0e-12))

    def believe(self, point, value):
        """Adds an observation of value at point, with the same hyperparameters."""
        self.inputs = np.vstack([self.inputs, point])
        self.targets = np.append(self.targets, (value - self.mean) / self.scale)
        self._set(np.r_[np.log(self.lengths), np.log(self.variance), np.log(self.noise)])


def pareto_front(objectives):
    """Mask of the points that no other point beats in both objectives (both maximized)."""
    objectives = np.asarray(objectives)
    return np.array([not np.any(np.all(objectives >= point, axis=1) & np.any(objectives > point, axis=1))
                     for point in objectives], dtype=bool)


def hypervolume(objectives, reference):
    """Area that the points of a two-objective maximization dominate above the reference point."""
    points = sorted((tuple(point) for point in objectives if np.all(np.asarray(point) > reference)), reverse=True)
    area, highest = 0.0, reference[1]
    for first, second in points:
        if second > highest:
            area += (first - reference[0]) * (second - highest)
            highest = second
    return area


def propose_batch(inputs, objectives, constraint, target, batch_size, rng):
    """The next batch of setups on the unit cube.

    inputs: (n, d) setups of the finished runs on the unit cube; objectives: (n, 2), both maximized; constraint: (n,),
    which must reach target.
    """
    inputs = np.asarray(inputs, dtype=float)
    objectives = np.asarray(objectives, dtype=float)
    constraint = np.asarray(constraint, dtype=float)
    objective_models = [GaussianProcess(inputs, objectives[:, k]) for k in range(objectives.shape[1])]
    constraint_model = GaussianProcess(inputs, constraint)
    candidates = qmc.Sobol(inputs.shape[1], seed=rng).random(4096)
    batch = []
    for _ in range(batch_size):
        feasibility = lambda points: ndtr((lambda m, s: (m - target) / s)(*constraint_model.predict(points)))
        feasible = constraint >= target
        if feasible.any():
            low, high = objectives.min(axis=0), objectives.max(axis=0)
            shortfall = 1.0 - (objectives - low) / np.where(high > low, high - low, 1.0)
            weight = rng.uniform()
            weights = np.array([weight, 1.0 - weight])
            # ParEGO's augmented Chebyshev scalarization, negated so that larger is better.
            scalar = -(np.max(weights * shortfall, axis=1) + 0.05 * shortfall @ weights)
            model = GaussianProcess(inputs, scalar)
            incumbent = scalar[feasible].max()

            def acquisition(points):
                mean, deviation = model.predict(points)
                z = (mean - incumbent) / deviation
                expected = deviation * (z * ndtr(z) + np.exp(-0.5 * z * z) / np.sqrt(2.0 * np.pi))
                return expected * feasibility(points)
        else:
            acquisition = feasibility

        scores = acquisition(candidates)
        best, best_score = candidates[np.argmax(scores)], scores.max()
        for start in candidates[np.argsort(-scores)[:5]]:
            result = minimize(lambda x: -acquisition(x[None, :])[0], start, method="L-BFGS-B",
                              bounds=[(0.0, 1.0)] * inputs.shape[1])
            if -result.fun > best_score:
                best, best_score = result.x, -result.fun
        batch.append(best)

        believed = np.array([model_k.predict(best[None, :])[0][0] for model_k in objective_models])
        believed_constraint = constraint_model.predict(best[None, :])[0][0]
        for model_k, value in zip(objective_models, believed):
            model_k.believe(best, value)
        constraint_model.believe(best, believed_constraint)
        inputs = np.vstack([inputs, best])
        objectives = np.vstack([objectives, believed])
        constraint = np.append(constraint, believed_constraint)
    return np.array(batch)
