"""A self-driving laboratory on a computer: closed-loop optimization of adsorption cycles with RUPTURA.

As in an autonomous laboratory, a robot runs the experiments, inline analytics measure them, and a brain chooses the
next batch from all results so far:

- psa.PressureSwingCycle is the robot and the analytics: a five-step pressure-swing cycle with product backfill on the
  column of a base input, run through ruptura.load_simulation(...).run() until the cyclic steady state, with purity,
  recovery, and productivity from SwingAdsorption.phase_records;
- optimizer.propose_batch is the brain: Gaussian processes, ParEGO, expected improvement times the probability of
  feasibility, and batches by Kriging believer;
- campaign.Campaign closes the loop, runs the experiments of a batch in parallel processes, and keeps the lab notebook.

    from ruptura.selfdriving import Campaign, PressureSwingCycle

    cycle = PressureSwingCycle("simulation.json", product="H2", purity_target=0.999)
    campaign = Campaign(cycle, "campaign", rounds=15, batch_size=8, workers=8)
    rows = campaign.run()

Needs scipy (pip install "ruptura[selfdriving]").
"""

from .campaign import Campaign
from .optimizer import GaussianProcess, hypervolume, pareto_front, propose_batch
from .psa import BASELINE, VARIABLES, PressureSwingCycle, key_performance, run_experiment

__all__ = ["BASELINE", "Campaign", "GaussianProcess", "PressureSwingCycle", "VARIABLES", "hypervolume",
           "key_performance", "pareto_front", "propose_batch", "run_experiment"]
