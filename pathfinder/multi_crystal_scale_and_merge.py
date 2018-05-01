#!/usr/bin/env dials.python
from __future__ import absolute_import, division, print_function

from collections import OrderedDict
import copy
import logging
import math
import os

import libtbx.load_env
from libtbx.utils import Sorry
import iotbx.phil
from iotbx.reflection_file_reader import any_reflection_file
from cctbx import crystal
from cctbx import sgtbx
from dxtbx.serialize import dump, load
from dxtbx.model import Crystal, ExperimentList

from dials.array_family import flex
from dials.algorithms.symmetry.cosym import analyse_datasets as cosym_analyse_datasets

from dials.command_line.export import phil_scope as export_phil_scope
from dials.util import log
from dials.util.options import OptionParser
from dials.util.options import flatten_experiments, flatten_reflections
from dials.util.export_mtz import export_mtz

from dials_research.multi_crystal_analysis import multi_crystal_analysis
from dials_research.multi_crystal_analysis import separate_unmerged

from xia2.lib.bits import auto_logfiler
from xia2.Handlers.Phil import PhilIndex
import xia2.Modules.Scaler.tools as tools
from xia2.Wrappers.CCP4.Aimless import Aimless
from xia2.Wrappers.CCP4.Pointless import Pointless
from xia2.Wrappers.Dials.Refine import Refine
from xia2.Wrappers.Dials.TwoThetaRefine import TwoThetaRefine

logger = logging.getLogger('dials.multi_crystal_scale_and_merge')

help_message = '''
'''

# The phil scope
phil_scope = iotbx.phil.parse('''
unit_cell_clustering {
  threshold = 5000
    .type = float(value_min=0)
    .help = 'Threshold value for the clustering'
  log = False
    .type = bool
    .help = 'Display the dendrogram with a log scale'
}

scaling
  .short_caption = "aimless"
{
  #intensities = summation profile *combine
    #.type = choice
  surface_tie = 0.001
    .type = float
    .short_caption = "Surface tie"
  surface_link = True
    .type = bool
    .short_caption = "Surface link"
  rotation.spacing = 2
    .type = int
    .expert_level = 2
    .short_caption = "Interval (in degrees) between scale factors on rotation axis"
  brotation.spacing = None
    .type = int
    .expert_level = 2
    .short_caption = "Interval (in degrees) between B-factors on rotation axis"
  secondary {
    frame = camera *crystal
      .type = choice
      .help = "Whether to do the secondary beam correction in the camera spindle"
              "frame or the crystal frame"
    lmax = 0
      .type = int
      .expert_level = 2
      .short_caption = "Aimless # secondary harmonics"
  }
}

symmetry {
  cosym {
    include scope dials.algorithms.symmetry.cosym.phil_scope
  }
}

resolution
  .short_caption = "Resolution"
{
  d_max = None
    .type = float(value_min=0.0)
    .help = "Low resolution cutoff."
    .short_caption = "Low resolution cutoff"
  d_min = None
    .type = float(value_min=0.0)
    .help = "High resolution cutoff."
    .short_caption = "High resolution cutoff"
  include scope dials.util.Resolutionizer.phil_str
}

multi_crystal_analysis {
  include scope dials_research.multi_crystal_analysis.master_phil_scope
}

identifiers = None
  .type = strings

''', process_includes=True)

# override default parameters
phil_scope = phil_scope.fetch(source=iotbx.phil.parse(
  """\
resolution {
  cc_half_method = sigma_tau
  cc_half_fit = tanh
}
"""))

def run():

  # The script usage
  usage  = "usage: %s [options] [param.phil] " \
           "experiments1.json experiments2.json reflections1.pickle " \
           "reflections2.pickle..." \
           % libtbx.env.dispatcher_name

  # Create the parser
  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_reflections=True,
    read_experiments=True,
    check_format=False,
    epilog=help_message)

  # Parse the command line
  params, options = parser.parse_args(show_diff_phil=True)

  # Configure the logging
  log.config(info='multi_crystal_scale_and_merge.log',
             #debug=params.output.debug_log
             )
  from dials.util.version import dials_version
  logger.info(dials_version())

  # Try to load the models and data
  if len(params.input.experiments) == 0:
    logger.info("No Experiments found in the input")
    parser.print_help()
    return
  if len(params.input.reflections) == 0:
    logger.info("No reflection data found in the input")
    parser.print_help()
    return
  try:
    assert len(params.input.reflections) == len(params.input.experiments)
  except AssertionError:
    raise Sorry("The number of input reflections files does not match the "
      "number of input experiments")

  expt_filenames = OrderedDict((e.filename, e.data) for e in params.input.experiments)
  refl_filenames = OrderedDict((r.filename, r.data) for r in params.input.reflections)

  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  reflections_all = flex.reflection_table()
  assert len(reflections) == len(experiments)
  for i, (expt, refl) in enumerate(zip(experiments, reflections)):
    expt.identifier = '%i' % i
    refl['identifier'] = flex.std_string(refl.size(), expt.identifier)
    refl['id'] = flex.int(refl.size(), i)
    #refl.experiment_identifiers()[i] = expt.identifier
    reflections_all.extend(refl)
    reflections_all.experiment_identifiers()[i] = expt.identifier

  assert reflections_all.are_experiment_identifiers_consistent(experiments)

  if params.identifiers is not None:
    identifiers = []
    for identifier in params.identifiers:
      identifiers.extend(identifier.split(','))
    params.identifiers = identifiers
  scaled = Scale(experiments, reflections_all, params)


class DataManager(object):
  def __init__(self, experiments, reflections):
    self._input_experiments = experiments
    self._input_reflections = reflections

    self._experiments = copy.deepcopy(experiments)
    self._reflections = copy.deepcopy(reflections)

    self._set_batches()

  def _set_batches(self):
    max_batches = max(e.scan.get_image_range()[1] for e in self._experiments)
    max_batches += 10 # allow some head room

    n = int(math.ceil(math.log10(max_batches)))

    for i, expt in enumerate(self._experiments):
      expt.scan.set_batch_offset(i * 10**n)
      logger.info("%s %s" % (expt.scan.get_batch_offset(), expt.scan.get_batch_range()))

  @property
  def experiments(self):
    return self._experiments

  @experiments.setter
  def experiments(self, experiments):
    self._experiments = experiments

  @property
  def reflections(self):
    return self._reflections

  @reflections.setter
  def reflections(self, reflections):
    self._reflections = reflections

  def select(self, experiment_identifiers):
    self._experiments = ExperimentList(
      [expt for expt in self._experiments
       if expt.identifier in experiment_identifiers])
    experiment_identifiers = self._experiments.identifiers()
    sel = flex.bool(len(self._reflections), False)
    for i_expt, identifier in enumerate(experiment_identifiers):
      sel_expt = self._reflections['identifier'] == identifier
      sel.set_selected(sel_expt, True)
      self._reflections['id'].set_selected(sel_expt, i_expt)
    self._reflections = self._reflections.select(sel)
    assert self.reflections.are_experiment_identifiers_consistent(
      self._experiments)

  def reflections_as_miller_arrays(self, intensity_key='intensity.sum.value'):
    from cctbx import crystal, miller
    variance_key = intensity_key.replace('.value', '.variance')
    assert intensity_key in self._reflections
    assert variance_key in self._reflections

    miller_arrays = []
    for expt in self._experiments:
      crystal_symmetry = crystal.symmetry(
        unit_cell=expt.crystal.get_unit_cell(),
        space_group=expt.crystal.get_space_group())
      sel = ((self._reflections.get_flags(self._reflections.flags.integrated_sum)
              & (self._reflections['identifier'] == expt.identifier)))
      assert sel.count(True) > 0
      refl = self._reflections.select(sel)
      data = refl[intensity_key]
      variances = refl[variance_key]
      # FIXME probably need to do some filtering of intensities similar to that
      # done in export_mtz
      miller_indices = refl['miller_index']
      assert variances.all_gt(0)
      sigmas = flex.sqrt(variances)

      miller_set = miller.set(crystal_symmetry, miller_indices, anomalous_flag=False)
      intensities = miller.array(miller_set, data=data, sigmas=sigmas)
      intensities.set_observation_type_xray_intensity()
      intensities.set_info(miller.array_info(
        source='DIALS',
        source_type='pickle'
      ))
      miller_arrays.append(intensities)
    return miller_arrays

  def reindex(self, cb_op=None, cb_ops=None, space_group=None):
    assert [cb_op, cb_ops].count(None) == 1

    if cb_op is not None:
      logger.info('Reindexing: %s' % cb_op)
      self._reflections['miller_index'] = cb_op.apply(self._reflections['miller_index'])

      for expt in self._experiments:
        cryst_reindexed = expt.crystal.change_basis(cb_op)
        if space_group is not None:
          cryst_reindexed.set_space_group(space_group)
        expt.crystal.update(cryst_reindexed)

    elif isinstance(cb_ops, dict):
      for cb_op, dataset_ids in cb_ops.iteritems():
        cb_op = sgtbx.change_of_basis_op(cb_op)

        for dataset_id in dataset_ids:
          expt = self._experiments[dataset_id]
          logger.info('Reindexing experiment %s: %s' % (
            expt.identifier, cb_op.as_xyz()))
          cryst_reindexed = expt.crystal.change_basis(cb_op)
          if space_group is not None:
            cryst_reindexed.set_space_group(space_group)
          expt.crystal.update(cryst_reindexed)
          sel = self._reflections['identifier'] == expt.identifier
          self._reflections['miller_index'].set_selected(sel, cb_op.apply(
            self._reflections['miller_index'].select(sel)))

    else:
      assert len(cb_ops) == len(self._experiments)
      for cb_op, expt in zip(cb_ops, self._experiments):
        logger.info('Reindexing experiment %s: %s' % (
          expt.identifier, cb_op.as_xyz()))
        cryst_reindexed = expt.crystal.change_basis(cb_op)
        if space_group is not None:
          cryst_reindexed.set_space_group(space_group)
        expt.crystal.update(cryst_reindexed)
        sel = self._reflections['identifier'] == expt.identifier
        self._reflections['miller_index'].set_selected(sel, cb_op.apply(
          self._reflections['miller_index'].select(sel)))

  def export_reflections(self, filename):
    self._reflections.as_pickle(filename)

  def export_experiments(self, filename):
    dump.experiment_list(self._experiments, filename)

  def export_mtz(self, filename=None, params=None):
    if params is None:
      params = export_phil_scope.extract()
    if filename is not None:
      params.mtz.hklout = filename

    m = export_mtz(
      self._reflections,
      self._experiments,
      params.mtz.hklout,
      include_partials=params.mtz.include_partials,
      keep_partials=params.mtz.keep_partials,
      scale_partials=params.mtz.scale_partials,
      min_isigi=params.mtz.min_isigi,
      force_static_model=params.mtz.force_static_model,
      filter_ice_rings=params.mtz.filter_ice_rings,
      ignore_profile_fitting=params.mtz.ignore_profile_fitting,
      apply_scales=params.mtz.apply_scales)
    m.show_summary()

    b1 = set(b.num() for b in m.batches())
    b2 = set(m.get_column('BATCH').extract_values().as_double().iround())
    assert len(b2.difference(b1)) == 0

    return params.mtz.hklout


class Scale(object):
  def __init__(self, experiments, reflections, params):

    self._data_manager = DataManager(experiments, reflections)
    self._params = params

    if self._params.identifiers is not None:
      self._data_manager.select(self._params.identifiers)

    experiments = self._data_manager.experiments
    reflections = self._data_manager.reflections

    self.unit_cell_clustering(plot_name='cluster_unit_cell_p1.png')

    self.cosym()

    # export reflections
    self._integrated_combined_mtz = self._data_manager.export_mtz(
      filename='integrated_combined.mtz')

    self.decide_space_group()

    self.two_theta_refine()

    self.unit_cell_clustering(plot_name='cluster_unit_cell_sg.png')

    self.scale()

    d_min, reason = self.estimate_resolution_limit()

    logger.info('Resolution limit: %.2f (%s)' % (d_min, reason))

    self.scale(d_min=d_min)

    self.multi_crystal_analysis()

  def unit_cell_clustering(self, plot_name=None):
    crystal_symmetries = []
    for expt in self._data_manager.experiments:
      crystal_symmetry = expt.crystal.get_crystal_symmetry(
        assert_is_compatible_unit_cell=False)
      crystal_symmetries.append(crystal_symmetry.niggli_cell())
    lattice_ids = [expt.identifier for expt in self._data_manager.experiments]
    from xfel.clustering.cluster import Cluster
    from xfel.clustering.cluster_groups import unit_cell_info
    ucs = Cluster.from_crystal_symmetries(crystal_symmetries, lattice_ids=lattice_ids)
    threshold = 1000
    if plot_name is not None:
      from matplotlib import pyplot as plt
      plt.figure("Andrews-Bernstein distance dendogram", figsize=(12, 8))
      ax = plt.gca()
    else:
      ax = None
    clusters, _ = ucs.ab_cluster(
      self._params.unit_cell_clustering.threshold,
      log=self._params.unit_cell_clustering.log,
      write_file_lists=False,
      schnell=False,
      doplot=(plot_name is not None),
      ax=ax
    )
    if plot_name is not None:
      plt.tight_layout()
      plt.savefig(plot_name)
      plt.clf()
    logger.info(unit_cell_info(clusters))
    largest_cluster = None
    largest_cluster_lattice_ids = None
    for cluster in clusters:
      cluster_lattice_ids = [m.lattice_id for m in cluster.members]
      if largest_cluster_lattice_ids is None:
        largest_cluster_lattice_ids = cluster_lattice_ids
      elif len(cluster_lattice_ids) > len(largest_cluster_lattice_ids):
        largest_cluster_lattice_ids = cluster_lattice_ids

    if len(largest_cluster_lattice_ids) < len(crystal_symmetries):
      logger.info(
        'Selecting subset of data sets for subsequent analysis: %s' %str(largest_cluster_lattice_ids))
      self._data_manager.select(largest_cluster_lattice_ids)
    else:
      logger.info('Using all data sets for subsequent analysis')

  def cosym(self):

    # per-dataset change of basis operator to ensure all consistent

    cb_op_best_min = None
    change_of_basis_ops = []
    for i, expt in enumerate(self._data_manager.experiments):
      crystal_symmetry = expt.crystal.get_crystal_symmetry()
      metric_subgroups = sgtbx.lattice_symmetry.metric_subgroups(
        crystal_symmetry, max_delta=5)
      subgroup = metric_subgroups.result_groups[0]
      cb_op_inp_best = subgroup['cb_op_inp_best']
      crystal_symmetry_best = crystal_symmetry.change_basis(cb_op_inp_best)
      if cb_op_best_min is None:
        cb_op_best_min = crystal_symmetry_best.change_of_basis_op_to_niggli_cell()
      cb_op_inp_min = cb_op_best_min * cb_op_inp_best
      change_of_basis_ops.append(cb_op_inp_min)
    self._data_manager.reindex(cb_ops=change_of_basis_ops)

    miller_arrays = self._data_manager.reflections_as_miller_arrays(
      intensity_key='intensity.sum.value')

    miller_arrays_p1 = []
    for ma in miller_arrays:
      cb_op_to_primitive = ma.change_of_basis_op_to_primitive_setting()
      ma = ma.change_basis(cb_op_to_primitive)
      space_group_info = sgtbx.space_group_info('P1')
      ma = ma.customized_copy(space_group_info=space_group_info)
      ma = ma.merge_equivalents().array()
      miller_arrays_p1.append(ma)

    params = self._params.symmetry.cosym

    result = cosym_analyse_datasets(miller_arrays_p1, params)

    space_groups = {}
    reindexing_ops = {}
    for dataset_id in result.reindexing_ops.iterkeys():
      if 0 in result.reindexing_ops[dataset_id]:
        cb_op = result.reindexing_ops[dataset_id][0]
        reindexing_ops.setdefault(cb_op, [])
        reindexing_ops[cb_op].append(dataset_id)
      if dataset_id in result.space_groups:
        space_groups.setdefault(result.space_groups[dataset_id], [])
        space_groups[result.space_groups[dataset_id]].append(dataset_id)

    logger.info('Space groups:')
    for sg, datasets in space_groups.iteritems():
      logger.info(str(sg.info().reference_setting()))
      logger.info(datasets)

    logger.info('Reindexing operators:')
    for cb_op, datasets in reindexing_ops.iteritems():
      logger.info(cb_op)
      logger.info(datasets)

    self._data_manager.reindex(cb_ops=reindexing_ops)

    return

  def decide_space_group(self):
    # decide space group
    self._sorted_mtz = 'sorted.mtz'
    space_group, reindex_op = self._decide_space_group_pointless(
      self._integrated_combined_mtz, self._sorted_mtz)

    # reindex to correct bravais setting
    self._data_manager.reindex(cb_op=reindex_op, space_group=space_group)
    self._experiments_filename = 'experiments_reindexed.json'
    self._reflections_filename = 'reflections_reindexed.pickle'
    self._data_manager.export_experiments(self._experiments_filename)
    self._data_manager.export_reflections(self._reflections_filename)

  def refine(self):
    # refine in correct bravais setting
    self._experiments_filename, self._reflections_filename = self._dials_refine(
      self._experiments_filename, self._reflections_filename)
    self._data_manager.experiments = load.experiment_list(
      self._experiments_filename, check_format=False)
    self._data_manager.reflections = flex.reflection_table.from_pickle(
      self._reflections_filename)

  def two_theta_refine(self):
    # two-theta refinement to get best estimate of unit cell
    self.best_unit_cell, self.best_unit_cell_esd = self._dials_two_theta_refine(
      self._experiments_filename, self._reflections_filename)
    tools.patch_mtz_unit_cell(self._sorted_mtz, self.best_unit_cell)

  def scale(self, d_min=None):

    # scale data with aimless
    aimless = self._aimless_scale(
      self._sorted_mtz, self._params.scaling, d_min=d_min)
    self._scaled_mtz = aimless.get_hklout()
    self._scaled_unmerged_mtz \
      = os.path.splitext(self._scaled_mtz)[0] + '_unmerged.mtz'

  @staticmethod
  def _decide_space_group_pointless(hklin, hklout):
    pointless = Pointless()
    auto_logfiler(pointless)
    pointless.set_hklin(hklin)
    pointless.set_hklout(hklout)
    pointless.set_allow_out_of_sequence_files(allow=True)
    pointless.decide_pointgroup()
    possible = pointless.get_possible_lattices()
    space_group = sgtbx.space_group_info(
      symbol=str(pointless.get_pointgroup())).group()
    cb_op =  sgtbx.change_of_basis_op(pointless.get_reindex_operator())
    #probably_twinned = pointless.get_probably_twinned()
    return space_group, cb_op

  @staticmethod
  def _dials_refine(experiments_filename, reflections_filename):
    refiner = Refine()
    auto_logfiler(refiner)
    refiner.set_experiments_filename(experiments_filename)
    refiner.set_indexed_filename(reflections_filename)
    refiner.run()
    return refiner.get_refined_experiments_filename(), refiner.get_refined_filename()

  @staticmethod
  def _dials_two_theta_refine(experiments_filename, reflections_filename):
    tt_refiner = TwoThetaRefine()
    auto_logfiler(tt_refiner)
    tt_refiner.set_experiments([experiments_filename])
    tt_refiner.set_pickles([reflections_filename])
    tt_refiner.run()
    unit_cell = tt_refiner.get_unit_cell()
    unit_cell_esd = tt_refiner.get_unit_cell_esd()
    return unit_cell, unit_cell_esd

  @staticmethod
  def _aimless_scale(hklin, params, d_min=None):
    PhilIndex.params.xia2.settings.multiprocessing.nproc = 1
    aimless = Aimless()
    auto_logfiler(aimless)
    hklout = '%i_scaled.mtz' % aimless.get_xpid()
    aimless.set_surface_link(False) # multi-crystal
    aimless.set_hklin(hklin)
    aimless.set_hklout(hklout)
    aimless.set_surface_tie(params.surface_tie)
    if params.secondary.frame == 'camera':
      secondary = 'secondary'
    else:
      secondary = 'absorption'
    lmax = params.secondary.lmax
    aimless.set_secondary(mode=secondary, lmax=lmax)
    aimless.set_spacing(params.rotation.spacing)
    if d_min is not None:
      aimless.set_resolution(d_min)
    aimless.scale()
    return aimless

  def estimate_resolution_limit(self):
    # see also xia2/Modules/Scaler/CommonScaler.py: CommonScaler._estimate_resolution_limit()
    from xia2.Wrappers.XIA.Merger import Merger
    params = self._params.resolution
    m = Merger()
    auto_logfiler(m)
    m.set_hklin(self._scaled_unmerged_mtz)
    m.set_limit_rmerge(params.rmerge)
    m.set_limit_completeness(params.completeness)
    m.set_limit_cc_half(params.cc_half)
    m.set_cc_half_significance_level(params.cc_half_significance_level)
    m.set_limit_isigma(params.isigma)
    m.set_limit_misigma(params.misigma)
    #if batch_range is not None:
      #start, end = batch_range
      #m.set_batch_range(start, end)
    m.run()

    resolution_limits = []
    reasoning = []

    if params.completeness is not None:
      r_comp = m.get_resolution_completeness()
      resolution_limits.append(r_comp)
      reasoning.append('completeness > %s' % params.completeness)

    if params.cc_half is not None:
      r_cc_half = m.get_resolution_cc_half()
      resolution_limits.append(r_cc_half)
      reasoning.append('cc_half > %s' % params.cc_half)

    if params.rmerge is not None:
      r_rm = m.get_resolution_rmerge()
      resolution_limits.append(r_rm)
      reasoning.append('rmerge > %s' % params.rmerge)

    if params.isigma is not None:
      r_uis = m.get_resolution_isigma()
      resolution_limits.append(r_uis)
      reasoning.append('unmerged <I/sigI> > %s' % params.isigma)

    if params.misigma is not None:
      r_mis = m.get_resolution_misigma()
      resolution_limits.append(r_mis)
      reasoning.append('merged <I/sigI> > %s' % params.misigma)

    if len(resolution_limits):
      resolution = max(resolution_limits)
      reasoning = [
          reason for limit, reason in zip(resolution_limits, reasoning)
          if limit >= resolution
      ]
      reasoning = ', '.join(reasoning)
    else:
      resolution = 0.0
      reasoning = None

    return resolution, reasoning

  def multi_crystal_analysis(self):

    result = any_reflection_file(self._scaled_unmerged_mtz)
    intensities = None
    batches = None

    for ma in result.as_miller_arrays(
      merge_equivalents=False, crystal_symmetry=None):
      if ma.info().labels == ['I(+)', 'SIGI(+)', 'I(-)', 'SIGI(-)']:
        assert ma.anomalous_flag()
        intensities = ma
      elif ma.info().labels == ['I', 'SIGI']:
        assert not ma.anomalous_flag()
        intensities = ma
      elif ma.info().labels == ['BATCH']:
        batches = ma

    assert batches is not None
    assert intensities is not None

    separate = separate_unmerged(
      intensities, batches)

    self._params.multi_crystal_analysis
    mca = multi_crystal_analysis(
      separate.intensities.values(),
      labels=separate.intensities.keys(),
    )

    logger.info('\nIntensity correlation clustering summary:')
    for cluster in mca.cc_clusters:
      logger.info(cluster)

    logger.info('\nCos angle clustering summary:')
    for cluster in mca.cos_angle_clusters:
      logger.info(cluster)

    return mca


if __name__ == "__main__":
  run()
