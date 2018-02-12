#!/usr/bin/env dials.python
from __future__ import absolute_import, division

import copy
import math

import libtbx.load_env
from libtbx.phil import parse
from libtbx.utils import Sorry
from cctbx import sgtbx
from dxtbx.serialize import dump, load

from dials.array_family import flex
from dials.util.options import OptionParser
from dials.util.options import flatten_experiments, flatten_reflections

help_message = '''
'''

# The phil scope
phil_scope = parse('''
''', process_includes=True)

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

  # Try to load the models and data
  if len(params.input.experiments) == 0:
    print "No Experiments found in the input"
    parser.print_help()
    return
  if len(params.input.reflections) == 0:
    print "No reflection data found in the input"
    parser.print_help()
    return
  try:
    assert len(params.input.reflections) == len(params.input.experiments)
  except AssertionError:
    raise Sorry("The number of input reflections files does not match the "
      "number of input experiments")

  from collections import OrderedDict
  expt_filenames = OrderedDict((e.filename, e.data) for e in params.input.experiments)
  refl_filenames = OrderedDict((r.filename, r.data) for r in params.input.reflections)

  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  reflections = reflections[0]

  data_manager = DataManager(experiments, reflections)

  experiments = data_manager.experiments
  reflections = data_manager.reflections

  # export reflections
  hklin = data_manager.export_mtz(filename='combined.mtz')

  # decide pointgroup
  hklout = 'sorted.mtz'
  pointgroup, reindex_op = decide_pointgroup(hklin, hklout)
  print pointgroup.info(), reindex_op

  # reindex to correct bravais setting
  data_manager.reindex(reindex_op, pointgroup)
  reindexed_expts = 'experiments_reindexed.json'
  reindexed_refl = 'reflections_reindexed.pickle'
  data_manager.export_experiments(reindexed_expts)
  data_manager.export_reflections(reindexed_refl)

  # refine in correct bravais setting
  refined_expts, refined_refl = refine(
    reindexed_expts, reindexed_refl)
  data_manager.experiments = load.experiment_list(
    refined_expts, check_format=False)
  data_manager.reflections = flex.reflection_table.from_pickle(refined_refl)

  #export_mtz(data_manager.experiments, data_manager.reflections)
  hklin = data_manager.export_mtz(filename='combined.mtz')

  # re-run pointless using above refined reflections
  hklout = 'sorted.mtz'
  pointgroup, reindex_op = decide_pointgroup(hklin, hklout)
  print pointgroup.info(), reindex_op

  assert data_manager.experiments[0].crystal.get_space_group() == pointgroup
  assert reindex_op.is_identity_op()

  # two-theta refinement to get best estimate of unit cell
  import xia2.Modules.Scaler.tools as tools
  best_unit_cell, best_unit_cell_esd = two_theta_refine(
    refined_expts, refined_refl)
  tools.patch_mtz_unit_cell(hklout, best_unit_cell)

  # scale data with aimless
  hklin = hklout
  hklout = 'scaled.mtz'
  scale(hklin, hklout)


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
      print expt.scan.get_batch_offset(), expt.scan.get_batch_range()

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

  def reindex(self, cb_op, space_group=None):
    import os
    from dxtbx.model import Crystal
    self._reflections['miller_index'] = cb_op.apply(self._reflections['miller_index'])

    if space_group is None:
      space_group = self._experiments[0].get_space_group()

    for expt in self._experiments:
      cryst_orig = copy.deepcopy(expt.crystal)
      cryst_reindexed = cryst_orig.change_basis(cb_op)
      a, b, c = cryst_reindexed.get_real_space_vectors()
      cryst_reindexed = Crystal(a, b, c, space_group=space_group)
      expt.crystal.update(cryst_reindexed)

  def export_reflections(self, filename):
    self._reflections.as_pickle(filename)

  def export_experiments(self, filename):
    dump.experiment_list(self._experiments, filename)

  def export_mtz(self, filename=None, params=None):
    if params is None:
      from dials.command_line.export import phil_scope as export_phil_scope
      params = export_phil_scope.extract()
    if filename is not None:
      params.mtz.hklout = filename

    from dials.util.export_mtz import export_mtz
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

def refine(experiments_filename, reflections_filename):
  from xia2.Wrappers.Dials.Refine import Refine
  from xia2.lib.bits import auto_logfiler
  refiner = Refine()
  auto_logfiler(refiner)
  refiner.set_experiments_filename(experiments_filename)
  refiner.set_indexed_filename(reflections_filename)
  refiner.run()
  return refiner.get_refined_experiments_filename(), refiner.get_refined_filename()

def two_theta_refine(experiments_filename, reflections_filename):
  from xia2.Wrappers.Dials.TwoThetaRefine import TwoThetaRefine
  from xia2.lib.bits import auto_logfiler
  tt_refiner = TwoThetaRefine()
  auto_logfiler(tt_refiner)
  tt_refiner.set_experiments([experiments_filename])
  tt_refiner.set_pickles([reflections_filename])
  tt_refiner.run()
  unit_cell = tt_refiner.get_unit_cell()
  unit_cell_esd = tt_refiner.get_unit_cell_esd()
  return unit_cell, unit_cell_esd

def decide_pointgroup(hklin, hklout):
  from xia2.Wrappers.CCP4.Pointless import Pointless
  from xia2.lib.bits import auto_logfiler
  pointless = Pointless()
  auto_logfiler(pointless)
  pointless.set_hklin(hklin)
  pointless.set_hklout(hklout)
  pointless.set_allow_out_of_sequence_files(allow=True)
  pointless.decide_pointgroup()
  possible = pointless.get_possible_lattices()
  pointgroup = pointless.get_pointgroup()
  reindex_op =  pointless.get_reindex_operator()
  probably_twinned = pointless.get_probably_twinned()

  space_group = sgtbx.space_group_info(symbol=str(pointgroup)).group()
  cb_op = sgtbx.change_of_basis_op(reindex_op)
  return space_group, cb_op

def scale(hklin, hklout):
  from xia2.Wrappers.CCP4.Aimless import Aimless
  from xia2.lib.bits import auto_logfiler
  from xia2.Handlers.Phil import PhilIndex
  PhilIndex.params.xia2.settings.multiprocessing.nproc = 1
  PhilIndex.params.ccp4.aimless.secondary.lmax = 0
  aimless = Aimless()
  auto_logfiler(aimless)
  aimless.set_surface_link(False) # multi-crystal
  aimless.set_hklin(hklin)
  aimless.set_hklout(hklout)
  aimless.set_surface_tie(PhilIndex.params.ccp4.aimless.surface_tie)
  spacing = 2
  secondary = 'secondary'
  lmax = PhilIndex.params.ccp4.aimless.secondary.lmax
  aimless.set_secondary(mode=secondary, lmax=lmax)
  aimless.set_spacing(spacing)
  aimless.scale()


if __name__ == "__main__":
  run()
