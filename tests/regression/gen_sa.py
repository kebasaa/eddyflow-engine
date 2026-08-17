"""Author a spectral assessment file covering every configured gas.

The engine can only *fit* transfer functions where the input data supports a
fit, and the regression fixture is two half-hours - too short for any gas,
including CO2. So the fitted path is exercised from the other end: a file that
declares fitted parameters is fed back in (sa_mode=0), which is the same reader
the on-the-fly assessment writes for. Gases 5+ are given a cut-off far from the
analytic one so a fitted correction cannot be mistaken for the analytic one.
"""
import os

MONTHS = ['January', 'February', 'March', 'April', 'May', 'June', 'July',
          'August', 'September', 'October', 'November', 'December']
RH = ['  5 - 15', ' 15 - 25', ' 25 - 35', ' 35 - 45', ' 45 - 55', ' 55 - 65',
      ' 65 - 75', ' 75 - 85', ' 85 - 95']
# Slot order for base_n_gas, with EVERY water slot excluded - they are
# classed by RH, and OutputSpectralAssessmentResults skips them by species.
#
# These names must be the ones SpectralGasNames generates, because the reader
# matches blocks by name. They were `CO2`, `N2O` and a stray `H2O_2` block:
# once repeated species started being numbered the first CO2 and the first N2O
# became CO2_1 and N2O_1, so those two blocks matched nothing, both gases were
# left unfitted, and the file read as short - which sent the whole run to the
# analytic fallback. Every fixture in this family was therefore exercising
# Moncrieff, not the fitted path it exists to cover.
#
# Each block also states the species and analyser it was fitted for. The name
# is an ordinal over repeats - CO2_1 is whichever CO2 record came first - so it
# cannot survive a project that lists its records in another order, and handing
# CO2_1's transfer function to the other analyser is not a small error when the
# two cells are 943 hPa and 70 hPa apart. The stamp goes past `fc`, where the
# reader stops, so an older build ignores it.
#
# N2O_1 and N2O_2 deliberately carry the SAME stamp: base_n_gas declares n2o on
# the MIRO twice. A stamp matching two records identifies neither, so those two
# fall back to the name match - which is the case the fallback exists for.
BLOCKS = [('CO2_1', 0.20, 1.00, 'co2', 'miro_mga4_6_2'),
          ('CH4', 0.20, 1.00, 'ch4', 'none'),
          ('COS', 0.20, 1.00, 'cos', 'miro_mga4_6_2'),
          ('N2O_1', 0.20, 0.05, 'n2o', 'miro_mga4_6_2'),
          ('CO2_2', 0.20, 0.05, 'co2', 'li7200_1'),
          ('N2O_2', 0.20, 0.05, 'n2o', 'miro_mga4_6_2')]


def stamp(var, instr):
    return '   var=%s instr=%s' % (var, instr)

L = []
L.append('Transfer_function_parameters_(TFP)_for_IIR-shaped_filter_'
         '(see_Ibrom_et_al._2007_AFM).')
L.append('fc:_IIR_cut-off_frequency')
L.append('Fn:_normalization_parameter')
L.append('Water_vapour_TFP_are_calculated_for_9_RH_classes.')
L.append('Other_gases_TFP_are_calculated_on_a_monthly_base.')
L.append('-' * 84)
L.append('Water vapour TFP              Fn          fc    numerosity'
         + stamp('h2o', 'miro_mga4_6_2'))
for r in RH:
    L.append('RH class %s%% = %11.5f %11.5f %12d' % (r, 0.20, 0.80, 500))
L.append('')
for name, fn, fc, var, instr in BLOCKS:
    L.append('%s            TFP            Fn          fc%s'
             % (name, stamp(var, instr)))
    for m in MONTHS:
        L.append('%-18s = %11.5f %11.5f ' % (m, fn, fc))
    L.append(' ')
# The second hygrometer's own RH table.
#
# base_n_gas declares two H2O records. The primary's table is the "Water
# vapour TFP" block at the top, in the fixed position the reader's seven-line
# skip expects; the second had nowhere to go, so it was fitted on every run
# and discarded. It is a named block like a gas's, and keeps `numerosity` in
# the header - that word is what tells the reader these are nine RH rows and
# not twelve monthly ones.
#
# The name is the bare tag. It used to read `H2O_2 vapour TFP`, and the reader
# takes everything before the word TFP as the name - so the block called itself
# `H2O_2 VAPOUR`, matched no tag, and was consumed and thrown away on every
# read. The block existed, the round trip did not.
#
# Given a cut-off far from the primary's 0.80 so that a second hygrometer
# reading its own row cannot be mistaken for one falling back to the first's.
#
# It also carries its own RH/cut-off coefficients. The standalone
# `RH/fc_exponential_fit_parameters` section below is the PRIMARY's, and the
# iir correction evaluates exp(A*RH^2 + B*RH + C) - so without a per-hygrometer
# triple the second hygrometer takes the primary's curve however far apart the
# two read, which is what made its nine fitted cut-offs decorative.
#
# -2,-1,-2 against the primary's 0.5,0.5,0.5 puts the two cut-offs about two
# orders of magnitude apart, so a hygrometer reading its own row cannot be
# mistaken for one falling back.
L.append('H2O_2            TFP              Fn          fc    numerosity'
         + stamp('h2o', 'li7200_1') + '   exp=-2.0,-1.0,-2.0')
for r in RH:
    L.append('RH class %s%% = %11.5f %11.5f %12d' % (r, 0.20, 0.05, 500))
L.append('')
L.append('RH/fc_exponential_fit_parameters_for_water_vapour_spectral_corrections')
L.append('-' * 35)
L.append('         exp1         exp2         exp3')
L.append('%13.6f%13.6f%13.6f' % (0.5, 0.5, 0.5))
L.append('')
L.append('')
L.append('High-pass_correction_factor_model_parameters')
L.append('Model: CF = [c1 * u / (c2 + f_co) + 1] after_Ibrom_et_al_(2007_AFM)')
L.append('-' * 76)
L.append('                   c1          c2')
L.append('unstable = %11.7f %11.7f ' % (0.5, 0.5))
L.append('stable   = %11.7f %11.7f ' % (0.5, 0.5))

here = os.path.dirname(os.path.abspath(__file__))
sa_path = os.path.join(here, 'sa_n_gas_fitted.txt')
with open(sa_path, 'w') as fh:
    fh.write('\n'.join(L) + '\n')

# A short file: the historical three blocks only, to prove the reader notices
# the count mismatch instead of consuming the next section as parameters.
cut = next(i for i, ln in enumerate(L) if ln.startswith('COS            TFP'))
short = L[:cut] + L[L.index(
    'RH/fc_exponential_fit_parameters_for_water_vapour_spectral_corrections'):]
with open(os.path.join(here, 'sa_n_gas_short.txt'), 'w') as fh:
    fh.write('\n'.join(short) + '\n')


# A file fitted under a two-group month grouping, and the project that
# declares it.
#
# The assessment file is keyed by month and RegPar by class, and the two only
# coincide while there is a single all-months group - which every other file
# here has. `1-2,3-12` breaks the coincidence for this fixture's own run
# month: June is in group 2, so the lookup asks for class 2, and reading the
# file straight back by month would answer with FEBRUARY's row - which is in
# group 1. The two groups therefore carry cut-offs an order of magnitude
# apart, so taking the wrong one cannot be mistaken for rounding.
GROUPING = '1-2,3-12'
FC_GROUP1 = 1.00
FC_GROUP2 = 0.05

G = list(L[:L.index('')])          # header + RH table, up to the blank line
G.append('')
for name, fn, _, var, instr in BLOCKS:
    G.append('%s            TFP            Fn          fc   groups=%s%s'
             % (name, GROUPING, stamp(var, instr)))
    for i, m in enumerate(MONTHS, start=1):
        fc = FC_GROUP1 if i <= 2 else FC_GROUP2
        G.append('%-18s = %11.5f %11.5f ' % (m, fn, fc))
    G.append(' ')
G += L[L.index('RH/fc_exponential_fit_parameters_for_water_vapour_'
               'spectral_corrections'):]
grp_path = os.path.join(here, 'sa_n_gas_2grp.txt')
with open(grp_path, 'w') as fh:
    fh.write('\n'.join(G) + '\n')


def variant(name, sa_file, months=None, only=None):
    src = open(os.path.join(here, 'base_n_gas.eddyflow')).read().splitlines()
    ngas = 0
    for ln in src:
        if ln.startswith('gas_num='):
            ngas = int(ln.split('=', 1)[1])
    out = []
    for ln in src:
        if ln.startswith('sa_mode='):
            out.append('sa_mode=0')
        elif ln.startswith('sa_file='):
            out.append('sa_file=' + sa_file.replace(os.sep, '/'))
        elif ln.startswith('hf_meth='):
            out.append('hf_meth=3')       # ibrom_07: consumes RegPar directly
        else:
            out.append(ln)
        # Every configured gas, not only the four that happen to carry an
        # sa_fmin record - the gases past the fourth are precisely the ones
        # the old three tables could not reach.
        if months and ln == '[FluxCorrection_SpectralAnalysis_General]':
            out += ['gas_%d_sa_months=%s' % (i, months)
                    for i in (only or range(1, ngas + 1))]
    with open(os.path.join(here, name), 'w') as fh:
        fh.write('\n'.join(out) + '\n')


variant('base_n_gas_sa.eddyflow', sa_path)
variant('base_n_gas_sa_short.eddyflow',
        os.path.join(here, 'sa_n_gas_short.txt'))
variant('base_n_gas_sa_2grp.eddyflow', grp_path, months=GROUPING)

# Twelve groups of one month, written in the bare-month form.
#
# Reads the same two-group file, so it is a cross-check rather than a case of
# its own: under `1,2,...,12` June is class 6 and takes month six's row, and
# under `1-2,3-12` June is class 2 and takes group two's - which the file
# gives the same cut-off. The two fixtures must therefore agree, by different
# routes through the class map.
#
# It is also the only fixture that reaches MaxGasClasses, and the only one
# that spells a group without a dash - a form the parser accepts and the
# interface never writes.
variant('base_n_gas_sa_permonth.eddyflow', grp_path,
        months=','.join(str(m) for m in range(1, 13)))

# A grouping the parser refuses, on record one only.
#
# Must be byte-identical to base_n_gas_sa: malformed is treated exactly as
# absent, which is one group over the calendar. `1-13` is out of range and
# `junk` is not a number, so both refusal paths are covered, and the run
# continues for the other seven gases rather than stopping.
variant('base_n_gas_sa_bad.eddyflow', sa_path, months='1-13,junk', only=[1])


# The same file, read by a project that lists its two CO2 records the other way
# round.
#
# This is the case the stamps exist for. Block names are ordinals over repeats
# of a species - CO2_1 is whichever CO2 record came first - so with names alone
# the MIRO's transfer function lands on the LI-7200 and vice versa the moment
# the records are re-ordered, which an ordinary re-save in the interface is
# enough to do. The two analysers' cells are 943 hPa and 70 hPa apart, so this
# is not a rounding-level error.
#
# The two CO2 blocks carry cut-offs an order of magnitude apart (1.00 and 0.05),
# so co2_1_scf and co2_2_scf must simply *swap* between this fixture and
# base_n_gas_sa. If they stay put, the block followed the position.
swapped = []
for ln in open(os.path.join(here, 'base_n_gas_sa.eddyflow')).read().splitlines():
    if ln.startswith('gas_1_instr='):
        swapped.append('gas_1_instr=li7200_1')
    elif ln.startswith('gas_1_col='):
        swapped.append('gas_1_col=21')
    elif ln.startswith('gas_6_instr='):
        swapped.append('gas_6_instr=miro_mga4_6_2')
    elif ln.startswith('gas_6_col='):
        swapped.append('gas_6_col=7')
    else:
        swapped.append(ln)
with open(os.path.join(here, 'base_n_gas_sa_swapped.eddyflow'), 'w') as fh:
    fh.write('\n'.join(swapped) + '\n')

print('wrote %d-line assessment file, %d-line short file, %d-line two-group '
      'file, 6 projects' % (len(L), len(short), len(G)))
