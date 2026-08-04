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
BLOCKS = [('CO2_1', 0.20, 1.00), ('CH4', 0.20, 1.00), ('COS', 0.20, 1.00),
          ('N2O_1', 0.20, 0.05), ('CO2_2', 0.20, 0.05), ('N2O_2', 0.20, 0.05)]

L = []
L.append('Transfer_function_parameters_(TFP)_for_IIR-shaped_filter_'
         '(see_Ibrom_et_al._2007_AFM).')
L.append('fc:_IIR_cut-off_frequency')
L.append('Fn:_normalization_parameter')
L.append('Water_vapour_TFP_are_calculated_for_9_RH_classes.')
L.append('Other_gases_TFP_are_calculated_on_a_monthly_base.')
L.append('-' * 84)
L.append('Water vapour TFP              Fn          fc    numerosity')
for r in RH:
    L.append('RH class %s%% = %11.5f %11.5f %12d' % (r, 0.20, 0.80, 500))
L.append('')
for name, fn, fc in BLOCKS:
    L.append('%s            TFP            Fn          fc' % name)
    for m in MONTHS:
        L.append('%-18s = %11.5f %11.5f ' % (m, fn, fc))
    L.append(' ')
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
cut = L.index('COS            TFP            Fn          fc')
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
for name, fn, _ in BLOCKS:
    G.append('%s            TFP            Fn          fc   groups=%s'
             % (name, GROUPING))
    for i, m in enumerate(MONTHS, start=1):
        fc = FC_GROUP1 if i <= 2 else FC_GROUP2
        G.append('%-18s = %11.5f %11.5f ' % (m, fn, fc))
    G.append(' ')
G += L[L.index('RH/fc_exponential_fit_parameters_for_water_vapour_'
               'spectral_corrections'):]
grp_path = os.path.join(here, 'sa_n_gas_2grp.txt')
with open(grp_path, 'w') as fh:
    fh.write('\n'.join(G) + '\n')


def variant(name, sa_file, months=None):
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
                    for i in range(1, ngas + 1)]
    with open(os.path.join(here, name), 'w') as fh:
        fh.write('\n'.join(out) + '\n')


variant('base_n_gas_sa.eddyflow', sa_path)
variant('base_n_gas_sa_short.eddyflow',
        os.path.join(here, 'sa_n_gas_short.txt'))
variant('base_n_gas_sa_2grp.eddyflow', grp_path, months=GROUPING)
print('wrote %d-line assessment file, %d-line short file, %d-line two-group '
      'file, 3 projects' % (len(L), len(short), len(G)))
