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
# Slot order for base_n_gas, water's own slot excluded (it has the RH table).
BLOCKS = [('CO2', 0.20, 1.00), ('CH4', 0.20, 1.00), ('COS', 0.20, 1.00),
          ('N2O', 0.20, 0.05), ('CO2_2', 0.20, 0.05), ('H2O_2', 0.20, 0.05),
          ('N2O_2', 0.20, 0.05)]

L = []
L.append('Transfer_function_parameters_(TFP)_for_IIR-shaped_filter_'
         '(see_Ibrom_et_al._2007_AFM).')
L.append('fc:_IIR_cut-off_frequency')
L.append('Fn:_normalization_parameter')
L.append('Water_vapour_TFP_are_calculated_for_9_RH_classes.')
L.append('Other_gases_TFP_are_calculated_on_a_monthly_base_'
         '(currently_all_months_together_).')
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


def variant(name, sa_file):
    src = open(os.path.join(here, 'base_n_gas.eddyflow')).read().splitlines()
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
    with open(os.path.join(here, name), 'w') as fh:
        fh.write('\n'.join(out) + '\n')


variant('base_n_gas_sa.eddyflow', sa_path)
variant('base_n_gas_sa_short.eddyflow',
        os.path.join(here, 'sa_n_gas_short.txt'))
print('wrote %d-line assessment file, %d-line short file, 2 projects'
      % (len(L), len(short)))
