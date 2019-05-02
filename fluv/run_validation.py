import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
from validation import historical_validation

V = historical_validation()
V.plot()