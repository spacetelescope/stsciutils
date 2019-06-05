from pkg_resources import get_distribution, DistributionNotFound
__version__ = 'UNKNOWN'
__version_date__ = '^__version_date__ is deprecated^'
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    pass
