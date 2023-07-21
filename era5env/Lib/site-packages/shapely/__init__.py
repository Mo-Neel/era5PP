

""""""# start delvewheel patch
def _delvewheel_init_patch_1_2_0():
    import ctypes
    import os
    import platform
    import sys
    libs_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, 'shapely.libs'))
    is_pyinstaller = getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS')
    is_conda_cpython = platform.python_implementation() == 'CPython' and (hasattr(ctypes.pythonapi, 'Anaconda_GetVersion') or 'packaged by conda-forge' in sys.version)
    if sys.version_info[:2] >= (3, 8) and not is_conda_cpython or sys.version_info[:2] >= (3, 10):
        if not is_pyinstaller or os.path.isdir(libs_dir):
            os.add_dll_directory(libs_dir)
    else:
        load_order_filepath = os.path.join(libs_dir, '.load-order-shapely-2.0.1')
        if not is_pyinstaller or os.path.isfile(load_order_filepath):
            with open(os.path.join(libs_dir, '.load-order-shapely-2.0.1')) as file:
                load_order = file.read().split()
            for lib in load_order:
                lib_path = os.path.join(os.path.join(libs_dir, lib))
                if not is_pyinstaller or os.path.isfile(lib_path):
                    ctypes.WinDLL(lib_path)


_delvewheel_init_patch_1_2_0()
del _delvewheel_init_patch_1_2_0
# end delvewheel patch

from .lib import GEOSException  # NOQA
from .lib import Geometry  # NOQA
from .lib import geos_version, geos_version_string  # NOQA
from .lib import geos_capi_version, geos_capi_version_string  # NOQA
from .errors import setup_signal_checks  # NOQA
from ._geometry import *  # NOQA
from .creation import *  # NOQA
from .constructive import *  # NOQA
from .predicates import *  # NOQA
from .measurement import *  # NOQA
from .set_operations import *  # NOQA
from .linear import *  # NOQA
from .coordinates import *  # NOQA
from .strtree import *  # NOQA
from .io import *  # NOQA

# Submodule always needs to be imported to ensure Geometry subclasses are registered
from shapely.geometry import (  # NOQA
    Point,
    LineString,
    Polygon,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
    GeometryCollection,
    LinearRing,
)

from . import _version

__version__ = _version.get_versions()["version"]

setup_signal_checks()
