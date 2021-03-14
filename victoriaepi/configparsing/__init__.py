"""
This package enables the use of configuration files
"""
from pkg_resources import resource_filename

filepath = resource_filename(__name__, 'resources/thefile')
calls_dir =  resource_filename(__name__, 'calls/')
models_dir =  resource_filename(__name__, 'models/')
buildingblocks_dir =  resource_filename(__name__, 'models/buildingblocks/')
indentStr="    "
