from distutils.core import setup
import py2exe


# h5py example
# NOTE:  h5py needs numpy so this builds on the numpy example
#options = {}
#options["py2exe"] = {"includes": ["h5py.defs", "h5py.utils", "h5py._proxy", "h5py.h5p"],
#                  "excludes": ["h5py.ipy_completer", "IPython", "Tkinter", "tcl"]}
#
#setup(console=['0minup.v0.43m.py'], options=options)

setup(console=['ampbalance.py'],options={"py2exe":{"includes":["h5py.*","cython.*","scipy.special.*","scipy.linalg.*","scipy.sparse.csgraph._validation","mlpy.*"], "dll_excludes": ["MSVCP90.dll"] ,"excludes": ["IPython", "Tkinter", "tcl"]}})
