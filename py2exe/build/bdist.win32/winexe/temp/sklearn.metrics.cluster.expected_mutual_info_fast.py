
def __load():
    import imp, os, sys
    try:
        dirname = os.path.dirname(__loader__.archive)
    except NameError:
        dirname = sys.prefix
    path = os.path.join(dirname, 'sklearn.metrics.cluster.expected_mutual_info_fast.pyd')
    #print "py2exe extension module", __name__, "->", path
    mod = imp.load_dynamic(__name__, path)
##    mod.frozen = 1
__load()
del __load
