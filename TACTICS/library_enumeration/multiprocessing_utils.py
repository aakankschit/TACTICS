import dill as pickle

def initializer():
    """Initializer function for the multiprocessing pool."""
    pickle.settings.update({'recurse': True}) 