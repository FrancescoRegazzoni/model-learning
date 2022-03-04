import os
from setuptools import setup
from setuptools import find_packages

try:
    model_learning_path = os.path.dirname(os.path.realpath(__file__))
    import json
    print('model-learning path: %s' % model_learning_path)
    config_file = os.path.expanduser('~') + '/.py-model-learning'
    config = {'model-learning_path' : model_learning_path}
    with open(config_file, 'w') as f:
        json.dump(config, f, indent = 4)
    print('Created configuration file %s.' % config_file)
except:
    print('Unable to create configuration file.')
    pass

with open("requirements.txt", "r") as f:
    requirements = [line.strip() for line in f.readlines()]

exec(open('pyModelLearning/_version.py').read())

setup(
    name = "pyModelLearning",
    version = __version__,
    description = "Python wrapper of the model-learning library.",
    author="Francesco Regazzoni",
    author_email = "francesco.regazzoni@polimi.it",
    license = "MIT",
    install_requires = requirements,
    packages = find_packages(),
)