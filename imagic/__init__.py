# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import pyworkflow.em
from pyworkflow.utils import Environ, join

_logo = "imagic_logo.png"
_references = ['vanHeel1981', 'vanHeel1996', 'vanHeel2012']

IMAGIC_HOME_VAR = 'IMAGIC_HOME'


# The following class is required for Scipion to detect this Python module
# as a Scipion Plugin. It needs to specify the PluginMeta __metaclass__
# Some function related to the underlying package binaries need to be
# implemented
class Plugin:
    #__metaclass__ = pyworkflow.em.PluginMeta

    @classmethod
    def getEnviron(cls):
        """ Load the environment variables needed for Imagic.
    IMAGIC_ROOT is set to IMAGIC_HOME MPI-related vars
    are set if IMAGIC_HOME/openmpi path exists
    IMAGIC_BATCH is needed for batch files to work.
    """
        env = Environ(os.environ)
        IMAGIC_HOME = os.environ[('%s' % IMAGIC_HOME_VAR)]

        env.update({'IMAGIC_ROOT': IMAGIC_HOME,
                    'IMAGIC_BATCH': "1",
                    'FFTW_IPATH': IMAGIC_HOME + '/fftw/include',
                    'FFTW_LPATH': IMAGIC_HOME + '/fftw/lib',
                    'FFTW_LIB': 'fftw3f'
                    })
        env.set('LD_LIBRARY_PATH', IMAGIC_HOME + '/fftw/lib', env.BEGIN)
        env.set('LD_LIBRARY_PATH', IMAGIC_HOME + '/lib', env.BEGIN)

        mpidir = IMAGIC_HOME + '/openmpi'

        if os.path.exists(mpidir):
            env.update({'MPIHOME': mpidir,
                        'MPIBIN': mpidir + '/bin',
                        'OPAL_PREFIX': mpidir
                        })
            env.set('PATH', mpidir + '/bin', env.BEGIN)

        else:
            print "Warning: IMAGIC_ROOT directory (", IMAGIC_HOME, ") does not contain openmpi folder.\n", \
                "No MPI support will be enabled."

        return env

    @classmethod
    def getVersion(cls):
        path = os.environ[IMAGIC_HOME_VAR]
        for v in cls.getSupportedVersions():
            versionFile = join(path, 'version_' + v)
            if os.path.exists(versionFile):
                return v
        return ''

    @classmethod
    def getSupportedVersions(cls):
        """ Return the list of supported binary versions. """
        return ['110308', '160418']

    @classmethod
    def validateInstallation(cls):
        """ This function will be used to check if package is properly installed. """
        environ = cls.getEnviron()
        missingPaths = ["%s: %s" % (var, environ[var])
                        for var in [IMAGIC_HOME_VAR]
                        if not os.path.exists(environ[var])]

        return (["Missing variables:"] + missingPaths) if missingPaths else []

    @classmethod
    def getScript(cls, *paths):
        """ Return the script that will be used. """
        cmd = join(__path__[0], 'scripts', cls.getVersion(), *paths)
        return str(cmd)


pyworkflow.em.Domain.registerPlugin(__name__)