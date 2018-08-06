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
from pyworkflow.utils import Environ

_logo = "imagic_logo.png"
_references = ['vanHeel1981', 'vanHeel1996', 'vanHeel2012']

IMAGIC_HOME = 'IMAGIC_HOME'


class Plugin(pyworkflow.em.Plugin):
    _homeVar = IMAGIC_HOME
    _pathVars = [IMAGIC_HOME]
    _supportedVersions = ['110308', '160418']

    @classmethod
    def getEnviron(cls):
        """ Load the environment variables needed for Imagic.
        IMAGIC_ROOT is set to IMAGIC_HOME MPI-related vars
        are set if IMAGIC_HOME/openmpi path exists
        IMAGIC_BATCH is needed for batch files to work.
        """
        env = Environ(os.environ)
        env.update({'IMAGIC_ROOT': cls.getHome(),
                    'IMAGIC_BATCH': "1",
                    'FFTW_IPATH': cls.getHome('fftw', 'include'),
                    'FFTW_LPATH': cls.getHome('fftw', 'lib'),
                    'FFTW_LIB': 'fftw3f'
                    })
        env.set('LD_LIBRARY_PATH', cls.getHome('fftw', 'lib'), env.BEGIN)
        env.set('LD_LIBRARY_PATH', cls.getHome('lib'), env.BEGIN)

        mpidir = cls.getHome('openmpi')

        if os.path.exists(mpidir):
            env.update({'MPIHOME': mpidir,
                        'MPIBIN': mpidir + '/bin',
                        'OPAL_PREFIX': mpidir
                        })
            env.set('PATH', mpidir + '/bin', env.BEGIN)

        else:
            print "Warning: IMAGIC_ROOT directory (", cls.getHome(), ") does not contain openmpi folder.\n", \
                "No MPI support will be enabled."

        return env

    @classmethod
    def getVersion(cls):
        for v in cls.getSupportedVersions():
            if os.path.exists(cls.getHome('version_' + v)):
                return v
        return ''

    @classmethod
    def getScript(cls, *paths):
        """ Return the script that will be used. """
        cmd = os.path.join(__path__[0], 'scripts', cls.getVersion(), *paths)
        return str(cmd)


pyworkflow.em.Domain.registerPlugin(__name__)