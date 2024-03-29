# **************************************************************************
# *
# * Authors:     Grigory Sharov (gsharov@mrc-lmb.cam.ac.uk)
# *              J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pyworkflow.protocol.params import LabelParam, LEVEL_ADVANCED
from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER, WEB_DJANGO
from pwem.viewers import DataView, EmPlotter

from .protocols.protocol_msa import ImagicProtMSA


class ImagicViewerMSA(ProtocolViewer):
    """ Visualization of MSA results. """

    _environments = [DESKTOP_TKINTER, WEB_DJANGO]
    _targets = [ImagicProtMSA]
    _label = 'viewer MSA'

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        form.addParam('doShowEigenImages', LabelParam,
                      label="Show eigenimages?", default=True)
        form.addParam('doShowEigPixImages', LabelParam,
                      expertLevel=LEVEL_ADVANCED,
                      label="Show eigenvectors in pixel space?",
                      default=True)
        form.addParam('doShowPixelVecCoordImages',
                      LabelParam, expertLevel=LEVEL_ADVANCED,
                      label="Show eigenvectors in image space?",
                      default=True)
        form.addParam('doShowHistogram', LabelParam,
                      label="Show eigenvalue histogram?",
                      default=True)
        form.addParam('doShowMsaLisFile', LabelParam,
                      expertLevel=LEVEL_ADVANCED,
                      label="Show MSA lis file?", default=True)
        form.addParam('doShowMsaPltFile', LabelParam,
                      expertLevel=LEVEL_ADVANCED,
                      label="Show MSA plt file?", default=True)

    def _getVisualizeDict(self):
        return {'doShowEigenImages': self._viewParam,
                'doShowEigPixImages': self._viewParam,
                'doShowPixelVecCoordImages': self._viewParam,
                'doShowHistogram': self._plotHistogram,
                'doShowMsaLisFile': self._viewParam,
                'doShowMsaPltFile': self._viewParam,
                }

    def _viewParam(self, param=None):
        view = None
        if param == 'doShowEigenImages':
            view = DataView(self.protocol.getOutputEigenImages())
        elif param == 'doShowEigPixImages':
            view = DataView(self.protocol._getFileName('msa_eigen_pixel'))
        elif param == 'doShowPixelVecCoordImages':
            view = DataView(self.protocol._getFileName('msa_pixvec_coord'))
        elif param == 'doShowMsaLisFile':
            view = self.textView([self.protocol.getOutputLis()], "MSA lis file")
        elif param == 'doShowMsaPltFile':
            view = self.textView([self.protocol.getOutputPlt()], "MSA plt file")

        return [view]

    def _plotHistogram(self, param=None):
        """ First we parse the MSA plt:
        first column: cumulative percent.
        second column: iteration number.
        """

        iters = []
        cumPercents = []
        fn = self.protocol.getOutputPlt()
        with open(fn) as f:
            lines_after_2 = f.readlines()[2:]
            for line in lines_after_2:
                values = line.split()
                cumPercents.append(float(values[0]))
                iters.append(int(float(values[1])))
        f.close()

        width = 0.85
        xplotter = EmPlotter()
        a = xplotter.createSubPlot('Behaviour of sum of eigenvalues during analysis',
                                   'Iteration number', '%')
        a.bar(iters, cumPercents, width, color='b')

        return [xplotter]
