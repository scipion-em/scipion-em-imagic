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

from pwem.protocols import EMProtocol

from ..utils import runTemplate


class ImagicProtocol(EMProtocol):
    """
    Provides the common execution framework required for IMAGIC-based cryo-EM
    image processing workflows. The protocol acts as a foundational interface
    between Scipion and the IMAGIC environment, enabling derived protocols to
    execute external IMAGIC operations in a controlled and reproducible manner.

    AI Generated:

    IMAGIC Base Protocol (ImagicProtocol) - User Manual
        Overview

        The IMAGIC Base Protocol provides the infrastructure needed to integrate
        IMAGIC image processing utilities into cryo-EM workflows managed within
        Scipion. Rather than representing a standalone biological analysis step,
        this protocol serves as the operational backbone for protocols that rely
        on IMAGIC programs for tasks such as multivariate statistical analysis,
        classification, averaging, and other image-processing procedures.

        In practical terms, the protocol standardizes how IMAGIC jobs are
        prepared, executed, and monitored. This allows higher-level workflows to
        focus on biological interpretation while ensuring that the underlying
        execution environment remains consistent and reliable. Users typically
        interact with this protocol indirectly through more specialized IMAGIC
        tools built on top of it.

        Role Within Cryo-EM Workflows

        IMAGIC has historically been an important software package for single-
        particle cryo-EM analysis, particularly in workflows involving image
        classification and multivariate statistical approaches. The base protocol
        provides the shared operational layer that enables those analyses to run
        smoothly inside Scipion.

        From the perspective of a biological user, this protocol ensures that
        image-processing tasks are executed within the correct working
        environment and that intermediate and final outputs are consistently
        organized. This is especially important in large cryo-EM projects where
        multiple processing stages and external software packages must coexist in
        a reproducible manner.

        Execution Environment and Workflow Integration

        The protocol manages the preparation of IMAGIC execution environments,
        including the organization of temporary files, working directories, and
        image stacks required by downstream procedures. By standardizing these
        operations, derived protocols can focus on scientific tasks such as
        dimensionality reduction, particle classification, or averaging without
        requiring users to manually manage external scripts.

        The protocol also supports automated execution of IMAGIC command
        templates. This simplifies the interaction between Scipion and IMAGIC
        while reducing the risk of user configuration mistakes. In facility or
        high-throughput environments, this level of automation is particularly
        valuable because it improves reproducibility across many independent
        datasets.

        Error Detection and Reliability

        One important aspect of the protocol is its ability to monitor IMAGIC
        execution logs and detect failures automatically. In complex cryo-EM
        workflows, silent execution failures can propagate downstream and lead to
        biologically misleading interpretations. Automated validation of
        execution status therefore contributes to workflow reliability and helps
        users identify problematic processing stages early.

        Biological and Practical Importance

        Although the protocol itself does not perform direct biological analysis,
        it enables many biologically meaningful operations implemented in IMAGIC-
        based workflows. Stable execution infrastructure is essential for
        analyses such as identifying structural heterogeneity, separating
        conformational states, or improving signal quality through particle
        classification and averaging.

        In practice, users benefit from this abstraction because they can access
        sophisticated IMAGIC functionality through consistent interfaces without
        needing detailed knowledge of the original command-line environment.

        Final Perspective

        The IMAGIC Base Protocol serves as the foundational execution layer for
        IMAGIC-related cryo-EM workflows inside Scipion. By standardizing job
        preparation, execution management, and error monitoring, it enables
        higher-level protocols to focus on biologically meaningful analyses
        while maintaining reproducibility and operational stability across
        complex image-processing pipelines.
    """
    _label = None

    def _getFileName(self, key, **kwargs):
        """ Give a key, append the img extension
        and prefix the protocol working dir.
        """
        template = '%(' + key + ')s' + '.img'

        return self._getPath(template % self._params)

    def runTemplate(self, inputScript, paramsDict):
        """ This function will create a valid Imagic script
        by copying the template and replacing the values in dictionary.
        After the new file is read, the Imagic interpreter is invoked.
        """
        self._enterWorkingDir()

        log = getattr(self, '_log', None)
        runTemplate(inputScript, paramsDict, log)

        self._leaveWorkingDir()

        with open(self.getLogPaths()[0], 'r') as f:
            for line in f:
                if '**ERROR' in line:
                    raise RuntimeError('IMAGIC script error!')
