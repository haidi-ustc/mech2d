#/usr/bin/env python
#-code- utf-8
#    mech2d is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    mech2d is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public Licensea
#    along with mech2d.  If not, see <http://www.gnu.org/licenses/>.
from abc import ABCMeta, abstractmethod

class Calculator:
    __metaclass__ = ABCMeta
    def __init__(self):
        """
        CodeRun is the superclass defining the operations for running codes either directly or asynchronously via a
        submission script on a cluster.

        :param executable:   Name or Path to the executable if can be found on the $PATH or path to the executable
        :param workdir:  Path to a folder where input and output will be located (Default: '.')
        :param use_mpi:  True if code relies on MPI for execution (Default: False)

        """
        pass

    @abstractmethod
    def set_inputs(self):
        """
        This method must be implemented by child classes, it should write all the input files and prepare environment
        for execution

        :return: None
        """
        pass

    @abstractmethod
    def get_energy(self):
        """
        This method must be implemented by child classes, it parse of output data.

        :return: None
        """
        pass

