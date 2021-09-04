#/usr/bin/env python
#-code- utf-8
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

