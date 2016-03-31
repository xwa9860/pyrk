import math


class PrecursorData(object):

    def __init__(self, nuc, e, n):
        """initializes the precursor group data for the fissioning nuclide.
        :param e: should be 'thermal' or 'fast' to indicate the energy spectrum
        or multipt.
        """
        self._betas = self._get_betas(nuc, e)
        self._lambdas = self._get_lambdas(nuc, e)
        self._Lambda = self._get_Lambda(nuc, e)
        self._nuc = nuc
        self._e = e

    def beta(self):
        """
        The Big Beta is the fraction of all fission neutrons that are delayed.
        """
        return sum(self._betas)

    def betas(self):
        """
        The betas are populations of neutron precursors.
        """
        return self._betas

    def lambdas(self):
        """
        The lambdas are the decay constants for the neutron precursor groups.
        """
        return self._lambdas

    def Lambda(self):
        """

        """
        return self._Lambda

    def v_d(self, nuc, e):
        """
        TODO: figure out why you felt you needed this
        """
        if nuc == "u235" and e == "thermal":
            return 0.01668

    def _get_betas(self, nuc, e):
        """
        Retrieves the values for beta_i
        Data for u235 was obtained from http://arxiv.org/pdf/1001.4100.pdf
        """
        beta_dict = {}
        beta_dict["u235"] = {}
        beta_dict["pu239"] = {}
        beta_dict["sfr"] = {}
        beta_dict["fhr"] = {}
        beta_dict["u235"]["thermal"] = [0.00247, 0.0013845, 0.001222,
                                        0.0026455, 0.000832, 0.000169]
        beta_dict["u235"]["fast"] = [0.000266, 0.001491, 0.001316, 0.002849,
                                     0.000896, 0.000182]
        beta_dict["pu239"]["thermal"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        beta_dict["pu239"]["fast"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        beta_dict["sfr"]["fast"] = [0.009, 0.087, 0.070, 0.0014, 0.0060, 0.0055]
        beta_dict["fhr"]["thermal"] = [1.48756E-04, 9.45436E-04, 8.29928E-04,
                                       2.21997E-03, 6.90778E-04, 2.31801E-04]
        beta_dict["fhr"]["multipt"] = [1.48756E-04, 9.45436E-04, 8.29928E-04,
                                       2.21997E-03, 6.90778E-04, 2.31801E-04,
                                       0.084349, 0.168983]
        return beta_dict[nuc][e]

    def _get_lambdas(self, nuc, e):
        """
        Retrieves the values for lambda_i
        Data for u235 was obtained from http://arxiv.org/pdf/1001.4100.pdf
        """
        lambda_dict = {}
        lambda_dict["u235"] = {}
        lambda_dict["pu239"] = {}
        lambda_dict["sfr"] = {}
        lambda_dict["fhr"] = {}
        lambda_dict["u235"]["thermal"] = [math.log(2)/x for x in
                                          [54.51, 21.84, 6.00, 2.23, 0.496,
                                           0.179]]
        lambda_dict["u235"]["fast"] = [0.0127, 0.0317, 0.155, 0.311, 1.4, 3.87]
        lambda_dict["pu239"]["thermal"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        lambda_dict["pu239"]["fast"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        lambda_dict["sfr"]["fast"] = [0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01]
        lambda_dict["fhr"]["thermal"] = [1.25723E-02, 3.11643E-02, 1.09837E-01,
                                         3.16355E-01, 1.27045E+00, 7.83939E+00]
        lambda_dict["fhr"]["multipt"] = [1.25723E-02, 3.11643E-02, 1.09837E-01,
                                         3.16355E-01, 1.27045E+00, 7.83939E+00,
                                         786.3172199, 1209.079474]
        return lambda_dict[nuc][e]

    def _get_Lambda(self, nuc, e):
        """
        Mean generation time of the neutrons in this reactor type. [s]
        """
        Lambda_dict = {}
        Lambda_dict["u235"] = {}
        Lambda_dict["pu239"] = {}
        Lambda_dict["sfr"] = {}
        Lambda_dict["fhr"] = {}
        Lambda_dict["u235"]["thermal"] = 1.08e-5
        Lambda_dict["u235"]["fast"] = 0
        Lambda_dict["pu239"]["thermal"] = 0
        Lambda_dict["pu239"]["fast"] = 0
        Lambda_dict["sfr"]["fast"] = 1.0e-5
        Lambda_dict["fhr"]["thermal"] = 5.35878E-04 #ADJ_NAUCHI_LIFETIME from serpent
        Lambda_dict["fhr"]["multipt"] = 0.000226807
        return Lambda_dict[nuc][e]
