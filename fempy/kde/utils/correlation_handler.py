'''
Module to compute the correlation function from mixed-event and mixed-event
probability density functions.
'''

from ROOT import TF1


class CorrelationHandler:
    '''
    Class to handle the correlation function.
    '''

    def __init__(self, same_event, mixed_event, lazy=False, name='cf', lambda_param=1) -> None:
        """Default constructor.

        Parameters
        ----------
        same_event : TF1
            same_event distribution
        mixed_event : TF1
            mixed_event distribution
        lazy : bool, default = False
            if true the correlation function is not computed during the
            construction of the object
        name : string, default = 'cf'
            name of the CorrelationHandler object
        lambda_param : float, default = 1
            lambda parameter associated to the correlation function. Used when
            the CF is combined with other CorrelationHandler bjects
        
        Returns
        -------
        string
            a value in a string

        Raises
        ------
        KeyError
            when a key error
        OtherError
            when an other error
        """
        if not isinstance(same_event, TF1):
            print("Error")
        if not isinstance(mixed_event, TF1):
            print("Error")

        self.same_event = same_event.Clone('se')
        self.mixed_event = mixed_event.Clone('me')
        self.name = name
        self.lambda_param = lambda_param

        # if not lazy:
        #     self.compute_cf()

        print('expformula: ', self.same_event.GetExpFormula())

    def cf_formula(self, x, par):
        return self.same_event.Eval(x[0]) / self.mixed_event.Eval(x[0])


    # def compute_cf(self):
    #     '''Computes the correlation function'''
    #     x_min = self.same_event.GetMinimumX()
    #     x_max = self.same_event.GetMaximumX()

    #     lambdaex = lambda x, par : self.same_event.Eval(x[0]) / self.mixed_event.Eval(x[0])
    #     self.corr_func = TF1(self.name, lambdaex, x_min, x_max)
    #     self.corr_func.Draw()

    # def eval(self, x, par):
    #     return self.same_event.Eval(x[0]) / self.mixed_event.Eval(x[0])

    # def __Add__(self, cf2):
    #     return TF1('sum', f'{self.lambda_param}*{self.name}+{cf2.lambda_param}*{cf2.name}', self.)
