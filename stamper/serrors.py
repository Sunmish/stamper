'''Custom errors for stamper.'''


class DeclinationError(ValueError):
    '''Raise when a specific subset of values in context of app is wrong'''
    def __init__(self, message, declination, *args):
        self.message = message # without this you may get DeprecationWarning
        # Special attribute you desire with your Error, 
        # perhaps the value that caused the error?:
        self.declination = declination         
        # allow users initialize misc. arguments as any other builtin Error
        super(MyValueError, self).__init__(message, foo, *args)