"""
Just an assortement of useful functions
"""

def roundbase(x,base):
    """
    Function for e.g. rounding simulation steps in postprocessing.
    
    Args:
        x : (int), value to round

        base : (int), base for rounding

    Return:
          : (int), rounded value

    """

    
    base = int(base)
        
    return int(base * round( float(x) / base ))
