import functools
import multiprocessing

def _instance_method_alias(obj, arg):
    """
    Alias for instance method that allows the method to be called in a 
    multiprocessing pool
    """
    obj.instance_method(arg)
    return

class MyClass(object):
    """
    Our custom class whose instance methods we want to be able to use in a 
    multiprocessing pool
    """

    def __init__(self):
        self.my_string = "From MyClass: {0}, {1}"


    def multi_fib(self, args):
        _bound_instance_method_alias = functools.partial(_instance_method_alias, self)
        pool = multiprocessing.Pool(processes=4)

        pool.map(_bound_instance_method_alias, args)
        pool.close()
        pool.join()

    def single_fib(self, args):
        for a in args:
            self.instance_method(a)

    def F(self, n):
        if n == 0: return 0
        elif n == 1: return 1
        else: return self.F(n-1)+self.F(n-2)


    def instance_method(self, arg):
        """
        Some arbitrary instance method
        """

        print(self.my_string.format(arg, self.F(arg)))
        return

# create an object of MyClass
obj = MyClass()
# obj.single_fib(range(35))
obj.multi_fib(range(35))

# use functools.partial to create a new method that always has the 
# MyClass object passed as its first argument


# create our list of things we will use the pool to map



# create the pool of workers


# call pool.map, passing it the newly created function


# cleanup
