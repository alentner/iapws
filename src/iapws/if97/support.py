""" Provides supporting routines and decorators for IF97 """

# type annotations
from __future__ import annotations
from typing import cast, TYPE_CHECKING

# system libraries
from functools import partial, reduce, wraps

# static analysis
if TYPE_CHECKING:
    from typing import Any, Callable, TypeVar
    F = TypeVar('F', bound = Callable[..., Any])
    D = Callable[[F], F]

# deal w/ runtime cast
else:
    F = None

def _apply(value: Any, function: F, **options) -> Any:
    """Apply the function, with options, to a value."""
    return function(value, **options)

def _conversion(constant: float) -> D:
    """Usefull decorator factory to implement unit conversions."""
    def decorator(function: F) -> F:
        @wraps(function)
        def wrapper(unit: float, /, *, english: bool = True) -> float:
            if english:
                return unit / constant
            else:
                return unit * constant
        return cast(F, wrapper)
    return decorator

def _english(inputs: tuple[F], outputs: tuple[F]) -> D:
    """Usefull decorator factory to implement flexible units."""
    def decorator(function: F) -> F:
        @wraps(function)
        def wrapper(*args, english: bool = False, **kwargs):
            if english:
                args = tuple(method(unit) for method, unit in zip(inputs, args)) + args[len(inputs):]
                return reduce(partial(_apply, english=False), outputs, function(*args, **kwargs))
            else:
                return function(*args, **kwargs)
        return cast(F, wrapper)
    return decorator

def _first(function: F) -> F:
    """Usefull decorator to pass only first positional argument."""
    @wraps(function)
    def wrapper(_, *args, **kwargs):
        return function(_)
    return cast(F, wrapper)

def _one(*args, **kwargs):
    """Usefull dummy function for zero result."""
    return 1.000

def _output(_, **kwargs):
    """Usefull dummy function for output functions."""
    return _

def _region(methods: dict[int, F], find: F) -> D:
    """Usefull decorator factory to implement dispatch to region."""
    def decorator(function: F) -> F:
        @wraps(function)
        def wrapper(*args, region: int = 0, **kwargs):
            if region == 0:
                region = find(*args)
            assert(region in methods), "Water properties not available!"
            return methods[region](*args)
        return cast(F, wrapper)
    return decorator

def _zero(*args, **kwargs):
    """Usefull dummy function for zero result."""
    return 0.000
