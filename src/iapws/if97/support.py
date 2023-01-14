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

def apply(value: Any, function: F, **options) -> Any:
    """Apply the function, with options, to a value."""
    return function(value, **options)

def conversion(constant: float) -> D:
    """Usefull decorator factory to implement unit conversions."""
    def decorator(function: F) -> F:
        @wraps(function)
        def wrapper(unit: float, english: bool = True):
            if english:
                unit = unit / constant
            else:
                unit = unit * constant
            return function(unit) 
        return cast(F, wrapper)
    return decorator

def english_support(inputs: tuple[F], outputs: tuple[F]) -> D:
    """Usefull decorator factory to implement flexible units."""
    def decorator(function: F) -> F:
        @wraps(function)
        def wrapper(*args, **kwargs):
            english = kwargs.pop('english', False)
            if english:
                args = tuple(method(unit) for method, unit in zip(inputs, args)) + args[len(inputs):]
                return reduce(partial(apply, english=False), outputs, function(*args, **kwargs))
            else:
                return function(*args, **kwargs)
        return cast(F, wrapper)
    return decorator
