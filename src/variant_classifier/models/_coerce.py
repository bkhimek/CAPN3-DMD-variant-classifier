"""Small shared helpers used by every model's from_dict() classmethod.

Kept separate from the models themselves so the coercion/validation style
stays identical across all six models rather than drifting file to file.
"""

from typing import Any, Optional, Type

from ..errors import SchemaValidationError


def require_dict(data: Any, context: str) -> dict:
    if not isinstance(data, dict):
        raise SchemaValidationError(f"{context}: expected a mapping, got {type(data).__name__}")
    return data


def require_str(data: dict, field: str, context: str) -> str:
    value = data.get(field)
    if not isinstance(value, str) or not value.strip():
        raise SchemaValidationError(f"{context}: '{field}' must be a non-empty string, got {value!r}")
    return value


def optional_str(data: dict, field: str) -> Optional[str]:
    value = data.get(field)
    if value is None:
        return None
    if not isinstance(value, str):
        raise SchemaValidationError(f"'{field}' must be a string if provided, got {type(value).__name__}")
    return value


def require_bool(data: dict, field: str, context: str) -> bool:
    value = data.get(field)
    if not isinstance(value, bool):
        raise SchemaValidationError(f"{context}: '{field}' must be true/false, got {value!r}")
    return value


def optional_int(data: dict, field: str, context: str, minimum: Optional[int] = None) -> Optional[int]:
    value = data.get(field)
    if value is None:
        return None
    if not isinstance(value, int) or isinstance(value, bool):
        raise SchemaValidationError(f"{context}: '{field}' must be an integer if provided, got {value!r}")
    if minimum is not None and value < minimum:
        raise SchemaValidationError(f"{context}: '{field}'={value} must be >= {minimum}")
    return value


def optional_float(data: dict, field: str, context: str, minimum: Optional[float] = None, maximum: Optional[float] = None) -> Optional[float]:
    value = data.get(field)
    if value is None:
        return None
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise SchemaValidationError(f"{context}: '{field}' must be a number if provided, got {value!r}")
    value = float(value)
    if minimum is not None and value < minimum:
        raise SchemaValidationError(f"{context}: '{field}'={value} must be >= {minimum}")
    if maximum is not None and value > maximum:
        raise SchemaValidationError(f"{context}: '{field}'={value} must be <= {maximum}")
    return value


def coerce_enum(enum_cls: Type, value: Any, field: str, context: str):
    if isinstance(value, enum_cls):
        return value
    try:
        return enum_cls(value)
    except ValueError:
        valid = ", ".join(e.value for e in enum_cls)
        raise SchemaValidationError(
            f"{context}: '{field}'={value!r} is not one of the allowed values [{valid}]"
        )


def require_list(data: dict, field: str, context: str) -> list:
    value = data.get(field, [])
    if not isinstance(value, list):
        raise SchemaValidationError(f"{context}: '{field}' must be a list, got {type(value).__name__}")
    return value
