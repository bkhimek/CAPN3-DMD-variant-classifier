"""Shared exception types for the variant_classifier package."""


class SchemaValidationError(ValueError):
    """Raised when a curated record fails schema validation.

    Every raise site includes which record (variant_id / gene / file) and
    which field failed, so a malformed fixture is rejected with a message a
    curator can act on directly, rather than a bare traceback deep inside a
    dataclass constructor.
    """
