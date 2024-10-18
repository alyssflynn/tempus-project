def cast_float(value: str):
    """Attempts to cast a value as a float, otherwise returns it unchanged."""
    try:
        return float(value)
    except (ValueError, TypeError):
        return value
