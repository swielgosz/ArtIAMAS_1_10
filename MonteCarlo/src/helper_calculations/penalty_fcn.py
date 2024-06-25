# This script should accept some map and return the penalty associated with a constraint violation

# List of constraints we need to implement (? denotes a maybe):
    # WSN enforcement (equality)
    # Minimum sensor-target distance (inequality)
    # Minimum sensor-sensor spacing (inequality)
    # Valid placement (equality)
    # Anything else???

# Equality constraint penalties are usually framed as a P = f(x)
# Inequality constraints are usually framed as P = if{x>limit}{0, g(x)} 
    # where it's zero if a constraint is satisfied and g(x) otherwise