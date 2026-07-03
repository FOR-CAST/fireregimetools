# summarize_fire_sizes errors on a missing size column

    Code
      summarize_fire_sizes(root, size_col = "nope")
    Condition
      Error:
      ! summarize_fire_sizes(): size column 'nope' not in dataset (have: scenario, size, replicate)

