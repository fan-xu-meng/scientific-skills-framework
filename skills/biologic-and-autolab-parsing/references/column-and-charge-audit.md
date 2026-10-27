# Column and charge audit

Record instrument, acquisition program, source filename, header units, chosen time/current/potential columns, and any sign conversion. BioLogic and Autolab exports can contain multiple current, potential, or charge channels; select by header meaning rather than column position. Preserve the raw headers in the provenance table.

For each continuous interval, compare recorded charge with trapezoidal integration of current over time. A large disagreement can indicate wrong units, reset charge counters, gaps, or sign convention. Never join CV segments across a reset without preserving segment IDs and open-circuit endpoints.
