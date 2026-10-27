# Contributing

Contributions should keep each skill portable, reviewable, and safe to load on demand.

## Skill requirements

1. Place the skill in `skills/<skill-name>/SKILL.md`.
2. Use a lowercase, hyphenated `name` that exactly matches the directory.
3. Explain both what the skill does and when an agent should use it in `description`.
4. Keep detailed material in `references/`, deterministic helpers in `scripts/`, and reusable resources in `assets/`.
5. Do not commit raw research data, unpublished results, credentials, personal filesystem paths, caches, or generated output directories.
6. Document required software and any network access in the skill itself.
7. Add or update validation when introducing executable behavior.

## Before opening a pull request

Run:

```bash
for skill in skills/*; do
  uvx --from skills-ref agentskills validate "$skill"
done
python tests/validate_repository.py
```

Review all scientific assumptions, units, calibrations, and decision thresholds with a domain expert. Keep generated examples synthetic or explicitly licensed.

