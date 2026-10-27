# Scientific Skills Framework

Agent skills for scientific data analysis, organized as small, task-specific folders. The collection covers electrochemistry, chromatography, spectroscopy, microscopy, and materials characterization.

Each folder in [`skills/`](skills/) has a `SKILL.md` entry point and, where useful, a script. This follows the open [Agent Skills specification](https://agentskills.io/specification): an agent can discover a skill from its name and description, then load its instructions when needed.

## Skills

| Area | Skills |
| --- | --- |
| Electrochemistry | [`biologic-and-autolab-parsing`](skills/biologic-and-autolab-parsing/), [`co2rr-fe-analysis`](skills/co2rr-fe-analysis/), [`corr-fe-analysis`](skills/corr-fe-analysis/) |
| Chromatography and mass spectrometry | [`gc-fe-analysis`](skills/gc-fe-analysis/), [`gcms-analysis`](skills/gcms-analysis/), [`hplc-analysis`](skills/hplc-analysis/), [`ic-hcoo-co3-analysis`](skills/ic-hcoo-co3-analysis/) |
| NMR and optical spectroscopy | [`nmr-hcook-analysis`](skills/nmr-hcook-analysis/), [`uv-vis-analysis`](skills/uv-vis-analysis/) |
| Materials characterization | [`tem-analysis`](skills/tem-analysis/), [`tristar-bet-analysis`](skills/tristar-bet-analysis/), [`xas-analysis`](skills/xas-analysis/), [`xps-analysis`](skills/xps-analysis/), [`xrd-analysis`](skills/xrd-analysis/) |

## Use

Install the collection with a compatible Agent Skills client:

```bash
npx skills add fan-xu-meng/scientific-skills-framework
```

You can also copy one folder from `skills/` into your agent's skills directory. Read its `SKILL.md` for input requirements, checks, and any bundled script. Scripts operate on data paths supplied by the user; instrument calibration and scientific interpretation still need review against the original records.

## Validate

```bash
for skill in skills/*; do uvx --from skills-ref agentskills validate "$skill"; done
python tests/validate_repository.py
```

The validator checks skill structure, local links, Python syntax, and repository hygiene. GitHub Actions runs the same checks on each pull request.

## References

The organization is inspired by [Google DeepMind's Science Skills](https://github.com/google-deepmind/science-skills) and [K-Dense Scientific Agent Skills](https://github.com/K-Dense-AI/scientific-agent-skills). For the research context, see [*Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*](https://arxiv.org/abs/2609.00065). This is an independent project.

Licensed under [Apache-2.0](LICENSE).

