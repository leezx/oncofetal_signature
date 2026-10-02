# Repository working instructions

- Keep raw and processed biological data under `/Volumes/Stelligen_SSD/Stelligen/DATA`; never commit them.
- Keep executable analysis code under the relevant step directory, currently `step2_fetal/scripts/`.
- Keep stable tables, figure source data, and figures under the matching `results/` subdirectories.
- Update `WORKLOG.md` at every substantial checkpoint; do not rewrite historical entries except to correct facts.
- Preserve frozen analysis gates unless a documented review decision explicitly changes them.
- Do not use downstream CRC/TWEAKR results to tune the Step 2 developmental definition.
- Re-run the relevant assertions and update the rendered report after changing analysis code or results.
