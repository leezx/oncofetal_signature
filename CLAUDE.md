# Claude Code instructions: oncofetal_signature

@AGENTS.md

## Start of every session

1. Read `docs/HANDOFF.md`. It covers status, the four frozen signature
   objects, gate rules, script run order, data paths and open items.
2. Read `core_oncofetal/HANDOFF_RESTRICTED.md`. It is local only and
   git-ignored, and holds memberships and key values.
3. Read the latest entries at the end of `WORKLOG.md`.
4. The working branch is `analysis/core-oncofetal` (PR #6). Check
   `git status` and `git log --oneline -5`.

## Non-negotiable rules

**Restricted data.** Joanito (Synapse syn26844071) derived material must
never be committed: C-gate values, any membership list (CIOC, 338, Extended
CIOC, ECOS-38), `core_oncofetal/Methods.md`, manuscript text and funnel
workbooks.
- Keep such files git-ignored.
- Mirror them to DATA
  `2.PROJECTS/1.TWEAKR-oncoFetal/results/2026-10-02_step3_cancer_axis_v0.1/restricted_joanito/core_oncofetal/`.
- Commit everything else normally; PRs are used for review.

**Freeze protocol.** Rule changes are written into the docs and committed
before computing, applied once, and reported without tuning, rescue or
manual gene removal. Signature construction is complete; later work is
validation and application only and must not alter membership.

**Honesty.** Report inconvenient outcomes plainly. Use these terms:
"tumour-derived epithelial cells", not malignant; ECOS-38 is
"epithelial-compatible", not specific.

**Data** stays under `/Volumes/Stelligen_SSD/Stelligen/DATA` and is never
committed. Update `WORKLOG.md` at each checkpoint.
