# Autonomy protocol for long unattended runs

Shared across projects. The session loads this at the start of a run, together
with the project's `MISSION.md`. It applies for the whole run.

Nobody is watching this run. Every decision that would normally go to a person
goes to the `phase-reviewer` subagent instead, and every decision gets written
down so it can be audited afterwards.

---

## 1. Before phase 1

1. Read `MISSION.md`, the project's `STATUS.md` / `DECISIONS.md` if present, and
   the last `runs/*/RUN_SUMMARY.md` if one exists.
2. Create the run folder: `runs/<YYYY-MM-DD>_<slug>/`. Everything this run
   writes that is not code goes here.
3. **Git checkpoint.** Create and switch to a branch `auto/<YYYY-MM-DD>_<slug>`.
   If the working tree is dirty, commit it on that branch as
   `checkpoint: state before autonomous run`. Never commit to, merge into, or
   rebase `main`. Never push.
4. Write `runs/<run>/PLAN.md`: the phases you intend, each with a one-line
   done-condition and the outputs it will produce. Keep it to one screen.
5. Send `PLAN.md` to `phase-reviewer` as **phase 0**. Log the verdict (section 4)
   before starting phase 1.

## 2. The phase loop

For each phase:

1. **Do the work.** Stay inside the phase's scope. If you discover something
   out of scope, note it for the report; don't chase it.
2. **Write `runs/<run>/phase_NN_report.md`:**
   - what the phase set out to do and the done-condition
   - what was actually run (scripts, key parameters, inputs)
   - outputs produced, with paths
   - key numbers, each with the file it can be checked against
   - anything that went wrong, was skipped, or was blocked
   - what you propose next and the alternatives you considered
3. **Commit** on the run branch: `phase NN: <one line>`.
4. **Call `phase-reviewer`.** Give it the paths to `MISSION.md`, `PLAN.md`, the
   phase report, `REVIEW_LOG.md`, and the output directory, plus the output of
   `git diff --stat HEAD~1`. Do not summarise the results for it; let it read them.
5. **Log the verdict** verbatim in `REVIEW_LOG.md` (section 4).
6. **Act on it:**
   - `PROCEED`: start the next phase using the reviewer's next-phase spec.
   - `REVISE`: redo the phase with the reviewer's fixes, then re-review.
   - `STOP`: go to section 6.

The reviewer's next-phase spec replaces your own proposal when they differ.
If you disagree, say so in the next report; don't override it silently.

## 3. Stop conditions (hard)

Stop the run and go to section 6 when any of these happens:

- the mission's done-condition is met
- the phase limit in `MISSION.md` is reached (default 8, not counting phase 0)
- the same phase gets `REVISE` twice without measurable progress
- three blocked actions in a row (section 5)
- the reviewer returns `STOP`
- continuing would need a decision tagged `NEEDS-JACOB` and no independent
  work is left to do

## 4. REVIEW_LOG.md format

One entry per review, appended, never edited afterwards:

```
## Phase NN · <title> · <HH:MM>
**Verdict:** PROCEED | REVISE | STOP
**Checked independently:** <what the reviewer verified itself, with file paths>
**Findings:** <numbered, most important first>
**Decisions made:**
- <decision> · alternatives: <...> · why: <...> · confidence: high/med/low
**NEEDS-JACOB:** <items, or "none">
**Next-phase spec:** <scope, done-condition, outputs>
```

## 5. When an action is blocked

Guardrails deny some commands outright (installs, pushes, remote access,
destructive git). A block is expected, not an error.

- Don't retry the same thing with a different spelling. Workarounds for a
  blocked action count as a breach of this protocol.
- Record it under `## Blocked actions` in `REVIEW_LOG.md`: what you tried, why,
  and what it would unlock. Tag it `NEEDS-JACOB`.
- Continue with work that doesn't depend on it, or stop (section 3).

## 6. Closing the run

1. Write `runs/<run>/RUN_SUMMARY.md`, in this order:
   1. **NEEDS-JACOB**: every item from the log, one line each, most important first
   2. **What changed**: phases done, verdicts, key outputs with paths
   3. **Low-confidence decisions** the reviewer made, so they get looked at first
   4. **Suggested next run**
2. Final commit on the run branch: `run summary`.
3. Stop. Don't start new work after the summary.

## 7. Standing rules

- **Nothing existing gets overwritten or deleted.** New outputs go in new
  versioned files or folders (`figures/pub_v1/`, `results/<run>/`). Old versions
  stay where they are.
- **Nothing gets invented.** If a marker, parameter, threshold, sample ID or
  reference isn't in the project files, don't guess it. Log `NEEDS-JACOB`.
- **Compute budget.** Several sessions share this machine. Use at most the
  thread count in `MISSION.md` (default 3) for any tool, and run nothing
  expected to take longer than 2 hours without splitting it.
- **No software installs.** Use the environments that already exist. If one is
  missing a package, that's a blocked action (section 5).
- **Local only.** No remote hosts, HPC, or pushes.
- **Numbers in reports must be traceable** to a file on disk. A number without
  a source doesn't count as a result.
- No reference to how these files were authored goes into code, commit
  messages, reports or figures.
