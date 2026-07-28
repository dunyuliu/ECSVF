# Project Rules for Claude

## Release Workflow

When the user says "release", execute this workflow end to end — do not skip any step:

1. Run `git diff` and inspect the full current change set.
2. Find the latest existing release note named `release_note_vA.B.C.md`.
3. Create the next release note by incrementing only C by 1. Keep A and B unchanged unless explicitly asked otherwise.
4. Before creating the new release note, move all older release notes into `docs/` if not already there. Create `docs/` if it does not exist.
5. Audit the project:
   - detect new or unprocessed files
   - detect naming-rule violations
   - detect duplicate or outdated files
   - detect inconsistencies across READMEs, scripts, and tracked artifacts
   - update master documents as needed
   - rename files when needed
6. Apply all fixes required by the audit before writing the release note.
7. Write the release note for the final post-audit state. Must include:
   - version and date
   - summary of scope
   - files added, removed, renamed, or cleaned up
   - content updates to master documents
   - audit findings and fixes
   - remaining open issues or pending items
8. Re-check the project after all edits to confirm the release note matches actual final state.
9. Commit all changes with a concise commit message that matches the release.

Rules:
- Do not write release notes from `git diff` alone — reconcile against the final project state.
- Do not delete old release notes — move them to `docs/`.
- Release note filename convention: `release_note_vA.B.C.md` (singular "note").

## General Rules

- **Be concise.** Cut unnecessary communications. No summaries of what was just done, no filler.
- **Docs locked with developments.** After any substantive change: audit, fix, commit, and push. Do not let docs drift from code.
