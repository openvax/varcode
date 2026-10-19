# Varcode #604: cancel superseded PR test matrices

Scope the Tests workflow concurrency group to its workflow name, event type
and PR number. Non-PR runs use their unique run ID, so separate PRs and every
main push remain independent. Cancel in-progress runs only on PR events.

Run lint and the full suite, then push two successive PR commits while the
first matrix is active. Verify GitHub cancels the first run automatically and
the latest completes every matrix job. Record the run links in the PR, merge
only after green CI, and deploy the patch release from clean main.

Reference: GitHub's workflow syntax, concurrency and cancel-in-progress:
https://docs.github.com/en/actions/reference/workflows-and-actions/workflow-syntax#concurrency
