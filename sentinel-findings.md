# sentinel-findings.md

*Written by Sentinel — WARN-tier findings.*

## [WARN] P1-empty-dataframe-access
- **Location:** `run_analysis.py:255`
- **Detail:** pattern matched: q90_p = q90_row['P'].values[0]
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-19T20:15:00.401347+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_simulation.py:83`
- **Detail:** pattern matched: q90_p = q90_row['P'].values[0]
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-19T20:15:00.408328+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_simulation.py:83`
- **Detail:** pattern matched: q90_p = q90_row['P'].values[0]
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T11:57:01.238717+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_analysis.py:255`
- **Detail:** pattern matched: q90_p = q90_row['P'].values[0]
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T11:57:01.241514+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_analysis.py:255`
- **Detail:** pattern matched: q90_p = q90_row['P'].values[0]
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T11:57:41.117744+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_simulation.py:83`
- **Detail:** pattern matched: q90_p = q90_row['P'].values[0]
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T11:57:41.121806+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_analysis.py:257`
- **Detail:** pattern matched: q90_p = q90_row['P'].values[0]
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:00:30.011464+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_simulation.py:85`
- **Detail:** pattern matched: q90_p = q90_row['P'].values[0]
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:00:30.012528+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_simulation.py:85`
- **Detail:** pattern matched: q90_p = q90_row['P'].values[0]
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:06:27.831173+00:00

## [WARN] P1-empty-dataframe-access
- **Location:** `run_analysis.py:257`
- **Detail:** pattern matched: q90_p = q90_row['P'].values[0]
- **Fix hint:** Guard with `if df.empty: return <sentinel>` or `if len(df) == 0: ...` immediately before the positional access, OR use `.iat[0]` / `.at[...]` with a prior existence check. If the file is a generator where this access is provably safe repo-wide, add `# sentinel:skip-file` near the top of the file.

- **Source:** MEMORY.md#top-5-cross-project-defects
- **When:** 2026-04-23T12:06:27.836908+00:00
