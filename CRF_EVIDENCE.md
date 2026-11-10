# Reproducing the P301 case-report-form audit

The participant CRFs are FDA-released exports of electronic forms. Download the originals from [PHMPT](https://phmpt.org/multiple-file-downloads/) and extract them using the project’s existing scripts. No raw participant PDFs, local paths or staff-name histories are required in Git. [corpus_manifest.csv](evidence/corpus_manifest.csv) gives each archive, member, release, size and ZIP CRC32, including placeholders. Use the archive/member when there is no direct PDF URL. Page numbers are one-based PDF pages; Bates locators identify the source even when two exports differ in pagination.

## Find the cited data

| Report question | Evidence | How to check |
|---|---|---|
| What did the unblinding form actually say? | [arm_checkbox.csv](evidence/arm_checkbox.csv), [dose2_flags.csv](evidence/dose2_flags.csv), [blinded EX rows in the flagged cases](evidence/blinded_exposure_in_flagged_cases.csv) | Read the printed checkbox mark on the source PDF; missing marks stay unreadable. Join by subject to original assignment and blinded Day 1/29 administrations. A dose named at unblinding is not itself an administration. |
| What explains the second-dose entries? | [171 eligible participants](evidence/crf_source/Moderna_171_dose2_visit_comparison.csv), [50 discrepancies](evidence/crf_source/Moderna_50_dose2_discrepancies.csv) | Compare visit location/date, dated field entries, actual open-label vaccination and per-protocol flag. [Script](crf_audit/dose2_visit_history.py). US331 comments establish a workaround in three participants; a common timestamp alone does not identify the source of a field. |
| Which illnesses and source comments were reviewed? | [reviewed_cases.csv](evidence/reviewed_cases.csv), [candidate pages](evidence/illness_visit_search_hits.csv), [case_PCR_records.csv](evidence/case_PCR_records.csv), [illness_DV.csv](evidence/illness_DV.csv), [dv_links.csv](evidence/dv_links.csv) | The extractor searches for explicit language about a missed illness visit; a hit is **not** a confirmed failure. Read the site response in context. The 19 reviewed example rows are a human selection. The 55 PCR rows are a focused exhibit, not all tests for these subjects. All 24 DV rows for the three specified subjects and the five linked decision rows for US3932018 can be regenerated with [illness_context.py](crf_audit/illness_context.py). |
| When was contact first entered? | [contact_lags.csv](evidence/contact_lags.csv) | Compare the recorded contact date with the dated initial entry; retain corrections and exclusions. The 642 eligible fields are version-specific observations, not 642 different people. Entry delay does not prove delayed care. |
| Did a completed AE form or diagnostic result reach its corresponding submission? | [AE_term_onset_checks.csv](evidence/AE_term_onset_checks.csv), [pcr_comparison.csv](evidence/pcr_comparison.csv) | Rebuild the form-to-table joins, keeping both mixed and unmatched PCR days. 1,118 EUA and 4,362 licensing form observations, the latter checked against **both** licensing tables. The 9,842 rows are comparisons, not 9,842 distinct events. |
| What happened to the previously flagged serious rows? | [previous_23_cases_CRF_crosswalk.csv](evidence/previous_23_cases_CRF_crosswalk.csv), [subsequent three-case follow-up](evidence/crf_source/serious_case_followup/) | Identify each original form, term, onset and seriousness field by PDF page. The blank same-number licensing field does not exclude a different form or diagnosis. The later US3172069 follow-up locates prostate cancer in licensing medical history; the earlier ledger is preserved as originally produced. |
| What happened to 566 earlier-backed form numbers? | [form_dispositions_566.csv](evidence/form_dispositions_566.csv), [serious_12.csv](evidence/serious_12.csv), [case decisions and candidate paths](evidence/ae_form_review/) | Reconstruct eligibility from the complete AE snapshots; read the evidence supporting each interpretation. A shared number, symptom, date or case ID is not a signed withdrawal. The 12 are a subset of the 566 and overlap the earlier 23, not additional serious events. |

## Run in RStudio

Set the **repository root** as your working directory. Run the entire [analyse_crf_evidence.R](analyse_crf_evidence.R) script. The default mode requires only base R, replays the archived counts and checks the subsets, and writes to a new folder under `analysis/`. **Replay does not read PDFs.**

To re-extract source fields, install Python 3.10+ in your chosen environment (`python -m pip install -r crf_audit/requirements.txt`). Place the extracted CRFs and the P301 S3 EX and S5/S3/S27 AE XPTs under one data root, keeping versions separate. Then run:

```r
options(phmpt.root = getwd(), phmpt.crf.mode = "source",
        phmpt.crf.data_root = "D:/Projects/moderna_trial_analysis",
        phmpt.crf.python = "D:/Python/python.exe", phmpt.crf.jobs = 2L)
source("analyse_crf_evidence.R")
```

Full mode requires the 1,576 substantive PDFs in the manifest. To intentionally run a subset, also set `options(phmpt.crf.allow_partial = TRUE)`. A partial run compares **only** scanned exports and records all absent files in `scan/coverage.csv`. It never treats an absent PDF as a negative finding. Use fresh output directories. The source stage rejects bad PDF sizes/CRC32s, verifies table hashes, and requires one matching XPT per snapshot. For multiple local copies, use the explicit-path commands below. Increase workers only as memory permits. Using `--keep-text` can require significant disk space; source extracts and candidate text stay local.

### Equivalent explicit command line

From the repository root, substitute actual local paths for the four XPTs. These commands also work outside RStudio:

```sh
python crf_audit/extract_fields.py --root DATA --out RUN/scan --jobs 2
python crf_audit/reconcile_fields.py --scan RUN/scan --out RUN/joins --ex S3_EX.xpt --ae-s5 S5_AE.xpt --ae-s3 S3_AE.xpt --ae-s27 S27_AE.xpt
python crf_audit/verify_extract.py --scan RUN/scan --joins RUN/joins
python crf_audit/dose2_visit_history.py --root DATA --joins RUN/joins --ex S3_EX.xpt --out RUN/dose2_workflow
python crf_audit/check_case_pages.py --root DATA --out RUN/named_cases
python crf_audit/check_release_tables.py --ex S3_EX.xpt --ae-s5 S5_AE.xpt --ae-s3 S3_AE.xpt --ae-s27 S27_AE.xpt
```

The first two stages create fresh observations without reading the published finding rows. `verify_extract.py` then compares them with the originals as multisets, independently of sorting or integer formatting. Both missing and extra rows cause failure. The visit-workflow check additionally reopens the dose and decision-visit fields for eligible participants. The targeted page checker verifies cited Bates pages and the 23 EUA AE terms/onsets/seriousness fields; it does **not** automatically adjudicate the reviewed narratives. The full-table checker independently reconstructs the 566 keys from **all three** AE tables, rather than starting with the 287 flagged people.

### Source rules and intermediate tables

* **Dose:** retain checked product, “Not Administered”, and unreadable as distinct states. Use the original randomized arm and blinded, positive-dose EX administrations at visits 1 and 2. Later open-label doses are kept separate. The Day 29 treatment-given audit field is checked separately; an absent extract is not “No”.
* **AE:** subject + AESPID numeric suffix, never AESEQ across releases. Within a snapshot compare normalized reported term and onset calendar day. Through-cutoff checks require a complete onset date (S5: 25 Nov 2020; S3/S27: 26 Mar 2021). S6 participant PDFs are compared with S5; S3/S10 licensing PDFs with S3 and S27. Unmatched rows remain visible. A gap with no surviving forms requires direct cross-release comparison.
* **PCR:** use final diagnostic test date/result in the versioned form audit, then date-matched usable binary results from the appropriate EUA or production PCR ledger. Emit positive, negative, mixed or no usable same-day binary result. Matching a calendar date does not prove specimen identity.
* **Contact:** enter the initial dated, nonempty contact observation and its entry timestamp; flag later revisions. Eligibility requires an unchanged date and nonnegative lag. A retrospective entry is not a missing swab.
* **Missed visits:** candidate pages, Bates locators and short excerpts are exported locally. The site response, exposure window and PCR/deviation records must be read before calling an illness visit missed. The regex is a targeted search, not proof that all complaints in free text were found.
* The compact [comparator inputs](evidence/crf_source/comparator_rows/) include the relevant ADSL columns and **all** records in each project PCR ledger for participants with substantive released CRFs. They are not the final CRF comparison results. To regenerate them, export the subject columns from the BLA ADSL JMP (example below), produce fresh `EUA/pcr_record_audit.csv` and `Prod/pcr_record_audit.csv` with [analyse_pcr_testing_by_lab.R](analyse_pcr_testing_by_lab.R), then run [prepare_comparators.py](crf_audit/prepare_comparators.py) with `--adsl FILE --pcr-eua FILE --pcr-prod FILE --out NEW_DIR`. The reference comparator files and source digests are documented in [comparators.csv](evidence/crf_source/comparators.csv). The raw data are not bundled or silently downloaded.
* To regenerate the deviation and rationale exhibits, export S3 SUPPDV decision fields as follows:
```r
source("crf_audit/prepare_source_columns.R")
crf_source_columns(getwd(), adsl_file="PATH/TO/S3_ADSL.jmp",
                   suppdv_file="PATH/TO/S3_SUPPDV.jmp", out="RUN/columns")
```
Then run `python crf_audit/illness_context.py --dv S3_DV.xpt --addv S3_ADDV.xpt --suppdv-csv RUN/columns/case_SUPPDV.csv --out RUN/illness_deviations`. The case set is deliberately fixed, not a trial-wide deviation count.
* [Original ledger hashes and row counts](evidence/crf_source/cited_files_manifest.csv) identify the twelve archived files byte-for-byte. The reports’ numeric keys are not fabricated for unmatched records. Human selections, alternative explanations and signed source actions are labelled separately. The tests concern recorded data; a complaint that never appeared on a form is outside the form-to-table check.

## Validation of this release

1. Offline R replay reproduced the archived report totals and the 566-form disposition breakdown.
2. Six independently obtained, CRC-verified original PDFs exercised both EUA and licensing, dose discrepancies, the US331 workaround, retrospective contact, missed illness visits, diagnostic forms and AE forms. The R source launcher regenerated all six comparison ledgers for those exports and matched every row/field, normalizing only numeric and boolean formatting. See [validation](evidence/crf_source/validation/).
3. The complete XPT-side check verified all 9,842 AE comparison rows and all 566 keys; the independent PCR comparator verified 640 comparisons and 55 focused records. The 24 DV and five linked decision rows were independently reproduced against the study tables. The compact ADSL input was compared with an export from the original licensing JMP (1,180 subjects, ten fields).
4. An additional targeted check on the 36 PDFs cited for named examples and the 23 serious EUA events verified all 44 page references in the 19 reviewed example rows and all 23 original serious forms.
5. The six-PDF source test is **not** a new full 1,576-PDF extraction. The full-corpus reference ledgers come from the archived audit. The script records missing sources and any disagreement so that others can run a genuine independent full-corpus check.
