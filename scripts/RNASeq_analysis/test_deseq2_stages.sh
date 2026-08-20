#!/bin/bash
# Tests for the 4-stage DESeq2 pipeline.
#   - unit tests on deseq2_common.R helpers
#   - end-to-end runs of stages 1 and 2 against synthetic counts
#   - the stage-signature guard that makes splitting safe
#   - submitter --stage parsing
# Stages 3-4 need DESeq2; their ARGUMENT GUARDS are tested here (they run before
# library(DESeq2)), the fitting itself is exercised on the cluster.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"

W=/tmp/dq; rm -rf $W; mkdir -p $W/counts
LAST=$W/last.out
PASS=0; FAIL=0
ok()  { if eval "$2"; then PASS=$((PASS+1)); echo "  PASS   $1"; else FAIL=$((FAIL+1)); echo "  FAIL   $1"; fi; }
cap() { "$@" > "$LAST" 2>&1; echo $?; }
has() { grep -qF -- "$1" "$LAST"; }

# ---- synthetic counts ------------------------------------------------------
mk() { # mk <name> <ambiguous>
  mkdir -p "$W/counts/$1"
  { for i in $(seq 1 200); do printf "ENSG%03d\t%d\n" "$i" $(( (RANDOM % 500) + 10 )); done
    printf "__no_feature\t5000\n__ambiguous\t%d\n__too_low_aQual\t0\n__not_aligned\t0\n__alignment_not_unique\t3000\n" "$2"
  } > "$W/counts/$1/$1_genecounts.tsv"
  printf "sample\tstrand_verdict\tstranded_used\n%s\treverse\treverse\n" "$1" > "$W/counts/$1/$1_library_qc.tsv"
}
for s in Parental_rep1 Parental_rep2 Parental_rep3 Bottom10_rep1 Bottom10_rep2 Bottom10_rep3; do mk "$s" 4000; done
mk UNC_rep1 4000; mk UNC_rep2 4000

S1="Rscript deseq2_01_setup.R --counts-dir $W/counts"
S2="Rscript deseq2_02_counts.R --counts-dir $W/counts"

echo
echo "== unit: deseq2_common.R =="
Rscript - <<'EOF' > $W/unit.out 2>&1
source("deseq2_common.R")
p <- 0; f <- 0
t <- function(l, c) if (isTRUE(c)) { p <<- p+1; cat("  PASS  ",l,"\n") } else { f <<- f+1; cat("  FAIL  ",l,"\n") }
x <- c(-3,-1,-.2,0,.2,1,3)
t("greaterAbs", identical(make_lfc_pass("greaterAbs",1)(x), c(TRUE,TRUE,FALSE,FALSE,FALSE,TRUE,TRUE)))
t("lessAbs is non-empty", sum(make_lfc_pass("lessAbs",1)(x)) == 3)
t("greater excludes wrong sign", !any(make_lfc_pass("greater",1)(c(-3,-2))))
t("less excludes wrong sign", !any(make_lfc_pass("less",1)(c(3,2))))
t("thr=0 passes all", all(make_lfc_pass("greaterAbs",0)(x)))
t("NA never passes", identical(make_lfc_pass("greaterAbs",1)(c(NA_real_,2)), c(FALSE,TRUE)))
t("log2(1.5): 1.6-fold in", make_lfc_pass("greaterAbs",0.585)(log2(1.6)))
t("log2(1.5): 1.4-fold out", !make_lfc_pass("greaterAbs",0.585)(log2(1.4)))
bal <- data.frame(condition=factor(c("P","P","B","B")), batch=factor(c(1,2,1,2)))
con <- data.frame(condition=factor(c("P","P","B","B")), batch=factor(c(1,1,2,2)))
t("balanced full rank", isTRUE(design_rank_report(bal, ~batch+condition)$ok))
t("confounded flagged", !isTRUE(design_rank_report(con, ~batch+condition)$ok))
t("dropping batch rescues", isTRUE(design_rank_report(con, ~condition)$ok))
t("sanitize IFN-gamma", sanitize_level("IFN-gamma") == "IFN_gamma")
t("sanitize leading digit", sanitize_level("48h/72h") == "X48h_72h")
t("sanitize keeps dots", sanitize_level("A.1_b") == "A.1_b")
o1 <- list(alpha="0.05", ref="grch38"); o2 <- list(ref="grch38", alpha="0.05")
t("signature ignores flag order", identical(signature_for(o1,4), signature_for(o2,4)))
t("alpha NOT in stage-3 signature", identical(signature_for(list(ref="a",alpha="0.05"),3),
                                              signature_for(list(ref="a",alpha="0.9"),3)))
t("alpha IS in stage-4 signature", !identical(signature_for(list(ref="a",alpha="0.05"),4),
                                              signature_for(list(ref="a",alpha="0.9"),4)))
t("ref IS in stage-1 signature", !identical(signature_for(list(ref="hg19"),1),
                                            signature_for(list(ref="grch38"),1)))
a <- parse_args(c("--ref","grch38","--alpha=0.01","--dry-run","--group1","B10"))
t("parse --key value", identical(a$ref,"grch38"))
t("parse --key=value", identical(a$alpha,"0.01"))
t("parse boolean", isTRUE(a[["dry-run"]]))
t("boolean consumes nothing", identical(a$group1,"B10"))
cat("UNIT:",p,"pass",f,"fail\n"); if (f>0) quit(status=1)
EOF
sed -n '/PASS\|FAIL/p' $W/unit.out
ok "unit tests all pass" 'grep -q "0 fail" $W/unit.out'

echo
echo "== stage 1: sample assignment =="
rc=$(cap $S1 --group1 Bottom10 --group2 Parental --control Parental --outdir $W/o1 --dry-run)
ok "keyword scan works"                '[[ $rc -eq 0 ]] && has "DRY RUN OK"'
ok "log tagged with stage"             'grep -q "\[1/setup\]" $W/o1/run.log'
ok "config table written"              '[[ -f $W/o1/tables/00_run_config.csv ]]'
ok "sample table written"              '[[ -f $W/o1/tables/01_sample_table.csv ]]'
ok "TEMPLATE written (no CSV ritual)"  '[[ -f $W/o1/tables/01c_samplesheet_template.csv ]]'
ok "template has real sample names"    'grep -q "Parental_rep1" $W/o1/tables/01c_samplesheet_template.csv'
ok "template batch column EMPTY"       '[[ $(tail -n +2 $W/o1/tables/01c_samplesheet_template.csv | cut -d, -f3 | tr -cd "A-Za-z0-9" | wc -c) -eq 0 ]]' 
ok "template carries path"             'grep -q "_genecounts.tsv" $W/o1/tables/01c_samplesheet_template.csv'

cat > $W/sheet.csv <<EOF
sample,condition,batch
Parental_rep1,Parental,A
Parental_rep2,Parental,B
Parental_rep3,Parental,C
Bottom10_rep1,Bottom10,A
Bottom10_rep2,Bottom10,B
Bottom10_rep3,Bottom10,C
EOF
rc=$(cap $S1 --samplesheet $W/sheet.csv --control Parental --batch batch --outdir $W/o2 --dry-run)
ok "sample sheet + balanced batch"     '[[ $rc -eq 0 ]] && has "FULL RANK"'
ok "batch x condition table written"   '[[ -f $W/o2/tables/01b_covariate_by_condition.csv ]]'
ok "covariate typing logged"           'has "batch=factor(3"'

cat > $W/conf.csv <<EOF
sample,condition,batch
Parental_rep1,Parental,kitA
Parental_rep2,Parental,kitA
Parental_rep3,Parental,kitA
Bottom10_rep1,Bottom10,kitB
Bottom10_rep2,Bottom10,kitB
Bottom10_rep3,Bottom10,kitB
EOF
rc=$(cap $S1 --samplesheet $W/conf.csv --control Parental --batch batch --outdir $W/o3 --dry-run)
ok "confounded batch REFUSED"          '[[ $rc -eq 2 ]] && has "NOT full rank"'
ok "error offers collapse-technical"   'has "collapse-technical"'
rc=$(cap $S1 --samplesheet $W/conf.csv --control Parental --outdir $W/o3b --dry-run)
ok "same data without --batch is fine" '[[ $rc -eq 0 ]]'

echo
echo "== stage 1: --batch-map (no CSV needed) =="
rc=$(cap $S1 --group1 Bottom10 --group2 Parental --control Parental \
       --batch-map "rep1=A,rep2=B,rep3=C" --outdir $W/o4 --dry-run)
ok "--batch-map creates the column"    '[[ $rc -eq 0 ]] && has "batch-map"'
ok "--batch-map design is full rank"   'has "FULL RANK"'
ok "--batch-map counts each pattern"   'has "rep1=A(n=2)"'
rc=$(cap $S1 --group1 Bottom10 --group2 Parental --control Parental \
       --batch-map "rep1=A" --outdir $W/o5 --dry-run)
ok "--batch-map leaving gaps REFUSED"  '[[ $rc -eq 2 ]] && has "unassigned"'
rc=$(cap $S1 --group1 Bottom10 --group2 Parental --control Parental \
       --batch-map "nomatch=A,rep1=B,rep2=C,rep3=D" --outdir $W/o6 --dry-run)
ok "--batch-map dead pattern REFUSED"  '[[ $rc -eq 2 ]] && has "matched no sample"'
rc=$(cap $S1 --group1 Bottom10 --group2 Parental --control Parental \
       --batch-map "repA" --outdir $W/o7 --dry-run)
ok "--batch-map bad syntax REFUSED"    '[[ $rc -eq 2 ]] && has "PATTERN=LABEL"'

echo
echo "== stage 1: design guards =="
rc=$(cap $S1 --samplesheet $W/sheet.csv --control Parental --design "~ batch*condition" --outdir $W/o8 --dry-run)
ok "INTERACTION design REFUSED"        '[[ $rc -eq 2 ]] && has "interaction term"'
ok "interaction error gives the fix"   'has "combine the factors into ONE"'
rc=$(cap $S1 --samplesheet $W/sheet.csv --control Parental --design "~ batch:condition" --outdir $W/o8b --dry-run)
ok "colon interaction also refused"    '[[ $rc -eq 2 ]]'
rc=$(cap $S1 --samplesheet $W/sheet.csv --control Parental --design "~ batch" --outdir $W/o9 --dry-run)
ok "design without condition refused"  '[[ $rc -eq 2 ]] && has "term named"'
rc=$(cap $S1 --samplesheet $W/sheet.csv --control Parental --design "~ rin + condition" --outdir $W/o10 --dry-run)
ok "unknown design term refused"       '[[ $rc -eq 2 ]] && has "not in the sample table"'

cat > $W/num.csv <<EOF
sample,condition,rin
Parental_rep1,Parental,7.1
Parental_rep2,Parental,8.2
Parental_rep3,Parental,9.0
Bottom10_rep1,Bottom10,7.4
Bottom10_rep2,Bottom10,8.8
Bottom10_rep3,Bottom10,9.3
EOF
rc=$(cap $S1 --samplesheet $W/num.csv --control Parental --design "~ rin + condition" --outdir $W/o11 --dry-run)
ok "covariate factored by DEFAULT"     'has "rin=factor"' 
ok "continuous-as-factor unfittable"   '[[ $rc -eq 2 ]] && has "NOT full rank"' 
rc=$(cap $S1 --samplesheet $W/num.csv --control Parental --design "~ rin + condition" \
       --numeric-covars rin --outdir $W/o12 --dry-run)
ok "--numeric-covars keeps numeric"    '[[ $rc -eq 0 ]] && has "rin=numeric"'
ok "continuous -> outlier NOTE"        'has "does NOT perform automatic outlier"'

cat > $W/nab.csv <<EOF
sample,condition,batch
Parental_rep1,Parental,A
Parental_rep2,Parental,
Parental_rep3,Parental,C
Bottom10_rep1,Bottom10,A
Bottom10_rep2,Bottom10,B
Bottom10_rep3,Bottom10,C
EOF
rc=$(cap $S1 --samplesheet $W/nab.csv --control Parental --batch batch --outdir $W/o13 --dry-run)
ok "blank batch cell REFUSED by name"  '[[ $rc -eq 2 ]] && has "Parental_rep2"'

echo
echo "== stage 1: contrast + replicates + numeric validation =="
cat > $W/three.csv <<EOF
sample,condition
Parental_rep1,Parental
Parental_rep2,Parental
Bottom10_rep1,Bottom10
Bottom10_rep2,Bottom10
UNC_rep1,UNC
UNC_rep2,UNC
EOF
rc=$(cap $S1 --samplesheet $W/three.csv --control Parental --outdir $W/o14 --dry-run)
ok ">2 conditions REQUIRE --treatment" '[[ $rc -eq 2 ]] && has "--treatment is required"'
rc=$(cap $S1 --samplesheet $W/three.csv --control Parental --treatment UNC --outdir $W/o15 --dry-run)
ok ">2 conditions with --treatment"    '[[ $rc -eq 0 ]] && has "contrast: UNC vs Parental"'
ok "all 6 samples kept by default"     'has "resolved 6 samples"'
rc=$(cap $S1 --samplesheet $W/three.csv --control Parental --treatment UNC \
       --subset-to-contrast --outdir $W/o16 --dry-run)
ok "--subset-to-contrast drops others" '[[ $rc -eq 0 ]] && has "resolved 4 samples"'
ok "subset reports what it dropped"    'has "dropping 2 sample(s)"'

cat > $W/one.csv <<EOF
sample,condition
Parental_rep1,Parental
Bottom10_rep1,Bottom10
EOF
rc=$(cap $S1 --samplesheet $W/one.csv --control Parental --outdir $W/o17 --dry-run)
ok "n=1 per group refused"             '[[ $rc -eq 2 ]] && has ">= 2 replicates"'

cat > $W/dash.csv <<EOF
sample,condition
Parental_rep1,IFN-gamma
Parental_rep2,IFN-gamma
Bottom10_rep1,siRNA #2
Bottom10_rep2,siRNA #2
EOF
rc=$(cap $S1 --samplesheet $W/dash.csv --control "IFN-gamma" --outdir $W/o18 --dry-run)
ok "non-conforming levels sanitised"   '[[ $rc -eq 0 ]] && has "IFN_gamma"'
ok "--control sanitised to match"      'has "contrast: siRNA_2 vs IFN_gamma"'
ok "rename is reported"                'has "DESeq2 factor levels allow only"'

rc=$(cap $S1 --samplesheet $W/sheet.csv --control Parental --min-count ten --outdir $W/o19 --dry-run)
ok "garbage --min-count refused"       '[[ $rc -eq 2 ]] && has "whole number"'
rc=$(cap $S1 --samplesheet $W/sheet.csv --control Parental --threads 8x --outdir $W/o20 --dry-run)
ok "garbage --threads refused"         '[[ $rc -eq 2 ]]'
rc=$(cap $S1 --samplesheet $W/sheet.csv --control Parental --bogus 1 --outdir $W/o21 --dry-run)
ok "unknown flag refused"              '[[ $rc -eq 2 ]] && has "unknown option"'

echo
echo "== stage 2: counts + QC =="
rc=$(cap $S1 --samplesheet $W/sheet.csv --control Parental --batch batch --outdir $W/run)
ok "stage 1 full run"                  '[[ $rc -eq 0 ]] && has "STAGE 1 COMPLETE"'
rc=$(cap $S2 --samplesheet $W/sheet.csv --control Parental --batch batch --outdir $W/run --no-annotation)
ok "stage 2 runs"                      '[[ $rc -eq 0 ]] && has "STAGE 2 COMPLETE"'
ok "log APPENDED not overwritten"      'grep -q "\[1/setup\]" $W/run/run.log && grep -q "\[2/counts\]" $W/run/run.log'
ok "run.log survives GC (no early close)" '[[ $(grep -c "\[2/counts\]" $W/run/run.log) -ge 5 ]]' 
ok "run.log has stage-4-style late lines"  '[[ $(wc -l < $W/run/run.log) -ge 20 ]]' 
ok "assignment QC table"               '[[ -f $W/run/tables/02_count_assignment_qc.csv ]]'
ok "strandedness read back"            '[[ -f $W/run/tables/02b_strandedness_qc.csv ]]'
ok "stranded_used surfaced"            'has "counted upstream with -s: reverse"'
ok "count matrix built"                'has "200 features x 6 samples"'
ok "ambiguity threshold reported"      'has "ambiguity"'
ok "denominator caveat in table"       'grep -q "NOT fragments sequenced" $W/run/tables/02_count_assignment_qc.csv'

for s in Parental_rep1 Parental_rep2 Parental_rep3 Bottom10_rep1 Bottom10_rep2 Bottom10_rep3; do mk "$s" 90000; done
rc=$(cap $S2 --samplesheet $W/sheet.csv --control Parental --batch batch --outdir $W/run --no-annotation)
ok "high ambiguity WARNS (>15%)"       'has "Unassigned_Ambiguity"'
ok "warning explains hg19 gene_id"     'has "TRANSCRIPT accession"'
for s in Parental_rep1 Parental_rep2 Parental_rep3 Bottom10_rep1 Bottom10_rep2 Bottom10_rep3; do mk "$s" 4000; done

echo
echo "== the stage-signature guard (the risk of splitting) =="
rc=$(cap $S1 --samplesheet $W/sheet.csv --control Parental --batch batch --outdir $W/run)
rc=$(cap $S2 --samplesheet $W/sheet.csv --control Parental --batch batch --outdir $W/run --no-annotation)
ok "matching signature -> stage 2 OK"  '[[ $rc -eq 0 ]]'
rc=$(cap $S2 --samplesheet $W/sheet.csv --control Bottom10 --batch batch --outdir $W/run --no-annotation)
ok "changed --control -> stage 2 STOPS" '[[ $rc -eq 2 ]] && has "different arguments"'
ok "guard names the changed flag"      'has "control="'
rc=$(cap $S2 --samplesheet $W/sheet.csv --control Bottom10 --batch batch --outdir $W/run --no-annotation --force)
ok "--force overrides the guard"       '[[ $rc -eq 0 ]] && has "WARNING"'
rc=$(cap Rscript deseq2_03_model.R --samplesheet $W/sheet.csv --control Parental --outdir $W/nostate)
ok "missing upstream state explained"  '[[ $rc -eq 2 ]] && has "Run the earlier stage"'

echo
echo "== stage 3/4 argument guards (run before library(DESeq2)) =="
rc=$(cap Rscript deseq2_03_model.R --samplesheet $W/sheet.csv --control Parental --batch batch \
       --outdir $W/run --no-annotation --test bogus)
ok "bad --test refused"                '[[ $rc -eq 2 ]] && has "wald or lrt"'
rc=$(cap Rscript deseq2_03_model.R --samplesheet $W/sheet.csv --control Parental --batch batch \
       --outdir $W/run --no-annotation --fit-type bogus)
ok "bad --fit-type refused"            '[[ $rc -eq 2 ]] && has "fit-type invalid"'
ok "glmGamPoi is accepted"             'grep -q "glmGamPoi" deseq2_03_model.R'

rc=$(cap $S1 --samplesheet $W/sheet.csv --control Parental --batch batch --outdir $W/run)
rc=$(cap $S2 --samplesheet $W/sheet.csv --control Parental --batch batch --outdir $W/run --no-annotation)
G4="Rscript deseq2_04_results.R --samplesheet $W/sheet.csv --control Parental --batch batch --outdir $W/run --no-annotation"
rc=$(cap $G4 --alt-hypothesis lessAbs)
ok "lessAbs w/o threshold HARD-STOPS"  '[[ $rc -eq 2 ]] && has "EMPTY SET"'
ok "lessAbs error explains equivalence" 'has "equivalence test"'
rc=$(cap $G4 --alt-hypothesis lessAbs --lfc-threshold 0.585)
ok "lessAbs WITH threshold passes args" '[[ $rc -ne 2 ]] || ! has "EMPTY SET"'
rc=$(cap $G4 --rank-on shrunk --shrink none)
ok "--rank-on shrunk + none refused"   '[[ $rc -eq 2 ]] && has "requires a --shrink"'
rc=$(cap $G4 --alpha 1.5)
ok "alpha out of range refused"        '[[ $rc -eq 2 ]] && has "between"'
rc=$(cap $G4 --shrink bogus)
ok "bad --shrink refused"              '[[ $rc -eq 2 ]] && has "shrink invalid"'
ok "ashr excluded when thresholded"    'grep -q "ashr" deseq2_04_results.R && grep -q "setdiff" deseq2_04_results.R' 
ok "lfcShrink gets res= (alpha carry)" 'grep -q "res = res" deseq2_04_results.R'
ok "svalue requested when thresholded" 'grep -q "a\$svalue <- TRUE" deseq2_04_results.R'
ok "svalue carried into results table" 'grep -q "add\$svalue <- sh\$svalue" deseq2_04_results.R'
ok "plotMA gets alpha explicitly"      '[[ $(grep -c "plotMA(.*alpha=alpha" deseq2_04_results.R) -eq 2 ]]' 
ok "altHypothesis forwarded one-sided" 'grep -q "one-sided test requested" deseq2_04_results.R'

echo
echo "== submitter =="
SUB=./submit_rnaseq_deseq2_dge.sh
export DESEQ2_SKIP_ENV=1
rc=$(cap bash $SUB --help)
ok "--help works"                      '[[ $rc -eq 0 ]]'
rc=$(cap bash $SUB --stage 9 --counts-dir $W/counts)
ok "bad --stage refused"               '[[ $rc -ne 0 ]] && has "must be 1|2|3|4"'
rc=$(cap bash $SUB --stage 3-2 --counts-dir $W/counts)
ok "backwards range refused"           '[[ $rc -ne 0 ]] && has "backwards"'
rc=$(cap bash $SUB --stage 3 --dry-run --counts-dir $W/counts)
ok "--dry-run from stage 3 refused"    '[[ $rc -ne 0 ]] && has "only makes sense"'
rc=$(cap bash $SUB --samplesheet $W/sheet.csv --counts-dir $W/counts --control Parental \
       --batch batch --outdir $W/sub --dry-run)
ok "dry-run runs stage 1 only"         '[[ $rc -eq 0 ]] && has "STAGE 1/setup"'
ok "dry-run does NOT reach stage 2"    '! has "STAGE 2/counts"'
ok "submitter forwards flags"          'has "DRY RUN OK"'
rc=$(cap bash $SUB --stage 1-2 --samplesheet $W/sheet.csv --counts-dir $W/counts \
       --control Parental --batch batch --outdir $W/sub2 --no-annotation)
ok "--stage 1-2 runs both"             '[[ $rc -eq 0 ]] && has "STAGE 1/setup" && has "STAGE 2/counts"'
ok "stage range stops before 3"        '! has "STAGE 3/model"'
ok "reports per-stage timing"          'has "finished in"'
rc=$(cap bash $SUB --stage 1 --samplesheet $W/sheet.csv --counts-dir /nonexistent \
       --control Parental --outdir $W/sub3)
ok "failure names the stage"           '[[ $rc -ne 0 ]] && has "FAILED at stage 1/setup"' 
ok "failure suggests how to resume"    'has "--stage 1-"' 

echo
echo "----------------------------------------"
echo "PASS: $PASS   FAIL: $FAIL"
[[ $FAIL -eq 0 ]]