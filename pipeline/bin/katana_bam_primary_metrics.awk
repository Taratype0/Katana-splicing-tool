# Exact aggregate metrics over samtools view -F 2304; no sequence or QNAME stored.
# Portable POSIX awk: arithmetic flag decoding instead of implementation bit ops.
BEGIN {
    FS = "\t"; OFS = "\t"
    nmetrics = split("primary_records paired_records unpaired_records read1_records read2_records invalid_pair_flag_records mapped_records unmapped_records missing_seq_records invalid_sequence_records malformed_cigar_records cigar_query_length_mismatch_records mapped_without_cigar_records soft_clipped_records hard_clipped_records duplicate_records qcfail_records spliced_records", metric_names, " ")
    for (m = 0; m < 2; m++) for (m2 = 0; m2 < 2; m2++) for (side = 1; side <= 2; side++)
        metric_names[++nmetrics] = (m ? "unmapped" : "mapped") "_mate_" (m2 ? "unmapped" : "mapped") "_read" side "_records"
    for (i = 1; i <= nmetrics; i++) metrics[metric_names[i]] = 0
}
function hasflag(flag, bit) { return int(flag / bit) % 2 }
function fail(message) { print message > "/dev/stderr"; failed = 1; exit 1 }
{
    if (substr($0, 1, 1) == "@" || NF == 0) next
    if (NF < 11 || $2 !~ /^[0-9]+$/ || $2 + 0 > 65535) fail("Malformed SAM record/FLAG at primary stream line " NR)
    flag = $2 + 0
    if (hasflag(flag, 256) || hasflag(flag, 2048)) next
    metrics["primary_records"]++
    paired = hasflag(flag, 1); first = hasflag(flag, 64); second = hasflag(flag, 128); unmapped = hasflag(flag, 4); mate_unmapped = hasflag(flag, 8)
    if (paired) metrics["paired_records"]++; else metrics["unpaired_records"]++
    metrics["read1_records"] += first; metrics["read2_records"] += second
    if ((paired && first == second) || (!paired && (hasflag(flag, 2) || mate_unmapped || hasflag(flag, 32) || first || second)) || (hasflag(flag, 2) && (unmapped || mate_unmapped))) metrics["invalid_pair_flag_records"]++
    if (paired && first != second) metrics[(unmapped ? "unmapped" : "mapped") "_mate_" (mate_unmapped ? "unmapped" : "mapped") "_read" (first ? "1" : "2") "_records"]++
    if (unmapped) metrics["unmapped_records"]++; else { metrics["mapped_records"]++; references[$3]++ }
    metrics["duplicate_records"] += hasflag(flag, 1024); metrics["qcfail_records"] += hasflag(flag, 512)
    sequence = $10; cigar = $6
    if (sequence == "*") metrics["missing_seq_records"]++
    else {
        lengths[length(sequence)]++
        if (sequence !~ /^[ACMGRSVTWYHKDBNacmgrsvtwyhkdbn]+$/) metrics["invalid_sequence_records"]++
    }
    if (cigar == "*") { if (!unmapped) metrics["mapped_without_cigar_records"]++; next }
    rest = cigar; query = 0; soft = 0; hard = 0; spliced = 0; malformed = 0; num_ops = 0
    while (length(rest) > 0) {
        if (!match(rest, /^[1-9][0-9]*[MIDNSHP=X]/)) { malformed = 1; break }
        token = substr(rest, RSTART, RLENGTH); op = substr(token, length(token), 1); count = substr(token, 1, length(token)-1) + 0
        ops[++num_ops] = op
        if (op ~ /^[MIS=X]$/) query += count
        if (op == "S") soft = 1; if (op == "H") hard = 1; if (op == "N") spliced = 1
        rest = substr(rest, RLENGTH + 1)
    }
    for (k = 1; k <= num_ops && !malformed; k++) {
        if (ops[k] == "H" && k != 1 && k != num_ops) malformed = 1
        if (ops[k] == "S") {
            before_h_only = 1; after_h_only = 1
            for (j = 1; j < k; j++) if (ops[j] != "H") before_h_only = 0
            for (j = k + 1; j <= num_ops; j++) if (ops[j] != "H") after_h_only = 0
            if (!before_h_only && !after_h_only) malformed = 1
        }
    }
    if (malformed) { metrics["malformed_cigar_records"]++; next }
    metrics["soft_clipped_records"] += soft; metrics["hard_clipped_records"] += hard; metrics["spliced_records"] += spliced
    if (sequence != "*" && query != length(sequence)) metrics["cigar_query_length_mismatch_records"]++
}
END {
    if (failed) exit 1
    print "kind", "name", "count"
    # Explicit integer formatting avoids scientific notation for large cohorts.
    for (i = 1; i <= nmetrics; i++) printf "metric\t%s\t%.0f\n", metric_names[i], metrics[metric_names[i]]
    for (i in lengths) printf "read_length\t%s\t%.0f\n", i, lengths[i]
    for (i in references) printf "mapped_reference\t%s\t%.0f\n", i, references[i]
}
