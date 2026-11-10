process PAIRS2BED {
    tag "$meta.id"
    label 'process_low'

    // same gawk environment and images as the nf-core gawk module
    conda "conda-forge::gawk=5.3.1"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a1/a125c778baf3865331101a104b60d249ee15fe1dca13bdafd888926cc5490a34/data' :
        'community.wave.seqera.io/library/gawk:5.3.1--e09efb5dfc4b8156' }"

    input:
    tuple val(meta), path(pairs), path(index), path(ligation_pairs)
    val read_length

    output:
    tuple val(meta), path("*.reads.bed.gz"), emit: bed
    path("versions.yml"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def lig = ligation_pairs instanceof List ? ligation_pairs : [ligation_pairs]
    // The reads FitHiChIP's PeakInferHiChIP.sh (Imp_Scripts/, FitHiChIP 11.0) gives MACS to infer the 1D ChIP peaks of a HiChIP library:
    //   - the HiC-Pro dangling-end, self-circle and re-ligation pairs (*.DEPairs, *.SCPairs, *.REPairs; HiC-Pro processing), and
    //   - the valid pairs whose two ends are on one chromosome less than 1 kb apart.
    // Without HiC-Pro's ligation classes (pairtools processing, or --no_digestion), those pairs are among the valid pairs and the
    // short-range valid pairs are the whole input. Every end becomes a read_length interval from its 5' end in the read direction,
    // without a strand column, as in PeakInferHiChIP; intervals starting before the chromosome start are dropped.
    // HiC-Pro pairs columns: readID chr1 pos1 strand1 chr2 pos2 strand2 ...; .pairs: readID chrom1 pos1 chrom2 pos2 strand1 strand2 ...
    def ends = 'if (s1 == "+") print c1, p1, p1 + R - 1; else print c1, p1 - R + 1, p1; if (s2 == "+") print c2, p2, p2 + R - 1; else print c2, p2 - R + 1, p2'
    """
    {
        for f in ${lig.join(' ')}; do
            awk -v R=${read_length} 'BEGIN{OFS="\\t"} !/^#/ { c1 = \$2; p1 = \$3; s1 = \$4; c2 = \$5; p2 = \$6; s2 = \$7; ${ends} }' "\$f"
        done
        gzip -dc ${pairs} \\
            | awk -v R=${read_length} 'BEGIN{OFS="\\t"} !/^#/ && \$2 == \$4 && (\$5 > \$3 ? \$5 - \$3 : \$3 - \$5) < 1000 {
                c1 = \$2; p1 = \$3; c2 = \$4; p2 = \$5; s1 = \$6; s2 = \$7; ${ends} }'
    } | awk '\$2 >= 0' | gzip -c > ${prefix}.reads.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(awk -Wversion | sed '1!d; s/.*Awk //; s/,.*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.reads.bed.gz
    touch versions.yml
    """
}
