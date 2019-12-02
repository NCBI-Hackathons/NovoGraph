#!/usr/bin/env nextflow

if (!(workflow.profile in ['local', 'hilbert'])) {
  println "ERROR: You need to set -profile (values: local, hilbert)"
  exit 1
}
println params.genome

//params.basedir='/home/houwaart/Repositories/NovoGraph/pipeline'
params.mafft_dir="$baseDir/forMAFFT/"
params.aligner='mafft'

Channel
    .fromPath(params.mafft_dir + "*.fa")
    .ifEmpty { error "Cannot find any data -- Check the path specified: `${params.mafft_dir}`" }
    .set { file_names }

process index_ref {

    input:
    file x from Channel.fromPath(params.genome)

    output:
    file "${x}.fai" into index_ch

    when:
    !(file(params.genome + ".fai")).exists()

    script:
    if (0) // don't need this ATM
        """
        samtools faidx ${x}
        """
    else
        """
        echo "all good"
        """
}


process callmafft {
    storeDir params.mafft_dir
    input:
        file(in_file) from file_names
    output:
        file "*.msa" into msa_files mode flatten
    shell:
    '''
    tmp_name=`basename !{in_file} | awk -F. '{print $1}'`
    out_name=${tmp_name}".msa"
    mafft !{in_file} > $out_name
    '''
}
