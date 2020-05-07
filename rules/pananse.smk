
rule sort_bdg:
    """
    Sort bdg file from macs2.
    """
    input:
        bdgfile=expand("{result_dir}/macs2/{{assembly}}-{{sample}}_treat_pileup.bdg", **config)
    output:
        sortbdgfile=expand("{result_dir}/sort_bdg/{{assembly}}-{{sample}}_sort.bdg", **config)
    shell:
        "sort -k1,1 -k2,2n {input.bdgfile} > {output.sortbdgfile}"


rule bdg2wig:
    """
    Switch bdg file to wig file.
    """
    input:
        sortbdgfile=expand("{result_dir}/sort_bdg/{{assembly}}-{{sample}}_sort.bdg", **config),
        genome_size=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config)
    output:
        wigfile=expand("{result_dir}/macs2_wig/{{assembly}}-{{sample}}.wig", **config)
    shell:
        "/home/qxu/bin/ucsc/bedGraphToBigWig {input.sortbdgfile} {input.genome_size} {output.wigfile}"


def read_chrsize(sizefile):
    sdic={}
    with open (sizefile) as sizes:
        for i in sizes:
            sdic[i.split()[0]]=int(i.split()[1])
    return sdic


rule make_enhancer:
    input:
        narrowpeak=expand("{result_dir}/macs2/{{assembly}}-{{sample}}_peaks.narrowPeak", **config),
        wigfile=expand("{result_dir}/macs2_wig/{{assembly}}-{{sample}}.wig", **config),
        genome_size=expand("{genome_dir}/{{assembly}}/{{assembly}}.fa.sizes", **config)
    output:
        enhancerbed=expand("{result_dir}/enhancer/{{assembly}}-{{sample}}.bed", **config)
    run:
        chrsizedic=read_chrsize(input.genome_size[0])
        with open(str(input.narrowpeak)) as bdgf, open(str(output.enhancerbed),"w") as ebed:
            for line in bdgf:
                a=line.split()
                p=int(a[9])+int(a[1])
                if int(p-100)>0 and chrsizedic[a[0]]>p+100:
                    commd1=os.popen("/home/qxu/bin/ucsc/bigWigSummary -type=max " + str(input.wigfile) + " " + a[0] + " " + str(p-100) + " " + str(p+100) + " 1")
                    r=commd1.read()
                    if r != "":
                        ebed.write(a[0]+"\t"+str(p-100)+"\t"+str(p+100)+"\t"+str(r))

