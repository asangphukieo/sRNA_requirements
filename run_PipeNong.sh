#!/bin/bash


INPUT_FOLDER="./Sa/Sa_HG001_sRNAs_for_GLASSgo/" 


# -----------
# run GLASSgo
# -----------

GLASSgo_OUT_FOLDER='01_GLASSgo_Results'

:<<"COMMENT"
mkdir -p $GLASSgo_OUT_FOLDER

source ./activate apiwatenv
#source ./activate aligners

for f in `ls $INPUT_FOLDER`; do
    echo "Running GLASSgo for sRNA: " $f
    fullpath=${INPUT_FOLDER}/$f
    output_sRNAs=GLASSgo_output_${f}
    output_json=${output_sRNAs%.*}.json
    ./GLASSgo_1_5_1.py -d /home/asangphukieo/Downloads/workflow-1/2_GLASSgo/GLASSgo_GI-Version-master/glassgo/GLASSGO/BLAST_NT/blast_nt_db/nt -e 1 -i $fullpath -t 5 -o $output_sRNAs 
    ./createJSON-1_0.py -g $output_sRNAs -o $output_json
    ./generate_GLASSgo_html.py $output_json
done
mv GLASSgo_output_* $GLASSgo_OUT_FOLDER
conda deactivate

COMMENT


# --------------------
# run GLASSgo2CopraRNA 
# --------------------

GLASSgo2CopraRNA="02_GLASSgo2CopraRNA"
:<<"COMMENT"
mkdir -p $GLASSgo2CopraRNA

source activate r342

for f in `ls $GLASSgo_OUT_FOLDER`; do
    if [[ $f =~ .*\.(fasta|fa|fas) ]]; then 
        echo "Running GLASSgo2CopraRNA for sRNA: " $f
        fullpath=$GLASSgo_OUT_FOLDER/$f
        filestem=${f%.*}
        ID=${filestem#GLASSgo_output_}
        # fasta generation
        R --slave -f  GLASSgo2CopraRNA_fasta_generation_by_accession.r --args input_file=$fullpath refpath=taxid_to_refseq cop_path=CopraRNA_available_organisms.txt output_file=${ID}_coprarna_candidates.txt
        # exclusion 
        R --slave -f  GLASSgo2CopraRNA_exclusion_script.r --args datapath=full_GLASSgo_table.Rdata wildcard=NC_007795,NC_004461,NC_013893,NC_020164,NZ_CP007601,NC_014925,NC_007350 ooi=NZ_CP018205
        # balanced selection
        R --slave -f  GLASSgo2CopraRNA_balanced_ooi_selection.r --args wildcard=NC_007795,NC_004461,NC_013893,NC_020164,NZ_CP007601,NC_014925,NC_007350 max_number=20 outfile_prefix=${ID}_sRNA ooi=NZ_CP018205
        # move results to its folder
        mkdir -p $GLASSgo2CopraRNA/$ID
        mv ${ID}_* $GLASSgo2CopraRNA/$ID
    fi
done
conda deactivate
#source deactivate

COMMENT


# ------------
# run CopraRNA XXXXXXXXXXXXXXXXXXXXXXX
# ------------

CopraRNA_OUT_FOLDER="03_CopraRNA_Results"
mkdir -p $CopraRNA_OUT_FOLDER

:<<"COMMENT" This approach work for me
# parallel run
source activate /home/asangphukieo/Documents/anaconda2/envs/JensCopraRNA2.yml
./parallel_CopraRNA_prediction.py $GLASSgo2CopraRNA --suffix .fasta --cores 6 --batch 20 --out_folder $CopraRNA_OUT_FOLDER \
              --ntup 200 --ntdown 100 --region 5utr --enrich 200 --topcount 200 --websrv --noclean 2>&1 | tee ./parallel_CopraRNA_log_file.txt
conda deactivate
COMMENT


# simplify CopraRNA_results
:<<"COMMENT"
python ./simplify_CopraRNA_results.py $CopraRNA_OUT_FOLDER ./Sa/Refseq2HG001.tsv
source ./activate r342
simplified_folder=03_CopraRNA_Results_simplified
for f in `ls $simplified_folder/*_locustag.txt`; do
    # GO enrichment
    R --slave -f GOstats_analysis.R --args workdir=. organism=Staphylococcus_HG001 target_gene_file=$f ipr2go_file=NZ_CP018205_ipr2go.tsv 
    # fetch revigo csv
    filename=`basename $f`
    filestem=${filename%.txt}
    ./fetch_revigo_csv.py $simplified_folder/${filestem}_GO_enrichment_BP.tsv -o $simplified_folder
    ./fetch_revigo_csv.py $simplified_folder/${filestem}_GO_enrichment_CC.tsv -o $simplified_folder
    ./fetch_revigo_csv.py $simplified_folder/${filestem}_GO_enrichment_MF.tsv -o $simplified_folder
    # revigo visualization
    R --slave -f REVIGO_plotter.R --args workdir=. REVIGO_BP=$simplified_folder/${filestem}_GO_enrichment_BP_revigo.csv REVIGO_CC=$simplified_folder/${filestem}_GO_enrichment_CC_revigo.csv REVIGO_MF=$simplified_folder/${filestem}_GO_enrichment_MF_revigo.csv out_pdf=$simplified_folder/${filestem}_GO_enrichment.pdf
done
conda deactivate
COMMENT

# -----------
# run mafft
# -----------

MAFFT_OUT_FOLDER="04_Mafft_Results"
:<<"COMMENT"
mkdir -p $MAFFT_OUT_FOLDER

source activate apiwatenv
#source activate aligners
unset MAFFT_BINARIES
MVIEW="/home/asangphukieo/Downloads/ShengWei_work/mview-1.64/bin/mview"
for ID_Folder in `ls $GLASSgo2CopraRNA`; do
    for f in `ls $GLASSgo2CopraRNA/$ID_Folder`; do
        if [[ $f =~ .*\.(fasta|fa|fas) ]]; then 
            echo "Running mafft-linsi for sRNA: " $f
            fullpath=$GLASSgo2CopraRNA/$ID_Folder/$f
            filestem=${f%.*}
            ID=${filestem#GLASSgo_output_}
            # run mafft-linsi
            mafft --localpair --maxiterate 1000 --reorder --thread 20 $fullpath > ${ID}_mafft.fasta
            # run mview
            $MVIEW -in fasta -html head -css on -coloring id -colormap clustal ${ID}_mafft.fasta > ${ID}_mafft.html
            # run fasttree
            cat ${ID}_mafft.fasta | tr : _ > ${ID}_modified.fa       
            fasttree -quiet -gtr -gamma -nt ${ID}_modified.fa > ${ID}_fasttree.nwk
            rm ${ID}_modified.fa
            #figtree  -graphic  PDF ${ID}_fasttree.nwk ${ID}_fasttree.pdf

            # move results to MAFFT folder
            mv ${ID}_mafft.* ${ID}_fasttree.* $MAFFT_OUT_FOLDER
        fi
    done
done
conda deactivate
COMMENT


# ---------------
# run RNAalifold
# ---------------

RNAalifold_OUT_FOLDER="05_RNAalifold_Results"


:<<"COMMENT"
mkdir -p $RNAalifold_OUT_FOLDER
source ./activate py27
for f in `ls $MAFFT_OUT_FOLDER`; do
    if [[ $f =~ .*\.(fasta|fa|fas) ]]; then 
        echo "Running RNAalifold for sRNA: " $f
        fullpath=$MAFFT_OUT_FOLDER/$f
        filestem=${f%.*}
        ID=${filestem%_mafft.fasta}
        # shorten header
        ./shorten_GLASSgo_fasta.py --force $fullpath -p $ID 
        # convert to sto and aln format
        seqmagick convert ${ID}_shorten.fa ${ID}_shorten.sto
        seqmagick convert ${ID}_shorten.fa ${ID}_shorten.aln
        # run RNAalifold and tidy output files
        RNAalifold --color --aln --sci ${ID}_shorten.aln > ${ID}_RNAalifold_output.txt
        mv aln.ps ${ID}_aln.ps
        mv alirna.ps ${ID}_alirna.ps
        ps2pdf -dEPSCrop ${ID}_aln.ps ${ID}_aln.pdf 
        ps2pdf -dEPSCrop ${ID}_alirna.ps ${ID}_alirna.pdf
        # prepare Rscape input 
        ./modify_sto_file_for_Rscape.py ${ID}_shorten.sto ${ID}_RNAalifold_output.txt
        mv ${ID}_* $RNAalifold_OUT_FOLDER
    fi
done
conda deactivate
COMMENT

# ------------
# run R-scape 
# ------------
Rscape_OUT_FOLDER="06_Rscape_Results"
mkdir -p $Rscape_OUT_FOLDER

RSCAPE="/home/asangphukieo/Downloads/ShengWei_work/rscape_v1.2.3/bin/R-scape"

:<<"COMMENT"
for f in `ls $RNAalifold_OUT_FOLDER`; do
    if [[ $f =~ .*\.sto ]]; then 
        echo "Running R-scape for sRNA: " $f
        fullpath=$RNAalifold_OUT_FOLDER/$f
        filestem=${f%.*}
        ID=${filestem%_mafft_shorten_for_Rscape.sto}
        # run Rscape
        $RSCAPE -E 0.1 --outmsa $Rscape_OUT_FOLDER/${ID}_Rscape_msa.fa -o $Rscape_OUT_FOLDER/${ID}_Rscape_output.txt --r2rall --outdir $Rscape_OUT_FOLDER $fullpath
    fi
done
COMMENT

# ------------
# run RNAcode 
# ------------
RNAcode_OUT_FOLDER="07_RNAcode_Results"

#RNACODE="/mnt/data/software/Module_ncRNA/RNAcode/bin/bin/RNAcode"

:<<"COMMENT"
source activate RNAtools
mkdir -p $RNAcode_OUT_FOLDER

for f in `ls $RNAalifold_OUT_FOLDER`; do
    if [[ $f =~ .*\.aln ]]; then 
        echo "Running RNAcode for sRNA: " $f
        fullpath=$RNAalifold_OUT_FOLDER/$f
        filestem=${f%.*}
        ID=${filestem%_sRNA_CopraRNA_input_balanced*.aln}
        # run RNAcode
        mkdir -p $ID
        RNAcode --gtf --tabular --eps --eps-dir $ID --eps-cutoff 0.1 --cutoff 0.1 --outfile ${ID}_RNAcode.tab $fullpath

        if [ -z "$(ls -A $ID)" ]; then
            echo "Protein coding regions were not found!"
            rm ${ID}_*
            rm -rf $ID
        else
            echo "Coding regions detected!"
            cd $ID
            for f in `ls .`; do
                if [[ $f =~ .*\.(eps|ps) ]]; then
                    ps2pdf -dEPSCrop $f
                fi
            done
            cd ../
            mv ${ID}_* $ID
            mv $ID $RNAcode_OUT_FOLDER
        fi

    fi
done
conda deactivate
COMMENT


# ------------
# run synteny 
# ------------

Synteny_OUT_FOLDER="08_Synteny_Results"

:<<"COMMENT"
source ./activate r342
mkdir -p $Synteny_OUT_FOLDER

# copy sRNA to Synteny_OUT_FOLDER
echo "Copying GLASSgo output sRNAs to $Synteny_OUT_FOLDER ..."
cp $GLASSgo_OUT_FOLDER/*.{fa,fas,fasta} $Synteny_OUT_FOLDER
#rm $Synteny_OUT_FOLDER/*.fa_old
cd $Synteny_OUT_FOLDER
cp ../synteny_pdf.r .
for f in `ls .`; do
    if [[ $f =~ .*\.(fasta|fa|fas) ]]; then 
        echo "Running synteny for sRNA: " $f
        filestem=${f%.*}
        ID=${filestem#GLASSgo_output_}
        if [ ! -f ${ID}_synteny.pdf ]; then
            R --slave -f synteny_pdf.r --args input_sRNA=$f output_prefix=${ID} 2>&1 | tee ./${ID}_synteny_log.txt
        fi
    fi
done
conda deactivate
rm synteny_pdf.r
cd ../
COMMENT


# -------------------------------------
# run motif_finding in promoter regions 
# -------------------------------------

Promoter_OUT_FOLDER="09_Promoter_Results"

:<<"COMMENT"
mkdir -p $Promoter_OUT_FOLDER
source ./activate r342 
for f in `ls $GLASSgo_OUT_FOLDER`; do
    

    if [[ $f =~ .*\.(fasta|fa|fas) ]]; then 
        echo "Fetching promoter sequences for sRNA $f ... "
        fullpath=$GLASSgo_OUT_FOLDER/$f
        filestem=${f%.*}
        ID=${filestem#GLASSgo_output_}
        # fasta generation
        R --slave -f  GLASSgo2CopraRNA_fasta_generation_by_accession.r --args input_file=$fullpath refpath=taxid_to_refseq cop_path=CopraRNA_available_organisms.txt output_file=${ID}_coprarna_candidates.txt
        # exclusion 
        R --slave -f  GLASSgo2CopraRNA_exclusion_script.r --args datapath=full_GLASSgo_table.Rdata wildcard=NC_007795,NC_004461,NC_013893,NC_020164,NZ_CP007601,NC_014925,NC_007350 ooi=NZ_CP018205
        # fetch promoters
        R --slave -f promoter_sequence_fetcing_from_coor.r --args datapath=refined_GLASSgo_table.Rdata output_prefix=${ID}_promoter
        # move results to its folder
        mkdir -p $Promoter_OUT_FOLDER/$ID
        mv ${ID}_* $Promoter_OUT_FOLDER/$ID

        # do alignment and mview
        
        # do motif finding and weblogo 
        #source ./activate        ################### do someting
        #meme $Promoter_OUT_FOLDER/$ID/${ID}_promoter.fasta -oc $Promoter_OUT_FOLDER/$ID/MEME_result -dna -mod anr -nmotifs 50 -revcomp -evt 0.05 -minw 4 -maxw 200 -maxsize 10000000000000

    fi
done
conda deactivate
COMMENT

#:<<"COMMENT"
source ./activate Genomics
for ID_FOLDER in `ls $Promoter_OUT_FOLDER`; do
    meme $Promoter_OUT_FOLDER/$ID_FOLDER/*_promoter.fasta -oc $Promoter_OUT_FOLDER/$ID_FOLDER/MEME_results -dna -mod anr -nmotifs 50 -revcomp -evt 0.05 -minw 4 -maxw 200 -maxsize 10000000000000
done
conda deactivate
#COMMENT