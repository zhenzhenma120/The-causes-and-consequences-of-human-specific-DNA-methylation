#!/bin/bash

# submit_all_celltypes.sh

# Change coverage and minimum read cutoffs in script 06 (e.g. 5bp cutoff for RRBS and 10bp for WGBS)
CELLTYPES=("IPSC" "SKM" "HEP" "DA" "CNCC")

# script dir
#SCRIPT_DIR=""

echo "Submitting cis-trans analysis jobs for all cell types..."
echo "Script directory: $SCRIPT_DIR"
echo

# Loop through each cell type and submit a job
for CELLTYPE in "${CELLTYPES[@]}"; do
    echo "Submitting job for cell type: $CELLTYPE"
    
    # Create individual SLURM script for this cell type
    SLURM_SCRIPT="${SCRIPT_DIR}/run_${CELLTYPE}_analysis.sh"
    
    cat > "$SLURM_SCRIPT" << EOF
#!/bin/bash
echo "Starting cis-trans analysis for cell type: ${CELLTYPE}"
echo "Job ID: \$SLURM_JOB_ID"
echo "Node: \$SLURM_NODELIST"
echo "Start time: \$(date)"
echo

cd ${SCRIPT_DIR}
module load R/4.3.2
echo "Running R script for ${CELLTYPE}..."
Rscript 06_RUN_ANALYSIS.R ${CELLTYPE}
if [ \$? -eq 0 ]; then
    echo "Analysis completed successfully for ${CELLTYPE}"
    echo "End time: \$(date)"
else
    echo "Analysis failed for ${CELLTYPE}"
    exit 1
fi
EOF

    chmod +x "$SLURM_SCRIPT"
    
    JOB_ID=$(sbatch "$SLURM_SCRIPT" | awk '{print $4}')
    echo "  Job submitted with ID: $JOB_ID"
    echo "  Log file: cistrans_${CELLTYPE}_${JOB_ID}.log"
    echo "  Error file: cistrans_${CELLTYPE}_${JOB_ID}.err"
    echo
    sleep 2
done

echo "All jobs submitted!"
echo "Expected output directories:"
for CELLTYPE in "${CELLTYPES[@]}"; do
    echo "  cistrans_results_${CELLTYPE}/"
done
