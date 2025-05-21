options(timeout = 300)  # Set timeout to 5 minutes
# Vectors
patient_ids <- c("GSM6592048","GSM6592049","GSM6592050", "GSM6592051", "GSM6592052", "GSM6592053",
                 "GSM6592054", "GSM6592055", "GSM6592056", "GSM6592057", "GSM6592058",
                 "GSM6592059", "GSM6592060", "GSM6592061", "GSM6592062")

sample_ids <- c("M1_CL-like3", "M2_CL-like1", "M3_Unstable1", "M4_Unstable2", "M5_Unstable3",
                "M6_Unstable4", "M7_CL-like4", "M8_CL-like2", "M9_CL-like3", "M10_CL-like4",
                "M11_CL-like3", "M13_MBC_spindle1", "M14_MBC_spindle2", "M15_MBC_chondroid", "M16_MBC_NST")
suffixes <- c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10", "M11", "M13", "M14", "M15","M16")

new_patient_id <- c("P1_M1", "P2_M2", "P3_M3", "P4_M4", "P5_M5", "P6_M6", "P7_M7", "P8_M8", "P9_M9", "P10_M10",
                    "P11_M11", "P12_M13", "P13_M14", "P14_M15","P15_M16")

patient_condition <- c("normal_P1",  # M1
                       "cancer_P2",  # M2
                       "cancer_P3",  # M3
                       "cancer_P4",  # M4
                       "cancer_P5",  # M5
                       "cancer_P6",  # M6
                       "normal_P7",  # M7
                       "cancer_P8",  # M8
                       "normal_P9",  # M9
                       "normal_P10", # M10
                       "cancer_P11", # M11
                       "cancer_P12", # M13
                       "cancer_P13", # M14
                       "cancer_P14",
                       "cancer_P15") # M15
# Make sure sample_ids length matches folder count
stopifnot(length(patient_ids) == 15)

# Create folders
sapply(patient_ids, function(x) {
  dir.create(x, showWarnings = TRUE, recursive = TRUE)
})
# Vectors

file_types <- list(
  filtered_matrix = "filtered_feature_bc_matrix.h5",
  annotations     = "pathologist_annotations.csv.gz",
  scalefactors    = "scalefactors_json.json.gz",
  hires_image     = "tissue_hires_image.png.gz",
  positions_list  = "tissue_positions_list.csv.gz"
)

# Loop through each patient and download files into their folder
for (i in seq_along(patient_ids)) {
  pid <- patient_ids[i]
  sid <- suffixes[i]
  
  # Create the directory if it doesn't exist
  if (!dir.exists(pid)) {
    dir.create(pid)
  }
  
  # Loop through each file type
  for (ftype in names(file_types)) {
    fname <- file_types[[ftype]]
    encoded_fname <- URLencode(paste0(pid, "_", sid, "_", fname))
    
    url <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?acc=", pid, "&format=file&file=", encoded_fname)
    destfile <- file.path(pid, paste0(pid, "_", sid, "_", fname))  # Save inside the folder
    
    # Download the file
    download.file(url = url, destfile = destfile, mode = "wb")
  }
}
# Load required libraries
library(R.utils)
library(Seurat)
library(dplyr)
# Working directory
wd <- getwd()

for (i in seq_along(patient_ids)) {
  pid    <- patient_ids[i]
  suf    <- suffixes[i]
  sampid <- sample_ids[i]
  
  combo <- paste0(pid, "_", suf)
  patient_dir <- file.path(wd, pid)
  spatial_dir <- file.path(patient_dir, "spatial")
  
  # Create spatial folder if it doesn't exist
  if (!dir.exists(spatial_dir)) dir.create(spatial_dir)
  
  # Step 1: Move and decompress spatial files
  spatial_files <- c("tissue_positions_list.csv.gz",
                     "tissue_hires_image.png.gz",
                     "scalefactors_json.json.gz")
  
  for (file in spatial_files) {
    src <- file.path(patient_dir, paste0(combo, "_", file))
    dst <- file.path(spatial_dir, paste0(combo, "_", file))
    if (file.exists(src)) {
      file.rename(src, dst)
    } else {
      message("Missing file: ", src)
    }
  }
  
  # Step 2: Decompress files
  gz_files <- list.files(spatial_dir, pattern="\\.gz$", full.names=TRUE)
  for (gzf in gz_files) {
    gunzip(gzf, destname=sub("\\.gz$", "", gzf), remove=TRUE)
  }
  
  # Step 3: Rename files to match standard names
  rename_map <- c("tissue_positions_list.csv",
                  "tissue_hires_image.png",
                  "scalefactors_json.json")
  
  for (file in rename_map) {
    oldname <- file.path(spatial_dir, paste0(combo, "_", file))
    newname <- file.path(spatial_dir, file)
    if (file.exists(oldname)) {
      file.rename(oldname, newname)
    } else {
      message("File not found for renaming: ", oldname)
    }
  }
  
  # Step 4: Load matrix and image
  h5f <- file.path(patient_dir, paste0(combo, "_filtered_feature_bc_matrix.h5"))
  if (!file.exists(h5f)) {
    message("Missing H5: ", h5f, " — skipping.")
    next
  }
  matrix <- Read10X_h5(h5f)
  
  imgf <- file.path(spatial_dir, "tissue_hires_image.png")
  if (!file.exists(imgf)) {
    message("Missing image: ", imgf, " — skipping.")
    next
  }
  img <- Read10X_Image(spatial_dir, image.name="tissue_hires_image.png")
  
  # Step 5: Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = matrix, assay = "Spatial")
  seurat_obj@images$spatial <- img
  
  # Step 6: Add annotation if available
  annf <- file.path(patient_dir, paste0(combo, "_pathologist_annotations.csv.gz"))
  if (file.exists(annf)) {
    gunzip(annf, remove=FALSE)
    ann_txt <- sub("\\.gz$", "", annf)
    if (file.exists(ann_txt)) {
      annotations <- read.csv(ann_txt, stringsAsFactors = FALSE)
      colnames(annotations) <- c("barcode", "annotations")
      
      seurat_obj$barcode <- rownames(seurat_obj@meta.data)
      metadata <- left_join(seurat_obj@meta.data, annotations, by = "barcode")
      rownames(metadata) <- metadata$barcode
      metadata$barcode <- NULL
      seurat_obj@meta.data <- metadata
    }
  }
  
  # Step 7: Add metadata columns
  seurat_obj[["sampleid"]]   <- sampid
  seurat_obj[["patient_id"]] <- pid
  seurat_obj[["suffix"]]     <- suf
  seurat_obj[["new_patient_id"]] <- new_patient_id[i]
  seurat_obj[["condition"]] <- patient_condition[i]
  
  # Step 8: QC & processing
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern="^MT-")
  seurat_obj <- subset(seurat_obj, subset = nCount_Spatial > 200 & nFeature_Spatial < 7500 & percent.mt < 10)
  seurat_obj <- SCTransform(seurat_obj, assay="Spatial", verbose=FALSE)
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method="vst", nfeatures=2000)
  seurat_obj <- RunPCA(seurat_obj, verbose=FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, reduction="pca", dims=1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution=0.5)
  seurat_obj <- RunUMAP(seurat_obj, reduction="pca", dims=1:30)
  
  # Step 9: Save
  saveRDS(seurat_obj, file=file.path(patient_dir, paste0(combo, "_seurat.rds")))
}
#####################################integration######################
GSM6592048_M1_seurat <- readRDS("~/integration with more metadata/GSM6592048/GSM6592048_M1_seurat.rds")
GSM6592049_M2_seurat <- readRDS("~/integration with more metadata/GSM6592049/GSM6592049_M2_seurat.rds")
GSM6592050_M3_seurat <- readRDS("~/integration with more metadata/GSM6592050/GSM6592050_M3_seurat.rds")
GSM6592051_M4_seurat <- readRDS("~/integration with more metadata/GSM6592051/GSM6592051_M4_seurat.rds")
GSM6592052_M5_seurat <- readRDS("~/integration with more metadata/GSM6592052/GSM6592052_M5_seurat.rds")
GSM6592053_M6_seurat <- readRDS("~/integration with more metadata/GSM6592053/GSM6592053_M6_seurat.rds")
GSM6592054_M7_seurat <- readRDS("~/integration with more metadata/GSM6592054/GSM6592054_M7_seurat.rds")
GSM6592055_M8_seurat <- readRDS("~/integration with more metadata/GSM6592055/GSM6592055_M8_seurat.rds")
GSM6592056_M9_seurat <- readRDS("~/integration with more metadata/GSM6592056/GSM6592056_M9_seurat.rds")
GSM6592057_M10_seurat <- readRDS("~/integration with more metadata/GSM6592057/GSM6592057_M10_seurat.rds")
GSM6592058_M11_seurat <- readRDS("~/integration with more metadata/GSM6592058/GSM6592058_M11_seurat.rds")
GSM6592059_M13_seurat <- readRDS("~/integration with more metadata/GSM6592059/GSM6592059_M13_seurat.rds")
GSM6592060_M14_seurat <- readRDS("~/integration with more metadata/GSM6592060/GSM6592060_M14_seurat.rds")
GSM6592061_M15_seurat <- readRDS("~/integration with more metadata/GSM6592061/GSM6592061_M15_seurat.rds")
GSM6592062_M16_seurat <- readRDS("~/integration with more metadata/GSM6592062/GSM6592062_M16_seurat.rds")
library(Seurat)
library(dplyr)
library(harmony)
library(patchwork)
GSM6592061_M15_seurat@meta.data
# 1. Read in all Seurat objects into a list
seurat_list <- list(
  GSM6592048_M1 = GSM6592048_M1_seurat,
  GSM6592049_M2 = GSM6592049_M2_seurat,
  GSM6592050_M3 = GSM6592050_M3_seurat,
  GSM6592051_M4 = GSM6592051_M4_seurat,
  GSM6592052_M5 = GSM6592052_M5_seurat,
  GSM6592053_M6 = GSM6592053_M6_seurat,
  GSM6592054_M7 = GSM6592054_M7_seurat,
  GSM6592055_M8 = GSM6592055_M8_seurat,
  GSM6592056_M9 = GSM6592056_M9_seurat,
  GSM6592057_M10 = GSM6592057_M10_seurat,
  GSM6592058_M11 = GSM6592058_M11_seurat,
  GSM6592059_M13 = GSM6592059_M13_seurat,
  GSM6592060_M14 = GSM6592060_M14_seurat,
  GSM6592061_M15 = GSM6592061_M15_seurat,
  GSM6592062_M16 = GSM6592062_M16_seurat)


  
  # Merge all Seurat objects while ensuring unique cell names
  merged <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = sample_ids,  # exclude first ID because it's already in x
    project = "HarmonyIntegration"
  )
library(stringr)
  
merged@meta.data <- merged@meta.data %>%
  mutate(annotations = annotations %>%
           str_replace_all("Artefacts|Artifacts", "Artifacts") %>%
           str_replace_all("Tumor Cells|Tumour Cells|Tumor cells|Tumour cells|\\?", "Tumor cells") %>%
           str_replace_all("Tumor Stroma|Tumour Stroma", "Tumor Stroma"))
unique(merged@meta.data$annotations)
#merged2 <- merge(
  #x = seurat_list[[1]],
  #y = seurat_list[-1],
  #add.cell.ids = patient_condition,  # exclude first ID because it's already in x
  #project = "HarmonyIntegration"
#)

# 3. Standard preprocessing on merged object
# Basic normalization
merged <- NormalizeData(merged, verbose = FALSE)

# Find highly variable genes
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)

# Scale data
merged <- ScaleData(merged, verbose = FALSE)

# Run PCA
merged <- RunPCA(merged, verbose = FALSE)

# 4. PCA visualization before Harmony
p1 <- DimPlot(merged, reduction = "pca", group.by = "sampleid") + ggtitle("PCA before Harmony")
p11 <- DimPlot(merged, reduction = "pca", group.by = "annotations") + ggtitle("PCA before Harmony")
p111 <- DimPlot(merged, reduction = "pca", group.by = "condition") + ggtitle("PCA before Harmony")
p1111 <- DimPlot(merged, reduction = "pca", group.by = "new_patient_id") + ggtitle("PCA before Harmony")


# 5. Run Harmony integration on patient_id
combined1 <- RunHarmony(merged, group.by.vars = "sampleid", verbose = TRUE)

combined11 <- RunPCA(combined1, verbose = FALSE)



# Extract PCA embeddings
pca_coords <- Embeddings(combined11, "pca")  # or "harmony" if post-Harmony
# Get unique sample IDs
samples <- unique(combined11$sampleid)
par(mfrow = c(1, 2))

# Assign distinct colors
cols <- setNames(rainbow(length(samples)), samples)

# Plot using base R
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols[combined11$sampleid],
  main = "PCA",
  pch = 20,
  xlab = "PC1",
  ylab = "PC2"
)

# Get unique sample IDs
samples <- unique(merged$annotations)

# Assign distinct colors
cols <- setNames(rainbow(length(samples)), samples)

# Plot using base R
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols[merged$annotations],
  main = "PCA",
  pch = 20,
  xlab = "PC1",
  ylab = "PC2"
)

# Plot using base R
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols[combined11$annotations],
  main = "PCA",
  pch = 20,
  xlab = "PC1",
  ylab = "PC2"
)

# Get unique sample IDs
samples <- unique(merged$annotations)

# Assign distinct colors
cols <- setNames(rainbow(length(samples)), samples)

# Plot using base R
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols[merged$annotations],
  main = "PCA",
  pch = 20,
  xlab = "PC1",
  ylab = "PC2"
)

par(mfrow = c(1, 2))

samples <- unique(merged$sampleid)
par(mfrow = c(1, 2))

# Assign distinct colors
cols <- setNames(rainbow(length(samples)), samples)

# Plot using base R
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols[merged$sampleid],
  main = "PCA",
  pch = 20,
  xlab = "PC1",
  ylab = "PC2"
)
####################
library(Seurat)
library(harmony)
library(irlba)

# Step 1: Scale data
merged <- ScaleData(merged)

# Step 2: Run PCA
combined11 <- RunPCA(merged, npcs = 50, verbose = FALSE)
pca_coords <- Embeddings(merged, "pca")

# Step 3: Harmony integration
harmony_coords <- Embeddings(combined11, "harmony")

# Step 4: Color settings
sample_ids <- unique(combined11$sampleid)
annotation_groups <- unique(combined11$annotations)

cols_sampleid <- setNames(rainbow(length(sample_ids)), sample_ids)
cols_annotations <- setNames(rainbow(length(annotation_groups)), annotation_groups)

# Step 5: Side-by-side PCA colored by sampleid
par(mfrow = c(1, 2))
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols_sampleid[combined11$sampleid],
  main = "PCA (Before Harmony)\nColored by SampleID",
  pch = 20, xlab = "PC1", ylab = "PC2"
)

plot(
  harmony_coords[, 1], harmony_coords[, 2],
  col = cols_sampleid[combined11$sampleid],
  main = "PCA (After Harmony)\nColored by SampleID",
  pch = 20, xlab = "PC1", ylab = "PC2"
)

# Step 6: Side-by-side PCA colored by annotations
par(mfrow = c(1, 2))
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols_annotations[combined11$annotations],
  main = "PCA (Before Harmony)\nColored by Annotations",
  pch = 20, xlab = "PC1", ylab = "PC2"
)

plot(
  harmony_coords[, 1], harmony_coords[, 2],
  col = cols_annotations[combined11$annotations],
  main = "PCA (After Harmony)\nColored by Annotations",
  pch = 20, xlab = "PC1", ylab = "PC2"
)
#################the true that we will include
# Split sample IDs into 3 groups for 3 lines
sample_ids <- names(cols_sampleid)
n <- length(sample_ids)
split_ids <- split(sample_ids, ceiling(seq_along(sample_ids) / ceiling(n / 3)))
dev.off()
# Setup layout: 2 plots on top row, legend box full width below
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(4, 1))  # 4 units tall plots, 1 unit tall legend

# Plot margins for PCA plots
par(mar = c(5, 5, 4, 2))  # leave space for axis labels and titles

# Plot 1
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols_sampleid[combined11$sampleid],
  main = "PCA (Before Harmony)\nColored by SampleID",
  pch = 20, xlab = "PC1", ylab = "PC2", bty = "o",xaxt = "s",    # show x-axis
  yaxt = "s"     # show y-axis
)

# Plot 2
plot(
  harmony_coords[, 1], harmony_coords[, 2],
  col = cols_sampleid[combined11$sampleid],
  main = "PCA (After Harmony)\nColored by SampleID",
  pch = 20, xlab = "PC1", ylab = "PC2", bty = "o",xaxt = "s",    # show x-axis
  yaxt = "s"     # show y-axis
)

# Now plot legend in bottom row
par(mar = c(0, 0, 0, 0), xpd = TRUE, xaxt = 'n', yaxt = 'n', bty = 'n')
plot.new()

# Legend plotting area dimensions [0,1] in both x and y

# Number of lines
lines_count <- length(split_ids)

# Vertical positions for lines (spread nicely)
y_pos <- seq(0.75, 0.25, length.out = lines_count)

# Horizontal spacing - spread across full width, adjust multiplier for spacing
max_per_line <- max(sapply(split_ids, length))
x_spacing <- 1 / (max_per_line + 1)

# Split sample IDs into 3 groups for 3 lines
sample_ids <- names(cols_sampleid)
n <- length(sample_ids)
split_ids <- split(sample_ids, ceiling(seq_along(sample_ids) / ceiling(n / 3)))

# Setup layout: 2 plots on top row, legend box full width below 
#that is code with the plots with legends below 
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(4, 1))  # 4 units tall plots, 1 unit tall legend

# Plot margins for PCA plots
par(mar = c(5, 5, 4, 2))  # leave space for axis labels and titles

# Plot 1
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols_sampleid[combined11$sampleid],
  main = "PCA (Before Harmony)\nColored by SampleID",
  pch = 20, xlab = "PC1", ylab = "PC2", bty = "o",xaxt = "s",    # show x-axis
  yaxt = "s"     # show y-axis
)

# Plot 2
plot(
  harmony_coords[, 1], harmony_coords[, 2],
  col = cols_sampleid[combined11$sampleid],
  main = "PCA (After Harmony)\nColored by SampleID",
  pch = 20, xlab = "PC1", ylab = "PC2", bty = "o",xaxt = "s",    # show x-axis
  yaxt = "s"     # show y-axis
)

# Now plot legend in bottom row
par(mar = c(0, 0, 0, 0), xpd = TRUE, xaxt = 'n', yaxt = 'n', bty = 'n')
plot.new()

# Legend plotting area dimensions [0,1] in both x and y

# Number of lines
lines_count <- length(split_ids)

# Vertical positions for lines (spread nicely)
y_pos <- seq(0.75, 0.25, length.out = lines_count)

# Horizontal spacing - spread across full width, adjust multiplier for spacing
max_per_line <- max(sapply(split_ids, length))
x_spacing <- 1 / (max_per_line + 1)

# Plot colored points and labels line by line
for (i in seq_along(split_ids)) {
  ids <- split_ids[[i]]
  n_ids <- length(ids)
  
  # Horizontal positions spaced evenly
  x_positions <- seq(from = x_spacing, by = x_spacing, length.out = n_ids)
  y <- y_pos[i]
  
  # Plot colored points
  points(x = x_positions, y = rep(y, n_ids),
         col = cols_sampleid[ids],
         pch = 20, cex = 2)
  
  # Add sample ID labels slightly above points
  text(x = x_positions, y = rep(y + 0.08, n_ids),
       labels = ids, cex = 1.4, adj = 0.5)
}

#########################
# Split annotation groups (or sample IDs) into 3 roughly equal groups
annotation_groups <- unique(combined11$annotations)  # or replace with your vector of IDs
n <- length(annotation_groups)
split_ids <- split(annotation_groups, ceiling(seq_along(annotation_groups) / ceiling(n / 3)))
dev.off()  # close all open graphics devices (run multiple times if needed)

# Layout: two plots on top row, one legend area full width below
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(4, 1))  # Adjust heights as needed

# Set margins for PCA plots
par(mar = c(5, 5, 4, 2))  # bottom, left, top, right margins

# Plot 1: PCA before Harmony
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols_annotations[combined11$annotations],
  main = "PCA (Before Harmony)\nColored by Annotations",
  pch = 20, xlab = "PC1", ylab = "PC2",
  bty = "o", xaxt = "s",   # show x-axis
  yaxt = "s"    # show y-axis
)

# Plot 2: PCA after Harmony
plot(
  harmony_coords[, 1], harmony_coords[, 2],
  col = cols_annotations[combined11$annotations],
  main = "PCA (After Harmony)\nColored by Annotations",
  pch = 20, xlab = "PC1", ylab = "PC2",
  bty = "o", xaxt = "s",   # show x-axis
  yaxt = "s"    # show y-axis
)

# Legend area: no margins, no axes, no box
par(mar = c(0, 0, 0, 0), xpd = TRUE, xaxt = 'n', yaxt = 'n', bty = 'n')
plot.new()

# Legend plotting coordinates [0,1] in x and y
lines_count <- length(split_ids)

# Vertical positions for 3 lines evenly spaced between 0.75 and 0.25
y_pos <- seq(0.75, 0.1, length.out = lines_count)

# Horizontal spacing - max number of items in any line + 1 buffer space
max_per_line <- max(sapply(split_ids, length))
x_spacing <- 1 / (max_per_line + 1)

# Loop through each line to draw colored points and labels
for (i in seq_along(split_ids)) {
  ids <- split_ids[[i]]
  n_ids <- length(ids)
  
  # Base horizontal spacing
  base_spacing <- 1 / (max_per_line + 1)
  
  # Increase spacing by multiplying (e.g., 1.5 times larger spacing)
  spacing_multiplier <- 3
  adjusted_spacing <- base_spacing * spacing_multiplier
  
  # Calculate x positions with bigger spacing
  # Start a bit more left so points fit inside [0,1]
  start_x <- 0.1
  x_positions <- start_x + (0:(n_ids - 1)) * adjusted_spacing
  
  # If points go beyond 1, rescale to fit inside [0,1]
  if (max(x_positions) > 1) {
    x_positions <- seq(from = 0.05, to = 0.95, length.out = n_ids)
  }
  
  y <- y_pos[i]
  
  # Plot colored points
  points(x = x_positions, y = rep(y, n_ids),
         col = cols_annotations[ids],
         pch = 20, cex = 2)
  
  # Add labels slightly above points, with center alignment (adj=0.5)
  text(x = x_positions, y = rep(y + 0.08, n_ids),
       labels = ids, cex = 1.4, adj = 0.5)
}




# Assign colors for new_patient_id
unique_patients <- unique(combined11$new_patient_id)
cols_patient <- setNames(rainbow(length(unique_patients)), unique_patients)

# Assign colors for condition
unique_conditions <- unique(combined11$condition)
cols_condition <- setNames(rainbow(length(unique_conditions)), unique_conditions)


### ---- PCA Before & After Harmony Colored by New Patient ID ----

# Layout: 2 plots (top) + 1 legend (bottom)
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(4, 1))

# Margins for PCA plots
par(mar = c(5, 5, 4, 2))

# Plot 1: PCA Before Harmony colored by new_patient_id
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols_patient[combined11$new_patient_id],
  main = "PCA (Before Harmony)\nColored by New Patient ID",
  pch = 20, xlab = "PC1", ylab = "PC2",
  bty = "o", xaxt = "s",   # show x-axis
  yaxt = "s"    # show y-axis
)

# Plot 2: PCA After Harmony colored by new_patient_id
plot(
  harmony_coords[, 1], harmony_coords[, 2],
  col = cols_patient[combined11$new_patient_id],
  main = "PCA (After Harmony)\nColored by New Patient ID",
  pch = 20, xlab = "PC1", ylab = "PC2",
  bty = "o", xaxt = "s",   # show x-axis
  yaxt = "s"    # show y-axis
)

# Legend for new_patient_id
par(mar = c(0, 0, 0, 0), xpd = TRUE, xaxt = 'n', yaxt = 'n', bty = 'n')
plot.new()

annotation_groups <- unique_patients
n <- length(annotation_groups)
split_ids <- split(annotation_groups, ceiling(seq_along(annotation_groups) / ceiling(n / 3)))

lines_count <- length(split_ids)
y_pos <- seq(0.75, 0.1, length.out = lines_count)
max_per_line <- max(sapply(split_ids, length))

for (i in seq_along(split_ids)) {
  ids <- split_ids[[i]]
  n_ids <- length(ids)
  base_spacing <- 1 / (max_per_line + 1)
  spacing_multiplier <- 3
  adjusted_spacing <- base_spacing * spacing_multiplier
  start_x <- 0.1
  x_positions <- start_x + (0:(n_ids - 1)) * adjusted_spacing
  
  if (max(x_positions) > 1) {
    x_positions <- seq(from = 0.05, to = 0.95, length.out = n_ids)
  }
  
  y <- y_pos[i]
  points(x = x_positions, y = rep(y, n_ids),
         col = cols_patient[ids],
         pch = 20, cex = 2)
  text(x = x_positions, y = rep(y + 0.08, n_ids),
       labels = ids, cex = 1.4, adj = 0.5)
}

dev.off()


### ---- PCA Before & After Harmony Colored by Condition ----

# Layout: 2 plots (top) + 1 legend (bottom)
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(4, 1))

# Margins for PCA plots
par(mar = c(5, 5, 4, 2))

# Plot 3: PCA Before Harmony colored by condition
plot(
  pca_coords[, 1], pca_coords[, 2],
  col = cols_condition[combined11$condition],
  main = "PCA (Before Harmony)\nColored by Condition",
  pch = 20, xlab = "PC1", ylab = "PC2",
  bty = "o"
)

# Plot 4: PCA After Harmony colored by condition
plot(
  harmony_coords[, 1], harmony_coords[, 2],
  col = cols_condition[combined11$condition],
  main = "PCA (After Harmony)\nColored by Condition",
  pch = 20, xlab = "PC1", ylab = "PC2",
  bty = "o"
)

# Legend for condition
par(mar = c(0, 0, 0, 0), xpd = TRUE, xaxt = 'n', yaxt = 'n', bty = 'n')
plot.new()

annotation_groups <- unique_conditions
n <- length(annotation_groups)
split_ids <- split(annotation_groups, ceiling(seq_along(annotation_groups) / ceiling(n / 3)))

lines_count <- length(split_ids)
y_pos <- seq(0.75, 0.1, length.out = lines_count)
max_per_line <- max(sapply(split_ids, length))

for (i in seq_along(split_ids)) {
  ids <- split_ids[[i]]
  n_ids <- length(ids)
  base_spacing <- 1 / (max_per_line + 1)
  spacing_multiplier <- 3
  adjusted_spacing <- base_spacing * spacing_multiplier
  start_x <- 0.1
  x_positions <- start_x + (0:(n_ids - 1)) * adjusted_spacing
  
  if (max(x_positions) > 1) {
    x_positions <- seq(from = 0.05, to = 0.95, length.out = n_ids)
  }
  
  y <- y_pos[i]
  points(x = x_positions, y = rep(y, n_ids),
         col = cols_condition[ids],
         pch = 20, cex = 2)
  text(x = x_positions, y = rep(y + 0.08, n_ids),
       labels = ids, cex = 1.4, adj = 0.5)
}

dev.off()

