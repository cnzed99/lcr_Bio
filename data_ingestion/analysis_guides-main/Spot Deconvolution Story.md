**1. 定义函数 (Define functions):**

*   **`marker_gene_list()`:** 根据给定的主题 (topic) 和基因表达阈值，从最佳 LDA 模型的基因表达矩阵中提取高表达基因及其 log2 倍数变化 (log2FC)。
    ```R
    marker_gene_list<-function(topic,exp_value,Gexp){
      ## highly expressed in cell-type of interest
      highgexp <- names(which(Gexp[topic,] > exp_value))
      ## high log2(fold-change) compared to other deconvolved cell-types
      log2fc <- sort(log2(Gexp[topic,highgexp]/colMeans(Gexp[-topic,highgexp])), decreasing=TRUE)
      return(tibble(Gene=names(log2fc),log2fc=log2fc))
    }
    ```

*   **`plot_spatial()`:** 绘制特定主题在空间上的分布图，颜色深浅代表该主题的比例。
    ```R
    plot_spatial<-function(plot_data=plot_data,suffix1='_prop.jpg',dir=dir,i){
      p1<-ggplot(plot_data, aes(x, y,fill = prop)) +
        geom_point(shape=21,size=4) +
        guides(size="none")+
        labs(title=str_c("topic ",i))+
        scale_fill_viridis_c() +
        theme_void()+
        coord_equal()
      ggsave(plot = p1, paste(dir,"topic_",i,suffix1, sep=''),
             height=5, width=5, units='in', dpi=300)
    }
    ```

*   **`run_me_QC()`:** 绘制模型的质量控制 (QC) 图，包括 perplexity、稀有细胞类型数量和 alpha 值随主题数 (K) 的变化趋势。
    ```R
    run_me_QC<-function(ldas,starting_k=2,ending_k=22,dir){
      if(!dir.exists(dir)){dir.create(dir)}
      alpha<-map(.x = (starting_k-1):(ending_k-1),~ldas$model[[.x]]@alpha)%>%unlist
      plot_df<-tibble(K=starting_k:ending_k,alpha=alpha,perplexities=ldas$perplexities,rare=ldas$numRare)
      p1<-ggplot(data = plot_df) +
        geom_line(mapping = aes(x = K,y = perplexities), color="red3",size=2) +
        geom_point(mapping = aes(x = K,y = perplexities),shape=21, color="black", fill=ifelse(alpha > 1, "white", "red3"), size=6)+
        theme_linedraw(base_size = 16,base_rect_size =2,base_line_size = 2)+ylab("perplexity")+
        ylim(min(plot_df$perplexities)-10,10+max(plot_df$perplexities))
      p2<-ggplot(data = plot_df) +
        geom_point(mapping = aes(x = K,y = rare),shape=21, color="black", fill="blue", size=4)+
        geom_line(mapping = aes(x = K,y = rare), color="blue",size=2)+
        theme_linedraw(base_size = 16,base_rect_size =2,base_line_size = 2)+ylab("cell−types with mean proportion < 5%")
      p3<-ggplot(data = plot_df) +
        geom_line(mapping = aes(x = K,y = alpha), color="darkgreen",size=2) +ylim(c(0,1))+
        geom_point(mapping = aes(x = K,y = alpha),shape=21, color="black", fill=ifelse(alpha > 1, "white", "darkgreen"), size=6)+
        theme_linedraw(base_size = 16,base_rect_size =2,base_line_size = 2)+ylab("alpha")
      ggsave(plot = p1+p2+p3,filename = str_c(dir,"merged_QC_plot.jpg"),height=5, width=12, units='in', dpi=300)
    }
    ```

*   **`run_me_results()`:** 对最佳 LDA 模型进行结果分析和可视化。它会生成每个主题的空间分布图，并输出每个主题的高表达基因列表及其 log2FC 值。
    ```R
    run_me_results<-function(opt,
                             dir,ldas ){
      optLDA <- optimalModel(models = ldas, opt = opt)
      results <- getBetaTheta(optLDA,
                              perc.filt = 0.05,
                              betaScale = 1000)
      deconProp <- results$theta
      deconGexp <- results$beta
      if(!dir.exists(dir)){dir.create(dir)}
      for(i in 1:dim(deconProp)[2]){
        plot_data<-merge(pos,deconProp[,i],by = 0)
        names(plot_data)<-c("barcode","x","y","prop")
        plot_spatial(plot_data=plot_data,suffix1='_prop.jpg',dir=dir,i=i)
      }
      marker_gene_output<-map(.x = 1:dim(deconGexp)[1],
                              ~marker_gene_list(topic = .x,exp_value = 2,Gexp = deconGexp))
      names(marker_gene_output)<-str_c("topic_genes_exp2.",1:dim(deconGexp)[1],".csv")
      map2(.x = names(marker_gene_output),
           .y = marker_gene_output,
           ~write_csv(x = .y,file = paste(dir,.x)))
    }
    ```

**2. 加载和预处理数据 (Load and preprocess data):**

*   `counts <- Read10X_h5(filename = "VisiumFFPE_Mouse_Brain_Transgenic_Age_17p9_Rep_1.h5")`: 使用 `Read10X_h5()` 函数读取 HDF5 格式的 10x Visium 数据，得到基因表达计数矩阵。

*   `spatial_barcodes <- read_csv("spatial_cord_subset_17p9_rep1.csv")`: 读取空间坐标子集文件。

*   `counts_subset <- counts[, colnames(counts) %in% spatial_barcodes$barcode]`: 根据空间坐标子集文件，筛选出特定区域的 spots 的基因表达数据。

*   `pos <- as.data.frame(spatial_barcodes)`: 将空间坐标子集文件转换为数据框，并提取 x 和 y 坐标。

*   `counts_subset_clean <- cleanCounts(...)`: 使用 `cleanCounts()` 函数对计数矩阵进行质量控制，过滤掉低质量的 spots 和低表达的基因。

*   `odGenes <- getOverdispersedGenes(...)`: 使用 `getOverdispersedGenes()` 函数识别过表达的基因。

*   `astro <- read_csv(file = "astro_markers.csv")`: 读取星形胶质细胞标记基因列表。

*   `gene_astro <- c(genes, astro_overlap) %>% unique()`: 将过表达基因和星形胶质细胞标记基因合并，作为模型的输入基因。

*   `corpus <- preprocess(...)`: 使用 `preprocess()` 函数对数据进行预处理，构建 LDA 模型的语料库。

**3. 拟合 LDA 模型 (Fit LDA model for a range of topics (K values)):**

*   `ldas <- fitLDA(corpus$corpus, Ks = seq(2, 22, by = 1), ...)`: 使用 `fitLDA()` 函数对不同主题数 (K) 的 LDA 模型进行拟合。这里 K 的范围是 2 到 22，步长为 1。

*   `run_me_QC(ldas, starting_k=2, ending_k=22, dir="output_18/")`: 绘制模型的 QC 图，并保存到 "output_18/" 文件夹中。

**4. 结果分析和可视化 (Exporting spatial plots and topics for the optimal model):**

*   `run_me_results(opt=18, dir = "output_18/", ldas=ldas)`: 选择主题数为 18 的最佳模型 (根据 QC 图选择最佳 K 值，这里假设为 18)，生成每个主题的空间分布图和高表达基因列表，并保存到 "output_18/" 文件夹中。

**输出文件详细解释:**

`output_18/` 文件夹中包含以下文件：

*   **merged_QC_plot.jpg:** 模型的 QC 图，包括 perplexity、稀有细胞类型数量和 alpha 值随主题数 (K) 的变化趋势。该图用于评估不同 K 值模型的拟合优度，并选择最佳的 K 值。

*   **topic__i__prop.jpg (i = 1, 2, ..., 18):** 每个主题 (topic) 的空间分布图，其中 i 代表主题编号。图片中每个点代表一个 spot，颜色深浅表示该主题在该 spot 中的比例。这些图可以直观地展示每个主题在组织中的空间分布模式。

*   **topic_genes_exp2._i_.csv (i = 1, 2, ..., 18):** 每个主题的高表达基因列表及其 log2FC 值，其中 i 代表主题编号。这些文件包含以下两列：

    *   **Gene:** 高表达基因的名称。
    *   **log2fc:** 该基因在该主题中相对于其他主题的 log2 倍数变化 (log2FC)。正值表示该基因在该主题中高表达，负值表示低表达。这些文件可以帮助我们了解每个主题的生物学意义，并识别潜在的细胞类型标记基因。

**后续步骤:**

1. **确定最佳主题数 (K):** 仔细检查 `merged_QC_plot.jpg`，综合考虑 perplexity、稀有细胞类型数量和 alpha 值，选择最佳的 K 值。一般来说，选择 perplexity 较低、稀有细胞类型数量较少且 alpha 值小于 1 的 K 值。

2. **解释主题的生物学意义:**

    *   查看每个主题的空间分布图 (`topic_[i]_prop.jpg`)，观察其在组织中的分布模式。
    *   分析每个主题的高表达基因列表 (`topic_genes_exp2.[i].csv`)，结合已知的细胞类型标记基因，推断每个主题可能代表的细胞类型。例如，如果某个主题富集星形胶质细胞标记基因，则该主题可能代表星形胶质细胞。可以使用 `enrichR` 等工具进行基因集富集分析。
    *   可以将基因和主题关联性进行可视化，例如可以利用 `view_cell_marker_correlations` 函数。
    *   使用 `cell_topic_heatmap` 函数可视化细胞类型和主题之间的关系。
    *   使用 `view_cell_subtype_marker_correlations` 函数可视化细胞亚型标记物之间的关系。

3. **与其他数据整合分析:**

    *   可以将反卷积结果与其他的空间转录组学数据或单细胞 RNA 测序数据进行整合分析，以更深入地了解组织中细胞类型的组成和空间分布。

4. **下游分析:** 基于上述分析，可以探究不同细胞类型在空间上的相互作用、信号通路变化、以及与疾病发生发展的关系等。