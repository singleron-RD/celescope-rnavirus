import numpy as np
import pandas as pd
import scanpy as sc

from celescope.__init__ import HELP_DICT
from celescope.rna.mkref import Mkref_rna
from celescope.tools import utils
from celescope.tools.plotly_plot import Tsne_plot
from celescope.tools.step import Step, s_common

# markers adjust p_value
PVAL_CUTOFF = 0.05
# scanpy mitochonrial variable name
MITO_VAR = "mito"
NORMALIZED_LAYER = "normalised"
RESOLUTION = 1.2
N_PCS = 25
MITO_GENE_PERCENT_LIST = [5, 10, 15, 20, 50]
# output marker top n in html
MARKER_TOP_N = 100


def read_tsne(tsne_file):
    df = pd.read_csv(tsne_file, sep="\t", index_col=0)
    df.index.names = ["barcode"]
    return df


def format_df_marker(df_marker, feature="gene"):
    avg_logfc_col = "avg_log2FC"  # seurat 4
    df_marker = df_marker.loc[
        :, ["cluster", feature, avg_logfc_col, "pct.1", "pct.2", "p_val_adj"]
    ]
    df_marker["cluster"] = df_marker["cluster"].apply(lambda x: f"cluster {x}")
    df_marker = df_marker[df_marker["p_val_adj"] < PVAL_CUTOFF]
    return df_marker


def get_opts_analysis(parser, sub_program):
    parser.add_argument("--genomeDir", help=HELP_DICT["genomeDir"], required=True)
    parser.add_argument(
        "--kbDir",
        help="kb reference directory path. Must contain index.idx and t2g.txt.",
    )
    if sub_program:
        parser.add_argument(
            "--gene_matrix_file",
            required=True,
        )
        parser.add_argument(
            "--transcript_matrix_file",
            required=True,
        )
        parser = s_common(parser)


class Analysis(Step):
    def __init__(self, args, display_title=None):
        super().__init__(args, display_title=display_title)

        self.adata = dict()

        # data
        self.adata["gene"] = sc.read_10x_mtx(
            args.gene_matrix_file,
            var_names="gene_symbols",
        )
        self.adata["transcript"] = sc.read_10x_mtx(
            args.transcript_matrix_file,
            var_names="gene_ids",
        )
        self.adata["gene"].layers["raw"] = self.adata["gene"].X.copy()
        self.mt_gene_list = Mkref_rna.get_config(args.genomeDir)["files"][
            "mt_gene_list"
        ]

        # out
        self.marker_file_dict = {}
        self.marker_raw_file_dict = {}
        self.h5ad_dict = {}
        for feature in ("gene", "transcript"):
            self.marker_file_dict[feature] = f"{self.outdir}/{feature}_markers.tsv"
            self.marker_raw_file_dict[feature] = (
                f"{self.outdir}/{feature}_markers_raw.tsv"
            )
            self.h5ad_dict[feature] = f"{self.outdir}/{feature}.h5ad"
            self.outs += [
                self.marker_file_dict[feature],
                self.marker_raw_file_dict[feature],
                self.h5ad_dict[feature],
            ]
        self.df_tsne_file = f"{self.outdir}/tsne_coord.tsv"
        self.outs += [self.df_tsne_file]

    @utils.add_log
    def calculate_qc_metrics(self, feature="gene"):
        if self.mt_gene_list:
            mito_genes, _ = utils.read_one_col(self.mt_gene_list)
            self.adata[feature].var[MITO_VAR] = self.adata[feature].var_names.map(
                lambda x: True if x in mito_genes else False
            )
            # if not astype(bool), it will be type object and raise an error
            # https://github.com/theislab/anndata/issues/504
            self.adata[feature].var[MITO_VAR] = (
                self.adata[feature].var[MITO_VAR].astype(bool)
            )
        else:
            self.adata[feature].var[MITO_VAR] = (
                self.adata[feature].var_names.str.upper().str.startswith("MT-")
            )

        sc.pp.calculate_qc_metrics(
            self.adata[feature],
            qc_vars=[MITO_VAR],
            percent_top=None,
            use_raw=False,
            log1p=False,
            inplace=True,
        )

    @utils.add_log
    def write_mito_stats(self):
        mt_pct_var = f"pct_counts_{MITO_VAR}"
        total_cell_number = self.adata["gene"].n_obs

        for mito_gene_percent in MITO_GENE_PERCENT_LIST:
            cell_number = sum(self.adata["gene"].obs[mt_pct_var] > mito_gene_percent)
            fraction = round(cell_number / total_cell_number * 100, 2)
            self.add_metric(
                name=f"Fraction of cells have mito gene percent>{mito_gene_percent}%",
                value=f"{fraction}%",
            )

    @utils.add_log
    def normalize(self, feature="gene"):
        """
        sc.pp.normalize_per_cell() and sc.pp.log1p()
        """

        sc.pp.normalize_total(
            self.adata[feature],
            target_sum=1e4,
            inplace=True,
        )
        sc.pp.log1p(
            self.adata[feature],
        )
        self.adata[feature].layers[NORMALIZED_LAYER] = self.adata[feature].X

    @utils.add_log
    def hvg(self):
        """
        Wrapper function for sc.highly_variable_genes()
        """
        sc.pp.highly_variable_genes(
            self.adata["gene"],
            layer=None,
            n_top_genes=None,
            min_disp=0.5,
            max_disp=np.inf,
            min_mean=0.0125,
            max_mean=3,
            span=0.3,
            n_bins=20,
            flavor="seurat",
            subset=False,
            inplace=True,
            batch_key=None,
            check_values=True,
        )

    @utils.add_log
    def scale(self):
        """
        Wrapper function for sc.pp.scale
        """
        sc.pp.scale(
            self.adata["gene"],
            zero_center=True,
            max_value=10,
            copy=False,
            layer=None,
            obsm=None,
        )

    @utils.add_log
    def pca(self):
        """
        Wrapper function for sc.pp.pca
        """
        sc.pp.pca(
            self.adata["gene"],
            n_comps=50,
            zero_center=True,
            svd_solver="auto",
            random_state=0,
            return_info=False,
            use_highly_variable=True,
            dtype="float32",
            copy=False,
            chunked=False,
            chunk_size=None,
        )

    @utils.add_log
    def neighbors(
        self,
    ):
        """
        Wrapper function for sc.pp.neighbors(), for supporting multiple n_neighbors
        """
        sc.pp.neighbors(
            self.adata["gene"],
            n_neighbors=15,
            n_pcs=N_PCS,
            use_rep=None,
            knn=True,
            random_state=0,
            method="umap",
            metric="euclidean",
            key_added=None,
            copy=False,
        )

    @utils.add_log
    def tsne(self):
        """
        Wrapper function for sc.tl.tsne, for supporting named slot of tsne embeddings
        """
        sc.tl.tsne(
            self.adata["gene"],
            n_pcs=N_PCS,
            copy=False,
        )

    @utils.add_log
    def umap(
        self,
    ):
        """
        Wrapper function for sc.tl.umap, for supporting named slot of umap embeddings
        """
        sc.tl.umap(
            self.adata["gene"],
            min_dist=0.5,
            spread=1.0,
            n_components=2,
            maxiter=None,
            alpha=1.0,
            gamma=1.0,
            negative_sample_rate=5,
            init_pos="spectral",
            random_state=0,
            a=None,
            b=None,
            copy=False,
            method="umap",
            neighbors_key=None,
        )

    @utils.add_log
    def leiden(self):
        """
        Wrapper function for sc.tl.leiden
        """
        sc.tl.leiden(
            self.adata["gene"],
            resolution=RESOLUTION,
            restrict_to=None,
            random_state=0,
            key_added="cluster",
            adjacency=None,
            directed=True,
            use_weights=True,
            n_iterations=-1,
            partition_type=None,
            neighbors_key=None,
            obsp=None,
            copy=False,
        )

    @utils.add_log
    def transfer_cluster(self):
        self.adata["transcript"].obs["cluster"] = (
            self.adata["gene"]
            .obs["cluster"]
            .reindex(self.adata["transcript"].obs_names)
        )

    @utils.add_log
    def find_marker_genes(self, feature="gene"):
        """
        Wrapper function for sc.tl.rank_genes_groups
        """
        sc.tl.rank_genes_groups(
            self.adata[feature],
            "cluster",
            reference="rest",
            pts=True,
            method="wilcoxon",
            use_raw=False,
            layer=NORMALIZED_LAYER,
        )

    @utils.add_log
    def get_transcript_gene_dict(self):
        df = pd.read_csv(f"{self.args.kbDir}/t2g.txt", sep="\t", header=None)
        transcript_gene_dict = dict(zip(df[0], df[2]))
        return transcript_gene_dict

    @utils.add_log
    def write_markers(self, feature="gene"):
        """
        write only p_val_adj < PVAL_CUTOFF to avoid too many markers
        """
        df_markers = sc.get.rank_genes_groups_df(
            self.adata[feature], group=None, pval_cutoff=PVAL_CUTOFF
        )
        df_markers = df_markers[df_markers["logfoldchanges"].notna()]
        markers_name_dict = {
            "group": "cluster",
            "names": feature,
            "logfoldchanges": "avg_log2FC",
            "pvals": "p_val",
            "pvals_adj": "p_val_adj",
            "pct_nz_group": "pct.1",
            "pct_nz_reference": "pct.2",
        }
        df_markers = df_markers.rename(markers_name_dict, axis="columns")
        df_markers["cluster"] = df_markers["cluster"].map(
            lambda x: f"cluster {int(x)+1}"
        )
        df_markers = df_markers.loc[
            :, ["cluster", feature, "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]
        ]
        df_markers = df_markers.loc[df_markers["p_val_adj"] < PVAL_CUTOFF,]
        df_markers.to_csv(self.marker_raw_file_dict[feature], index=None, sep="\t")

        df_markers_filter = (
            df_markers.loc[df_markers["avg_log2FC"] > 0]
            .sort_values("p_val_adj")
            .groupby("cluster")
            .head(100)
        )
        df_markers_filter = df_markers_filter.round(
            {
                "avg_log2FC": 3,
                "pct.1": 3,
                "pct.2": 3,
            }
        )
        if feature == "transcript":
            transcript_gene_dict = self.get_transcript_gene_dict()
            df_markers_filter["gene"] = df_markers_filter["transcript"].map(
                lambda x: transcript_gene_dict[x]
            )
            df_markers_filter = df_markers_filter.loc[
                :,
                [
                    "cluster",
                    "gene",
                    "transcript",
                    "avg_log2FC",
                    "pct.1",
                    "pct.2",
                    "p_val_adj",
                ],
            ]
        df_markers_filter.to_csv(self.marker_file_dict[feature], index=None, sep="\t")
        return df_markers_filter

    @utils.add_log
    def write_tsne(self):
        df_tsne = self.adata["gene"].obsm.to_df()[["X_tsne1", "X_tsne2"]]
        df_tsne["cluster"] = self.adata["gene"].obs.cluster
        df_tsne["Gene_Counts"] = self.adata["gene"].obs.n_genes_by_counts
        df_tsne["Transcript_Counts"] = self.adata["transcript"].obs.n_genes_by_counts
        df_tsne["Transcript_Counts"] = (
            df_tsne["Transcript_Counts"].fillna(0).astype(int)
        )
        df_tsne.index.names = ["barcode"]
        tsne_name_dict = {"X_tsne1": "tSNE_1", "X_tsne2": "tSNE_2"}
        df_tsne = df_tsne.rename(tsne_name_dict, axis="columns")
        df_tsne["cluster"] = df_tsne["cluster"].map(lambda x: int(x) + 1)
        df_tsne.to_csv(self.df_tsne_file, sep="\t")
        return df_tsne

    @utils.add_log
    def write_h5ad(self, feature):
        self.adata[feature].write(self.h5ad_dict[feature])

    @utils.add_log
    def add_report_data(self, df_tsne, gene_markers, transcript_markers):
        tsne_cluster = Tsne_plot(df_tsne, "cluster").get_plotly_div()
        self.add_data(tsne_cluster=tsne_cluster)

        tsne_gene = Tsne_plot(df_tsne, "Gene_Counts", discrete=False).get_plotly_div()
        self.add_data(tsne_gene=tsne_gene)

        table_dict_gene = self.get_table_dict(
            title="Marker Genes by Cluster",
            table_id="marker_genes",
            df_table=gene_markers,
        )
        self.add_data(table_dict_gene=table_dict_gene)

        table_dict_transcript = self.get_table_dict(
            title="Marker Transcripts by Cluster",
            table_id="marker_transcripts",
            df_table=transcript_markers,
        )
        self.add_data(table_dict_transcript=table_dict_transcript)

    def add_marker_help(self):
        self.add_help_content(
            name="Marker by Cluster",
            content="differential expression analysis based on the non-parameteric Wilcoxon rank sum test",
        )
        self.add_help_content(
            name="avg_log2FC",
            content="log fold-change of the average expression between the cluster and the rest of the sample",
        )
        self.add_help_content(
            name="pct.1",
            content="The percentage of cells where the gene is detected in the cluster",
        )
        self.add_help_content(
            name="pct.2",
            content="The percentage of cells where the gene is detected in the rest of the sample",
        )
        self.add_help_content(
            name="p_val_adj",
            content="Adjusted p-value, based on bonferroni correction using all genes in the dataset",
        )

    @utils.add_log
    def run(self):
        self.calculate_qc_metrics(feature="gene")
        self.calculate_qc_metrics(feature="transcript")
        self.write_mito_stats()
        self.normalize(feature="gene")
        self.hvg()
        self.scale()
        self.pca()
        self.neighbors()
        self.tsne()
        self.umap()
        self.leiden()
        self.find_marker_genes(feature="gene")
        gene_markers = self.write_markers(feature="gene")

        self.normalize(feature="transcript")
        self.transfer_cluster()
        self.find_marker_genes(feature="transcript")
        transcript_markers = self.write_markers(feature="transcript")

        df_tsne = self.write_tsne()
        self.write_h5ad("gene")
        self.write_h5ad("transcript")
        self.add_report_data(df_tsne, gene_markers, transcript_markers)
        self.add_marker_help()


class Report_runner(Step):
    def __init__(self, args, display_title=None, feature="Gene"):
        super().__init__(args, display_title=display_title)
        self.feature = feature

    def add_marker_help(self):
        self.add_help_content(
            name=f"Marker {self.feature} by Cluster",
            content="differential expression analysis based on the non-parameteric Wilcoxon rank sum test",
        )
        self.add_help_content(
            name="avg_log2FC",
            content="log fold-change of the average expression between the cluster and the rest of the sample",
        )
        self.add_help_content(
            name="pct.1",
            content="The percentage of cells where the gene is detected in the cluster",
        )
        self.add_help_content(
            name="pct.2",
            content="The percentage of cells where the gene is detected in the rest of the sample",
        )
        self.add_help_content(
            name="p_val_adj",
            content="Adjusted p-value, based on bonferroni correction using all genes in the dataset",
        )

    def run(self):
        pass


@utils.add_log
def analysis(args):
    with Analysis(args, display_title="Analysis") as runner:
        runner.run()
