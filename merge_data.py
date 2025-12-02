import os
import argparse
import pandas as pd


def extract_msa_name_from_path(path_str: str) -> str:
    """
    Từ full path:
      D:\\test\\TreeBASEMirror-main\\output\\10000_0.phy\\msa.fasta
    -> '10000_0.phy' (tên folder chứa msa.fasta)
    """
    path_norm = os.path.normpath(str(path_str))
    parent_dir = os.path.dirname(path_norm)
    msa_name = os.path.basename(parent_dir)
    return msa_name


def load_feature_csv(csv_path: str) -> pd.DataFrame:
    """
    Đọc 1 file CSV feature (có cột 'filename'),
    thêm cột 'msa_name' = tên folder chứa msa.fasta.
    """
    df = pd.read_csv(csv_path)
    if "filename" not in df.columns:
        raise ValueError(f"{csv_path} không có cột 'filename'.")

    df["msa_name"] = df["filename"].apply(extract_msa_name_from_path)
    return df


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge 2 file CSV feature (msa_features_*.csv) vào training_data.parquet dựa trên tên folder msa."
    )
    parser.add_argument(
        "--train-parquet",
        required=True,
        help="Đường dẫn tới file parquet gốc (training_data.parquet).",
    )
    parser.add_argument(
        "--csv",
        required=True,
        action="append",
        help="Đường dẫn tới 1 file CSV feature (có thể truyền nhiều lần). "
             "Ví dụ: --csv msa_features_fixed.csv --csv msa_features_output.csv",
    )
    parser.add_argument(
        "--out-parquet",
        required=True,
        help="Đường dẫn file parquet sau khi merge (file mới).",
    )
    parser.add_argument(
        "--out-csv",
        default=None,
        help="(Optional) Nếu muốn xuất thêm file CSV sau khi merge.",
    )

    args = parser.parse_args(argv)

    print(f"Loading parquet gốc: {args.train_parquet}")
    df_train = pd.read_parquet(args.train_parquet)

    if "verbose_name" not in df_train.columns:
        raise SystemExit("File parquet không có cột 'verbose_name' để làm khóa join.")

    feature_dfs = []
    for csv_path in args.csv:
        print(f"Loading feature CSV: {csv_path}")
        df_feat = load_feature_csv(csv_path)
        feature_dfs.append(df_feat)

    df_features = pd.concat(feature_dfs, ignore_index=True)
    df_features = df_features.drop_duplicates(subset=["msa_name"], keep="last")

    if "filename" in df_features.columns:
        df_features = df_features.drop(columns=["filename"])

    print("Merging features vào parquet ...")
    df_merged = df_train.merge(
        df_features,
        left_on="verbose_name",
        right_on="msa_name",
        how="left",
    )

    if "msa_name" in df_merged.columns:
        df_merged = df_merged.drop(columns=["msa_name"])

    print(f"Saving merged parquet to: {args.out_parquet}")
    df_merged.to_parquet(args.out_parquet, index=False)

    if args.out_csv is not None:
        print(f"Saving merged csv to: {args.out_csv}")
        df_merged.to_csv(args.out_csv, index=False)

    print("Done.")


if __name__ == "__main__":
    raise SystemExit(main())
