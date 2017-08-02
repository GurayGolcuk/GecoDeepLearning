"""
This code is prepared assumning the kirc data from GECO/GMQL is used

"""

import pickle
import pandas as pd
import numpy as np
from collections import Counter

def importData(data_url = ".\gmql_data\kirc_data_df.p"):

    """
    Get and reshape the input dataframe
    :param data_url: string url of the input file
    :return: Dataframe
    """

    df = pickle.load(open(data_url, "rb"))
    return df

def changeStructure(df):
    """
    Change unnnecesery multiindex structure to standart dataframe structure
    :param df:
    :return:
    """

    df = df.normalized_count
    df.index = df.index.get_level_values(-1)

    return df

def getLabels(df):
    """
    #Colecting all labels from inout data
    :param regiondata: input dataframe
    :return: labels of the dataframe
    """

    lbls = df.index.get_level_values('biospecimen_sample__sample_type_id')

    #Since we have only 1 label as "05" and it has same meaning as "01" we convert "05" to "01"
    lbls = ["01" if x == "05" else x for x in lbls]
    for i in range(len(lbls)):
        lbls[i] = int(lbls[i])

    return lbls


def getFeatures(onogenes_url =".\onogenes\kirc.txt"):

    """
    :param onogenes_url: url of onogenes file
    :return: List of necessery features
    """
    df = pd.read_csv(onogenes_url, delimiter="\t")
    df = pd.DataFrame(df)
    features = list(df.Gene_Symbol)

    # delete the "RHEBP1" feature since our gene data don"t have this property
    features.remove("RHEBP1")
    return features

def fillEmptyCells(df):

    """
    Fill empty values as the mean value of the colomn so it won't change the output
    :param df: the input dataframe
    :return: the resulting dataframe
    """
    df = df.fillna(df.mean())

    return df

def to_zero_mean(df):
    """
    Standardizes the data by shifting the mean value to zero
    :param df: the input dataframe
    :return: the resulting dataframe
    """
    df_norm = (df - df.mean())
    return df_norm


def to_unit_variance(df):
    """
    Makes the variance of each gene equal to one
    :param df: the dataframe
    :return: the resulting dataframe
    """
    df_norm = df / df.std()
    return df_norm

def pickleNormalSave(df, output_url = ".\gmql_data\kirc_data_with_label_deneme.p"):

    """
    Save df to picke file
    :param output_url: Output path
    :param df: final dataframe
    """

    try:

        f = open(output_url, 'wb')
        pickle.dump(df, f, pickle.HIGHEST_PROTOCOL)
        f.close()
    except Exception as e:
        print('Unable to save data to', output_url, ':', e)
        raise

if __name__ == "__main__":

    data = importData()
    labels = getLabels(data)
    data = changeStructure(data)
    features = getFeatures()

    data = to_unit_variance(to_zero_mean(fillEmptyCells(data)))

    pickleNormalSave(data)


