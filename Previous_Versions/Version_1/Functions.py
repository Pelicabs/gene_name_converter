import csv
import os
import re
import pandas as pd
import time
import requests
import logging

#Figure out which columns need to be included in the API
def addColumns(df, names_col):
    def_cols = ["app_sym","app_name","prev_sym","aliases"]
    #Check for HGNC,Ensembl,and NCBI until all true or all done

    #Set all to False initially
    HGNC = False
    Ensembl = False
    NCBI = False
    
    for label in df[names_col]:
        if not HGNC and re.match(r"^HGNC:\d+$", label):
            HGNC = True
        elif not Ensembl and re.match(r"^ENSG\d+$", label):
            Ensembl = True
        elif not NCBI and re.match(r"^\d+(\.\d+)?$", label):
            NCBI = True

    if HGNC == True:
        def_cols.append("hgnc_id")
    if Ensembl == True:
        def_cols.append("pub_ensembl_id")
    if NCBI == True:
        def_cols.append("pub_eg_id")

    return def_cols

#Create the download link and download

def createDownloadURL(columns):
    #columns is a list of strings e.g. ["gd_hgnc_id","gd_app_name"]
    URL = ""
    for entry in columns:
        full_add = f"col=gd_{entry}&"
        URL = f"{URL}{full_add}"
    return URL

def makeAndFetchURL(columns):
    BASE_URL = "https://www.genenames.org/cgi-bin/download/custom?"
    REST = "status=Approved&hgnc_dbtag=off&order_by=gd_app_sym_sort&format=text&submit=submit" #Status, HGNC DB Tag, sorting, formatting, submit
    COLS = createDownloadURL(columns)
    FULL_URL = f"{BASE_URL}{COLS}{REST}"
    
    df = pd.read_csv(FULL_URL, delimiter="\t")
    temp_path = 'tempData.csv'
    df.to_csv('tempData.csv', index=False)

    return temp_path

def transform_string(input_string):
    #Check if the string starts with 'ENSG' followed by numbers
    if input_string.startswith("ENSG"):
        return re.sub(r"\.\d+$", "", input_string)
    
    #Check if the string starts with 'HGNC:' followed by numbers
    if input_string.startswith("HGNC:"):
        return input_string.split(":")[1]

    return input_string

def search_single_gene(database_path, gene_name):
    results = []  # To store the matching rows

    gene_name = transform_string(gene_name)
    
    with open(database_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        columns = next(reader) #Get the header row (column names)
    
    with open(database_path, mode='r', encoding='utf-8') as file:
        reader = csv.DictReader(file)  #Use DictReader for column access by name
        
        for row in reader:
            for column in columns:  # Dynamically iterate over all columns and rows
                if row[column] and gene_name in row[column].split(', '):
                    results.append(row)
                    break  #Stop checking other columns if a match is found
    
    result_df = pd.DataFrame(results)
    return result_df['Approved symbol'],result_df['Approved name'],result_df['Previous symbols'],result_df['Alias symbols']