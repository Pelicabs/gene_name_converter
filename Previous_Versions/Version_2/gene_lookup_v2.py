import csv
import os
import re
import pandas as pd
import time
import requests
import logging
from urllib.parse import quote

#Setup logging file
log_path = f"{os.getcwd()}/Logs/gene_lookup.log"
logging.basicConfig(
    filename=log_path,
    filemode='a',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

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
    Type = None
    if input_string.startswith("ENSG"):
        Type = "Ensembl gene ID"
        return re.sub(r"\.\d+$", "", input_string), Type
    
    #Check if the string starts with 'HGNC:' followed by numbers
    if input_string.startswith("HGNC:"):
        Type = "HGNC ID"
        return input_string.split(":")[1], Type

    tryNCBI = str(input_string)
    if '.' in tryNCBI:
        Type = "NCBI Gene ID"
        return input_string, Type
    else:
        return input_string, Type

def search_single_gene(database_path, gene_name):
    results = []  # To store the matching rows
    match_types = []

    gene_name, Type = transform_string(gene_name)
    
    with open(database_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        columns = next(reader) #Get the header row (column names)
    
    with open(database_path, mode='r', encoding='utf-8') as file:
        reader = csv.DictReader(file)  #Use DictReader for column access by name
        
        for row in reader:
            if Type:
                if row[Type] and gene_name in row[Type].split(', '):
                    results.append(row)
                    break
            else:
                for column in columns:  # Dynamically iterate over all columns and rows
                    if row[column] and gene_name in row[column].split(', '):
                        results.append(row)
                        match_types.append(column)
                        break  #Stop checking other columns if a match is found
                    
    result_df = pd.DataFrame(results)
    return result_df['Approved symbol'],result_df['Approved name'],result_df['Previous symbols'],result_df['Alias symbols'],match_types

def getData(URL, label):
    try:
        response = requests.get(URL, headers = headers)
        if response.status_code == 200:
            try:
                data = response.json()
            except Exception as json_err:
                logging.info(f"Error parsing JSON for {label}")

            if data.get("response", {}).get("numFound", 0) > 0:
                record = data["response"]["docs"][0]
                return record
        else:
            logging.info(f"API returned error status code: {response.status_code} for gene symbol '{label}'")
            
    except Exception as e:
        logging.info(f"Exception occurred while querying HGNC for {label}: {e}")

    time.sleep(0.1) #10 requests per second

def find_API(label):
    label, Type = transform_string(label)
    label = quote(label, safe='')

    headers = {"Accept": "application/json"}
    ENDPOINTS = [
    ("Approved symbol", "https://rest.genenames.org/fetch/symbol/"),
    ("Alias gene symbol", "https://rest.genenames.org/fetch/alias_symbol/"),
    ("Alias name", "https://rest.genenames.org/fetch/alias_name/"),
    ("Previous HGNC symbol", "https://rest.genenames.org/fetch/prev_symbol/")]


    if Type == "Ensembl gene ID":
        URL = f"https://rest.genenames.org/fetch/ensembl_gene_id/{label}"
        data = getData(URL, label)
    elif Type == "NCBI Gene ID":
        URL = f"https://rest.genenames.org/fetch/entrez_id/{label}"
        data = getData(URL, label)
    elif Type == "HGNC ID":
        URL = f"https://rest.genenames.org/fetch/hgnc_id/{label}"
        data = getData(URL, label)
    else:
        for endpoint_type, base_url in ENDPOINTS:
            URL = f"{base_url}{label}"
            data = getData(URL, label)
            if data:
                break

    a_sym = data['symbol']
    a_name = data['name']
    p_sym = data['prev_symbol']
    alias = data['alias_symbol']
    sleep(0.1)
    return a_sym, a_name, p_sym, alias

def convert_gene_names(df_original, name_col, to_return):    
    logging.info(f"Began new lookup using gene_lookup_v2.")
    #1. Find what info is needed and then download data, returning path
    columns_needed = addColumns(df_original, name_col)
    path = makeAndFetchURL(columns_needed)

    #2. Prep new df (add columns)
    df = df_original.copy()
    position = df.columns.get_loc(name_col)

    new_columns = pd.DataFrame({
    'Approved symbol': [None] * len(df),
    'Approved name': [None] * len(df),
    'Previous symbols': [None] * len(df),
    'Alias symbols': [None] * len(df)})

    for i, col in enumerate(new_columns.columns):
        df.insert(position + 1 + i, col, new_columns[col])

    #3. Iterate through each entry and add names
    for idx, row in df.iterrows():
        name = row[name_col]
        aSym, aName, pSym, alias, match_types = search_single_gene(path, name)

        if len(aSym.index) == 0:
            logging.info(f"Entry not in downloaded database, using API",name)
            #Use API to fetch in two rounds (first for approved, then using approved to find rest of data)
            aSym, aName, pSym, alias = find_API(name)
        else:
            # If multiple entries are returned, log them.
            if len(aSym.index) > 1:
                logging.info(f"  Multiple entries found for {name}")
                for i in range(len(aSym)):
                    if match_types[i] == 'Approved symbol':
                        match = aSym[i]
                    elif match_types[i] == 'Approved name':
                        match = aName[i]
                    elif match_types[i] == 'Previous symbols':
                        match = pSym[i]
                    elif match_types[i] == 'Alias symbols':
                        match = alias[i]
                    logging.info(f"    Match {i+1}: {name} is a {match_types[i]} for {aSym[i]}")
        df.at[idx, 'Approved symbol'] = aSym[0]
        df.at[idx, 'Approved name'] = aName[0]
        df.at[idx, 'Previous symbols'] = pSym[0]
        df.at[idx, 'Alias symbols'] = alias[0]

    os.remove(path)
    df = df.drop(columns=[name_col])
    output_path = f"{os.getcwd()}/Outputs/result.csv"
    df.to_csv(output_path, index = False)

    if to_return:
        return df