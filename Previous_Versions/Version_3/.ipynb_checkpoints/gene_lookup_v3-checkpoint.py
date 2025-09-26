import csv
import os
import re
import pandas as pd
import time
import requests
import logging
from urllib.parse import quote

def setup_logging(output_name):
    log_path = f"{os.getcwd()}/Logs/gene_lookup_{output_name}.log"
    # Remove all handlers associated with the root (old) logger
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(
        filename=log_path,
        filemode='a',
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s')
    logging.info(f"Began new lookup using gene_lookup_v3.")

#Figure out which columns need to be included in the API
def addColumns(df, names_col):
    default_columns = ["app_sym", "app_name", "prev_sym", "aliases"]
    found = {"HGNC": False, "Ensembl": False, "NCBI": False}
    
    for label in df[names_col]:
        if not found["HGNC"] and re.match(r"^HGNC:\d+$", label):
            found["HGNC"] = True
        elif not found["Ensembl"] and re.match(r"^ENSG\d+$", label):
            found["Ensembl"] = True
        elif not found["NCBI"] and re.match(r"^\d+(\.\d+)?$", label):
            found["NCBI"] = True
        
        if all(found.values()): # Stop early if all found
            break

    if found["HGNC"]:
        default_columns.append("hgnc_id")
    if found["Ensembl"]:
        default_columns.append("pub_ensembl_id")
    if found["NCBI"]:
        default_columns.append("pub_eg_id")
    
    return default_columns

#Create the download link and download
def createDownloadURL(columns):
    # columns: e.g. ["hgnc_id", "app_name"]
    url_parts = [f"col=gd_{entry}" for entry in columns]
    return "&".join(url_parts) + "&" if url_parts else ""
    
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
    if input_string.startswith("ENSG"): # Ensembl gene ID: starts with 'ENSG'
        return input_string.split('.', 1)[0], "Ensembl gene ID" # Remove version suffix after the last dot (if present)
    if input_string.startswith("HGNC:"): # HGNC ID: starts with 'HGNC:'
        return input_string.partition(":")[2], "HGNC ID" # Get the part after 'HGNC:'
    if '.' in input_string: # NCBI Gene ID: contains a dot (e.g. '1234.1')
        return input_string, "NCBI Gene ID"
    return input_string, None # Default: return as-is, with None type

def _contains_gene(value, gene_name): # Helper function to return True if gene_name is in value (comma-separated or exact)
    if not value:
        return False
    if gene_name == value:
        return True
    return gene_name in {item.strip() for item in value.split(',')}

def _extract_columns(result_df): # Helper function to extract columns or empty Series if result_df is empty
    cols = ['Approved symbol', 'Approved name', 'Previous symbols', 'Alias symbols']
    return tuple(result_df.get(col, pd.Series(dtype=str)) for col in cols)

def search_single_gene(database_path, gene_name):
    gene_name, gene_type = transform_string(gene_name)

    with open(database_path, mode='r', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        columns = reader.fieldnames

        for row in reader:
            if gene_type:
                if _contains_gene(row.get(gene_type), gene_name):
                    result_df = pd.DataFrame([row])
                    return (*_extract_columns(result_df), [gene_type])
            else:
                for column in columns:
                    if _contains_gene(row.get(column), gene_name):
                        result_df = pd.DataFrame([row])
                        return (*_extract_columns(result_df), [column])
    # If no match found
    empty = pd.Series(dtype=str)
    return empty, empty, empty, empty, []
    

def _parse_json(response, label): # Helper function to parse a json
    try:
        return response.json()
    except Exception as e:
        logging.info(f"Error parsing JSON for {label}: {e}")
        return None

def _extract_record(data): # Helper function to get records
    response = data.get("response", {})
    if response.get("numFound", 0) > 0:
        docs = response.get("docs", [])
        if docs:
            return docs[0]
    return None

def getData(URL, label):
    headers = {"Accept": "application/json"}
    try:
        response = requests.get(URL, headers=headers)
    except Exception as e:
        logging.info(f"Exception occurred while querying HGNC for {label}: {e}")
        time.sleep(0.1)
        return None

    if response.status_code != 200:
        logging.info(f"API returned error status code: {response.status_code} for gene symbol '{label}'")
        time.sleep(0.1)
        return None

    data = _parse_json(response, label)
    
    if data is None:
        time.sleep(0.1)
        return None
    record = _extract_record(data)
    if record:
        return record

    # If nothing is found
    return None

def find_API(label):
    label, Type = transform_string(label)
    label = quote(label, safe='')

    headers = {"Accept": "application/json"}
    ENDPOINTS = [
    ("Approved symbol", "https://rest.genenames.org/fetch/symbol/"),
    ("Alias gene symbol", "https://rest.genenames.org/fetch/alias_symbol/"),
    ("Alias name", "https://rest.genenames.org/fetch/alias_name/"),
    ("Previous HGNC symbol", "https://rest.genenames.org/fetch/prev_symbol/")]


    data = None
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
    
    if data:
        a_sym = [data.get('symbol')]
        a_name = [data.get('name')]
        p_sym = [data.get('prev_symbol')]
        alias = [data.get('alias_symbol')]
    else:
        a_sym = [None]
        a_name = [None]
        p_sym = [None]
        alias = [None]

    time.sleep(0.1)
    return a_sym, a_name, p_sym, alias

def _insert_result_columns(df, name_col): # Helper function to insert result columns into df
    position = df.columns.get_loc(name_col)
    new_columns = pd.DataFrame({
        'Approved symbol': [None] * len(df),
        'Approved name': [None] * len(df),
        'Previous symbols': [None] * len(df),
        'Alias symbols': [None] * len(df)})
    for i, col in enumerate(new_columns.columns):
        df.insert(position + 1 + i, col, new_columns[col])
    return df

def _assign_gene_names(df, idx, aSym, aName, pSym, alias): # Helper function to safely assign gene data to correct row
    if aSym is not None and len(aSym) > 0 and aSym[0]:
        df.at[idx, 'Approved symbol'] = aSym[0]
    if aName is not None and len(aName) > 0 and aName[0]:
        df.at[idx, 'Approved name'] = aName[0]
    if pSym is not None and len(pSym) > 0 and pSym[0]:
        df.at[idx, 'Previous symbols'] = pSym[0]
    if alias is not None and len(alias) > 0 and alias[0]:
        df.at[idx, 'Alias symbols'] = alias[0]

def _handle_multiple_matches(name, match_types, aSym, aName, pSym, alias): # Helper function to log details for multiple matches
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
        else:
            match = None
        logging.info(f"    Match {i+1}: {name} is a {match_types[i]} for {aSym[i]}")

def process_multiple_names(df, idx, name_str, path): # Helper function for if there are multiple entries seperated by a comma
    gene_names = [n.strip() for n in name_str.split(',') if n.strip()]
    aSym_set, aName_set, pSym_set, alias_set = set(), set(), set(), set()
    matched = False

    for gene in gene_names:
        aSym, aName, pSym, alias, _ = search_single_gene(path, gene)
        found = len(aSym.index) > 0 if hasattr(aSym, "index") else bool(aSym and aSym[0])
        
        if not found:
            aSym, aName, pSym, alias = find_API(gene)
            found = bool(aSym and aSym[0])

        if found:
            matched = True
            if aSym is not None and len(aSym) > 0 and aSym[0]:
                aSym_set.add(aSym[0])
            if aName is not None and len(aName) > 0 and aName[0]:
                aName_set.add(aName[0])
            if pSym is not None and len(pSym) > 0 and pSym[0]:
                pSym_set.add(pSym[0])
            if alias is not None and len(alias) > 0 and alias[0]:
                alias_set.add(alias[0])

    # Assign unique values (or empty string if none)
    df.at[idx, 'Approved symbol'] = "; ".join(sorted(aSym_set)) if aSym_set else ""
    df.at[idx, 'Approved name'] = "; ".join(sorted(aName_set)) if aName_set else ""
    df.at[idx, 'Previous symbols'] = "; ".join(sorted(pSym_set)) if pSym_set else ""
    df.at[idx, 'Alias symbols'] = "; ".join(sorted(alias_set)) if alias_set else ""

    df.at[idx, 'matching_status'] = "matched" if matched else "un-matched"
        
def process_row(df, idx, name, path): # Helper function to process a row by firstly using the downloaded database, and API otherwise
    aSym, aName, pSym, alias, match_types = search_single_gene(path, name)

    if len(aSym.index) == 0:
        logging.info(f"Entry not in downloaded database, using API: {name}")
        aSym, aName, pSym, alias = find_API(name)
        if aSym is not None and len(aSym) > 0 and aSym[0]:
            df.at[idx, 'matching_status'] = "matched"
        else:
            logging.info(f"Unmatched entry found for gene: {name}")
            df.at[idx, 'matching_status'] = "un-matched"
    else:
        df.at[idx, 'matching_status'] = "matched"
        if len(aSym.index) > 1:
            _handle_multiple_matches(name, match_types, aSym, aName, pSym, alias)
    _assign_gene_names(df, idx, aSym, aName, pSym, alias)

def convert_gene_names(df_original, name_col, to_return, output_name = 'results'):    
    setup_logging(output_name)
    columns_needed = addColumns(df_original, name_col)
    path = makeAndFetchURL(columns_needed)

    df = df_original.copy()
    df['matching_status'] = "un-matched"
    df = _insert_result_columns(df, name_col)

    for idx, row in df.iterrows():
        name = row[name_col]
        if ',' in str(name):
            process_multiple_names(df, idx, name, path)
        else:
            process_row(df, idx, name, path)

    os.remove(path)
    df.rename(columns={name_col: "user_input"}, inplace=True)
    output_path = os.path.join(os.getcwd(), "Outputs", f"{output_name}_results.csv")
    df.to_csv(output_path, index=False)

    if to_return:
        return df