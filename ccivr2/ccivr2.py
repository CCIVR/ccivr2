# -*- coding: utf-8 -*-

import pandas as pd
from ccivr2.preproc import get_prior_inf
import ccivr2.exception as e

def find_cisnats(df0,df1,filter):
    
    cols = df0.columns

    dict_pairs ={}
    counter = 0
 
    df1_chr_map = dict(list(df1.groupby('Chr')))

    for i in range(df0.shape[0]):

        # Chromosomal selection
        chr_i = df0.at[i,'Chr']
        min0 = df0.at[i,'Start']
        max0 = df0.at[i,'End']

        if chr_i not in df1_chr_map.keys():
            continue
        
        extrby_chr = df1_chr_map[chr_i]

        # Locational selection
        min1 = extrby_chr.Start
        max1 = extrby_chr.End

        extrby_st_ed = extrby_chr[eval(filter)]
        
        # Forming a pair
        if len(extrby_st_ed.index) != 0:

            gene_target = pd.DataFrame(df0.iloc[i]).T.reset_index(drop=True)
            gene_cisnats = extrby_st_ed.reset_index(drop=True).rename(columns=lambda x: '_' + x)

            pairs = pd.concat([gene_target, gene_cisnats], axis =1)
            pairs[cols] = pairs[cols].fillna(method = 'ffill')

            for j in range(pairs.shape[0]):
                dict_pairs[counter] = pairs.iloc[j]
                
                counter += 1

    result = pd.DataFrame.from_dict(dict_pairs, orient='index')
    return result


def extract_one_type_cisnats(method, genes_plus, genes_minus, tss_range, tes_range):
    
    # Criteria for locational extraction
    FILTER_PLUS_MINUS= {
        "EB":"(min1 >= min0) & (max1 <= max0)",
        "FO":"(min1 <= min0) & (max1 >= max0)",
        "HH":"(min1 < min0) & (max1 >= min0) & (max1 <= max0)",
        "TT":"(min1 >= min0) & (min1 <= max0) & (max1 > max0)",
        "TSS comparison":"(max1 >= min0 + {min}) & (max1 <= min0 + {max})".format(min=tss_range[0], max=tss_range[1]),
        "TES comparison":"(min1 >= max0 + {min}) & (min1 <= max0 + {max})".format(min=tes_range[0], max=tes_range[1])
        }

    FILTER_MINUS_PLUS= {
        "EB":"(min1 >= min0) & (max1 <= max0)",
        "FO":"(min1 <= min0) & (max1 >= max0)",
        "HH":"(min1 >= min0) & (min1 <= max0) & (max1 > max0)",
        "TT":"(min1 < min0) & (max1 >= min0) & (max1 <= max0)",
        "TSS comparison":"(min1 >= max0 - {max}) & (min1 <= max0 - {min})".format(min=tss_range[0], max=tss_range[1]),
        "TES comparison":"(max1 >= min0 - {max}) & (max1 <= min0 - {min})".format(min=tes_range[0], max=tes_range[1])
        }

    # funtion to add other structural information to each pair
    def adinf(row):
        if row['Strand'] == '+':
            A,B,a,b = row['Start'],row['End'],row['_End'],row['_Start']
        else:
            A,B,a,b = -row['End'],-row['Start'],-row['_Start'],-row['_End']

        # Add 'Overlap' and 'Type'
        if b > B+5000:
            row['Overlap'] = "Non_overlap"
            row['Type'] = "fTT"
        
        elif B < b <= B+5000:
            row['Overlap'] = "Non_overlap"
            row['Type'] = "nTT"

        elif (a > B) & (A <= b <= B):
            row['Overlap'] = "Overlap"
            row['Type'] = "TT"

        elif (a > B) & (b < A):
            row['Overlap'] = "Overlap"
            row['Type'] = "FO"

        elif (A <= a <= B) & (A <= b <= B):
            row['Overlap'] = "Overlap"
            row['Type'] = "EB"

        elif (A <= a <= B) & (b < A):
            row['Overlap'] = "Overlap"
            row['Type'] = "HH"

        elif A-5000 <= a < A:
            row['Overlap'] = "Non_overlap"
            row['Type'] = "nHH"

        elif a < A-5000:
            row['Overlap'] = "Non_overlap"
            row['Type'] = "fHH"

        # Add 'AS_TSS_position' in mode2
        if method == 'TSS comparison':
            row['AS_TSS_position'] = a-A

        # Add 'AS_TES_position' in mode3
        if method == 'TES comparison':
            row['AS_TES_position'] = b-B

        return row

    print(method)

    # Extraction
    # Plus-to-minus extracting
    print('Plus-to-minus extracting')
    cnpairs1 = find_cisnats(
        df0 = genes_plus,
        df1 = genes_minus,
        filter = FILTER_PLUS_MINUS[method]
        )

    # Minus-to-plus extracting
    print('Minus-to-plus extracting')
    cnpairs2 = find_cisnats(
        df0 = genes_minus,
        df1 = genes_plus,
        filter = FILTER_MINUS_PLUS[method]
        )

    cnpairs_combined = pd.concat([cnpairs1,cnpairs2], axis=0)
    cnpairs_combined = cnpairs_combined.apply(adinf, axis=1)

    return cnpairs_combined


def write_summary(df, total, setting):

    result = pd.DataFrame(columns = ['count','rate'])

    result.at['mode','count'] = setting.mode

    if setting.mode == '2':
        result.at['AS_TSS position [min]','count'] = setting.tss[0]
        result.at['AS_TSS position [max]','count'] = setting.tss[1]
    if setting.mode == '3':
        result.at['AS_TES position [min]','count'] = setting.tes[0]
        result.at['AS_TES position [max]','count'] = setting.tes[1]

    result.at['total genes','count'] = total

    # All type
    count_all = df['id'].nunique()
    rate_all = '{:.2%}'.format(count_all/total)

    result.at['extracted genes','count'] = count_all
    result.at['extracted genes','rate'] = rate_all

    def zero_or_count(what,keys:list,map:dict):
        if what not in keys:
            count = 0
            rate = '0.00%'
        else:
            what_df = map[what]

            count = what_df['id'].nunique()
            rate = '{:.2%}'.format(count/total)

        return count,rate

    # Overlap/Non_overlap count (in case mode 2 or 3)
    if setting.mode in ['2','3']:
        overlap_map = dict(list(df.groupby('Overlap')))
        overlap_keys = list(overlap_map.keys())
        result.at['Overlap','count'],result.at['Overlap','rate'] = zero_or_count('Overlap',overlap_keys,overlap_map)
        result.at['Non_overlap','count'],result.at['Non_overlap','rate'] = zero_or_count('Non_overlap',overlap_keys,overlap_map)

    # Each type count
    type_map = dict(list(df.groupby('Type')))
    type_keys = list(type_map.keys())

    if setting.mode == '1':
        type_list = ['EB','FO','HH','TT']
    elif setting.mode in ['2','3']:
        type_list = ['fTT','nTT','TT','FO','EB','HH','nHH','fHH']
    
    for t in type_list:
        result.at[t,'count'],result.at[t,'rate'] = zero_or_count(t,type_keys,type_map)

    # print summary
    def print_result_row(idx):
        print(idx.ljust(17) + str(result.at[idx,'count']).ljust(6) + '[' + str(result.at[idx,'rate']).center(7) + ']')
    
    print('\n< Result >')
    print('total genes'.ljust(17) + str(result.at['total genes','count']).ljust(6))
    print_result_row('extracted genes')
    print('')
    
    if setting.mode in ['2','3']:
        print_result_row('Overlap')
        print_result_row('Non_overlap')
        print('')

    for t in type_list:
        print_result_row(t)
    print('')

    return result
      

def main():

    inf = get_prior_inf()

    try:
        print('Reading ' + inf.paths.input)
        df = pd.read_csv(inf.paths.input, header=0)

        e.item_check(df.columns)

        genes_total = df['id'].nunique()

        # Dividing df into two groups, minus (-) or plus (+) strand
        df_plus=df.query('Strand == "+"').reset_index(drop=True)
        df_minus=df.query('Strand == "-"').reset_index(drop=True)

        # Extracting
        result_list = list()

        targets = {
            "1":["EB","FO","HH","TT"],
            "2":["TSS comparison"],
            "3":["TES comparison"]
            }

        for i in targets[inf.mode]:
            result = extract_one_type_cisnats(
                method = i,
                genes_plus = df_plus,
                genes_minus = df_minus,
                tss_range = inf.tss,
                tes_range = inf.tes
                )
            result_list.append(result)

        table = pd.concat(result_list, axis=0)

        if len(table.index) == 0:
            print("No genes extracted")
            from shutil import rmtree
            rmtree(inf.paths.output[0])

        else:
            summary = write_summary(table, genes_total, setting=inf)

            print('Writing the table to ' + inf.paths.output[1])
            table.to_csv(inf.paths.output[1], index=False)
            print('Writing the summary to ' + inf.paths.output[2])
            summary.to_csv(inf.paths.output[2], index=True, header=False)

    except:
        import traceback
        traceback.print_exc()
        from shutil import rmtree
        rmtree(inf.paths.output[0])


if __name__ == '__main__':
    main()
