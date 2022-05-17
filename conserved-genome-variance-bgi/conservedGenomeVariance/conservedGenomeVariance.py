import os
import sys
import pandas as pd
import numpy as np

def conglomerateVariance(SNV_INDEL_DIR):
    # Combine variance from all chromosomes into single dataframe

    for root, dirs, files in os.walk(SNV_INDEL_DIR):
        # Create dataframe to hold whole organism variance
        organismVarianceDf = pd.DataFrame(data=None, columns=['Label', 'Chr ID', 'Source', 'Variation type', 'Site', 'Refbase',
        'Seqbase', 'RefReadsNumber', 'SeqReadsNumber', 'Approximate read depth', 'SeqReadsRate', 'Quality value', 'GenoType',
        'Gene ID', 'Codon changes', 'Amino acid changes', 'Amino acid Substitution', 'Gene name', 'Alias', 'ORF classification',
        'Essential status', 'Functional description' ])
        # Add individual variance
        for file in files:
            if file.endswith('.xls'):
                chrVarianceFilePath = os.path.join(root, file)
                # File path is .xls but content is tab seperated
                chrVarianceDf = pd.read_csv(chrVarianceFilePath, sep=r'\t')
                # Join chromosome onto whole organism df
                organismVarianceDf = pd.concat([organismVarianceDf, chrVarianceDf])

        # Find name of strain and insert into dataframe
        strainDir = os.path.join(root, files[0])
        strainName = os.path.basename(os.path.dirname(os.path.dirname(strainDir)))
        organismVarianceDf.insert(1, 'Strain', strainName )
        # Reset index of combined table
        organismVarianceDf = organismVarianceDf.reset_index(drop=True)

    return organismVarianceDf

def removeVarianceOverlap(STRAIN_SNV_INDEL_DIR, REF_SNV_INDEL_DIR):
    # Removes overlap in variance from reference to strain

    refOrganismVarianceDf = conglomerateVariance(REF_SNV_INDEL_DIR)
    strainOrganismVarianceDf = conglomerateVariance(STRAIN_SNV_INDEL_DIR)

    # Go through ref organism and remove SNVs from strain which occur in ref
    cleanStrainOrganismVarianceDf = strainOrganismVarianceDf
    preCleanVariantCount = (strainOrganismVarianceDf.size / 23)

    for index, row in refOrganismVarianceDf.iterrows():
            cleanStrainOrganismVarianceDf = cleanStrainOrganismVarianceDf.drop(
            cleanStrainOrganismVarianceDf[(cleanStrainOrganismVarianceDf['Chr ID']==row['Chr ID']) &
            (cleanStrainOrganismVarianceDf['Site']==row['Site']) &
            (cleanStrainOrganismVarianceDf['Seqbase']==row['Seqbase'])].index)

    postCleanVariantCount = (cleanStrainOrganismVarianceDf.size / 23)

    return cleanStrainOrganismVarianceDf, preCleanVariantCount, postCleanVariantCount

def filterVariance(organismVarianceDf):
    # Remove variance with filter label (poor quality)

    preFilterVariantCount = (organismVarianceDf.size / 23)
    filteredOrganismVarianceDf = organismVarianceDf[organismVarianceDf['Label']!='filter']
    filteredOrganismVarianceDf = filteredOrganismVarianceDf.reset_index(drop=True)
    postFilterVariantCount = (filteredOrganismVarianceDf.size / 23)

    strainName = filteredOrganismVarianceDf.Strain[0]
    filteredOrganismVarianceDf.to_csv(strainName + '_filtered.csv') # Leave uncommented to generate indivudal .csv sheats for strain variants


    return filteredOrganismVarianceDf, preFilterVariantCount, postFilterVariantCount

def findVarianceConservation(varianceDfList):
    # Generate df with combined variance and showing relative and absolute counts of individual variants

    # Collate all dfs into super df
    collectedVarianceDf = pd.DataFrame(data=None, columns=['Label', 'Strain', 'Chr ID', 'Source', 'Variation type', 'Site', 'Refbase',
    'Seqbase', 'RefReadsNumber', 'SeqReadsNumber', 'Approximate read depth', 'SeqReadsRate', 'Quality value', 'GenoType',
    'Gene ID', 'Codon changes', 'Amino acid changes', 'Amino acid Substitution', 'Gene name', 'Alias', 'ORF classification',
    'Essential status', 'Functional description' ])
    for df in varianceDfList:
        collectedVarianceDf = pd.concat([collectedVarianceDf, df])
    collectedVarianceDf = collectedVarianceDf.reset_index(drop=True)

    # returns a dataframe with all variants that occur >1 times
    conservedVarianceDf = collectedVarianceDf[collectedVarianceDf.duplicated(subset=['Chr ID', 'Site', 'Seqbase'], keep=False)]
    conservedVarianceDf = conservedVarianceDf.reset_index(drop=True)

    # insert frequency count column, default=1 and absolute occurrence column, default =1
    conservedVarianceDf.insert(8, 'Variant Frequency %', 1)
    conservedVarianceDf.insert(9, 'Variant Occurrence', 1)

    # Get counts for each variants occurence and insert into dataframe
    count = conservedVarianceDf.value_counts(['Chr ID', 'Site', 'Seqbase'])
    variant = conservedVarianceDf.value_counts(['Chr ID', 'Site', 'Seqbase']).index
    for SNV, i in zip(variant, count):
        # Caculate relative frequency to all strains
        freq = (i / len(varianceDfList)) * 100
        conservedVarianceDf['Variant Frequency %'] = np.where(
        ((conservedVarianceDf['Chr ID']==SNV[0]) &
        (conservedVarianceDf['Site']==SNV[1]) & (conservedVarianceDf['Seqbase']==SNV[2]))
        , freq, conservedVarianceDf['Variant Frequency %'])

        # Insert absolute occurence
        conservedVarianceDf['Variant Occurrence'] = np.where(
        ((conservedVarianceDf['Chr ID']==SNV[0]) &
        (conservedVarianceDf['Site']==SNV[1]) & (conservedVarianceDf['Seqbase']==SNV[2]))
        , i, conservedVarianceDf['Variant Occurrence'])

    return conservedVarianceDf


def main(input):

    strainData = input[0]
    refData = input[1]

    varianceDfList = []
    for strain in strainData:
        varianceDfList.append(filterVariance(removeVarianceOverlap(strain, refData)[0])[0])

    return findVarianceConservation(varianceDfList)


if __name__ == '__main__':

    # Assign variables to individual strains to search conserved variance
    NGS_36 = r'C:/Users/James/Documents/conserved-genome-variance-bgi/input-example/JSy036/SNV_INDEL/'
    NGS_37 = r'C:/Users/James/Documents/conserved-genome-variance-bgi/input-example/JSy037/SNV_INDEL/'
    NGS_38 = r'C:/Users/James/Documents/conserved-genome-variance-bgi/input-example/JSy038/SNV_INDEL/'
    NGS_39 = r'C:/Users/James/Documents/conserved-genome-variance-bgi/input-example/JSy039/SNV_INDEL/'
    REF_SNV_INDEL_DIR = r'C:/Users/James/Documents/conserved-genome-variance-bgi/input-example/JSy001/SNV_INDEL/'

    # Strucutre into input for main fn
    input = ([NGS_36, NGS_37, NGS_38, NGS_39], REF_SNV_INDEL_DIR)

    conservedStrainDataFrame = main(input) # df holds collated strain information with conservation and filtering
    #conservedStrainDataFrame.to_csv('URA_escapees.csv') # output to .csv if required
