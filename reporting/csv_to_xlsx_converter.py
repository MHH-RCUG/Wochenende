import pandas
import numpy
import sys

# Fabian Charly Friedrich, April 2020
# Convert CSV from Wochenende pipeline or Nextflow blast pipeline to XLSX using pandas
# Use the runbatch_csv_to_xlsx.sh script
# Or python3 csv_to_xlsx_converter.py <csv file>


def read_annot(file):
    lines = open(file, 'r').readlines()
    data = []


    for line in lines:	
        split_line = line.strip().split(' ')

        tmp = []
        for i in range(3):
            try:
                tmp.append(split_line[i])
            except IndexError:
                tmp.append("")

        try:
            tmp.append(' '.join(split_line[3:]))
        except IndexError:
            tmp.append("")

        data.append(tmp)

    return pandas.DataFrame(numpy.array(data))

def read_krakenuniq(file):
    lines = open(file, 'r').readlines()
    data = []
    for line in lines[3:]:
        tmp = []
        split_line = line.replace('\n', '').split('\t')

        for i in range(9):
            try:
                tmp.append(split_line[i])
            except IndexError:
                tmp.append('')

        data.append(tmp)

    return pandas.DataFrame(numpy.array(data[1:]), columns=data[0])


def main():
    file = sys.argv[1]

    if ("annot" in str(file)):
        # used for nextflow_blast output
        df = read_annot(file)
        # convert to Excel
        df.to_excel(file.replace('.csv', '.xlsx'), index=None, header=False)

    elif (".rep." in str(file)):
        # used for Wochenende reporting output
        df = pandas.read_csv(file)
        # convert to Excel
        df.to_excel(file.replace('.csv', '.xlsx'), index=None, header=True)

    elif "se.report.txt" in str(file):
        # used for kraken2 output
        df = pandas.read_csv(file, sep='\t', header=None)
        # convert to excel
        df.to_excel(file.replace('.txt', '.xlsx'), index=False, header=False)

    elif "report.txt" in str(file):
        # used for krakenuniq output
        df = read_krakenuniq(file)
        # convert to excel
        df.to_excel(file.replace('.txt', '.xlsx'), index=False, header=True)

    else:
        print("Could not detect CSV / text type")


if __name__ == '__main__':
    main()
