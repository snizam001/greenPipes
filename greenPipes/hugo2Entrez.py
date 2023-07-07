import pandas as pd
import mygene
mg = mygene.MyGeneInfo()

def hugo2Entrez (infile,Species):
    try:
        d=pd.read_csv(infile,sep='\t',header=None)
    except IOError:
        print(colored("Hugo2Entrez: " + infile+' does not exist',
                      'red',
                      attrs=['bold']
                     )
             )
        exit()

    Species_mg=["human",
                "mouse",
                "rat",
                "fruitfly",
                "nematode",
                "zebrafish",
                "thale-cress",
                "frog",
                "pig"]

    if Species not in Species_mg:
        print(colored("Hugo2Entrez: " +
                      Species+
                      ' is not correct. Possible options are (Case sensitive):',
                      'red',
                      attrs=['bold']
                     )
             )
        print(Species_mg)
        exit()

    d.columns=['Gene']
    entrezIds=[]
    for j in range(0,d.shape[0]):
        if ";" in d.iloc[j,0]:
            print(d.iloc[j,0])
            print(colored('The file '+ infile +' contains ";" or any other special character.'+
                          ' Did you read the manual? it should not be present in the file',
                          'red',
                          attrs=['bold']
                         )
                 )
            exit()
        else:
            mgout=mg.query(d.iloc[j,0],
                           fields='entrezgene,taxid,symbol',
                           as_dataframe=True,
                           species=Species)
            if mgout.empty:
                print('Can not find entrez gene id for '+
                      d.iloc[j,0] +" for " +
                      Species)
            else:
                entrezIds.append(mgout['entrezgene'].tolist()[0])
    print ('entrez ids for genes present in MassSpectrometry datasets: ')
    print (entrezIds)
    return(entrezIds)
