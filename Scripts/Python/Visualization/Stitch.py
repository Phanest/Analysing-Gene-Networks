import os
from os import listdir
from os.path import join, isdir
from PIL import Image

fromFolder = './Images/AllImages/Stage T/Stage TIS'
images = ['clinical_data_date_of_initial_pathologic_diagnosis',
          'clinical_data_gender',
          'clinical_data_pathologic_stage',
          'clinical_data_pathology_M_stage',
          'clinical_data_pathology_N_stage',
          'clinical_data_pathology_T_stage',
          'ips_ctla4_neg_pd1_neg',
          'ips_ctla4_neg_pd1_pos',
          'ips_ctla4_pos_pd1_pos',
          'ips_ctla4_pos_pd1_neg'
          ]
toFolder = './Images/Stage TIS'

files = [f for f in listdir(fromFolder) if isdir( join(fromFolder, f) ) ]

print( files )

for i in images:

    print(i)

    imgs = []

    FWidth = [0]
    FHeight = [0]

    for f in files:

        path = join(fromFolder, f)
        path = join(path, i + '.eps')

        Open = Image.open(path)
        imgs.append( Open )

        (width, height) = Open.size

        FWidth.append(FWidth[len(FWidth)-1] + width)
        FHeight.append(height)



    NewImage = Image.new('CMYK', (FWidth[len(FWidth)-1], max(FHeight)) )

    for j in range(0, len(imgs)):
        NewImage.paste(imgs[j], box=(FWidth[j], 0) )

    if not os.path.isdir(toFolder):
        os.mkdir(toFolder)

    try:
        name = join(toFolder, i + '.eps')
        NewImage.save(name, "eps")
        NewImage.close()
    except IOError:
        print('Error: ' + i)