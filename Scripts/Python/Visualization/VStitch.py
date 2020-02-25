from os import listdir
from os.path import join, isdir
from PIL import Image

fromFolder = './Images'
images = ['Text',
          'clinical_data_date_of_initial_pathologic_diagnosis',
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
toFolder = './Images'

fromFiles = [f for f in listdir(fromFolder) if isdir( join(fromFolder, f) ) ]

print(fromFiles)

for folder in fromFiles:


    fromf = join(fromFolder, folder)
    tof = join(toFolder, folder)

    imgs = []

    FWidth = [0]
    FHeight = [0]

    print (folder)

    for i in images:

        path = join(fromf)
        path = join(path, i + '.eps')

        Open = Image.open(path)
        imgs.append( Open )

        (width, height) = Open.size

        FWidth.append(width)
        FHeight.append(FHeight[len(FHeight)-1] + height)


    NewImage = Image.new('CMYK', (max(FWidth), FHeight[len(FHeight)-1] + 550) )

    for j in range(0, len(imgs)):
        NewImage.paste(imgs[j], box=(0, FHeight[j] + 150) )

    try:
        name = join(tof, folder + '.eps')
        NewImage.save(name, "eps")
        NewImage.close()
    except IOError:
        print('Error: ' + folder)