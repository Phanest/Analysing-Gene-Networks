from PIL import ImageDraw, ImageFont
from os import listdir
from os.path import join, isdir
from PIL import Image

fromFolder = './Images/AllImages/Pathologic Stage/Stage I'
images = ['clinical_data_date_of_initial_pathologic_diagnosis'
          ]
toFolder = './Images'

files = [f for f in listdir(fromFolder) if isdir( join(fromFolder, f) ) ]

print(files)

imgs = []

FWidth = [0]
FHeight = [0]

for f in files:
    print(f)

    path = join(fromFolder, f)
    path = join(path, images[0] + '.eps')

    Open = Image.open(path)

    (width, height) = Open.size

    FWidth.append(FWidth[len(FWidth) - 1] + width)
    FHeight.append(height)

    font = ImageFont.truetype("arial.ttf", 70, encoding="unic")
    canvas = Image.new('CMYK', (width + 10, height + 10))
    draw = ImageDraw.Draw(canvas)
    draw.text((int(width / 2), int(height / 2)), f, 'black', font)

    imgs.append(canvas)

NewImage = Image.new('CMYK', (FWidth[len(FWidth) - 1], max(FHeight)))

for j in range(0, len(imgs)):
    # NewImage.paste(imgs[j], box=(FWidth[j], 0) )
    NewImage.paste(imgs[j], box=(FWidth[j], 0))

try:
    name = join(toFolder, 'Text' + '.eps')
    NewImage.save(name, "eps")
    NewImage.close()
except IOError:
    print('Error: ' + f)