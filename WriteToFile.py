import csv


def man():
    

    with open('pressure_drop_heat.txt', newline='') as f:
        
        for line in f:
            line = line.split()
            print(line[0])
            # line2 = line.split()
            # line3 =line2[0]
            xcalibor = int(line[0]) +3756
            heattransfer = abs(float(line[2]))
            with open('shape_{}.txt'.format(xcalibor), 'w') as ff:
                for i in range(159):
                    ff.write('{} {} \n'.format(i, heattransfer))



if __name__ == '__main__':
    man()
