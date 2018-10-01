#! /usr/bin/python3

def openVariantsPickle(fileName):
    if not fileName:
        raise RuntimeError("File name has to be set")
    import pickle
    file = open(fileName, 'rb')
    data = pickle.load(file)
    file.close()
    return data

if __name__ == "__main__":
    fileName = "LB2964tester.pkl"
    data = openVariantsPickle(fileName)
    keycount = 0
    for key in data.keys():
        print(key)
        keycount += 1
    print(keycount)
    print("Data headers: %s" %data.keys())