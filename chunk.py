import numpy as np


# Here we created a sequence chunk. The sequence chunk is the combination of several consecutive residue.
# We will compute the average interaction strength across the entire chunk
class Chunk():
    def __init__(self, residue, distance, name, chunk_size, contactmap, pairs):
        self.name = name
        self.residue = residue
        self.distance = distance
        self.chunk_size = chunk_size
        self.pair = self.search_pairs()
        self.left_residue = self.pair[0]
        self.right_residue = self.pair[1]
        self.left_residue_list = self.generate_residue_list('left')
        self.right_residue_list = self.generate_residue_list('right')
        self.pair_list = self.pair_list()

    def search_pairs(self):
        return [self.residue, self.residue + self.distance]

    def generate_residue_list(self, position):
        residue_list = []
        if position is 'left':
            i = self.left_residue
        if position is 'right':
            i = self.right_residue
        k = 0
        size = (self.chunk_size - 1) / 2
        while k < self.chunk_size:
            residue_list.append(i - size + k)
            k += 1
        return residue_list

    def pair_list(self):
        plist = []
        for index, i in enumerate(self.left_residue_list):
            pair = [i, self.right_residue_list[index]]
            plist.append(pair)
        return plist

    def pair_list(self):
        plist = []
        for index, i in enumerate(self.left_residue_list):
            pair = [i, self.right_residue_list[index]]
            plist.append(pair)
        return plist

    def calc_interact(self, contactmap, pairs):
        interaction_list = []
        temp_pairs = pairs.tolist()
        for i in self.pair_list:
            indexpairs = temp_pairs.index(i)
            interaction_value = contactmap[indexpairs]
            interaction_list.append(interaction_value)
        self.interact_list = interaction_list

    def sum_interact(self):
        average = np.average(self.interact_list)
        self.interact = average


#        print(average)


def mapping_chunk(name, distance, contact_map, pairs, chunk_size=5):
    location = []
    object = []
    value = []
    start = pairs[0][0]
    end = pairs[-1][-1]
    i = start
    while i < end:
        test = Chunk(i, distance, name + str(i), chunk_size, contact_map, pairs)
        if test.left_residue_list[0] >= start and test.right_residue_list[-1] <= end:
            chunk = test
            chunk.calc_interact(contact_map, pairs)
            chunk.sum_interact()
            object.append(chunk)
            value.append(chunk.interact)
            location.append(chunk.pair)
        i += 1
    location = np.array(location)
    value = np.array(value)
    # simple=dict(zip(location,value))
    return location, value, object


def full_mapping_chunk(name, contact_map, pairs, chunk_size=5):
    location_list = []
    value_list = []
    start = pairs[0][0]
    end = pairs[-1][-1]
    length = end - start
    max_distance = length - 5
    for j in range(5, max_distance + 1):
        location_fragment, value_fragment , a = mapping_chunk(name, j, contact_map, pairs)
        location_list.append(location_fragment)
        value_list.append(value_fragment)
    return (location_list, value_list)
# Input a list containing all special residues. Analysing list and plot all the interaction containing at least
# one special residue
def mapping_chunk_special_residue_list(special_list,name,residue,contact_map,pairs):
    sum_location_list=[]
    sum_value_list=[]
    for i in special_list:
        location,value=mapping_chunk_certain_residue_single(name,i,contact_map,pairs)
        sum_location_list.append(location)
        sum_value_list.append(value)
    return sum_location_list,sum_value_list

def mapping_chunk_certain_residue_single(name,residue,contact_map,pairs,chunk_size=5):
    value=[]
    location=[]
    start = pairs[0][0]
    end = pairs[-1][-1]
    index_list=[]
    for index,i in pairs:
        if i[0]==residue:
            index_list.append(index)
        elif i[1]==residue:
            index_list.append(index)
    for k in index_list:
        distance= pairs[index][1]-pairs[index][0]
        left_residue=pairs[index][0]
        chunk = Chunk(i, distance, name + str(left_residue), chunk_size, contact_map, pairs)
        chunk.calc_interact(contact_map, pairs)
        chunk.sum_interact()
        value.append(chunk.interact)
        location.append(chunk.pair)
    location = np.array(location)
    value = np.array(value)
    return location,value
