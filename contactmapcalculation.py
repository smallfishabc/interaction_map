import numpy as np
import interaction_strength

def interaction_plotting_calc(interaction, raw_value, pairs, intertype, strength, targetmap):
    # Remove unnecessary pairs
    inot = strength
    indexlist = []
    for index, i in enumerate(interaction):
        if i != inot:
            indexlist.append(index)
    pairsnew = np.delete(pairs, indexlist, axis=0)
    reaction_value = np.delete(interaction, indexlist, axis=0)
    raw_value_new = np.delete(raw_value, indexlist, axis=0)
    target_new = np.delete(targetmap, indexlist, axis=0)
    # Test new function 10/14/2021
    overall_strength=interaction_strength.calculate_overall_strength_multiple(pairsnew, raw_value_new,
                                                    reaction_value, target_new)
    # 10/7/2021 Bugfixed: wrong data is passed to following functions
    # full_strength, overall_strength = interaction_strength.calculate_overall_strength(pairsnew, raw_value_new, reaction_value, target_new)
    # else:
    #     full_strength = []
    #     overall_strength = -1
    np.savetxt(intertype + str(inot) + 'interaction_pair.csv', pairsnew, delimiter=',')
    # I save two text file for testing
    # np.savetxt(intertype + str(inot) + 'full_strength.csv', full_strength, delimiter=',')
    return overall_strength


def interaction_map_calc(seq, length, interaction, raw_value, pairs, figname, targetmap):
    att1 = interaction_plotting_calc(interaction, raw_value, pairs, 'att', 1, targetmap)
    att2 = interaction_plotting_calc(interaction, raw_value, pairs, 'att', 2, targetmap)
    rep1 = interaction_plotting_calc(interaction, raw_value, pairs, 'rep', -1, targetmap)
    rep2 = interaction_plotting_calc(interaction, raw_value, pairs, 'rep', -2, targetmap)
    # Save and demonstrate the plot
    return (att1, att2, rep1, rep2)

