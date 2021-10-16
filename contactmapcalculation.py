import numpy as np
import interaction_strength


def interaction_plotting_calc(interaction, raw_value, pairs, intertype, strength, targetmap):
    # Based on interaction type value to remove the unnecessary pairs
    # inot is used to identify interaction type
    inot = strength
    # Create an empty list to store the index of unnecessary pairs
    indexlist = []
    # Find the unnecessary paris based on the interaction type and the inot value we desired
    for index, i in enumerate(interaction):
        if i != inot:
            indexlist.append(index)
    # Remove unnecessary pairs
    pairsnew = np.delete(pairs, indexlist, axis=0)
    reaction_value = np.delete(interaction, indexlist, axis=0)
    raw_value_new = np.delete(raw_value, indexlist, axis=0)
    target_new = np.delete(targetmap, indexlist, axis=0)
    # Test new function 10/14/2021
    # Calculate the overall interaction strength
    overall_strength = interaction_strength.calculate_overall_strength_multiple(pairsnew, raw_value_new,
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
    # Calculate the interaction strength only without producing an interaction plot
    att1 = interaction_plotting_calc(interaction, raw_value, pairs, 'att', 1, targetmap)
    att2 = interaction_plotting_calc(interaction, raw_value, pairs, 'att', 2, targetmap)
    rep1 = interaction_plotting_calc(interaction, raw_value, pairs, 'rep', -1, targetmap)
    rep2 = interaction_plotting_calc(interaction, raw_value, pairs, 'rep', -2, targetmap)
    # Save and demonstrate the plot
    return att1, att2, rep1, rep2
