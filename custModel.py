# Cargar paquetes
import cobra
import cobra.core
import json
from cobra.core import Model, Reaction, Metabolite
from cobra.io import write_sbml_model
from cobra.io import read_sbml_model
from json import dump
from cobra.flux_analysis.parsimonious import optimize_minimal_flux
import xml.etree.ElementTree as ET
from cobra.flux_analysis.loopless import add_loopless, loopless_solution
from cobra.flux_analysis import pfba
from cobra.util import solvers

def convert_to_irreversible_with_indicators(cobra_model, mutually_exclusive_directionality_constraint = False):
    """Will break all of the reversible reactions into two separate irreversible
     reactions with different directions.  This function call modified from
     a version in the core cobra to facilitate the MILP formulation and
     include gene_reaction_rules with the reverse reaction

     Arguments:
      cobra_model: A model object which will be modified in place.
      mutually_exclusive_directionality_constraint: Boolean.  If True, turnover
       reactions are constructed to serve as MILP constraints to prevent loops.

     Returns:
      None, cobra_model is modified in place

    NOTE: This function has been modified from the manipulate module


    """
    reactions_to_add = []
    #from cobra.core.Reaction import Reaction
    #from cobra.core import Metabolite

    for reaction in cobra_model.reactions:
        # Potential artifact because a reaction might run backwards naturally
        # and this would result in adding an empty reaction to the
        # model in addition to the reverse reaction.
        if reaction.lower_bound < 0:
            #reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction = reaction.copy()
            reverse_reaction.id = reaction.id + "_reverse"
            reverse_reaction.lower_bound = 0
            reverse_reaction.upper_bound = reaction.lower_bound * -1
            reaction.lower_bound = 0
            # Make the directions aware of each other
            reaction.reflection = reverse_reaction
            reverse_reaction.reflection = reaction
            reaction.reversibility = 0
            reverse_reaction.reversibility = 0
            reaction_dict = {}
            current_metabolites = list(reaction.metabolites)
            for the_metabolite in current_metabolites:
                reaction_dict[the_metabolite] = -2 * reaction.get_coefficient(the_metabolite.id)
            reverse_reaction.add_metabolites(reaction_dict)
            reactions_to_add.append(reverse_reaction)
            # Also: GPRs should already copy
            # reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            # reverse_reaction._genes = reaction._genes

            if mutually_exclusive_directionality_constraint:
                # A continuous reaction bounded by 0., 1.
                # Serves as a source for the indicator metabolites
                tmp_source = Reaction('IRRMILP_direction_constraint_source_for_%s_and_%s'
                                                           %(reaction.id,
                                                             reverse_reaction.id))
                tmp_source.upper_bound = 1.
                tmp_source.lower_bound = 0.
                # The reverse indicator reaction is
                # an integer-valued reaction bounded by 0,1
                # that activates flux to the reverse reaction
                # and deactivates the forward reaction only when it is
                # turned on to 1
                tmp_indicator = Reaction('IRRMILP_reverse_indicator_for_%s_and_%s'
                                                           %(reaction.id,
                                                             reverse_reaction.id))
                tmp_indicator.upper_bound = 1
                tmp_indicator.lower_bound = 0
                tmp_indicator.variable_kind = 'integer'
                flux_constraint_forward = Metabolite(id =
                     'IRRMILP_direction_constraint_for_%s'%reaction.id)
                flux_constraint_reverse = Metabolite(id =
                     'IRRMILP_direction_constraint_for_%s'%reverse_reaction.id)
                flux_constraint_reverse._constraint_sense = 'G'
                flux_constraint_reverse._bound = 0.

                tmp_source.add_metabolites({flux_constraint_forward: 1})

                tmp_indicator.add_metabolites({flux_constraint_forward: -1,
                                      flux_constraint_reverse: 1})
                if reaction.upper_bound != 0:
                        reaction.add_metabolites({flux_constraint_forward: -1./reaction.upper_bound})
                else:
                    # could put 1.01 X the tolerance here,
                    # This is arbitrary.  Use 0.001
                    # since 1000 is a typical upper bound
                    reaction.add_metabolites({flux_constraint_forward: -0.001})
                if reverse_reaction.upper_bound != 0:
                    reverse_reaction.add_metabolites({flux_constraint_reverse: -1./reverse_reaction.upper_bound})
                else:
                    reverse_reaction.add_metabolites({flux_constraint_reverse: -0.001})
                reactions_to_add.append(tmp_indicator)
                reactions_to_add.append(tmp_source)
    cobra_model.add_reactions(reactions_to_add)

def add_turnover_metabolites(cobra_model, metabolite_id_list, epsilon):
    """ NOTE: Model must first be converted to irreversible!
    This entry creates a corresponding turnover metabolite
    that ensures flux through the metabolite of interest.

    Arguments:
     cobra_model: the model to be updated.
     metabolite_id_list: list of model metabolites for
                         which to add a turnover metabolite.
     epsilon: minimal flux to force through turnover metabolites.
      recommend 1.01 X solver_tolerance


    """
    from cobra.core import Reaction
    from cobra.core import Metabolite
    from numpy import abs

    # Set the minimum flux for metabolites equal to some factor larger than the solver's tolerance
    the_min_flux = epsilon

    turnover_metabolites = []
    sink_reactions = []
    for metabolite_id in metabolite_id_list:
        r_metabolite = cobra_model.metabolites.get_by_id(metabolite_id)
        v_metabolite = Metabolite("TM_" + metabolite_id,
                                  compartment = r_metabolite.compartment,
                                  name = ("TM_" + r_metabolite.name),
                                  charge = r_metabolite.charge,
                                  formula = r_metabolite.formula)
        # Now for reactions.  We include all reactions
        # that create or consume the real metabolite.
        # These reactions therefore also drive creation of the
        # turnover metabolite, and we need to add a reaction
        # that is the sink.  By constraining this single
        # sink reaction, we ensure flux through the real reactions
        # to create the real metabolite.
        sum_abs_source_reaction_bounds = 0.01

        the_reaction_id_list = [x.id for x in r_metabolite.reactions]
        for the_reaction_id in the_reaction_id_list:
            the_reaction = cobra_model.reactions.get_by_id(the_reaction_id)
            coefficient = abs(the_reaction.get_coefficient(r_metabolite))
            the_reaction._metabolites[v_metabolite] = coefficient
            v_metabolite._reaction.add(the_reaction)
            # The model should be irreversible, and the upper bound
            # should be greater than the lower bound and positive.
            # Use 1000 for each reaction to be conservative.
            # E.g. might flip reactions back on later, don't
            # use the current upper bound for the reactions.
            sum_abs_source_reaction_bounds += 100.

        # Add the sink reaction for the turnover metabolite
        # Since both creation and consumption of
        # the real metabolite drive the turnover,
        # we require 2 units of metabolite per
        # 1 unit of flux sink so this matches
        # the turnover through the real metabolite.
        sink_reaction = Reaction("TMS_" + metabolite_id)
        sink_reaction.add_metabolites({v_metabolite:-2})

        # Ensure a positive flux through the sink.
        # and maximum equal to maximal needed to
        # sustain reactions at their maximum.
        sink_reaction.lower_bound = the_min_flux
        sink_reaction.upper_bound = sum_abs_source_reaction_bounds / 2.
        v_metabolite._reaction.add(sink_reaction)
        turnover_metabolites.append(v_metabolite)
        sink_reactions.append(sink_reaction)
    cobra_model.add_metabolites(turnover_metabolites)
    cobra_model.add_reactions(sink_reactions)



#from gim3e.core import gim3e

print("Load initial model")
# Load sbml model
sbml_fname = input("Type path to metabolic model (sbml format): ")
model = cobra.io.read_sbml_model(sbml_fname)


print("Irreversible reactions")
## Transform reversible reactions
convert_to_irreversible_with_indicators(model)

print("Add Turnover metabolites and Sink reactions")
### Run the following piece of code (commented) if the identification of measured metabolites
### in the model is needed.

# kegg = ["C00033","C00670","C00158","C00055","C01697","C00469","C00037","C00051","C00191",
#         "C00116","C00025","C00082","C00079","C00041","C00148","C00152","C00407","C00135",
#         "C00047","C00065","C00186","C00049","C00137","C00378","C00042","C00029","C00043",
#         "C00221","C00002","C00064","C00123","C00183","C00003","C00078","C01083",
#         "C00354","C00144","C00267",
#         'C00031']
# metabolites = ["Acetic acid", "Glycerophosphocholine", "Citric acid", "Cytidine monophosphate",
#                "Galactitol", "Ethanol", "Glycine", "Glutathione", "D-Glucuronic acid", "Glycerol",
#                "L-Glutamic acid", "L-Tyrosine", "L-Phenylalanine", "L-Alanine", "L-Proline", "L-Asparagine",
#                "L-Isoleucine", "L-Histidine", "L-Lysine", "L-Serine", "L-Lactic acid", "L-Aspartic acid",
#                "myo-Inositol", "Thiamine", "Succinic acid", "Uridine diphosphate glucose",
#                "Uridine diphosphate-N-acetylglucosamine", "L-Arginine",
#                "Adenosine triphosphate", "L-Glutamine", "L-Leucine", "L-Valine", "NAD", "L-Tryptophan",
#                "Trehalose", "Fructose 1,6-bisphosphate", "Guanosine monophosphate", "Alpha-D-Glucose"]
#

### Get metabolites id in the model

# tree = ET.parse(sbml_fname)
# root = tree.getroot()
# 
# prefix = 'http://identifiers.org/kegg.compound/'
# kegg_match = []
# for k in kegg:
#   kegg_match.append(prefix + k)
# 
# mets = []
# for met in root[1][7]:
#   istrue = 0
# 	for i in met[0][0][0][0][0]:
#   	if i.attrib['{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource'] in kegg_match:
#   		istrue += 1
# 	if istrue > 0:
#   	name = met.attrib['name']
# 		aux = met.attrib['id'].split("_")
# 		body = "_".join(aux[1:(len(aux)-1)])
# 		id_new = ("_".join(aux[1:(len(aux)-1)]) + "_ms")
# 		idm = "_".join(aux[1:(len(aux))])
# 		if aux[(len(aux)-1)] == "c":
#   		mets.append(idm)

### For toy example
mets = ["ACETATE"]
            
add_turnover_metabolites(model, mets, 0)

print("Save model")
cobra.io.write_sbml_model(cobra_model=model, filename="toy/toy.xml")
cobra.io.save_matlab_model(model, "toy/toy.mat")
