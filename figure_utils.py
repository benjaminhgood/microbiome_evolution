
def get_pretty_species_name(species_name, include_number=False):
    
    items = species_name.split("_")
    
    pretty_name = "%s %s" % (items[0], items[1])
    
    if include_number:
        pretty_name += (" (%s)" % (items[2]))
        
    return pretty_name
    
def get_abbreviated_species_name(species_name):
    
    items = species_name.split("_")
    
    pretty_name = "%s. %s" % (items[0][0], items[1])
        
    return pretty_name