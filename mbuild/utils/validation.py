def assert_port_exists(port_name, compound):
    """Ensure that a Port label exists in a Compound.  """
    if port_name in compound.labels:
        return True
    else:
        from mbuild.port import Port
        available_ports = [name for name in compound.labels
                           if isinstance(compound.labels[name], Port)]
        compound_name = compound.__class__.__name__
        raise ValueError("No port named '{port_name}' in {compound_name}'s"
                         " labels. Labeled Ports in {compound_name} are:"
                         " {available_ports}".format(**locals()))
