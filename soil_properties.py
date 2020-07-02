def soil_props(soil_type, depth):
    """
    Parameters c, Ks, n, Beta, s_h, s_w, s_bal, s_fc, bulk_d:
    Laio et al., 2001, Plants in water-controlled ecosystems: active role in hydrologic processes and response to water stress: II. Probabilistic soil moisture dynamic

    Parameters p1 through p5:
    Ezlit et al., 2013, Modification of the McNeal Clay Swelling Model Improves Prediction of Saturated Hydraulic Conductivity as a Function of Applied Water Quality
    """
    # 5.7% clay
    class_1 = {'c':4.8, 'Ks':1000.0, 'n':0.42, 'Beta':12.7, 's_h':0.08,
               's_w':0.11, 's_bal':0.31, 's_fc':0.52, 'bulk_d':1.5,
               'p1':0.649, 'p2':0.003, 'p3':8.837, 'p4':4.046,
               'p7':0.008, 'p6':6.356, 'p5':30.818, 'CEC': 50}

    # 16.2% clay
    class_2 = {'c':6.5, 'Ks':800.0, 'n':0.43, 'Beta':13.8, 's_h':0.14,
               's_w':0.18, 's_bal':0.46, 's_fc':0.56, 'bulk_d':1.5,
               'p1':1.00, 'p2':0.912, 'p3':1.438, 'p4':7.29,
               'p7':0.204, 'p6':4.105, 'p5':-5.054, 'CEC': 150}

    # 48.5% clay
    class_3 = {'c':9.8, 'Ks':200.0, 'n':0.45, 'Beta':14.8, 's_h':0.19,
               's_w':0.24, 's_bal':0.57, 's_fc':0.65, 'bulk_d':1.2,
               'p1':0.449, 'p2':1.005, 'p3':0.846, 'p4':10.968,
               'p7':0.53, 'p6':4.0799, 'p5':-11.15, 'CEC': 300}

    if soil_type == "class_1":
        soil_dict = {**class_1}
    elif soil_type == "class_2":
        soil_dict = {**class_2}
    elif soil_type == "class_3":
        soil_dict = {**class_3}

    gapon = 0.01475
    mass = soil_dict['bulk_d']*depth
    soil_dict.update(Kg=gapon, Zr=depth, Msoil = mass)

    return soil_dict
