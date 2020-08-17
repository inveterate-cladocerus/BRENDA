from collections import namedtuple, defaultdict
import re
from pandas import DataFrame, ExcelWriter
from math import ceil

Enzyme = namedtuple("Enzyme", ["min_tn", "max_tn", "mean_tn", "subs", "min_ph", "max_ph", "mean_ph", "temp", "wt",
                               "mut", "descrp"])


class OrganismEnzymes:

    def __init__(self):
        self.specie = None
        self.org_num = None
        self.org_enzymes = list()

    def add_organism(self, data):
        pattern = r"PR\t\#(\d+)\# ([^(<]+)"
        search = re.match(pattern, data)
        specie = search.group(2).strip()
        org_num = search.group(1)
        self.specie = specie
        self.org_num = org_num

    @staticmethod
    def get_turnover_number(data):
        pattern = r"TN\t\#\d+(?:,\d+)*\# (?P<min_tn>\d+\.?\d*)\-?(?P<max_tn>\d+\.?\d*)?"
        search = re.match(pattern, data)
        tn = search.groupdict() if search else None
        min_tn = None
        max_tn = None
        mean_tn = None
        if tn:
            min_tn = tn["min_tn"]
            max_tn = tn["max_tn"] if tn["max_tn"] else tn["min_tn"]
            min_tn = float(min_tn) if "." in min_tn else int(min_tn)
            max_tn = float(max_tn) if "." in max_tn else int(max_tn)
            mean_tn = (min_tn + max_tn) / 2
        return min_tn, max_tn, mean_tn

    @staticmethod
    def get_substrate(data):
        pattern = r"TN\t\#\d+(?:,\d+)*\# \d+\.?\d*(?:\-\d+\.?\d*)? \{((?:(?!more).)+)\}"
        search = re.match(pattern, data, re.DOTALL)
        substrate = search.group(1) if search else None
        return substrate

    def get_ph(self, data):
        pattern_ph = r".+\#(?:\d+,)*" + self.org_num + r"(?:,\d+)*\#.+pH (?P<min_ph>\d{1,2}\.?\d*)\-?" \
                                                       r"(?P<max_ph>\d{1,2}\.?\d*)?"
        search_ph = re.match(pattern_ph, data, re.DOTALL)
        ph = search_ph.groupdict() if search_ph else None
        min_ph = None
        max_ph = None
        mean_ph = None
        if ph:
            min_ph = ph["min_ph"]
            max_ph = ph["max_ph"] if ph["max_ph"] else ph["min_ph"]
            min_ph = float(min_ph) if "." in min_ph else int(min_ph)
            max_ph = float(max_ph) if "." in max_ph else int(max_ph)
            mean_ph = (min_ph + max_ph) / 2
        return min_ph, max_ph, mean_ph

    def get_temperature(self, data):
        pattern_temp = r".+\#(?:\d+,)*" + self.org_num + r"(?:,\d+)*\#.+?(\d{1,3})\°C"
        search_temp = re.match(pattern_temp, data, re.DOTALL)
        temp = search_temp.group(1) if search_temp else None
        return temp

    def is_wild_type(self, data):
        pattern_wt = r".+\#(?:\d+,)*" + self.org_num + r"(?:,\d+)*\#(?!.+(?i:mutant|recombinant)).+(?i:wild|native)"
        search_wt = re.match(pattern_wt, data)
        result = True if search_wt else False
        return result

    def is_mutant(self, data):
        pattern_mut = r".+\#(?:\d+,)*" + self.org_num + r"(?:,\d+)*\#.+(?i:mutant|recombinant)"
        search_mut = re.match(pattern_mut, data)
        result = True if search_mut else False
        return result

    def add_enzyme(self, data):
        pattern = r"\{.+\}(?: )+"
        _, data2 = re.split(pattern, data)
        splitted_data2 = data2.split(";")
        min_tn, max_tn, mean_tn = self.get_turnover_number(data)
        substrate = self.get_substrate(data)
        for dat2 in splitted_data2:
            min_ph, max_ph, mean_ph = self.get_ph(dat2)
            temp = self.get_temperature(dat2)
            is_wild_type = self.is_wild_type(dat2)
            is_mutant = self.is_mutant(dat2)
            enzyme = Enzyme(min_tn, max_tn, mean_tn, substrate, min_ph, max_ph, mean_ph, temp, is_wild_type, is_mutant,
                            data2)
            self.org_enzymes.append(enzyme)


class EnzymesDatabase:

    def __init__(self):
        self.input_ec_numbers = set()
        self.reactions = dict()

    def process_input_ec_numbers(self, path_input_ec_numbers):
        with open(path_input_ec_numbers, 'r', encoding="utf-8") as file:
            for line in file:
                self.input_ec_numbers.add(line.strip())

    @staticmethod
    def get_organisms_number(data):
        org_nums = data.split()[1].strip('#')
        org_nums = org_nums.split(',')
        return org_nums

    @staticmethod
    def process_protein_data(data, enzymes):
        enzyme = OrganismEnzymes()
        enzyme.add_organism(data)
        enzymes[enzyme.org_num] = enzyme

    def process_turnover_number_data(self, data, enzymes):
        org_nums = self.get_organisms_number(data)
        for org_num in org_nums:
            enzymes[org_num].add_enzyme(data)

    def process_database(self, path_input_db):
        fields = {"AC", "AP", "CF", "CL", "CR", "EN", "EXP", "GI", "GS", "IC50", "ID", "IN", "KKM", "KI", "KM", "LO",
                  "ME", "MW", "NSP", "OS", "OSS", "PHO", "PHR", "PHS", "PI", "PM", "PR", "PU", "RE", "RF", "REN", "RN",
                  "RT", "SA", "SN", "SP", "SS", "ST", "SU", "SY", "TN", "TO", "TR", "TS"}
        fields_names = {"ACTIVATING_COMPOUND", "APPLICATION", "COFACTOR", "CLONED", "CRYSTALLIZATION", "ENGINEERING",
                        "EXPRESSION", "GENERAL_INFORMATION", "GENERAL_STABILITY", "IC50_VALUE", "INHIBITORS",
                        "KCAT_KM_VALUE", "KI_VALUE", "KM_VALUE", "LOCALIZATION", "METALS_IONS", "MOLECULAR_WEIGHT",
                        "NATURAL_SUBSTRATE_PRODUCT", "OXIDATION_STABILITY", "ORGANIC_SOLVENT_STABILITY", "PH_OPTIMUM",
                        "PH_RANGE", "PH_STABILITY", "PI_VALUE", "POSTTRANSLATIONAL_MODIFICATION", "PROTEIN",
                        "PURIFICATION", "REACTION", "REFERENCE", "RENATURED", "RECOMMENDED_NAME", "REACTION_TYPE",
                        "SPECIFIC_ACTIVITY", "SYNONYMS", "SUBSTRATE_PRODUCT", "STORAGE_STABILITY", "SOURCE_TISSUE",
                        "SUBUNITS", "SYSTEMATIC_NAME", "TURNOVER_NUMBER", "TEMPERATURE_OPTIMUM", "TEMPERATURE_RANGE",
                        "TEMPERATURE_STABILITY"}
        searched_fields = {"PR", "TN"}
        with open(path_input_db, 'r', encoding="utf-8") as file:
            ec_number = None
            current_field = None
            past_line = None
            past_field = None
            reaction = None
            for line in file:
                line = line.strip()
                if line:
                    if line.split('\t')[0] in fields:
                        current_field = line.split('\t')[0]
                        if current_field == "ID":
                            if past_field in searched_fields:
                                options = {"PR": self.process_protein_data, "TN": self.process_turnover_number_data}
                                options[past_field](past_line, enzymes)
                            past_field = None
                            ec_number = line.split()[1]
                            if reaction:
                                reaction["rxn_enzymes"] = enzymes
                                self.reactions[reaction["ec_number"]] = reaction
                            reaction = {"ec_number": ec_number}
                            enzymes = dict()
                    if ec_number in self.input_ec_numbers and current_field in searched_fields \
                            and line not in fields_names:
                        current_line = line
                        if current_line.split('\t')[0] != current_field:
                            current_line = past_line + ' ' + current_line
                        else:
                            if past_field in searched_fields:
                                options = {"PR": self.process_protein_data, "TN": self.process_turnover_number_data}
                                options[past_field](past_line, enzymes)
                        past_field = current_field
                        past_line = current_line
            if ec_number in self.input_ec_numbers and current_field in searched_fields:
                if past_field in searched_fields:
                    options = {"PR": self.process_protein_data, "TN": self.process_turnover_number_data}
                    options[past_field](past_line, enzymes)

    def generate_tables(self, tables_amount):
        input_size = len(self.input_ec_numbers)
        ecnum_per_table = ceil(input_size / tables_amount)
        ecnum_amount = 0
        i_table = 1
        writer = None
        table_data = None
        for i_ecnum, ec_number in enumerate(self.input_ec_numbers, 1):
            if ecnum_amount == 0:
                writer = ExcelWriter(f"output/EC_numbers{i_table}.xlsx", engine="xlsxwriter")
                table_data = defaultdict(list)
            reaction = self.reactions[ec_number]
            table_data["EC number"] += [ec_number for org_num in reaction["rxn_enzymes"] for _ in
                                        reaction["rxn_enzymes"][org_num].org_enzymes]
            table_data["Specie"] += [reaction["rxn_enzymes"][org_num].specie for org_num in
                                     reaction["rxn_enzymes"].keys() for _ in
                                     reaction["rxn_enzymes"][org_num].org_enzymes]
            table_data["Number"] += [org_num for org_num in reaction["rxn_enzymes"] for _ in
                                     reaction["rxn_enzymes"][org_num].org_enzymes]
            table_data["Tn [mmol subs/(mmol enz * s)]"] += [enz.min_tn if enz.min_tn == enz.max_tn
                                                            else f"{enz.min_tn}-{enz.max_tn}"
                                                            for org_num in reaction["rxn_enzymes"] for enz in
                                                            reaction["rxn_enzymes"][org_num].org_enzymes]
            table_data["Substrate"] += [enz.subs for org_num in reaction["rxn_enzymes"] for enz in
                                        reaction["rxn_enzymes"][org_num].org_enzymes]
            table_data["pH"] += [enz.min_ph if enz.min_ph == enz.max_ph else f"{enz.min_ph}-{enz.max_ph}" for org_num in
                                 reaction["rxn_enzymes"].keys() for enz in reaction["rxn_enzymes"][org_num].org_enzymes]
            table_data["Temperature [°C]"] += [enz.temp for org_num in reaction["rxn_enzymes"] for enz in
                                               reaction["rxn_enzymes"][org_num].org_enzymes]
            table_data["Wild-type"] += ["Yes" if enz.wt else None for org_num in
                                        reaction["rxn_enzymes"] for enz in
                                        reaction["rxn_enzymes"][org_num].org_enzymes]
            table_data["Mutant"] += ["Yes" if enz.mut else None for org_num in reaction["rxn_enzymes"] for enz in
                                     reaction["rxn_enzymes"][org_num].org_enzymes]
            table_data["Description"] += [enz.descrp for org_num in reaction["rxn_enzymes"] for enz in
                                          reaction["rxn_enzymes"][org_num].org_enzymes]
            ecnum_amount += 1
            if ecnum_amount == ecnum_per_table or i_ecnum == input_size:
                df = DataFrame(table_data)
                df.to_excel(writer, sheet_name="Enzymes", index=False)
                writer.save()
                ecnum_amount = 0
                i_table += 1

                
turnover_number_db = EnzymesDatabase()
turnover_number_db.process_input_ec_numbers("input/ecnumberCain.txt")
turnover_number_db.process_database("brenda_db/08_08_2020.txt")
turnover_number_db.generate_tables(1)
