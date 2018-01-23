"""MODULE FOR POPULATING MONGO DATABASE WITH VITEK REPORT OBJECTS OBTAINED FROM XML FILES"""

"""Import Dependencies"""
from sys import argv, exit
from bs4 import BeautifulSoup as Soup
from collections import defaultdict
import re
import PyPDF2 as pdfreader
import pymongo
import os

class BuildReportTree:
    """Generate a tree of hash tables to represent the reports extracted from XML file"""
    
    def __init__(self, path):
        """Instantiate BuildReportTree object using XML file path. Will generate report_array property, a list of string
        elements repesenting the list
        params:
        path -- binary string"""
        #There may be multiple reports in an xml file i.e. multiple isolates
        self.reports = []
        with open(path, "r") as f:
            handler = f.read()
            soup = Soup(handler, 'lxml')
            report_body = soup.findAll("report_body")
            for content in report_body:
                if str(content.contents).find("ReportData&gt") != -1:
                    report_strings = str(content.contents).replace("\n", "").split("&gt;&lt;")
                    report_array = list(map(lambda x: x.replace("/", ""), report_strings))
                    self.reports.append(report_array)

    def drop_not_body(self, tag):
        """Remove any elements that are not part of the report body
        params:
        tag -- element of report soup object"""
        
        return str(tag).find("ReportData&gt") != -1
    
    def build_tree(self):
        """Using current object property report_array, generate a tree structure to represent the report"""
        report_trees = []
        for report in self.reports:
            report_tree = dict()
            headings = {'ReportData': 0,
                'AstDetailedInfo': 0,
                'AstTestInfo':0}
            #Drop source_xmlstring
            def not_sourcexmlstring(x):
                return x.find("source_xmlstring") == -1
            report = list(filter(not_sourcexmlstring, report))
            #Find start index for each section
            for i, row in enumerate(report):
                if row in headings.keys():
                    headings[row] = i
            if not all(val == 0 for key, val in headings.items()):
                sections = dict()
                sections["ReportData"] = report[headings["ReportData"]:headings["AstDetailedInfo"]]
                sections["AstDetailedInfo"] = report[headings["AstDetailedInfo"]: headings['AstTestInfo']]
                sections["AstTestInfo"] = report[headings["AstTestInfo"]:len(report)]
                for header, section in sections.items():
                    report_tree = self.process_data(header, section, report_tree)
                report_trees.append(report_tree)
            else:
                return {"error":"Section index error, check report_array property for inconsistencies"}
        return report_trees
    
    def process_data(self, header, section, report_tree):
        """Create branch and leaves for passed section, add too tree and return structure
        params:
        header -- section header, as string, to be used as branch key
        section -- array of strings of section elements
        report_tree -- report tree structure"""
        
        section_data = dict()
        #remove any elements containing a single item
        section = list(filter(lambda x: len(x.split(" ")) > 1, section))
        #If this section is the AstTestInfo then pull out drug family names as keys, and assign value of phenotype
        #All other elements add as key, value pairs according to string value
        if header == "AstTestInfo":
            phenotype_data = []
            for row in section:
                if row.split(" ")[0] == "DrugFamily" or row.split(" ")[0] == "Phenotype":
                    phenotype_data.append(" ".join(row.split(" ")[1:]))
                else:
                    section_data[row.split(" ")[0]] = self.create_dict(" ".join(row.split(" ")[1:]))
            phenotype_data = list(self.split_list(phenotype_data, 2))
            phenotype_data = list(self.create_dict(" ".join(x)) for x in phenotype_data)
            section_data["phenotype_info"] = dict()
            for phenotype in phenotype_data:
                section_data["phenotype_info"].update({phenotype["familyName"]: phenotype["phenotypeName"]})
            report_tree[header] = section_data
            return report_tree
        #For all other sections, build dictionary using key value pairs according to string value
        for row in section:
            #If row is for Drug information, use drug name as key
            if row.split(" ")[0] == "AstDrugResultInfo":
                drug_key, values = self.get_drug_data(" ".join(row.split(" ")[1:]))
                section_data[drug_key] = values
            else:
                section_data[row.split(" ")[0]] = self.create_dict(" ".join(row.split(" ")[1:]))
            report_tree[header] = section_data
        return report_tree
    
    def get_drug_data(self, drug_info):
        """Take string of drug information, create dictionary with key as drug name, and value as dictionary of attributes"""
        
        drug_dict = self.create_dict(drug_info)
        drug_key = drug_dict["drugName"]
        values = {key: value for key, value in drug_dict.items() if key is not "drugName"}
        return drug_key, values

    def create_dict(self, string):
        """Take in string containing substrings of format *key*=*val*, seperate into key, value pairs and return as dictionary"""
        element_dict = dict()
        key_vals = list(map(lambda x: x.replace("\"", ""), string.split("\" ")))
        for key_val in key_vals:
            key, value = key_val[0:key_val.find("=")], key_val[key_val.find("=")+1:len(key_val)]
            if not self.confidential_data(key):
                element_dict[key] = self.format_val(value)
        return element_dict
    
    def split_list(self, l, n):
        """Split list into list of lists with length n. List length must equal n to yield
        params:
        l -- list to split
        n -- disired number of elements per list"""
        
        for i in range(0, len(l), n):
            if len(l[i:i+n]) == n:
                yield l[i:i+n]
                
    def format_val(self, string):
        """Check if value is interget or float"""
        
        if len(string) == 0:
            return string
        if all(x.isdigit() for x in list(string)):
            return int(string)
        else:
            try:
                return float(string)
            except:
                return string
        
    def confidential_data(self, string):
        """If key is a patient identifier return true"""
        
        if string.find('patient') != -1:
            return True
        else:
            return False

class BuildDatabase:
    """Using a supplied mongodb client, database name, and CD-ROM file pathway, this object attempts to populate the designated
    mongo database with report objects obtained from XML files on the target CD-ROM"""
    
    def __init__(self, mongoclient, dbname, dir_path):
        """Initislise object and set global variables"""
        
        self.db = mongoclient[dbname]
        self.file_path = os.fsencode(dir_path)
        
    def build(self):
        """Iterate over files in path specified, if they correspond to a report, add to database"""
        
        for file in os.listdir(self.file_path):
            filename = os.fsdecode(file)
            if 'reports_isolate' in filename:
                xml_obj = BuildReportTree(filename)
                report_trees = xml_obj.build_tree()
                for report_tree in report_trees:
                    if 'error' in report_tree.keys():
                        print('{}: {}'.format(filename, report_tree['error']))
                    else:
                        try:
                            insert_id = self.db.reports.insert_one(report_tree).inserted_id
                            print('{} inserted with id {}'.format(filename, insert_id))
                        except:
                            print('Failed to save {}'.format(filename))           

def getopts(argv):
    """Collect command-line options in a dictionary
    params:
    argv: list of command line arguments"""
    
    opts = {} 
    while argv:  
        if argv[0][0] == '-':
            opts[argv[0][1:]] = argv[1] 
        argv = argv[1:]
    return opts

if __name__ == '__main__':
    myargs = getopts(argv)
    if 'dbname' in myargs.keys():
        dbname = myargs['dbname']
    else:
        print("Please specify database name e.g '-dbname database1'")
        exit()
    if 'dir_path' in myargs.keys():
        dir_path = myargs['dir_path']
    else:
        print("Please specify target directory e.g '-dir_path C:\data\reports'")
        exit()
    client = pymongo.MongoClient()
    BuildDatabase(client, dbname, dir_path).build()
