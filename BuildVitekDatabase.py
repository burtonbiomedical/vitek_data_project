"""MODULE FOR POPULATING MONGO DATABASE WITH VITEK REPORT OBJECTS OBTAINED FROM XML FILES"""

"""Import Dependencies"""
from sys import argv, exit
from bs4 import BeautifulSoup as Soup
from collections import defaultdict
import re
import PyPDF2 as pdfreader
import pymongo
import os
from datetime import datetime

class BuildReportTree:
    """Generate a tree of hash tables to represent the reports extracted from XML file"""

    def __init__(self, path):
        """Instantiate BuildReportTree object using XML file path. Will generate report_array property, a list of string
        elements repesenting the list
        params:
        path -- binary string"""
        #There may be multiple reports in an xml file i.e. multiple isolates
        self.lab_reports = {}

        with open(path, "r") as f:
            handler = f.read()
            soup = Soup(handler, 'lxml')
            lab_reports_soup = soup.find_all('lab_report')
            for report in lab_reports_soup:
                id_ = re.compile(r'<lab_report id="([0-9]+)">').search(str(report)).group(1)
                report_array = str(report.find("source_xmlstring")).replace("\n", "").split("&gt;&lt;")
                report_array = list(map(lambda x: x.replace("/", ""), report_array))
                self.lab_reports[id_] = report_array

    def build_trees(self):
        """Using current object property report_array, generate a tree structure to represent the report"""
        try:
            document_tree = {}
            lab_reports = []
            for id_, report in self.lab_reports.items():
                isolate_branch = dict()
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
                        isolate_branch = self.init_document_tree(header, section, isolate_branch)
                    lab_reports.append({
                        'isolate_id': id_,
                        'isolate_data': isolate_branch
                    })
                else:
                    return {"error":"Section index error, check report_array property for inconsistencies"}
            document_tree['lab_reports'] = lab_reports
            try:
                document_tree['organism_summary'] = self.init_org_summary_tree(lab_reports)
            except:
                return {"error":"Failed to build organism_summary branch"}
            return document_tree
        except:
            return {"error":"Fatal error when building document tree"}

    def init_org_summary_tree(self, lab_reports):
        """Build organism summary - unique organism and MIC data from from report
        params:
        lab_reports -- tree structure for laboratory reports"""

        org_summary_array = []
        org_summary_tree = []
        iso_num = 0
        for isolate_branch in lab_reports:
            isolate_summary = {}
            id_ = isolate_branch['isolate_id']
            isolate_data = isolate_branch['isolate_data']
            #Get organism name
            org = isolate_data['AstTestInfo']['SelectedOrg']['orgFullName']
            isolate_summary['organism_name'] = org
            isolate_summary['mic_data'] = []
            #Get drug data
            drug_data = isolate_data['AstDetailedInfo']
            for drug in drug_data:
                if type(drug['details']['mic']) == str:
                    drug_result = drug['details']['interpretation']
                    key = 'interpretation'
                else:
                    drug_result = drug['details']['mic']   
                    key = 'mic'
                isolate_summary['mic_data'].append({'drug':drug['drug'], key: drug_result})
            #If this organism is not unique in MIC values/organism species, do not add to tree
            if isolate_summary not in org_summary_array:
                org_summary_tree.append({
                    'isolate_id': 'isolate_'+str(iso_num),
                    'isolate_data': isolate_summary
                })
                iso_num += 1
                org_summary_array.append(isolate_summary)
        return org_summary_tree

    def init_document_tree(self, header, section, document_tree):
        """Create branch and leaves for passed section, add too tree and return structure
        params:
        header -- section header, as string, to be used as branch key
        section -- array of strings of section elements
        document_tree -- report tree structure"""

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
            document_tree[header] = section_data
            return document_tree
        #If section is AstDetailedInfo, sort for Drug information
        elif header == 'AstDetailedInfo':
            document_tree[header] = []
            for row in section:
                drug_key, values = self.get_drug_data(" ".join(row.split(" ")[1:]))
                document_tree[header].append({
                    'drug': drug_key,
                    'details': values
                })
        #Report data save as just key value pairs
        else:
            for row in section:
                section_data[row.split(" ")[0]] = self.create_dict(" ".join(row.split(" ")[1:]))
            document_tree[header] = section_data
        return document_tree

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

    def __init__(self, mongoclient, dbname, dir_path, error_path):
        """Initislise object and set global variables"""

        self.db = mongoclient[dbname]
        self.file_path = dir_path
        self.error_path = error_path
        self.errors = []

    def build(self):
        """Iterate over files in path specified, if they correspond to a report, add to database"""

        for file in os.listdir(self.file_path):
            filename = os.fsdecode(file)
            if 'reports_isolate' in filename:
                xml_obj = BuildReportTree(str(self.file_path) + filename)
                try:
                    document_tree = xml_obj.build_trees()
                    if 'error' in document_tree.keys():
                         print('{}: {}'.format(filename, document_tree['error']))
                         self.errors.append("{} ERROR: {} FILENAME: {}".format(str(datetime.now()), document_tree['error'], filename))
                    else:
                         self.insert_report(document_tree, filename)
                    self.log_errors()
                except:
                    print("Fatal error on {}, failed to build document tree".format(filename))
                    self.errors.append("{} FATAL ERROR, UNABLE TO BUILD DOC TREE. FILENAME: {}".format(str(datetime.now()), filename))
                    self.log_errors()


    def insert_report(self, document_tree, filename):
        """Attempt to save report tree structure as new document in report collection
        params:
        document_tree -- nested hash tables representing the report
        filename -- string path of file currently being processed"""

        try:
            insert_id = self.db.reports.insert_one(document_tree).inserted_id
            print('{} inserted with id {}'.format(filename, insert_id))
            self.insert_org(document_tree, insert_id)
        except:
            print('Failed to save {}'.format(filename))
            self.errors.append("{} RECORD NOT SAVED. FILENAME: {}".format(str(datetime.now()), filename))

    def insert_org(self, document_tree, document_id):
        """Check if organism exists in organism collection, if not add new organism, else add report ID to list
        of report id's for this organism
        params:
        document_tree -- nested hash tables representing the report
        document_id -- mongo id for report document"""
        org_summary = document_tree['organism_summary']
        unique_orgs = set()
        for isolate in org_summary:
            unique_orgs.add(isolate['isolate_data']['organism_name'])
        for org_name in unique_orgs:
            if self.db['orgs'].find_one({org_name: {'$exists': True}}):
                try:
                    org_doc = self.db['orgs'].find_one({org_name: {'$exists': True}})
                    org_doc[org_name].append(document_id)
                    self.db.orgs.update_one({'_id': org_doc['_id']}, {'$set': org_doc}, upsert=False)
                    print("{} summary updated".format(org_name))
                except:
                    print("Failed to update org summary: {} with report id: {}".format(org_name, document_id))
                    self.errors.append("{} REPORT: {} NOT ADDED TO ORG SUMMARY FOR {}".format(str(datetime.now()), document_id, org_name))
            else:
                try:
                    new_org = {org_name:[document_id]}
                    insert_id = self.db.orgs.insert_one(new_org).inserted_id
                    print('Create new summary entry for organism {}, with id {}'.format(org_name, insert_id))
                except:
                    print("Failed to insert new organism summary for: {}".format(org_name))
                    self.errors.append("{} FAILED TO CREATE NEW ORG SUMMARY: {}".format(str(datetime.now()), orgname))

    def log_errors(self):
        """Insert errors into error log"""
        with open(self.error_path, 'w') as f:
            for error in self.errors:
                f.write(error+'\n')


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
    if 'error_path' in myargs.keys():
        error_path = myargs['error_path']
    else:
        print("Please specify target for error log e.g '-error_path C:\data\errors.txt'")
        exit()
    client = pymongo.MongoClient()
    BuildDatabase(client, dbname, dir_path, error_path).build()
