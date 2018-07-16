#!/usr/bin/env python3

from collections import OrderedDict
from string import Template
import argparse
import json

VERSION="0.1.0"

def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('json', help='Input JSON')
    parser.add_argument('xml_out', help='Output XML')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args

class JsonXml(object):
    """
      Under arguments, dict's look like this:

      "summary": "read one or more arguments files and add them to the command line",
      "name": "--arguments_file",
      "synonyms": "NA",
      "type": "List[File]",
      "required": "no",
      "fulltext": "",
      "defaultValue": "[]",
      "minValue": "NA",
      "maxValue": "NA",
      "minRecValue": "NA",
      "maxRecValue": "NA",
      "kind": "optional",
      "options": []
    """
    def __init__(self, blob):
        """
        :param blob:
        """
        self.blob = blob
        self.type = blob['type']
        self.xml_out = self.reblob()
        self.xml_write = self.template_xml()
        self.xml_param_out = []
        self.xml_param_out.append(self.xml_write)
        for entry in self.template_sel():
            self.xml_param_out.append(entry)
        if self.xml_out['type'] == 'select':
            self.xml_param_out.append('\t\t</param>\n')

    def template_xml(self):
        """
        Get those formatting braces on the param entry.
        :return:
        """

        templates = {'integer': Template('\t\t<param argument="$argument" type="$type" optional="$optional" value="$value" min="$min" max="$max" label="$label" help="$help" />\n'),
                     'float': Template('\t\t<param argument="$argument" type="$type" optional="$optional" value="$value" min="$min" max="$max" label="$label" help="$help" />\n'),
                     'text': Template('\t\t<param argument="$argument" type="$type" optional="$optional" value="$value" label="$label" help="$help" />\n'),
                     'data': Template('\t\t<param argument="$argument" type="$type" optional="$optional" format="" label="$label" help="$help" />\n'),
                     'select': Template('\t\t<param argument="$argument" type="$type" optional="$optional" label="$label" help="$help" >\n'),
                     'boolean': Template('\t\t<param argument="$argument" truevalue="$truevalue" falsevalue="$falsevalue" type="$type" optional="$optional" checked="$checked" label="$label" help="$help" />\n')}
        try:
            return templates[self.xml_out['type']].substitute(self.xml_out)
        except:
            raise Exception('Type ' + self.xml_out['type'] + " not recognized.")

    def template_sel(self):
        """
        Include templates for the output section as well.
        :return:
        """
        sel_out = Template('\t\t\t<option value="$value" selected="$selected">$value</option>\n')
        for entry in self.sel_blob():
            yield sel_out.substitute(entry)

    def sel_blob(self):
        """
        Define the select blob, and then we can make a template out of the dict.
        <option value="normal_yes" selected="true">Yes</option>
        :return:
        """
        for entry in self.get_select_opts():
            if self.blob['defaultValue'] == entry:
                sel_blob = {'value': entry, 'selected': 'true'}
            else:
                sel_blob = {'value': entry, 'selected': 'false'}
            yield sel_blob

    def reblob(self):
        """
        <param name="max_records_in_ram" type="integer" size="10" value="500000" min="1" label="Max Records in RAM" help="When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed." />
        :return:
        """
        xml_out = {'name': self.blob['name'].lstrip('-').lower(),
                   'argument': self.blob['name'],
                   'type': self.xml_json_type_map(self.blob['type']),
                   'label': self.blob['name'].lstrip('-').replace('_', ' ').title(),
                   'optional': self.xml_json_req_map(self.blob['required']),
                   'value': self.blob['defaultValue'],
                   'truevalue': self.blob['name'],
                   'falsevalue': '',
                   'checked': self.blob['defaultValue'],
                   'max': self.xml_json_num_map(self.blob['maxValue']),
                   'min': self.xml_json_num_map(self.blob['minValue']),
                   'help': self.blob['summary']}
        return xml_out

    def get_select_opts(self):
        """
        For select types, need to know what the legal options are.
        :return:
        """
        try:
            for entry in self.blob['options']:
                yield entry['name']
        except:
            raise Exception("This select statement does not provide any legal options.")

    def xml_json_type_map(self, json_type):
        """
        Transform type field in GATK json file to what would be recognized in xml.
        Types can be:

        Boolean
        File
        Integer
        List[File]
        LogLevel
        String
        ValidationStringency
        boolean
        int

        I think if it is capital Boolean, this is a picard native command.

        :param json_type:
        :return:
        """
        xml_json_type_map = {
                             'Boolean': 'boolean',
                             'boolean': 'boolean',
                             'Integer': 'integer',
                             'int': 'integer',
                             'byte': 'integer',
                             'Double': 'float',
                             'double': 'float',
                             'String': 'text',
                             'File': 'data',
                             'List[File]': 'data',
                             'List[Double]': 'data',
                             'List[Integer]': 'data',
                             'List[String]': 'data',
                             'FeatureInput[VariantContext]': 'data',
                             'LogLevel': 'select',
                             'Implementation': 'select',
                             'IntervalMergingRule': 'select',
                             'IntervalSetRule': 'select',
                             'ValidationStringency': 'select',
                             'GenotypingOutputMode': 'select',
                             'OutputMode': 'select',
                             'PCRErrorModel': 'select',
                             'WriterType': 'select'
                             }

        if json_type in xml_json_type_map:
            return xml_json_type_map[json_type]
        return json_type

    def xml_json_req_map(self, json_req):
        """
        Transform required field in GATK json file to what would be recognized in xml.
        :param json_req:
        :return:
        """
        xml_json_req_map = {'no': 'true', 'yes': 'false'}
        if json_req in xml_json_req_map:
            return xml_json_req_map[json_req]
        return json_req

    def xml_json_num_map(self, json_num):
        """
        Transform numerical fields in GATK json file to what would be recognized in xml.
        :param json_num:
        :return:
        """
        xml_json_num_map = {'Infinity': '', '-Infinity': '', 'NA': ''}
        if json_num in xml_json_num_map:
            return xml_json_num_map[json_num]
        return json_num


class JsonCheetah(JsonXml):
    """
    Manage Cheetah string creation, to complement XML section creation.
    """

    def cheetah_template(self):
        """
        Produce templates for how things looks in the Cheetah section.
        :return:
        """
        cht_tmpl = Template('$argument $${$name}')
        return cht_tmpl.substitute(self.xml_out)


class JsonShell(object):
    """
    if picard, don't provide macro for annotations, not necessary
    summary == description
    description == help
    name == name
    """
    def __init__(self, filename, profile='17.09'):
        """

        """
        self.common = ("--version", "--showHidden", "--help", "--arguments_file", "--VERBOSITY",
                       "--VALIDATION_STRINGENCY", "--USE_JDK_INFLATER", "--USE_JDK_DEFLATER", "--TMP_DIR", "--QUIET",
                       "--MAX_RECORDS_IN_RAM", "--GA4GH_CLIENT_SECRETS", "--CREATE_MD5_FILE", "--CREATE_INDEX",
                       "--COMPRESSION_LEVEL", "--REFERENCE_SEQUENCE", "--OUTPUT")
        self.profile = profile
        self.xml_params = []
        self.cheetah_params = []
        with open(filename, 'rU') as myfile:
            self.json_file = json.load(myfile)
            for entry in self.json_file['arguments']:
                if entry['name'] not in self.common:
                    self.xml_params.extend(JsonXml(entry).xml_param_out)
                    self.cheetah_params.append(JsonCheetah(entry).cheetah_template())

    def get_shell(self, outfile):
        """
        Return the xml shell.
        :return:
        """
        self.handle_out = open(outfile, 'w')
        for k, v in self.build_shell_tmpl().items():
            self.handle_out.write(v.substitute(self.build_shell_dict()))
            self.handle_out.write('\n')
            if k == 'inputs':
                self.write_params()
        self.handle_out.close()

    def write_params(self):
        """
        Write through the params section.
        :return:
        """
        for entry in self.xml_params:
            self.handle_out.write(entry)

    def build_shell_dict(self):
        """
        This will house all values the templates need.
        :return:
        """
        shell_dict = {'id': self.json_file['name'].lower().split(' ')[0],
                      'name': self.json_file['name'],
                      'short_name': self.json_file['name'].split(' ')[0],
                      'profile': self.profile,
                      'description': self.json_file['summary'].rstrip(' '),
                      'summary': self.json_file['description']}
        return shell_dict

    def build_shell_tmpl(self):
        """
        Provide templates for the shell of the XML file.
        :return:
        """
        shell_tmpl = OrderedDict({'tool': Template('<tool id="gatk4_$id" name="$name" version="@VERSION@.0" profile="$profile">\n'),
                                'description': Template('\t<description>- $description</description>\n'),
                                'macros': Template('\t<macros>\n\t\t<import>macros.xml</import>\n\t\t<import>picard_macros.xml</import>\n\t</macros>\n'),
                                'expand': Template('\t<expand macro="requirements"/>\n\t<expand macro="version_cmd"/>\n'),
                                'command': Template('\t<command detect_errors="exit_code"><![CDATA[\n\t\t@CMD_BEGIN@ $short_name\n\t\t#include source=$$picard_ref_opts#\n\t\t#include source=$$picard_opts#\n\t\t#include source=$$picard_output_opts#\n\t]]></command>\n'),
                                'inputs': Template('\t<inputs>\n\t\t<expand macro="gatk_req_params" />\n\t\t<expand macro="picard_params" />'),
                                'inputs_close': Template('\t</inputs>\n'),
                                'outputs': Template('\t<outputs>\n\t\t<expand macro="picard_output_params" />\n\t</outputs>\n'),
                                'tests': Template('\t<tests>\n\t</tests>\n'),
                                'help': Template('\t<help><![CDATA[\n\t$summary\n\t@PICARD_HELP@\n\t]]></help>\n'),
                                'citations': Template('\t<expand macro="citations"/>\n'),
                                'tool_close': Template('</tool>')})
        return shell_tmpl


def main():
    """
    Not including (Picard):
    CREATE_INDEX - Handled by Galaxy internally.
    help - duh
    QUIET - can't think of a good reason to include this, people can just ignore it
    version - Should be included, if anywhere, under the version tag
    :return:
    """
    args = supply_args()
    myshell = JsonShell(args.json)
    myshell.get_shell(args.xml_out)

if __name__ == "__main__":
    main()