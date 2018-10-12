#!/usr/bin/env python

# Potentially make List[String] types a repeat xml element.  Otherwise, it appears these are separated by commas on the command line.
# Annotation options need macros.
# Macros for read filter options?  These have their own json files it seems...
# gVCF not currently set to be different than VCF.

from lxml import etree
from string import Template
from xml.sax.saxutils import escape
import argparse
import json
import pypandoc

VERSION="0.2.0"

def supply_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('--json', help='Input JSON')
    parser.add_argument('--xml_out', help='Output Directory')
    parser.add_argument('--version', action='version', version='%(prog)s ' + VERSION)
    args = parser.parse_args()
    return args


class Mappings(object):
    """
    Hold all of the mappings we need to provide for parameters.
    """
    def __init__(self):
        # These are parameters that are already covered in the macros.xml file.  They also may be options that we
        # don't want included in any GATK wrapper.
        # REFERENCE_SEQUENCE: If reference is denoted optional then it should go in the optional section.  Build function to place it.
        # Should be able to utilize self.blob['required'] to determine section (either optional or main), then look at self.macro_to_chth to determine
        # the macro that REFERENCE_SEQUENCE is associated to.
        self.xml_json_type_map = {
                             'Boolean': 'boolean',
                             'boolean': 'boolean',
                             'Long': 'integer',
                             'long': 'integer',
                             'Integer': 'integer',
                             'int': 'integer',
                             'byte': 'integer',
                             'Double': 'float',
                             'double': 'float',
                             'Float': 'float',
                             'float': 'float',
                             'String': 'text',
                             'Set[String]': 'text',
                             'File': 'data',
                             'Set[File]': 'data',
                             'List[File]': 'data',
                             'List[Double]': 'text',
                             'List[Integer]': 'integer',
                             'List[String]': 'text',
                             'List[Type]': 'text',
                             'List[FeatureInput[VariantContext]]': 'data',
                             'ArrayList[String]': 'text',
                             'FeatureInput[VariantContext]': 'data',
                             'LogLevel': 'select',
                             'Implementation': 'select',
                             'IntervalMergingRule': 'select',
                             'IntervalSetRule': 'select',
                             'ReferenceConfidenceMode': 'select',
                             'ValidationStringency': 'select',
                             'GenotypingOutputMode': 'select',
                             'OutputMode': 'select',
                             'PCRErrorModel': 'select',
                             'WriterType': 'select',
                             'NumberAlleleRestriction': 'select'
                             }
        self.xml_json_num_map = {'Infinity': '', '-Infinity': '', 'NA': ''}
        self.xml_json_req_map = {'no': 'true', 'yes': 'false'}
        # If the parameter is in here, it won't be included in either the cheetah section or the param section.
        self.common_args = ("version", "showHidden", "help", "arguments_file", "VERBOSITY",
                            "VALIDATION_STRINGENCY", "USE_JDK_INFLATER", "USE_JDK_DEFLATER", "TMP_DIR", "QUIET",
                            "MAX_RECORDS_IN_RAM", "GA4GH_CLIENT_SECRETS", "CREATE_MD5_FILE", "CREATE_INDEX",
                            "COMPRESSION_LEVEL", "REFERENCE_SEQUENCE", "OUTPUT", "SEQUENCE_DICTIONARY", "INPUT",
                            "input", "reference", "output", "annotation", "annotation_group", "annotations_to_exclude",
                            "intervals", "exclude_intervals", "read_index", "interval_padding",
                            "interval_exclusion_padding")
        self.read_filter_args = ('')
        # There seems to be a significant amount of information that can't be retrieved from the provided json.
        # Make this in to a config file.
        # Might be able to classify a bunch of this based on input output files for a particular tool.  For instance, vcf_tabix
        # is applicable if the input is a VCF.
        # Optional status of a partiular parameters could be connected to macros as well.  For instance, ref_sel is the macro for the reference genome.
        # So, the 'ref_sel': 'reference' relation could be used along with 'required' status to place these.
        # One challenge with this is the macros would have to either be prebuilt with the necessary section variable, or
        # would have to be rewritten on the fly.
        self.tool_data = {'SortVcf':
                              {'output_fmt': {},
                               'pre_tmpls': ['vcf_tabix_multi'],
                               'post_tmpls': ['picard_opts', 'picard_vcf_output_opts', 'vcf_input_multi', 'picard_ref_opts_opt', 'picard_seqdict_opts'],
                               'pre_params': ['picard_vcf_params', 'vcf_input_params_multi'],
                               'opt_params': ['ref_sel', 'seq_dict_sel'],
                               'post_params': ['picard_params'],
                               'output_params': ['picard_vcf_output_params', 'picard_output_params']},
                          'HaplotypeCaller':
                              {'output_fmt': {},
                               'input_fmt': {},
                               'pre_tmpls': ['bam_index', 'pre_gatk_ints_chth', 'pre_gatk_excl_ints_chth'],
                               'post_tmpls': ['ref_opts', 'vcf_output_opts', 'gatk_bam_input', 'gatk_ints_chth', 'gatk_excl_ints_chth'],
                               'pre_params': ['ref_sel', 'gatk_req_params', 'picard_vcf_params', 'gatk_ints', 'gatk_excl_ints'],
                               'opt_params': [],
                               'adv_params': [],
                               'post_params': [],
                               'output_params': ['picard_vcf_output_params']},
                          'SelectVariants':
                              {'output_fmt': {},
                               'input_fmt': {},
                               'pre_tmpls': ['vcf_tabix', 'pre_gatk_ints_chth', 'pre_gatk_excl_ints_chth'],
                               'post_tmpls': ['ref_opts_opt', 'vcf_output_opts', 'vcf_input', 'gatk_ints_chth', 'gatk_excl_ints_chth'],
                               'pre_params': ['ref_sel', 'gatk_req_params', 'picard_vcf_params', 'gatk_ints',
                                              'gatk_excl_ints'],
                               'opt_params': [],
                               'adv_params': [],
                               'post_params': [],
                               'output_params': ['picard_vcf_output_params']}
                          }

        self.gen_out_fmt = {'activity_profile_out': 'tabular',
                            'assembly_region_out': 'tabular',
                            'graph_output': 'txt',
                            'bam_output': 'bam'}
        self.gen_in_fmt = {'alleles': 'vcf,vcf_bgzip',
                           'pedigree': 'tabular',
                           'population_callset': 'vcf,vcf_bgzip',
                           'contamination_fraction_per_sample_file': 'tabular',
                           'concordance': 'vcf,vcf_bgzip',
                           'comp': 'vcf,vcf_bgzip',
                           'discordance': 'vcf,vcf_bgzip',
                           'gatk_config_file': 'txt',
                           'dbsnp': 'vcf,vcf_bgzip'}

        self.param_tmpls = {'integer': ['name', 'argument', 'type', 'optional', 'value', 'min', 'max', 'label', 'help'],
                            'float': ['name', 'argument', 'type', 'optional', 'value', 'min', 'max', 'label', 'help'],
                            'text': ['name', 'argument', 'type', 'optional', 'value', 'label', 'help'],
                            'data': ['name', 'argument', 'type', 'optional', 'format', 'label', 'help'],
                            'select': ['name', 'argument', 'type', 'optional', 'label', 'help'],
                            'boolean': ['name', 'argument', 'type', 'truevalue', 'falsevalue', 'optional', 'checked', 'label', 'help'],
                            'output': ['format', 'name', 'label', 'help']}


class XmlTemplates(object):
    def __init__(self):
        self.xml_tmpl = Template('<expand macro="$macro_name" />')
        self.chth_tmpl = PercentTemplate('#include source=$%macro#')
        self.sel_out = Template('<option value="$value" selected="$selected">$value</option>')
        self.out_tmpl = Template('<data format="$format" name="$name" label="$${tool.name} on $${on_string}: $format" />')
        self.vcf_choose = PercentTemplate('#if $%section.%name:'
                                          '\n#if $%section.%name.is_of_type("vcf_bgzip"):'
                                          '\n%argument %name.vcf.gz'
                                          '\n#else'
                                          '\n%argument %name.vcf'
                                          '\n#end if'
                                          '\n#end if')
        self.vcf_tabix = PercentTemplate('#set datatype = $%section.%name.datatype'
                                         '\n#if $%section.%name:'
                                         '\n#if $%section.%name.is_of_type("vcf_bgzip"):'
                                         '\nln -s $%section.%name %name.vcf.gz &&'
                                         '\ntabix %name.vcf.gz &&'
                                         '\n#else'
                                         '\nln -s $%section.%name %name.vcf &&'
                                         '\n#end if'
                                         '\n#end if')
        self.file_chth = PercentTemplate('#if $output_opt.%out_sel_name:\n%argument $%name\n#end if')
        self.ext_arg = '#if $%section.%name:\n%argument $%section.%name\n#end if'
        self.reg_arg = '#if $%name:\n%argument $%name\n#end if'

class JsonXml(Mappings, XmlTemplates):
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
    def __init__(self, blob, tool_name):
        """
        :param blob:
        """
        Mappings.__init__(self)
        XmlTemplates.__init__(self)
        self.blob = blob
        self.tool_name = tool_name
        self.type = blob['type']
        self.section = blob['kind']
        self.pname = self.blob['name'].lstrip('-').replace('-', '_')
        self.is_output = (self.pname in self.tool_data[self.tool_name]['output_fmt']) or (self.pname in self.gen_out_fmt)
        self.is_input = (self.pname in self.tool_data[self.tool_name]['input_fmt']) or (self.pname in self.gen_in_fmt)
        self.is_input_vcf = self.is_input and self.gen_in_fmt[self.pname] == 'vcf,vcf_bgzip'
        if self.is_output:
            self.out_sel_name = self.pname + '_sel'
            self.out_sel_arg = self.blob['name'] + '_sel'
        self.xml_out = self.reblob()
        # Set to common is this argument is seen inside the known common arguments.
        self.common = self.xml_out['name'] in self.common_args
        #self.has_mcro_xml = self.xml_out['name'] in self.macro_xml
        #self.has_mcro_tmpl = self.xml_out['name'] in self.macro_tmpl
        # Check to see if the type, as defined in the json file, is recognized.
        if self.blob['type'] not in self.xml_json_type_map:
            print('Argument type %s not recognized, skipping.' % self.blob['type'])

        if self.blob['options']:
            self.is_select = True
            self.sel_blob = self.sel_prep()
        else:
            self.is_select = False
            self.sel_blob = None

        # Shouldn't need to call this twice, try to figure out how to synthesize this.
        if not self.common and self.section != 'common' and self.section != 'deprecated':
            if self.is_input_vcf:
                self.chth_pre = self.cheetah_template(pre=True)
            else:
                self.chth_pre = None
            self.chth = self.cheetah_template()
        else:
            self.chth = None
            self.chth_pre = None

#        self.xml_param = {label:self.xml_out[label] for label in self.param_tmpls[self.xml_out['type']]}

        if self.is_output:
            self.xml_param = {label:self.xml_out[label] for label in self.param_tmpls['output']}
        elif not self.common:
            self.xml_param = {label:self.xml_out[label] for label in self.param_tmpls[self.xml_out['type']]}
        else:
            self.xml_param = None

    def cheetah_template(self, pre=False):
        """
        Produce templates for how things looks in the Cheetah section.
        :return:
        """
        if self.is_output:
            xml_out = self.xml_out
            xml_out['out_sel_name'] = self.out_sel_name
            cht_tmpl = self.file_chth
            return cht_tmpl.substitute(self.xml_out)
        elif self.is_input and not pre:
            cht_tmpl = self.vcf_choose
            return cht_tmpl.substitute(self.xml_out)
        elif self.is_input and pre:
            cht_tmpl = self.vcf_tabix
            return cht_tmpl.substitute(self.xml_out)
        else:
            if self.xml_out['section'] not in ['required', 'common']:
                template_string = self.ext_arg
            else:
                template_string = self.reg_arg
            if self.xml_out['type'] == 'boolean':
                cht_tmpl = PercentTemplate(template_string.replace('%argument ', ''))
            else:
                cht_tmpl = PercentTemplate(template_string)
            return cht_tmpl.substitute(self.xml_out)

    def sel_prep(self):
        """
        Define the select blob, and then we can make a template out of the dict.
        <option value="normal_yes" selected="true">Yes</option>
        :return:
        """
        sel_blob = []
        for sel in self.blob['options']:
            if self.blob['defaultValue'] == sel['name']:
                sel_blob.append({'value': sel['name'], 'selected': 'true'})
            else:
                sel_blob.append({'value': sel['name'], 'selected': 'false'})

        return sel_blob

    def correct_sci(self, val):
        """
        If the defaultValue contains a scientific notation string, change it to something that looks like a float.
        In GATK json files, these are of the format "1.0E-6"
        Going to convert these instead of using xml validator since the smallest one is reasonable to include as float.
        :return:
        """
        try:
            dividend = float(val.split('E-')[0])
            divisor = float('1' + int(val.split('E-')[1]) * '0')
            return str('%f' % (dividend / divisor))
        except:
            return val

    def _value_correct(self, val):
        """
        Combine different operations necessary to provide a value that is acceptable for the wrapper.
        :return:
        """
        if ',' in val:
            return ''
        elif '[' in val and ']' in val:
            return val.replace('[', '').replace(']', '')
        elif val == 'null':
            return ''
        else:
            return self.correct_sci(val)

    def reblob(self):
        """
        Restructuring the information from the json input file to be used downstream.
        <param name="max_records_in_ram" type="integer" size="10" value="500000" min="1" label="Max Records in RAM" help="When writing files that need to be sorted, this will specify the number of records stored in RAM before spilling to disk. Increasing this number reduces the number of file handles needed to sort the file, and increases the amount of RAM needed." />
        :return:
        """
        xml_out = {'name': self.blob['name'].lstrip('-').replace('-', '_'),
                   'argument': self.blob['name'],
                   'type': self._assign_type(self.xml_json_type_map[self.type]),
                   'label': self.blob['name'].lstrip('-').replace('-', ' ').title(),
                   'optional': self.xml_json_req_map[self.blob['required']],
                   'value': self._value_correct(self.blob['defaultValue']),
                   'format': self.assign_format(),
                   'truevalue': self.blob['name'],
                   'falsevalue': '',
                   'checked': self.blob['defaultValue'],
                   'max': self.xml_json_num_map[self.blob['maxValue']],
                   'min': self.xml_json_num_map[self.blob['minValue']],
                   'help': self.blob['summary'],
                   'section': self.section}

        # if self.in_frmt:
        #     xml_out['format'] = self.in_frmt
        # if self.param_format_map(self.blob['name']):
        #     xml_out['format'] = self.param_format_map(self.blob['name'])
        for key in ['argument', 'checked', 'falsevalue', 'format', 'help', 'label', 'max', 'min', 'name', 'optional', 'truevalue', 'type', 'value']:
            xml_out[key] = escape(xml_out[key], {'"': '&quot;', "'": '&apos;'})
        return xml_out

    def _assign_type(self, type):
        """
        Make adjustments to the 'type' in reblob, for non-standard situations.
        :return:
        """
        if self.is_input:
            return 'data'
        else:
            return type

    def assign_format(self):
        """
        Assign the format where necessary, like if we are dealing with an output param.
        :return:
        """
        if self.is_output:
            if self.pname in self.tool_data[self.tool_name]['output_fmt']:
                return self.tool_data[self.tool_name]['output_fmt'][self.pname]
            elif self.pname in self.gen_out_fmt:
                return self.gen_out_fmt[self.pname]
        elif self.is_input:
            if self.pname in self.tool_data[self.tool_name]['input_fmt']:
                return self.tool_data[self.tool_name]['input_fmt'][self.pname]
            elif self.pname in self.gen_in_fmt:
                return self.gen_in_fmt[self.pname]
        else:
            # Not sure yet what this will be used for, but I think we need it.
            return ''


class JsonShell(object):
    """
    if picard, don't provide macro for annotations, not necessary
    summary == description
    description == help
    name == name
    """
    def __init__(self, args, profile='17.09'):
        self.profile = profile
        self.args = args
        self.pre_chth = []
        self.tool_chth = []
        self.tool_xml = []
        self.xml_opt = []
        self.xml_adv = []
        self.xml_out = []
        self.xml_comm = []
        self.sel_dict = {}
        with open(args.json, 'rU') as myfile:
            self.json_file = json.load(myfile)
            self.shell_dict = self.build_shell_dict()
            for entry in self.json_file['arguments']:
                self.my_xml = JsonXml(entry, self.shell_dict['short_name'])
                if self.my_xml.chth:
                    self.tool_chth.append(self.my_xml.chth)
                if self.my_xml.chth_pre:
                    self.pre_chth.append(self.my_xml.chth_pre)
                if self.my_xml.xml_param and not self.my_xml.is_output:
                    if self.my_xml.section == 'advanced':
                        self.xml_adv.append(self.my_xml.xml_param)
                    elif self.my_xml.section == 'optional':
                        self.xml_opt.append(self.my_xml.xml_param)
                    elif self.my_xml.section == 'common':
                        self.xml_comm.append(self.my_xml.xml_param)
                    elif self.my_xml.section == 'deprecated':
                        pass
                    else:
                        self.tool_xml.append(self.my_xml.xml_param)
                elif self.my_xml.is_output:
                    self.xml_out.append(self.my_xml.xml_param)
                else:
                    pass

                if self.my_xml.sel_blob:
                    self.sel_dict[self.my_xml.pname] = self.my_xml.sel_blob

    def build_shell_dict(self):
        """
        This will house all values the templates need.
        :return:
        """
        shell_dict = {'id': self.json_file['name'].lower().split(' ')[0],
                      'name': Template('GATK4 $name').substitute(self.json_file),
                      'short_name': self.json_file['name'].split(' ')[0],
                      'profile': self.profile,
                      'description': self.json_file['summary'].rstrip(' '),
                      'summary': pypandoc.convert_text(self.json_file['description'], 'rst', format='html')}
        return shell_dict

    def inputs_create(self):
        """
        Arrange the inputs section.
        :return:
        """
        inputs = []
        for macro in self.my_xml.tool_data[self.shell_dict['short_name']]['pre_params']:
            inputs.append(self.my_xml.xml_tmpl.substitute(macro_name=macro))
        for macro in self.my_xml.tool_data[self.shell_dict['short_name']]['post_params']:
            inputs.append(self.my_xml.xml_tmpl.substitute(macro_name=macro))

        return '\n'.join(inputs)


    def command_create(self):
        """
        Arrange the commands section.
        :return:
        """
        command = []
        for macro in self.my_xml.tool_data[self.shell_dict['short_name']]['pre_tmpls']:
            command.append(self.my_xml.chth_tmpl.substitute(macro=macro))
        command.extend(self.pre_chth)
        command.append(Template('@CMD_BEGIN@ $short_name').substitute(self.shell_dict))
        command.extend(self.tool_chth)
        for macro in self.my_xml.tool_data[self.shell_dict['short_name']]['post_tmpls']:
            command.append(self.my_xml.chth_tmpl.substitute(macro=macro))

        return '\n'.join(command)


    def xml_chth_expand(self, start_str, in_tmpl, *args, pre=False):
        """
        Based on the json_type, expand the template
        :return:
        """
        temp_str = start_str
        for macros in args:
            for entry in macros:
                if pre:
                    temp_str = in_tmpl.substitute(macro_name=entry)
                    temp_str += start_str
                else:
                    temp_str += in_tmpl.substitute(macro_name=entry)

        return Template(temp_str)


class XmlEtrees(JsonShell):
    """
    Hold all of the etrees we need for the structure.
    """
    def __init__(self, args):
        """
        Provide templates for the shell of the XML file.
        :return:
        """
        #        etree.write(stdout, xml_declaration=True, encoding='UTF-8')
        JsonShell.__init__(self, args)
        tool = etree.Element('tool', id='gatk4_' + self.shell_dict['id'], name=self.shell_dict['name'],
                             version="@WRAPPER_VERSION@0", profile=self.profile)
        description = etree.SubElement(tool, 'description')
        description.text = '- ' + self.shell_dict['description']
        macros = etree.SubElement(tool, 'macros')
        macros_imp = etree.SubElement(macros, 'import')
        macros_imp.text = 'macros.xml'
        exp_reqs = etree.SubElement(tool, 'expand', macro='requirements')
        exp_vers = etree.SubElement(tool, 'expand', macro='version_cmd')
        command = etree.SubElement(tool, 'command', detect_errors='exit_code')
        command.text = etree.CDATA(self.command_create())

        # INPUT section
        self.inputs = etree.SubElement(tool, 'inputs')
        for entry in self.my_xml.tool_data[self.shell_dict['short_name']]['pre_params']:
            etree.SubElement(self.inputs, 'expand', macro=entry)
        self.build_inputs(self.tool_xml, self.inputs, 'param')
        self.opt_sect = etree.SubElement(self.inputs, 'section', name='optional', title='Optional Parameters', expanded='False')
        for entry in self.my_xml.tool_data[self.shell_dict['short_name']]['opt_params']:
            etree.SubElement(self.opt_sect, 'expand', macro=entry)
        self.build_inputs(self.xml_opt, self.opt_sect, 'param')

        if self.xml_adv or self.my_xml.tool_data[self.shell_dict['short_name']]['adv_params']:
            self.adv_sect = etree.SubElement(self.inputs, 'section', name='advanced', title='Advanced Parameters', expanded='False')
            for entry in self.my_xml.tool_data[self.shell_dict['short_name']]['adv_params']:
                etree.SubElement(self.adv_sect, 'expand', macro=entry)
            self.build_inputs(self.xml_adv, self.adv_sect, 'param')
        for entry in self.my_xml.tool_data[self.shell_dict['short_name']]['post_params']:
            etree.SubElement(self.inputs, 'expand', macro=entry)

        # Temporary common section to work on macros.
        self.comm_sect = etree.SubElement(self.inputs, 'section', name='common', title='Common Parameters', expanded='False')
        self.build_inputs(self.xml_comm, self.comm_sect, 'param')


        # OUTPUT section
        if self.xml_out:
            self.output_sect = etree.SubElement(self.inputs, 'section', name='output_opt', title='Additional Output Parameters', expanded='False')
            self.build_inputs_out_sel(self.xml_out)
        self.outputs = etree.SubElement(tool, 'outputs')
        for entry in self.my_xml.tool_data[self.shell_dict['short_name']]['output_params']:
            etree.SubElement(self.outputs, 'expand', macro=entry)
        self.build_inputs(self.xml_out, self.outputs, 'data')

        tests = etree.SubElement(tool, 'tests')
        help = etree.SubElement(tool, 'help')
        help.text = etree.CDATA(self.shell_dict['summary'])
        citations = etree.SubElement(tool, 'citations')
        exp_cit = etree.SubElement(citations, 'expand', macro='citations')

        self.to_write = etree.tostring(tool, pretty_print=True, encoding="unicode")

    def build_inputs_out_sel(self, params):
        """
        Build extra options under the additional outputs section to provide booleans for choosing optional outputs.
        <param name="sites-only-vcf-output" argument="--sites-only-vcf-output" type="boolean" truevalue="--sites-only-vcf-output" falsevalue="" optional="true" checked="false" label="Sites Only Vcf Output" help="If true, don&amp;apos;t emit genotype fields when writing vcf file output."/>
        :return:
        """
        for param in params:
            new_name = param['name'] + '_sel'
            new_arg = '--' + new_name
            this_param = etree.SubElement(self.output_sect, 'param', name=new_name, argument=new_arg, type='boolean', truevalue=new_arg, falsevalue='', optional='true', checked='false', label=param['label'], help=param['help'])

    def create_output_loc(self):
        """
        Create the output file name for writing, based on input folder.
        :return:
        """
        self.output_name = [self.args.xml_out, 'gatk4_' + self.json_file['name'].lower().split(' ')[0] + '.xml']
        if not self.args.xml_out.endswith('/'):
            return '/'.join(self.output_name)
        else:
            return ''.join(self.output_name)

    def write_me(self):
        """
        Write to file.
        :return:
        """
        handle_out = open(self.create_output_loc(), 'w')
        handle_out.write(self.to_write)
        handle_out.close()

    def build_inputs(self, params, parent, elem):
        """

        :return:
        """
        for param in params:
#            if param['optional'] == 'true':
            this_param = etree.SubElement(parent, elem)
            for k, v in param.items():
                if (k == 'min') and (v == ''):
                    pass
                elif (k == 'max') and (v == ''):
                    pass
                else:
                    this_param.set(k, v)

            if param['name'] in self.sel_dict:
                self.build_sel_opt(this_param, self.sel_dict[param['name']])
#            yield this_param

            if elem == 'data':
                # <filter>not gzipped_output</filter>
                filt = etree.SubElement(this_param, 'filter')
                filt.text = 'output_opt[\'' + param['name'] + '_sel\']'


    def build_sel_opt(self, this_param, sel_blob):
        """
        Add select option statements if this is of select type.
        [{'summary': '', 'name': 'ALL'}, {'summary': '', 'name': 'OVERLAPPING_ONLY'}]
        :return:
        """
        for sel in sel_blob:
            this_sel = etree.SubElement(this_param, 'option', selected=sel['selected'], value=sel['value'])
            this_sel.text = sel['value']

class PercentTemplate(Template):
    delimiter = '%'


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
    myshell = XmlEtrees(args)
    myshell.write_me()

if __name__ == "__main__":
    main()
