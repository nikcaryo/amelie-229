"""
Extraction of positions of identified gene, phenotypes, and variants.
This uses the PDFTextToPos Java program to do the actual extraction and
communicates with messages sent in line-orient JSON, as define below.

"""
import subprocess
import json
import os.path
import sys

class PositionExtractor:

    def __init__(self):

        _ROOT = os.path.abspath(os.path.dirname(__file__))
        path_to_jar = os.path.join(_ROOT, 'avada_data/PDFTextToPos.jar')

        self.process = subprocess.Popen(
            # if java problems, try the XX option:
            # args=["java", "-Xmx1g", "-XX:MaxMetaspaceSize=128m", "-jar", path_to_jar],
            args=["java", "-Dpdfbox.fontcache=/tmp", "-Xmx1g", "-jar", path_to_jar],
            stdin=subprocess.PIPE, 
            stdout=subprocess.PIPE,
            universal_newlines=True,
            encoding='utf8',
        )

    def extract(self, pdf_path, strings_to_extract):
        """
        Returns the positions of the requested strings.
        Returns None if the pdf could not be processed.

        Example return:
            {
                "IVS27CAGT":[
                    {
                        "x":109.99825,
                        "y":579.6157,
                        "found":true,
                        "page":6,
                        "orientation":0.0,
                        "text":"IVS27CAGT",
                        "wordIndex":3026,
                        "numStopwordsLeft":7,
                        "numStopwordsRight":7,
                        "numLetterLeft":30,
                        "numLetterRight":18
                    }
                ]
            }
        """
        
        json_request = json.dumps({'pdfPath': pdf_path, 'searchStrings': strings_to_extract})
        self.process.stdin.write(json_request + '\n')
        self.process.stdin.flush()

        result = self.process.stdout.readline()

        response = json.loads(result)

        if 'errorMessage' in response:
            raise ValueError('Got error parsing file {} {} {}'.format(pdf_path, strings_to_extract, response['errorMessage']))
        else:
            return response['results']

    def close(self):
        try:
            self.process.stdin.close()
            self.process.wait(timeout=10)
        except:
            # This is in theory not optimal, but the program doesn't need to do any cleanup, so we should be ok.
            # This might need to be revisited
            self.process.kill()
            raise Exception("Forced to kill PositionExtractor subprocess. Something is wrong")

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        self.close()
        return False
