# find descriptions of variants in text

import logging, marshal, zlib, copy, struct, random, sqlite3, types, pkg_resources, csv, io
from collections import defaultdict, namedtuple
from os.path import join

from amelie.dataloaders import gene_data

import os
import sys

logging.basicConfig(stream=sys.stdout, level=logging.DEBUG)
logger = logging.getLogger('root')


# re2 is often faster than re
# we fallback to re just in case
try:
    import re2 as re
except ImportError:
    import re


from .psl import Psl
from . import maxbio, pslMapBed
from ..tuples import VariantMention, DNAModification
from .pubSeqTables import threeToOneLower, threeToOne, oneToThree, aaToDna, dnaToAa
# from pycbio.hgdata.Psl import Psl
# import pslMapBed, pubAlg, maxbio, pubConf, maxCommon, pubKeyVal

# from pygr.seqdb import SequenceFileDB

from pyfaidx import Fasta

pointRe = re.compile(r'[.] (?=[A-Z]|$)')


def findBestSnippet(text, start, end, minPos, maxPos, isLeft=False):
    " get end or start pos of best snippet for (start, end) in range (minPos, maxPos)"
    textLen = len(text)
        
    # make sure that (min,max) stays within string boundaries
    # and does not go into (start,end)
    if isLeft:
       minPos = max(0, minPos)
       maxPos = max(maxPos, 0)
       dotPos = text.rfind(". ", minPos, maxPos)
    else:
       minPos = max(0, minPos)
       maxPos = min(maxPos, textLen)
       #dotPos = text.find(". ", minPos, maxPos)
       # better: attempt to eliminate cases like E. coli 
       subText = text[minPos:minPos+250]
       match = None
       for m in pointRe.finditer(subText):
           match = m
           break
       if match!=None:
             dotPos = minPos+match.start()
       else:
             dotPos = -1

    if dotPos==-1:
        if isLeft:
            dotPos = maxPos
            if dotPos==start:
                dotPos=minPos
        else:
            dotPos = minPos
            if dotPos==end:
                dotPos=maxPos
    elif isLeft:
        dotPos+=2
    else:
        dotPos+=1

    return dotPos


def getSnippet(text, start, end, minContext=0, maxContext=150):
    start = int(start)
    end = int(end)
    if start==end==0:
        return ""

    rightDotPos = findBestSnippet(text, start, end, end+minContext, end+maxContext, isLeft=False)
    leftDotPos = findBestSnippet(text, start, end, start-maxContext, start-minContext, isLeft=True)

    leftSnip = text[leftDotPos:start]
    mainSnip = text[start:end]
    rightSnip = text[end:rightDotPos]
    snippet = leftSnip+"<<<"+mainSnip+">>>"+rightSnip
    snippet = snippet.replace("\n", " ")
    snippet = snippet.replace("\t", " ")
    return snippet

regexes = None

# this setting can be changed to allow protein variants
# that require a change of two base pairs. By default, it
# is off to reduce false positives
# Johannes: going for high recall, not high precision
# allowTwoBpVariants = False

# ===== DATA TYPES ========
Mention = namedtuple("Mention", "patName,start,end")

""" A mapped variant is a type-range-sequence combination from a text, 
    can be located on none, one or multiple types of sequences
    All positions are 0-based
"""
VariantFields = [
    "mutType",  # sub, del or ins
    "seqType",  # prot, dna or dbSnp
    "seqId",  # protein or nucleotide accession
    "geneId",  # entrez gene id, if a gene was found nearby in text
    "start",  # original position in text
    "end",  # end position in text
    "origSeq",  # wild type seq, used for sub and del, also used to store rsId for dbSnp variants
    "mutSeq",  # mutated seq, used for sub and ins
    "origStr",  # variant name as it occurs in the original paper
    "offset",  # offset for splicing variants, e.g., -3 for c.1184-3A>T
    "firstAa",  # first amino acid in a notation like p.A100_F102del2
    "secondAa",  # second amino acid in a notation like p.A100_F102del2
    "length",  # length of inserted or deleted string such as in p.A100_F102del2 (has length 2)
    "docId",  # document in which variant occurs (usually a pubmed ID)
    "variantType",  # missense, nonframeshift, frameshift, stopgain, stoploss, splicing
    "ivsNumber",  # IVS splicing notation number
    "vcfPos",
    "vcfRef",
    "vcfAlt",
    "correctProba"
    ]

# A completely resolved mutation
mutFields = \
    (
    "chrom",  # chromosome
    "start",  # on chrom
    "end",  # on chrom
    "varId",  # a unique id
    "inDb",  # list of db names where was found
    "patType",  # the type of the patterns (sub, del, ins)
    "hgvsProt",  # hgvs on protein, can be multiple, separated with |
    "hgvsCoding",  # hgvs on cdna, can be multiple, separated with |
    "hgvsRna",  # hgvs on refseq, separated by "|"
    "comment",  # comment on how mapping was done
    "rsIds",  # possible rsIds, separated by "|", obtained by mapping from hgvsRna
    "protId",  # the protein ID that was used for the first mapping
    "texts",  # mutation match in text
    # "mutSupport",  # prot, dna, protDna
    # "mutCount",    # how often was this mutation mentioned?
    "rsIdsMentioned",  # mentioned dbSnp IDs that support any of the hgvsRna mutations
    "dbSnpStarts" ,  # mentioned dbSnp IDs in text, start positions
    "dbSnpEnds",  # mentioned dbSNP Ids in text, end positions

    "geneSymbol",  # symbol of gene
    "geneType",  # why was this gene selected (entrez, symNearby, symInTitle, symInAbstract)
    "entrezId",  # entrez ID of gene
    "geneStarts",  # start positions of gene mentions in document
    "geneEnds",  # end positions of gene mentions in document

    "seqType",  # the seqType of the patterns, dna or protein
    "mutPatNames",  # the names of the patterns that matched, separated by |
    "mutStarts",  # start positions of mutation pattern matches in document
    "mutEnds",  # end positions of mutation pattern matches in document
    "mutSnippets",  # the phrases around the mutation mentions, separated by "|"
    "geneSnippets",  # the phrases around the gene mentions, separated by "|"
    "dbSnpSnippets"  # mentioned dbSNP Ids in text, snippets
    )

# fields of the output file
MutRec = namedtuple("mutation_desc", mutFields)

# ======= GLOBALS ===============
# this can be used to shuffle all protein sequences before matching
# to get a random background estimate

# ===== FUNCTIONS TO INIT THE GLOBALS =================


def openIndexedPsls(mutDataDir, fileBaseName):
    " return a dict-like object that returns psls given a transcript ID "
    liftFname = join(mutDataDir, fileBaseName)
    logger.debug("Opening %s" % liftFname)
    pslDict = SqliteKvDb(liftFname)
    return pslDict


def parseEntrez(fname):
    """ parse a tab-sep table with headers and return one dict with entrez to refprots
    and another dict with entrez to symbol
    """
    entrez2Sym = dict()
    entrez2RefseqProts = dict()
    entrez2RefseqCodingSeqs = dict()

    with open(fname) as file:
        for row in csv.DictReader(file, delimiter='\t'):
            entrez2Sym[int(row['entrezId'])] = row['sym']
            if row['refseqProtIds'] == "":
                refProts = None
            else:
                refProts = row['refseqProtIds'].split(",")
                # assert(len(refProts)==len(refseqs))
            if row['refseqIds'] == "":
                refseqIds = None
            else:
                refseqIds = row['refseqIds'].split(",")

            entrez2RefseqProts[int(row['entrezId'])] = refProts
            entrez2RefseqCodingSeqs[int(row['entrezId'])] = refseqIds
    return entrez2Sym, entrez2RefseqProts, entrez2RefseqCodingSeqs


def refSeqToNumber(refSeqId):
    return int(refSeqId.split('.')[0].split('_')[1])


class SqliteKvDb(object):
    """ wrapper around sqlite to create an on-disk key/value database (or just keys)
        Uses batches when writing.
        On ramdisk, this can write 50k pairs / sec, tested on 40M uniprot pairs
        set onlyUnique if you know that the keys are unique.
    """
    def __init__(self, fname, singleProcess=False, newDb=False, tmpDir=None, onlyKey=False, compress=False, keyIsInt=False, eightBit=False, onlyUnique=False):
        self.onlyUnique = onlyUnique
        self.compress = compress
        self.batchMaxSize = 100000
        self.batch = []
        self.finalDbName = None
        self.onlyKey = onlyKey
        self.dbName = "%s.sqlite" % fname
        if newDb and os.path.isfile(self.dbName):
            os.remove(self.dbName)
        isolLevel = None
        self.singleProcess = singleProcess
        if singleProcess:
            isolLevel = "exclusive"
        if not os.path.isfile(self.dbName) and tmpDir!=None:
            # create a new temp db on ramdisk
            self.finalDbName = self.dbName
            #self.dbName = join(pubConf.getFastTempDir(), basename(self.dbName))
            self.dbName = join(tmpDir, os.path.basename(self.dbName))
            logging.debug("Creating new temp db on ramdisk %s" % self.dbName)
            if os.path.isfile(self.dbName):
                os.remove(self.dbName)
            maxCommon.delOnExit(self.dbName)  # make sure this is deleted on exit
            self.con = sqlite3.connect(self.dbName, isolation_level=isolLevel)
        else:
            try:
                self.con = sqlite3.connect(self.dbName)
            except sqlite3.OperationalError:
                logging.warn("Could not open %s" % self.dbName)

        logging.debug("Opening sqlite DB %s" % self.dbName)

        keyType = "TEXT"
        if keyIsInt:
            keyType = "INT"
        if onlyKey:
            self.con.execute("CREATE TABLE IF NOT EXISTS data (key %s PRIMARY KEY)" % keyType)
        else:
            self.con.execute("CREATE TABLE IF NOT EXISTS data (key %s PRIMARY KEY,value BLOB)" % keyType)
        self.con.commit()

        self.cur = self.con
        if singleProcess:
            self.cur.execute("PRAGMA synchronous=OFF") # recommended by
            self.cur.execute("PRAGMA count_changes=OFF") # http://blog.quibb.org/2010/08/fast-bulk-inserts-into-sqlite/
            self.cur.execute("PRAGMA cache_size=8000000") # http://web.utk.edu/~jplyon/sqlite/SQLite_optimization_FAQ.html
            self.cur.execute("PRAGMA journal_mode=OFF") # http://www.sqlite.org/pragma.html#pragma_journal_mode
            self.cur.execute("PRAGMA temp_store=memory")
            self.con.commit()

        if eightBit:
            self.con.text_factory = str

    def get(self, key, default=None):
        try:
            val = self[key]
        except KeyError:
            val = default
        return val

    def __contains__(self, key):
        row = self.con.execute("select key from data where key=?",(key,)).fetchone()
        return row!=None

    def dispName(self):
        " return a name for log messages "
        if self.finalDbName is not None:
            return "sqlite:"+self.finalDbName
        else:
            return "sqlite:"+self.dbName

    def __getitem__(self, key):
        row = self.con.execute("select value from data where key=?",(key,)).fetchone()
        if not row: raise KeyError
        value = row[0]
        if self.compress:
            value = zlib.decompress(value)
        return value.decode('utf8')

    def __setitem__(self, key, value):
        if self.compress:
            value = zlib.compress(value)
        self.batch.append( (key, sqlite3.Binary(value)) )
        if len(self.batch)>self.batchMaxSize:
            self.update(self.batch)
            self.batch = []

    def __delitem__(self, key):
        if self.con.execute("select key from data where key=?",(key,)).fetchone():
            self.con.execute("delete from data where key=?",(key,))
            self.con.commit()
        else:
             raise KeyError

    def update(self, keyValPairs):
        " add many key,val pairs at once "
        logging.debug("Writing %d key-val pairs to db" % len(keyValPairs))
        if self.onlyUnique:
            sql = "INSERT INTO data (key, value) VALUES (?,?)"
        else:
            sql = "INSERT OR REPLACE INTO data (key, value) VALUES (?,?)"
        self.cur.executemany(sql, keyValPairs)
        self.cur.commit()

    def keys(self):
        return [row[0] for row in self.con.execute("select key from data").fetchall()]

    def close(self):
        if len(self.batch)>0:
            self.update(self.batch)
        self.con.commit()
        self.con.close()
        if self.finalDbName!=None:
            logging.info("Copying %s to %s" % (self.dbName, self.finalDbName))
            shutil.copy(self.dbName, self.finalDbName)
            os.remove(self.dbName)


# ===== CLASSES =================
class SeqData(object):
    """ functions to get sequences and map between identifiers for entrez genes,
    uniprot, refseq, etc """

    def __init__(self, taxId, pubmunch_data_dir):
        " open db files, compile patterns, parse input as far as possible "
        mutDataDir = pubmunch_data_dir + '/variants'
        geneDataDir = pubmunch_data_dir + '/genes'
        if mutDataDir == None:
            return
        self.mutDataDir = mutDataDir
        self.entrez2sym, self.entrez2refprots, self.entrez2refseqs = parseEntrez(join(geneDataDir, "entrez.tab"))
        self.symToEntrez = None  # lazy loading

        # refseq sequences
        fname = join(mutDataDir, "seqs")
        logger.info("opening %s" % fname)
        seqs = SqliteKvDb(fname)
        self.seqs = seqs

        # refprot to refseqId
        # refseq to CDS Start
        fname = join(mutDataDir, "refseqInfo.tab")
        logger.debug("Reading %s" % fname)
        self.refProtToRefSeq = {}
        self.refSeqCds = {}
        with open(fname) as file:
            for row in csv.DictReader(file, delimiter='\t'):
                self.refProtToRefSeq[row['refProt']] = row['refSeq']
                self.refSeqCds[row['refSeq']] = int(row['cdsStart']) - 1  # NCBI is 1-based

        # refseq to genome
        self.pslCache = {}
        self.refGenePsls = openIndexedPsls(mutDataDir, "refGenePsls.9606")
        self.intronChrStartEndStrands = {}
        with open(join(mutDataDir, "refGeneIntrons.hg19.bed")) as file:
            for row in csv.reader(file, delimiter='\t'):
                chrom = str(row[0])
                start = int(row[1])
                end = int(row[2])
                names = [str(x) for x in row[3].split('_')]
                refSeqId = names[0] + '_' + names[1]
                # logger.info("Adding refSeqId %s" % refSeqId)
                strand = row[5]
                if refSeqId not in self.intronChrStartEndStrands:
                    self.intronChrStartEndStrands[refSeqId] = []
                self.intronChrStartEndStrands[refSeqId].append((chrom, start, end, strand))
        for refSeqId in self.intronChrStartEndStrands:
            l = self.intronChrStartEndStrands[refSeqId]
            l = sorted(l, key=lambda x: x[1])
            if l[0][3] == "+":
                pass
            elif l[0][3] == "-":
                l = l[::-1]
            else:
                assert False, l[0][3]
            self.intronChrStartEndStrands[refSeqId] = l
        self.exonStartEnds = {}
        self.exonStrands = {}
        with open(join(mutDataDir, "refGeneExons.hg19.bed")) as file:
            for row in csv.reader(file, delimiter='\t'):
                start = int(row[1])
                end = int(row[2])
                names = [str(x) for x in row[3].split('_')]
                refSeqId = names[0] + '_' + names[1]
                # logger.info("Adding refSeqId %s" % refSeqId)
                strand = row[5]
                if refSeqId not in self.exonStartEnds:
                    self.exonStartEnds[refSeqId] = []
                self.exonStartEnds[refSeqId].append((start, end))
                self.exonStrands[refSeqId] = strand
        for refSeqId in self.exonStartEnds.keys():
            strand = self.exonStrands[refSeqId]
            l = self.exonStartEnds[refSeqId]
            l = sorted(l, key=lambda x: x[0])
            if strand == "+":
                pass
            elif strand == "-":
                l = l[::-1]
            else:
                assert False, strand
            self.exonStartEnds[refSeqId] = l

        # dbsnp db
        fname = join(self.mutDataDir, "dbSnp.sqlite")
        self.snpDb = sqlite3.connect(fname)

        logger.info("Reading of data finished")

    def incrementRefSeqVersionIfNotFound(self, seqId):
        if seqId in self.seqs:
            return seqId
        seqId_split = seqId.split('.')
        baseId = seqId_split[0]
        versionNum = int(seqId_split[1])
        # heuristic ... the largest majority of refseq IDs have version numbers less than 10
        for versionNum in range(versionNum + 1, 10):
            modifiedSeqId = baseId + "." + str(versionNum)
            if modifiedSeqId in self.seqs:
                logger.debug("Taking seqId %s instead of %s" % (modifiedSeqId, seqId))
                return modifiedSeqId
        return None

    def getSeq(self, seqId):
        " get seq from db , cache results "
        logger.log(5, "Looking up sequence for id %s" % seqId)
        seqId = str(seqId)  # no unicode
        if seqId in self.seqs:
            return self.seqs[seqId]
        return None

    def entrezToProtDbIds(self, entrezGene):
        " return protein accessions (list) in otherDb for entrezGene "
        entrezGene = int(entrezGene)
        protIds = self.mapEntrezToRefseqProts(entrezGene)
        return protIds

    def entrezToCodingSeqDbIds(self, entrezGene):
        " return protein accessions (list) in otherDb for entrezGene "
        entrezGene = int(entrezGene)
        seqIds = self.mapEntrezToRefseqCodingSeq(entrezGene)
        return seqIds

    def entrezToSym(self, entrezGene):
        entrezGene = str(entrezGene)
        if "/" in entrezGene:
            logger.debug("Got multiple entrez genes %s. Using only first to get symbol." % entrezGene)
        entrezGene = entrezGene.split("/")[0]

        entrezGene = int(entrezGene)
        if entrezGene in self.entrez2sym:
            geneSym = self.entrez2sym[entrezGene]
            logger.debug("Entrez gene %s = symbol %s" % (entrezGene, geneSym))
            return geneSym
        else:
            return None

    def mapSymToEntrez(self, sym):
        " return a list of entrez IDs for given symbol "
        if self.symToEntrez == None:
            self.symToEntrez = defaultdict(list)
            for e, s in self.entrez2sym.iteritems():
                self.symToEntrez[s].append(e)
        entrezIds = self.symToEntrez.get(sym)
        return entrezIds

    def mapEntrezToRefseqProts(self, entrezGene):
        " map entrez gene to refseq prots like NP_xxx "
        if entrezGene not in self.entrez2refprots:
            logger.debug("gene %s is not valid or not in selected species" % str(entrezGene))
            return []

        protIds = self.entrez2refprots[entrezGene]
        if protIds is None:
            logger.debug("gene %s is a non-coding gene, no protein seq available")
            return []

        protIds = sorted(protIds, key=refSeqToNumber)

        logger.debug("Entrez gene %d is mapped to proteins %s" % \
            (entrezGene, ",".join(protIds)))
        return protIds

    def mapEntrezToRefseqCodingSeq(self, entrezGene):
        if entrezGene not in self.entrez2refseqs:
            logger.debug("gene %s is not valid or not in selected species" % str(entrezGene))
            return []

        seqIds = self.entrez2refseqs[entrezGene]
        if seqIds is None:
            logger.debug("gene %s is a non-coding gene, no coding seq available")
            return []

        # lowest sequence first ... helps with mapping variant type
        seqIds = sorted(seqIds, key=refSeqToNumber)

        logger.debug("Entrez gene %d is mapped to coding sequence %s" % \
            (entrezGene, ",".join(seqIds)))
        return seqIds

    def getCdsStart(self, refseqId):
        " return refseq CDS start position "
        if refseqId not in self.refSeqCds:
            return None
        cdsStart = self.refSeqCds[refseqId]
        return cdsStart

    def getRefSeqId(self, refProtId):
        " resolve refprot -> refseq using refseq data "
        refseqId = self.refProtToRefSeq.get(refProtId, None)
        return refseqId

    def getRefseqPsls(self, refseqId):
        """ return psl objects for regseq Id
            as UCSC refseq track doesn't support version numbers, we're stripping those on input
        """
        psls = getPsls(refseqId, self.pslCache, self.refGenePsls, stripVersion=True)
        return psls

    # end of class seqData


class VariantDescription(object):
    """ A variant description fully describes a variant
        It has at least a type (sub, del, etc), a start-end position on a
        (potentially unknown) sequence a tuple (origSeq, mutSeq) that describes
        the mutation e.g. ("R", "S"). The class can generate a descriptive name
        for the mutation like "p.R123S"

        It can optionally include a sequence ID, when the sequence ID was part of the
        mutation description in the text, e.g. the HGVS "NM_000925:p.R123S"

    >>> VariantDescription("sub", "prot", 10, 11, "R", "S")
    VariantDescription(mutType=u'sub',seqType=u'prot',seqId=u'None',geneId=u'',start=u'10',end=u'11',origSeq=u'R',mutSeq=u'S')
    """
    __slots__ = VariantFields

    def __init__(self, mutType, seqType, start, end, origSeq, mutSeq, docId,
                 seqId=None, geneId="", offset=0, origStr="", firstAa=None, secondAa=None,
                 length=None, variantType=None, ivsNumber=None,
                 vcfPos=None, vcfRef=None, vcfAlt=None, correctProba=0.5):
        self.mutType = mutType
        self.seqType = seqType
        self.seqId = seqId
        self.geneId = geneId
        self.start = int(start) if start is not None else None
        self.end = int(end)  if end is not None else None
        self.origSeq = origSeq
        self.mutSeq = mutSeq
        self.origStr = origStr
        self.offset = offset
        self.firstAa = firstAa
        self.secondAa = secondAa
        if length == None and self.origSeq:
            self.length = len(self.origSeq)
        else:
            self.length = length
        self.docId = docId
        self.variantType = variantType
        self.ivsNumber = ivsNumber
        self.vcfPos = vcfPos
        self.vcfRef = vcfRef
        self.vcfAlt = vcfAlt
        self.correctProba = correctProba

    def getName(self):
        # used to be HGVS text for variant; now just unique identifier
        return str(self)

    def asRow(self):
        row = []
        for i in self.__slots__:
            row.append(getattr(self, i))
        return row

    def __repr__(self):
        # return ",".join(self.asRow())
        parts = []
        for field in self.__slots__:
            parts.append(field + "=" + repr(getattr(self, field)))
        return "VariantDescription(%s)" % ",".join(parts)


class SeqVariantData(object):
    """ the full information about variant located on a sequence, with mentions from the text that support it
        This is the final output of this module, including all information about mapped variants and their genes.
    """
    __slots__ = mutFields

    def __init__(self, varId="", protVars=[], rnaVars=[], \
                 comment="", beds=[], entrezGene="", geneSym="", mentions=[],
                 text="", seqType="prot", patType="sub"):
        self.varId = varId
        self.inDb = ""
        self.patType = patType
        self.seqType = seqType
        self.chrom = ""
        self.start = ""
        self.end = ""
        self.geneSymbol = geneSym
        self.entrezId = entrezGene
        self.comment = comment
        self.protId = ""
        self.geneType = "entrez"
        self.geneStarts = ""
        self.geneEnds = ""
        self.geneSnippets = ""

        mutStarts, mutEnds, patNames, snippets, texts = mentionsFields(mentions, text)
        self.mutStarts = ",".join(mutStarts)
        self.mutEnds = ",".join(mutEnds)
        self.mutPatNames = "|".join(patNames)
        self.mutSnippets = "|".join(snippets)
        self.texts = "|".join(set(texts))

    def asRow(self, rawStr=False):
        row = []
        for i in self.__slots__:
            s = getattr(self, i)
            if rawStr:
                s = str(s)
            else:
                s = unicode(s)
            row.append(s)
        return row

    def __repr__(self):
        # return ",".join(self.asRow())
        parts = []
        for field in self.__slots__:
            parts.append(field + "=" + repr(unicode(getattr(self, field))))
        return "SeqVariantData(%s)" % ",".join(parts)

# ===== FUNCTIONS =================
# helper methods for SeqData


def getPsls(qId, cache, dbm, stripVersion=False):
    """ load psls from compressed dbm, create Psl objects, use a cache
    reverse complement if on negative strand
    """
    qId = str(qId)
    if stripVersion:
        qId = str(qId).split(".")[0]
    logger.debug("Getting mapping psl for %s" % qId)
    if qId in cache:
        psls = cache[qId]
    else:
        if not qId in dbm:
            logger.error("Could not find PSL for %s" % qId)
            return []
        pslLines = dbm[qId]
        psls = []
        for line in pslLines.split("\n"):
            psl = Psl(line.split("\t"))
            psls.append(psl)
        cache[qId] = psls
    logger.debug("Got mapping psl %s" % str(psls[0]))

    corrPsls = []
    for p in psls:
        if p.strand == "-":
            p2 = p.reverseComplement()
        else:
            p2 = p
        corrPsls.append(p2)

    return corrPsls


def makeMention(match, patName):
    start = match.start()
    end = match.start() + len(match.groups(1))
    mention = Mention(patName, start, end)
    return mention


def parseRegex(pubmunch_data_dir):
    """ parse and compile regexes to list (seqType, mutType, patName, pat) """
    # read regexes, translate placeholders to long form and compile
    replDict = {
    "sep"         : r"""(?:^|[:;\s\(\[\'\"/,])""",
    "fromPos"     : r'(?P<fromPos>[1-9][0-9]*)',
    "toPos"       : r'(?P<toPos>[1-9][0-9]*)',
    "pos"         : r'(?P<pos>[1-9][0-9]*)',
    "offset"         : r'(?P<offset>[1-9][0-9]*)',
    "fromPoss"     : r'(?P<fromPos>[1-9][0-9]+)',
    "toPoss"       : r'(?P<toPos>[1-9][0-9]+)',
    "poss"         : r'(?P<pos>[1-9][0-9]+)',
    "plusMinus"    : r'(?P<plusMinus>[+\^\-\u00FE\u2AF9\u2AFA])',  # plus often gets picked up as \xc3\xbe, whatever that is
    "arrow"        : r'(?P<arrow>[/r\>4\u2192\.!\u2B0E\u02DA\[])',  # u2192 is the right arrow in unicode; some people choose to write c.234G.A instead of c.234G>A, maybe their shift key is not working?
    "underscore"   : r'(?P<underscore>[_ ])',
    "origAaShort" : r'(?P<origAaShort>[CISQMNPKDTFAGHLRWVEYX])',
    "origAasShort" : r'(?P<origAasShort>[CISQMNPKDTFAGHLRWVEYX]+)',
    "origAaLong"  : r'(?P<origAaLong>(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR|TER|GLUTAMINE|GLUTAMIC ACID|LEUCINE|VALINE|ISOLEUCINE|LYSINE|ALANINE|GLYCINE|ASPARTATE|METHIONINE|THREONINE|HISTIDINE|ASPARTIC ACID|ARGININE|ASPARAGINE|TRYPTOPHAN|PROLINE|PHENYLALANINE|CYSTEINE|SERINE|GLUTAMATE|TYROSINE|STOP|X))',
    "origAasLong"  : r'(?P<origAasLong>(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR|TER|GLUTAMINE|GLUTAMIC ACID|LEUCINE|VALINE|ISOLEUCINE|LYSINE|ALANINE|GLYCINE|ASPARTATE|METHIONINE|THREONINE|HISTIDINE|ASPARTIC ACID|ARGININE|ASPARAGINE|TRYPTOPHAN|PROLINE|PHENYLALANINE|CYSTEINE|SERINE|GLUTAMATE|TYROSINE|STOP|X)+)',
    "mutAaShort"  : r'(?P<mutAaShort>[fCISQMNPKDTFAGHLRWVEYX*])',
    "mutAaLong"  : r'(?P<mutAaLong>(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR|TER|GLUTAMINE|GLUTAMIC ACID|LEUCINE|VALINE|ISOLEUCINE|LYSINE|ALANINE|GLYCINE|ASPARTATE|METHIONINE|THREONINE|HISTIDINE|ASPARTIC ACID|ARGININE|ASPARAGINE|TRYPTOPHAN|PROLINE|PHENYLALANINE|CYSTEINE|SERINE|GLUTAMATE|TYROSINE|STOP|X))',
    "firstAaShort"  : r'(?P<firstAaShort>[CISQMNPKDTFAGHLRWVEY])',
    "firstAaLong"  : r'(?P<firstAaLong>(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR|TER|GLUTAMINE|GLUTAMIC ACID|LEUCINE|VALINE|ISOLEUCINE|LYSINE|ALANINE|GLYCINE|ASPARTATE|METHIONINE|THREONINE|HISTIDINE|ASPARTIC ACID|ARGININE|ASPARAGINE|TRYPTOPHAN|PROLINE|PHENYLALANINE|CYSTEINE|SERINE|GLUTAMATE|TYROSINE))',
    "secondAaShort"  : r'(?P<secondAaShort>[CISQMNPKDTFAGHLRWVEY])',
    "secondAaLong"  : r'(?P<secondAaLong>(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR|TER|GLUTAMINE|GLUTAMIC ACID|LEUCINE|VALINE|ISOLEUCINE|LYSINE|ALANINE|GLYCINE|ASPARTATE|METHIONINE|THREONINE|HISTIDINE|ASPARTIC ACID|ARGININE|ASPARAGINE|TRYPTOPHAN|PROLINE|PHENYLALANINE|CYSTEINE|SERINE|GLUTAMATE|TYROSINE))',
    "mutAasShort"  : r'(?P<mutAasShort>[fCISQMNPKDTFAGHLRWVEYX*]+)',
    "mutAasLong"  : r'(?P<mutAasLong>(CYS|ILE|SER|GLN|MET|ASN|PRO|LYS|ASP|THR|PHE|ALA|GLY|HIS|LEU|ARG|TRP|VAL|GLU|TYR|TER|GLUTAMINE|GLUTAMIC ACID|LEUCINE|VALINE|ISOLEUCINE|LYSINE|ALANINE|GLYCINE|ASPARTATE|METHIONINE|THREONINE|HISTIDINE|ASPARTIC ACID|ARGININE|ASPARAGINE|TRYPTOPHAN|PROLINE|PHENYLALANINE|CYSTEINE|SERINE|GLUTAMATE|TYROSINE|STOP|X)+)',
    "dna"         : r'(?P<dna>[ACTG])',
    "dnas"         : r'(?P<dnas>[ACTG]+)',
    "noDnas"         : r'([ACTG]+)',
    "origDna"     : r'(?P<origDna>[ACTG])',
    "origDnas"     : r'(?P<origDnas>[ACTG]+)',
    "mutDna"      : r'(?P<mutDna>[ACTG])',
    "mutDnas"      : r'(?P<mutDnas>[ACTG]+)',
    "fs"          : r'(?P<fs>fs([\*\u2217\u204E\u26B9\u2731\u066D]?[0-9]*[*\u2217\u204E\u26B9\u2731\u066D]?))',
    "length"      : r'(?P<length>[1-9][0-9]*)',
    "space"       : r'([ ?]{1,3})',
    "ivsNumber"   : r'(?P<ivsNumber>([IV]+)|[0-9]+)'
    }
    regexList = []
    counts = defaultdict(int)
    with open(pubmunch_data_dir + '/variants/regex.txt') as regexData:
        csvReader = csv.DictReader(regexData, delimiter='\t')
        for row in csvReader:
            logger.info("Translating %s" % row['pat'])
            patName = row['patName']
            if row['seqType'].startswith('#'):
                continue
            if patName == "":
                patName = row['pat']
            patFull = '(?=(' + row['pat'].format(**replDict) + '))'
            logger.info("full pattern is %s" % patFull)
            flags = 0
            if "Long}" in row['pat']:
                flags = re.IGNORECASE
                logger.info("ignoring case for this pattern")
            logger.info("Full pattern: %s" % patFull)
            patComp = re.compile(patFull, flags=flags | re.UNICODE)
            regexList.append((row['seqType'], row['mutType'], patName, patComp))
            counts[(row['seqType'], row['mutType'])] += 1

    for regexType, count in counts.items():
            logger.info("regexType %s, found %d regexes" % (str(regexType), count))
    return regexList


def parsePlusMinus(plusMinus):
    if plusMinus == "+" or plusMinus == u"\u00FE" or plusMinus == u"\u2AF9":
        return 1
    elif plusMinus == "-" or plusMinus == "^" or plusMinus == u"\u2AFA":
        return -1
    else:
        assert False, str(plusMinus)


def parseMatchSplicing(match, patName, seqType, docId):
    # dna splicing        {sep}c\.{pos}{plusMinus}{offset}{origDna}>{mutDna}
    groups = match.groupdict()
    seqStart = int(groups["pos"])
    seqEnd = seqEnd = seqStart + 1
    plusMinus = parsePlusMinus(groups["plusMinus"])
    offset = int(groups["offset"])
    offset *= plusMinus
    logger.info("Match %s orig offset %d, plusMinus: %s" % (match.group(1), offset, plusMinus))
    logger.info("Match %s signed offset %d" % (match.group(1), offset))
    origSeq = groups["origDna"]
    mutSeq = groups["mutDna"]
    var = VariantDescription(mutType="splicing", seqType=seqType, start=seqStart, end=seqEnd, origSeq=origSeq.upper(),
                             mutSeq=mutSeq.upper(), origStr=match.group(1).strip(), offset=offset, docId=docId,
                             variantType="splicing", correctProba=0.5)
    return var


def convertIVSNumber(ivsNumber):
    if re.match(r'[IV]+', ivsNumber):
        if ivsNumber == "I":
            return 1
        elif ivsNumber == "II":
            return 2
        elif ivsNumber == "III":
            return 3
        elif ivsNumber == "IV":
            return 4
        elif ivsNumber == "V":
            return 5
        elif ivsNumber == "VI":
            return 6
        elif ivsNumber == "VII":
            return 7
        elif ivsNumber == "VIII":
            return 8
        elif ivsNumber == "IX":
            return 9
        else:
            # what is this BS ... do they think people want exercise in translating roman numerals
            logger.debug("too lazy to convert IVS roman numeral %s" % (ivsNumber))
            return None
    else:
        return int(ivsNumber)
    assert False


def parseMatchIVS(match, patName, seqType, docId):
    groups = match.groupdict()
    ivsNumber = convertIVSNumber(groups["ivsNumber"])
    if not ivsNumber:
        return None
    plusMinus = parsePlusMinus(groups["plusMinus"])
    offset = int(groups["offset"])
    offset *= plusMinus
    logger.info("Match %s orig offset %d, plusMinus: %s" % (match.group(1), offset, plusMinus))
    logger.info("Match %s signed offset %d" % (match.group(0), offset))
    origSeq = groups["origDna"]
    mutSeq = groups["mutDna"]
    var = VariantDescription(mutType="ivssub", seqType=seqType, start=None, end=None, origSeq=origSeq.upper(),
                             mutSeq=mutSeq.upper(), origStr=match.group(1).strip(), offset=offset, docId=docId,
                             variantType="splicing", ivsNumber=ivsNumber, correctProba=0.3)
    return var


def parseMatchFS(match, patName, seqType, docId):
    " given a regular expression match object, return mutation and mention objects "
    groups = match.groupdict()
    # grab long and short versions of amino acid
    correctProba = 0.4
    if "origAaShort" in groups:
        origSeq = groups["origAaShort"]
        correctProba = 0.9
    if "origAaLong" in groups:
        correctProba = 0.9
        origSeq = threeToOneLower[groups["origAaLong"].lower()]

    origSeq = origSeq.upper()
    pos = int(groups["pos"])
    # if isBlacklisted(origSeq, pos, mutSeq):
    #     return None
    seqStart = pos
    seqEnd = pos + 1

    variantType = "frameshift"

    var = VariantDescription(mutType="fs", seqType=seqType, start=seqStart, end=seqEnd, origSeq=origSeq.upper(),
                             mutSeq=None, origStr=match.group(1).strip(), length=None, docId=docId,
                             variantType=variantType, correctProba=correctProba)
    return var


def parseMatchSub(match, patName, seqType, docId):
    " given a regular expression match object, return mutation and mention objects "
    groups = match.groupdict()
    # grab long and short versions of amino acid
    correctProba = 0.4
    if "origAaShort" in groups and "mutAaShort" in groups:
        if groups["origAaShort"].isupper() and groups["mutAaShort"].islower() or \
           groups["origAaShort"].islower() and groups["mutAaShort"].isupper():
            # usually something like p.S234f , where f is actually part of "fs" (frameshift)
            return None
    if "origAaShort" in groups:
        origSeq = groups["origAaShort"]
        correctProba = 0.9
    if "origAaLong" in groups:
        origSeq = threeToOneLower[groups["origAaLong"].lower()]
        correctProba = 0.9

    if "mutAaShort" in groups:
        mutSeq = groups["mutAaShort"]
    if "mutAaLong" in groups:
        mutSeq = threeToOneLower[groups["mutAaLong"].lower()]

    if "origDna" in groups:
        origSeq = groups["origDna"]
        correctProba = 0.7
    if "mutDna" in groups:
        mutSeq = groups["mutDna"]

    if "origDna" in groups and "mutDna" in groups:
        if groups["origDna"].isupper() and groups["mutDna"].islower() or \
           groups["origDna"].islower() and groups["mutDna"].isupper():
            return None

    mutSeq = mutSeq.upper()
    origSeq = origSeq.upper()

    if "fromPos" in groups:
        pos = int(groups["fromPos"])
        seqStart = pos

    if "toPos" in groups:
        seqEnd = int(groups["toPos"])
    else:
        pos = int(groups["pos"])
        seqStart = pos
        seqEnd = pos + 1

    if "length" in groups:
        length = int(groups["length"])
    else:
        length = len(origSeq)

    variantType = None
    if "fs" in groups and groups["fs"]:
        variantType = "frameshift"
    elif origSeq and (("*" in origSeq.lower()) or ("x" in origSeq.lower())):
        variantType = "stoploss"
    elif mutSeq and (("*" in mutSeq.lower()) or ("x" in mutSeq.lower())):
        variantType = "stopgain"
    elif variantType is None and seqType == "prot":
        if origSeq == mutSeq:
            variantType = "synonymous"
        else:
            variantType = "missense"

    if variantType == "frameshift":
        var = VariantDescription(mutType="fs", seqType=seqType, start=seqStart, end=seqEnd, origSeq=origSeq.upper(),
                                 mutSeq=mutSeq.upper(), origStr=match.group(1).strip(), length=length, docId=docId,
                                 variantType=variantType, correctProba=correctProba)
    else:
        var = VariantDescription(mutType="sub", seqType=seqType, start=seqStart, end=seqEnd, origSeq=origSeq.upper(),
                                 mutSeq=mutSeq.upper(), origStr=match.group(1).strip(), length=length, docId=docId,
                                 variantType=variantType, correctProba=correctProba)
    return var


def parseMatchDel(match, patName, seqType, docId):
    " given a regular expression match object, return mutation and mention objects "
    groups = match.groupdict()
    logger.debug("Parsing match del %s into: %s" % (match.group(1), str(groups)))

    correctProba = 0.4
    haveStart = False
    if "fromPos" in groups:
        pos = int(groups["fromPos"])
        seqStart = pos
        haveStart = True
    firstAa = None
    if "firstAaShort" in groups:
        firstAa = groups["firstAaShort"].upper()
        correctProba = 0.9
    if "firstAaLong" in groups:
        firstAa = threeToOneLower[groups["firstAaLong"].lower()].upper()

    haveEnd = False
    if "toPos" in groups:
        seqEnd = int(groups["toPos"]) + 1
        haveEnd = True
    else:
        pos = int(groups["pos"])
        seqStart = pos
        seqEnd = pos + 1
    secondAa = None
    if "secondAaShort" in groups:
        secondAa = groups["secondAaShort"].upper()
        correctProba = 0.95
    if "secondAaLong" in groups:
        secondAa = threeToOneLower[groups["secondAaLong"].lower()].upper()
        correctProba = 0.95

    # these are heuristics ... we might still deal with a stoploss if origSeq is not given, for example ...
    origSeq = None
    if "origAasShort" in groups:
        origSeq = groups["origAasShort"]
        correctProba = 0.9
    if "origAasLong" in groups:
        origSeq = ''.join([threeToOneLower[groups["origAasLong"].lower()[i:3 + i]] for i in range(0, len(groups["origAasLong"]), 3)])
        correctProba = 0.9
    if "origDnas" in groups:
        origSeq = groups["origDnas"]
        correctProba = 0.7
    if "origDna" in groups:
        origSeq = groups["origDna"]
        correctProba = 0.7
    if "origAaShort" in groups:
        origSeq = groups["origAaShort"]
        correctProba = 0.9
    elif "origAaLong" in groups:
        origSeq = threeToOneLower[groups["origAaLong"].lower()]
        correctProba = 0.9
    if origSeq:
        origSeq = origSeq.upper()

    if "length" in groups and groups["length"]:
        length = int(groups["length"])
        if haveStart and haveEnd and (length != seqEnd - seqStart):
            logger.debug("Warning: length %d != seqStart %d - seqEnd %d; setting to seqEnd-seqStart (we might confuse reference numbers for length)" % (length, seqEnd, seqStart))
            length = seqEnd - seqStart  # toPos + 1 already up there ...
    elif origSeq:
        length = len(origSeq)
    else:
        length = seqEnd - seqStart

    variantType = None
    if "fs" in groups and groups["fs"]:
        variantType = "frameshift"
    elif firstAa and (("*" in firstAa.lower()) or ("x" in firstAa.lower())):
        variantType = "stoploss"
    elif secondAa and (("*" in secondAa.lower()) or ("x" in secondAa.lower())):
        variantType = "stoploss"
    elif origSeq and (("*" in origSeq.lower()) or ("x" in origSeq.lower())):
        variantType = "stoploss"
    elif seqType == "prot":
        variantType = "nonframeshift"
    elif seqType == "dna":
        if length % 3 == 0:
            variantType = "nonframeshift"
        else:
            variantType = "frameshift"

    var = VariantDescription(mutType="del", seqType=seqType, start=seqStart, end=seqEnd,
                             origSeq=origSeq, mutSeq=None, origStr=match.group(1).strip(), firstAa=firstAa,
                             secondAa=secondAa, length=length, docId=docId, variantType=variantType,
                             correctProba=correctProba)
    return var


def parseMatchIns(match, patName, seqType, docId):
    " given a regular expression match object, return mutation and mention objects "
    groups = match.groupdict()

    correctProba = 0.4
    if "fromPos" in groups and groups["fromPos"]:
        pos = int(groups["fromPos"])
        seqStart = pos
    firstAa = None
    if "firstAaShort" in groups:
        firstAa = groups["firstAaShort"]
        correctProba = 0.9
    if "firstAaLong" in groups:
        firstAa = threeToOneLower[groups["firstAaLong"].lower()]
        correctProba = 0.9

    if "toPos" in groups and groups["toPos"]:
        seqEnd = int(groups["toPos"]) + 1
    else:
        pos = int(groups["pos"])
        seqStart = pos
        seqEnd = pos + 2
    secondAa = None
    if "secondAaShort" in groups:
        secondAa = groups["secondAaShort"]
        correctProba = 0.95
    if "secondAaLong" in groups:
        secondAa = threeToOneLower[groups["secondAaLong"].lower()]
        correctProba = 0.95

    mutSeq = None
    if "mutAasShort" in groups:
        mutSeq = groups["mutAasShort"]
    if "mutAasLong" in groups:
        mutSeq = ''.join([threeToOneLower[groups["mutAasLong"].lower()[i:3 + i]] for i in range(0, len(groups["mutAasLong"]), 3)])
    if "dnas" in groups:
        mutSeq = groups["dnas"]
    logger.debug("match: %s" % match.group(1))
    logger.info(str([g for g in groups]))
    if mutSeq:
        mutSeq = mutSeq.upper()

    if "length" in groups:
        length = int(groups["length"])
    else:
        length = len(mutSeq)

    variantType = None
    if "fs" in groups and groups["fs"]:
        variantType = "frameshift"
    elif mutSeq and (("*" in mutSeq.lower()) or ("x" in mutSeq.lower())):
        variantType = "stopgain"
    elif seqType == "prot":
        variantType = "nonframeshift"
    elif seqType == "dna":
        if length % 3 == 0:
            variantType = "nonframeshift"
        else:
            variantType = "frameshift"

    var = VariantDescription(mutType="ins", seqType=seqType, start=seqStart, end=seqEnd,
                             origSeq=None, mutSeq=mutSeq, origStr=match.group(1).strip(),
                             firstAa=firstAa, secondAa=secondAa, length=length, docId=docId,
                             variantType=variantType, correctProba=correctProba)
    return var


def parseMatchDelIns(match, patName, seqType, docId):
    groups = match.groupdict()

    correctProba = 0.4
    if "fromPos" in groups:
        pos = int(groups["fromPos"])
        seqStart = pos
    firstAa = None
    if "firstAaShort" in groups:
        firstAa = groups["firstAaShort"]
        correctProba = 0.9
    if "firstAaLong" in groups:
        firstAa = threeToOneLower[groups["firstAaLong"].lower()]
        correctProba = 0.9

    if "toPos" in groups:
        seqEnd = int(groups["toPos"]) + 1
    else:
        pos = int(groups["pos"])
        seqStart = pos
        seqEnd = pos + 1
    secondAa = None
    if "secondAaShort" in groups:
        secondAa = groups["secondAaShort"]
        correctProba = 0.95
    if "secondAaLong" in groups:
        secondAa = threeToOneLower[groups["secondAaLong"].lower()]
        correctProba = 0.95

    mutSeq = None
    origSeq = None
    if "mutAasShort" in groups:
        mutSeq = groups["mutAasShort"]
    if "mutAasLong" in groups:
        mutSeq = ''.join([threeToOneLower[groups["mutAasLong"].lower()[i:3 + i]] for i in range(0, len(groups["mutAasLong"]), 3)])
    if "dnas" in groups:
        mutSeq = groups["dnas"]
    if "origDnas" in groups:
        origSeq = groups["origDnas"]
    if "origDna" in groups:
        origSeq = groups["origDna"]
    if "mutDnas" in groups:
        mutSeq = groups["mutDnas"]
    logger.debug("match: %s" % match.group(1))
    logger.info(str([g for g in groups]))
    if mutSeq:
        mutSeq = mutSeq.upper()

    if "length" in groups:
        length = int(groups["length"])
    else:
        if origSeq:
            length = abs(len(origSeq) - len(mutSeq))
        else:
            length = len(mutSeq)

    variantType = None
    if "fs" in groups and groups["fs"]:
        variantType = "frameshift"
    elif mutSeq and (("*" in mutSeq.lower()) or ("x" in mutSeq.lower())):
        variantType = "stopgain"
    elif seqType == "prot":
        variantType = "nonframeshift"
    elif seqType == "dna":
        if (seqEnd - seqStart + len(mutSeq)) % 3 == 0:
            variantType = "nonframeshift"
        else:
            variantType = "frameshift"

    var = VariantDescription(mutType="delins", seqType=seqType, start=seqStart, end=seqEnd,
                             origSeq=origSeq, mutSeq=mutSeq, origStr=match.group(1).strip(),
                             firstAa=firstAa, secondAa=secondAa, length=length, docId=docId,
                             variantType=variantType, correctProba=correctProba)
    return var


def parseMatchDup(match, patName, seqType, docId):
    groups = match.groupdict()

    correctProba = 0.4
    origSeq = None
    if "origDna" in groups:
        origSeq = groups["origDna"]
        correctProba = 0.9
    elif "origDnas" in groups:
        origSeq = groups["origDnas"]
        correctProba = 0.9

    if "pos" in groups:
        pos = int(groups["pos"])
        seqStart = pos
        seqEnd = seqStart + 1
        length = 1
    elif "fromPos" in groups:
        assert "toPos" in groups
        seqStart = int(groups["fromPos"])
        seqEnd = int(groups["toPos"]) + 1
        length = seqEnd - seqStart
    else:
        assert False

    if origSeq:
        origSeq = origSeq.upper()
        # duplication ... in the old Indian fashion:
        mutSeq = origSeq + origSeq
        # ! duplication ... in the Assyro-Babylonian fashion:
        # ! mutSeq = origSeq * 2
    else:
        mutSeq = None

    if origSeq:
        if len(origSeq) != length:
            # assume that we picked up a reference accidentally ...
            logger.debug(("Warning: length %d != length of origSeq %d; setting to seqEnd-seqStart " + \
                         "(we might confuse reference numbers for length)") % (length, len(origSeq)))
            length = len(origSeq)

    variantType = None
    if "fs" in groups and groups["fs"]:
        variantType = "frameshift"
    elif seqType == "dna" and length % 3 == 0:
        variantType = "nonframeshift"
    elif seqType == "dna" and length % 3 != 0:
        variantType = "frameshift"

    var = VariantDescription(mutType="dup", seqType=seqType, start=seqStart,
                             end=seqEnd, origSeq=origSeq, mutSeq=mutSeq,
                             origStr=match.group(1).strip(), length=length,
                             docId=docId, variantType=variantType, correctProba=correctProba)
    return var


def firstDiffNucl(str1, str2):
    """Return first pos and all letters where strings differ. Returns None if more than maxDiff chars are different"""
    assert(len(str1) == len(str2))
    if str1 == str2:
        return None
    diffCount = 0
    i = 0
    diffPos = []
    diffCh1 = []
    diffCh2 = []

    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffCount += 1
            diffCh1.append(ch1)
            diffCh2.append(ch2)
            diffPos.append(i)
        i += 1

    if diffCount == 1:
        return (diffPos[0], diffCh1[0], diffCh2[0])
    # return (diffPos[0], "".join(diffCh1), "".join(diffCh2))
    return None


def possibleDnaChanges(origAa, mutAa, origDna):
    """ figure out which nucleotides were possibly mutated by an amino acid change
    will only look for single-bp mutations
    returns list of: position of nucleic acid, original and new basepair
    >>> possibleDnaChanges("V", "V", "GTA")
    [(2, 'A', 'T'), (2, 'A', 'C'), (2, 'A', 'G')]
    >>> possibleDnaChanges("V", "I", "GTA")
    [(0, 'G', 'A')]
    >>> possibleDnaChanges("G", "G", "GGC")
    [(2, 'C', 'T'), (2, 'C', 'G'), (2, 'C', 'A')]
    """

    origDna = origDna.upper()
    ret = set()
    mutDnas = backTrans(mutAa)
    logger.debug("Looking for possible DNA change. Aa change %s -> %s, original dna %s" % (origAa, mutAa, origDna))
    for mutDna in mutDnas:
        diffTuple = firstDiffNucl(origDna, mutDna)
        if diffTuple != None:
            ret.add(diffTuple)
            logger.debug("found possible mutated DNA: %s" % (mutDna))

    return list(ret)


def backTrans(aa):
    """ back translate protein to all nucleotide strings
    Returns the back-translated nucleotide sequences for a protein and codon
    table combination.
    copied from http://www.biostars.org/p/3129/
    >>> protein = 'FVC'
    >>> len(backTrans(protein))
    16
    >>> backTrans('CD')
    ['TGTGAT', 'TGCGAT', 'TGTGAC', 'TGCGAC']
    """
    # create initial sequences == list of codons for the first amino acid
    sequences = [codon for codon in aaToDna[aa[0]]]
    for amino_acid in aa[1:]:
        # add each codon to each existing sequence replacing sequences
        # leaves (num_codons * num_sequences) for next amino acid
        to_extend = sequences
        sequences = []
        for codon in aaToDna[amino_acid]:
            for sequence in to_extend:
                sequence += codon
                sequences.append(sequence)
    return sequences


def translate(dna):
    " return the aa translation of a dna seq "
    aaSeq = []
    for i in range(0, len(dna), 3):
        codon = dna[i:i + 3].upper()
        if len(codon) != 3:
            # continue should be break here ...
            continue
        aa = dnaToAa[codon]
        aaSeq.append(aa)
    return "".join(aaSeq)


def mentionsFields(mentions, text):
    " convert the mention objects to something that fits into a tab-sep file "
    mutStarts = []
    mutEnds = []
    snippets = []
    patNames = []
    texts = []
    for m in mentions:
        mutStarts.append(str(m.start))
        mutEnds.append(str(m.end))
        snippets.append(getSnippet(text, m.start, m.end).replace("|", " "))
        patNames.append(m.patName)
        texts.append(text[m.start:m.end].strip("() -;,."))
    return mutStarts, mutEnds, patNames, snippets, texts


def rewriteToRefProt(variant, protId):
    " create new VariantDescriptions, one for each protein Id "
    varNew = copy.copy(variant)
    varNew.seqId = protId
    return varNew


def pslMapVariant(variant, psl):
    " map variant through psl on target, given query position, and create a new variant "
    maker = pslMapBed.PslMapBedMaker()
    maker.mapQuery(psl, variant.start, variant.end)
    bed = maker.getBed()
    if bed == None:
        return None

    varNew = copy.deepcopy(variant)
    varNew.seqId = bed[0]
    varNew.start = int(bed[1])
    varNew.end = int(bed[2])
    return varNew


class VariantExtractor:
    def __init__(self, pubmunch_data_dir):
        self.geneData = SeqData(9606, pubmunch_data_dir)
        self.regexes = parseRegex(pubmunch_data_dir)

        self.ensembl_to_entrez = gene_data.load_ensembl_to_entrez()
        self.hg19 = Fasta(join(pubmunch_data_dir, 'variants/hg19.fa'), rebuild=False)

    def extract(self, text, genes, ensembl_to_mentioned_refseq):
        variants = self.extract_raw_variants(text)
        entrez_with_ensembl, ensemble_to_source_text = self.get_entrez_with_ensembl(genes)

        all_mentions = set()

        resulting_variants = defaultdict(set)

        for variant, mentions in variants:
            print(variant, flush=True)
            for possible_variant in self.ground_variant(variant, mentions, entrez_with_ensembl, ensemble_to_source_text):
                # print("Possible variant: %s" % str(possible_variant))
                chrom, start, end, name_extra_info, score, strand, thick_start, thick_end, item_rgb, block_count, \
                block_sizes, block_starts, variant_text, vcf_pos, vcf_ref, vcf_alt, nm_info = possible_variant
                # print("Have vcf coords %s, %s, %s, %s" % (chrom, vcf_pos, vcf_ref, vcf_alt))

                nm_refseq, codon_num, region = nm_info.split(':')
                nm_refseq_without_dot = nm_refseq.split('.')[0]

                variant_text = re.sub(r'^\W+', '', variant_text)
                variant_text = re.sub(r'\W+$', '', variant_text)
                if 'hap' in chrom.lower() or 'random' in chrom.lower():
                    continue

                if vcf_pos is None or vcf_ref is None or vcf_alt is None:
                    continue

                for ensembl_id in name_extra_info['eids']:
                    if ensembl_id in ensembl_to_mentioned_refseq:
                        if nm_refseq_without_dot not in ensembl_to_mentioned_refseq[ensembl_id] and \
                                name_extra_info['grounded_refseq'] not in ensembl_to_mentioned_refseq[ensembl_id]:
                            continue
                    # if vcf_pos is None or vcf_pos == "None":
                    #     continue
                    # if vcf_ref is None or vcf_ref == "None":
                    #     continue
                    # if vcf_alt is None or vcf_alt == "None":
                    #     continue
                    resulting_variants[(ensembl_id, variant_text)].add(
                        DNAModification(chrom=chrom,
                                        bed_start=int(start),
                                        bed_end=int(end),
                                        gene_symbol=name_extra_info['gene_symbol'],
                                        entrez_id=name_extra_info['entrez_id'],
                                        checked=bool(name_extra_info['checked']),
                                        variant_type=name_extra_info['variant_type'],
                                        grounded_refseq=name_extra_info['grounded_refseq'],
                                        strand=strand,
                                        vcf_pos=int(vcf_pos) if vcf_pos != "None" else None,
                                        vcf_ref=vcf_ref if vcf_ref != "None" else None,
                                        vcf_alt=vcf_alt if vcf_alt != "None" else None,
                                        nm_refseq=nm_refseq,
                                        codon_num=int(codon_num),
                                        region=region))

        for (gene, variant_text), modifications in resulting_variants.items():
            if len(modifications) > 0:
                all_mentions.add(
                    VariantMention(
                        gene=gene,
                        variant_text=variant_text,
                        modifications=tuple(modifications),
                    )
                )

        return all_mentions

    def extract_raw_variants(self, text):
        """ put mutation mentions from document together into dicts indexed by normal form
            return dict of "prot"|"dna"|"dbSnp" -> list of (VariantDescription, list of Mention)
            uses global variable "regexes", see loadDb()

        >>> findVariantDescriptions("The R71G BRCA1 mutation is really a p.R71G mutation")
        {'prot': [(VariantDescription(mutType=u'sub',seqType=u'prot',seqId=u'None',geneId=u'',start=u'70',end=u'71',origSeq=u'R',mutSeq=u'G'), [Mention(patName=u'{sep}p\\\\.\\\\(?{origAaShort}{pos}{mutAaShort}{fs}', start=35, end=42), Mention(patName=u'{sep}{origAaShort}{pos}{mutAaShort}', start=3, end=8)])]}
        """
        varMentions = defaultdict(list)
        varDescObj = {}
        for seqType, mutType, patName, pat in self.regexes:
            for match in pat.finditer(text):
                print("Match: %s" % match.group(1), flush=True)
                logger.debug("Match: Pattern %s, text %s" % (patName, match.group(1)))
                if mutType == "sub":
                    variant = parseMatchSub(match, patName, seqType, docId=None)
                elif mutType == "del":
                    variant = parseMatchDel(match, patName, seqType, docId=None)
                elif mutType == "ins":
                    variant = parseMatchIns(match, patName, seqType, docId=None)
                elif mutType == "dup":
                    variant = parseMatchDup(match, patName, seqType, docId=None)
                elif mutType == "splicing":
                    variant = parseMatchSplicing(match, patName, seqType, docId=None)
                elif mutType == "ivssub":  # a different splicing notation
                    variant = parseMatchIVS(match, patName, seqType, docId=None)
                elif mutType == "delins":
                    variant = parseMatchDelIns(match, patName, seqType, docId=None)
                elif mutType == "fs":
                    variant = parseMatchFS(match, patName, seqType, docId=None)
                else:
                    logger.debug("Ignoring match %s; don't know how to handle" % str(match.group(1).encode('utf-8').strip()))
                    continue
                if variant == None:
                    continue

                mention = makeMention(match, patName)
                varDescObj[variant.getName()] = variant
                varMentions[variant.getName()].append(mention)
                debugSnip = getSnippet(text, mention.start, mention.end, maxContext=60)
                logger.debug("Found Variant: %s, snippet %s" % (str(variant), debugSnip))

        # convert to dict of "prot"|"dna"|"dbSnp" -> list (variant, mentions)
        variants = []

        for varName, mentions in varMentions.items():
            variant = varDescObj[varName]
            variants.append((variant, mentions))
        return variants

    def get_entrez_with_ensembl(self, genes_found):
        ensembl_to_source_text = defaultdict(set)

        entrez_to_ensembl = defaultdict(lambda: set())
        for eid in genes_found:
            if eid in self.ensembl_to_entrez:
                entrez_to_ensembl[self.ensembl_to_entrez[eid]].add(eid)

        return entrez_to_ensembl, ensembl_to_source_text

    def ground_variant(self, variant, mention, entrez_with_ensembl, ensembl_to_source_text):
        allBeds = []
        logger.debug("Grounding mutation %s onto genes %s" % (variant, str(entrez_with_ensembl)))
        for entrezGene, ensemblIds in entrez_with_ensembl.items():
            geneSym = self.geneData.entrezToSym(entrezGene)
            if not geneSym:
                logger.warn("No symbol for entrez gene %s. Skipping gene." % str(entrezGene))
                continue
            for (seqId, mappedVariant, checked) in self.checkVariantAgainstSequence(variant, entrezGene, geneSym):
                if not seqId:
                    logger.debug("Don't have a sequence ID for variant %s" % str(variant))
                    continue
                varId = {"eids": ensemblIds, "gene_symbol": geneSym,
                         "entrez_id": entrezGene, "checked": checked,
                         "variant_type": variant.variantType, "grounded_refseq": seqId}
                if variant.seqType == "prot":
                    logger.info("Rewriting variant %s to refProt %s" % (str(variant), str(seqId)))
                    protVar = rewriteToRefProt(variant, seqId)
                    logger.info("Got protVar %s, mapping to rnaVars" % str(protVar))
                    rnaVars = self.mapToRna(protVar)
                    logger.info("Got rnaVars %s, mapping to beds" % str(rnaVars))
                    beds = self.mapToGenome(rnaVars, varId)
                elif variant.seqType == "dna":
                    rnaVars = []
                    if mappedVariant:
                        assert variant.mutType == "splicing" or variant.mutType == "ivssub"
                        logger.info("Taking precomputed mappings for variant %s" % variant.origStr)
                        beds = [mappedVariant]
                        for bed in beds:
                            bed[3] = varId
                        logger.info("beds: %s" % str(beds))
                    else:
                        cdsStart = self.geneData.getCdsStart(seqId)
                        if cdsStart is None:
                            continue
                        rVariant = copy.copy(variant)
                        # Johannes: not sure why -1, but so be it ... these wacko off by one
                        # errors are not worth hunting down to the original source lol
                        rVariant.start = variant.start + cdsStart - 1
                        rVariant.end = variant.end + cdsStart - 1
                        rVariant.seqId = seqId
                        rnaVars.append(rVariant)
                        beds = self.mapToGenome(rnaVars, varId)
                else:
                    assert False, "can only ground prot and dna variants; have %s" % variant.seqType
                allBeds.extend(beds)
        return allBeds


    def checkVariantAgainstSequence(self, variant, entrezGene, sym):
        """ given a variant and a gene ID, 
        try to resolve gene to transcript sequence via  various protein databases 
        and check if they have matches for the wildtype aa at the right position 
        seqDbs can be any of "refseq", "oldRefseq", "uniprot", "genbank"
        - variant is a namedtuple with VariantFields defined above
        - entrezGene has to be a number as a string or a list of numbers separated by "/"
        - sym is only used for the logger system
        """
        entrezGene = str(entrezGene)
        for entrezGene in entrezGene.split("/"):
            entrezGene = int(entrezGene)
            logger.debug("Trying to ground %s to entrez gene %s / %s" % (str(variant), entrezGene, sym))
            if variant.seqType == "prot":
                seqIds = self.geneData.entrezToProtDbIds(entrezGene)
            elif variant.seqType == "dna":
                seqIds = self.geneData.entrezToCodingSeqDbIds(entrezGene)
            else:
                logger.debug("variant is neither DNA nor prot variant, but %s instead" % variant.seqType)
                continue
            if len(seqIds) == 0:
                continue
            seqIds = [self.geneData.incrementRefSeqVersionIfNotFound(x) for x in seqIds]
            seqIds = [x for x in seqIds if x is not None]
            seqsWithVariant = self.hasSeqAtPos(seqIds, variant)
            for seqId, mappedVariant, checked in seqsWithVariant:
                yield seqId, mappedVariant, checked

    def checkSeq(self, vStart, vEnd, seqType, seqId, expectedSeq, variant):
        seq = self.geneData.getSeq(seqId)

        logger.info("vStart: %d" % vStart)
        logger.info("vEnd: %d" % vEnd)
        if seq == None:
            # uniprot sometimes uses other species as support
            logger.debug("sequence %s is not human or not available" % seqId)
            return False
        logger.info("len(seq): %d" % len(seq))

        if not vEnd <= len(seq):
            logger.debug("sequence %s is too short (position %d, sequence length %d, sequence %s)" % (seqId, vEnd, len(seq), seq))
            return False

        if seqType == "dna":
            cdsStart = self.geneData.getCdsStart(seqId)
            if cdsStart is None:
                logger.debug("Didn't get cdsStart for sequence %s" % seqId)
                return False
        else:
            cdsStart = 0
        genomeSeq = seq[vStart + cdsStart:vEnd + cdsStart].upper()
        logger.debug('Have %s %s', repr(seq[0:70]), repr(seq[70:140]))
        if genomeSeq == expectedSeq:
            logger.debug("Seq match: Found %s at pos %d-%d in seq %s (surrounding sequence: %s)" % \
                (genomeSeq, vStart, vEnd, seqId, 
                 seq[max(0, vStart+cdsStart-15):vStart+cdsStart] + '<' + seq[vStart+cdsStart:vEnd+cdsStart] + '>' + \
                 seq[vEnd+cdsStart:min(len(seq), vEnd+cdsStart+15)]))
            if variant is not None and variant.seqType == "dna" and variant.variantType is None:
                startOfCodon = vStart
                inCodonOffset = 0
                while startOfCodon % 3 != 0:
                    startOfCodon -= 1
                    inCodonOffset += 1
                logger.debug("Start of codon: %d" % startOfCodon)
                endOfCodon = startOfCodon + 3
                logger.debug("End of codon: %d" % endOfCodon)
                codon = seq[cdsStart + startOfCodon:cdsStart + endOfCodon].upper()
                if len(codon) != 3:
                    return False
                logger.debug("Codon: %s" % codon)
                mutCodon = codon[0:inCodonOffset].upper() + variant.mutSeq.upper() + codon[inCodonOffset + 1:3].upper()
                logger.debug("Mutated codon: %s" % mutCodon)
                origAa = dnaToAa[codon]
                logger.debug("Original amino acid: %s" % origAa)
                mutAa = dnaToAa[mutCodon]
                logger.debug("Mutated amino acid: %s" % mutAa)
                if origAa == mutAa:
                    variant.variantType = "synonymous"
                elif origAa == "_":
                    variant.variantType = "stoploss"
                elif mutAa == "_":
                    variant.variantType = "stopgain"
                else:
                    variant.variantType = "missense"
                logger.debug("Inferred variant type: %s" % (variant.variantType))
            return True
        else:
            surroundingSeq = seq[vStart + cdsStart - 5:vEnd + 5 + cdsStart].upper()
            logger.debug("No seq match: Need %s, but found %s at pos %d-%d in seq %s (surrounding: %s) (vStart: %d, vEnd: %d, cdsStart: %d)" % \
                (expectedSeq, genomeSeq, vStart, vEnd, seqId, surroundingSeq, vStart, vEnd, cdsStart))
            return False

    def isSeqCorrect(self, seqId, variant):
        " check if wild type sequence in protein corresponds to mutation positions "
        if seqId.startswith("NR_"):
            logger.info("Skipping noncoding sequence ID %s" % seqId)
            return False, None, False
        if seqId.startswith("XM_"):
            logger.info("Skipping XM sequence ID %s" % seqId)
            return False, None, False
        if seqId.startswith("XR_"):
            logger.info("Skipping XR sequence ID %s" % seqId)
            return False, None, False
        seq = self.geneData.getSeq(seqId)
        if not seq:
            logger.info("Can't find seqId %s in my database" % seqId)
            return False, None, False

        if variant.mutType == "splicing":
            correct, mappedVariants = self.isSplicingSeqCorrect(seqId, variant)
            return correct, mappedVariants, True

        if variant.mutType == "ivssub":
            correct, mappedVariants = self.isIVSSubSeqCorrect(seqId, variant)
            return correct, mappedVariants, True

        if (variant.mutType == "ins" or variant.mutType == "del") and variant.firstAa and variant.secondAa:
            correct = self.checkSeq(vStart=variant.start - 1,
                            vEnd=variant.start,
                            seqType=variant.seqType,
                            seqId=seqId,
                            expectedSeq=variant.firstAa,
                            variant=variant) \
                    and \
                    self.checkSeq(vStart=variant.end - 1,
                             vEnd=variant.end,
                             seqType=variant.seqType,
                             seqId=seqId,
                             expectedSeq=variant.secondAa,
                             variant=variant)
            return correct, None, True

        if not variant.origSeq and not variant.firstAa and not variant.secondAa:
            return True, None, False

        correct = self.checkSeq(vStart=variant.start - 1,
                        vEnd=variant.end - 1,
                        seqType=variant.seqType,
                        seqId=seqId,
                        expectedSeq=variant.origSeq,
                        variant=variant)
        return correct, None, True

    def hasSeqAtPos(self, seqIds, variant):
        " check a list of IDs return those with a wild-type sequence at a position "
        if seqIds == None:
            return []
        rv = []
        for seqId in seqIds:
            logger.info("checking seqId %s for variant %s" % (seqId, str(variant)))
            seqCorrect, mappedVariants, c = self.isSeqCorrect(seqId, variant)
            logger.info("hasSeqAtPos got: %s" % (str((seqCorrect, mappedVariants))))
            if seqCorrect:
                if mappedVariants:
                    for mappedVariant in mappedVariants:
                        rv.append((seqId, mappedVariant, c))
                else:
                    rv.append((seqId, None, c))
        return rv


    def mapToRna(self, protVar):
        """ given ref protein positions and refseq proteinIds, try to figure out the nucleotide 
        changes on the refseq cdna sequence and add these to the variant object
        """
        rnaVars = []
        transId = self.geneData.getRefSeqId(protVar.seqId)
        if transId == None:
            logger.error("could not resolve refprot to refseq of protein %s. This is due to a difference between "
                    "UniProt and Refseq updates. Skipping this protein." % protVar.seqId)
            return []

        logger.info("Variant to map to coding and RNA: %s" % str(protVar))
        if protVar.mutType == "ins":
            origDnaSeq, cdnaStart, cdnaEnd = self.dnaAtCodingPos(transId, protVar.start - 1, \
                                                            protVar.end, protVar.origSeq)
            if origDnaSeq is None:
                return []
            rnaVar = VariantDescription(mutType=protVar.mutType, seqType="dna", start=cdnaStart, end=cdnaEnd, \
                                        origSeq="", mutSeq=None, seqId=transId, origStr=protVar.origStr,
                                        docId=protVar.docId)
            rnaVars.append(rnaVar)
        elif (protVar.mutType == "del" and protVar.mutSeq is None) or (protVar.mutType == "fs"):
            pos = protVar.start - 1
            origDnaSeq, cdnaStart, cdnaEnd = self.dnaAtCodingPos(transId, pos, \
                                                            pos + protVar.length, protVar.origSeq)
            if origDnaSeq is None:
                return []
            # the salomonian solution ... I don't exactly know where this deletion starts and ends,
            # at least some papers (25278557) do a poor job of telling you what's going on, so just settling
            # on this position for the deletion

            # also dumping frameshifts in here. These could be done more accurately, but that's maybe not even necessary
            if protVar.origSeq:
                cdnaEnd = cdnaStart + 3 * len(protVar.origSeq)
            else:
                cdnaEnd = cdnaStart + 3 * protVar.length
            rnaVar = VariantDescription(mutType=protVar.mutType, seqType="dna", start=cdnaStart, end=cdnaEnd, \
                origSeq=origDnaSeq.upper(), mutSeq=None, seqId=transId, origStr=protVar.origStr, docId=protVar.docId)
            rnaVars.append(rnaVar)
        else:
            pos = protVar.start - 1
            psls = self.geneData.getRefseqPsls(transId)
            strand = None
            if psls:
                strand = psls[0].getTStrand()
            origDnaSeq, cdnaStart, cdnaEnd = self.dnaAtCodingPos(transId, pos, \
                pos + len(protVar.origSeq), protVar.origSeq)
            if origDnaSeq is None:
                return []
            # mark the whole amino acid no matter what
            possChanges = possibleDnaChanges(protVar.origSeq, protVar.mutSeq, origDnaSeq)
            if possChanges:
                for relPos, oldNucl, newNucl in possChanges:
                    logger.info("Got relPos %d, oldNucl %s, newNucl %s, cdnaStart %d" % (relPos, oldNucl, newNucl, cdnaStart))
                    # cdnaNuclEnd = cdnaNuclStart + len(newNucl)
                    vcfPos = cdnaStart + relPos + 1
                    vcfRef = None
                    vcfAlt = None
                    if strand == "+":
                        vcfRef = oldNucl.upper()
                        vcfAlt = newNucl.upper()
                    elif strand == "-":
                        vcfRef = maxbio.revComp(oldNucl).upper()
                        vcfAlt = maxbio.revComp(newNucl).upper()
                    cdnaEnd = cdnaStart + len(origDnaSeq)
                    rnaVar = VariantDescription(mutType=protVar.mutType, seqType="dna", start=cdnaStart, end=cdnaEnd, \
                        origSeq=origDnaSeq.upper(), mutSeq=None, seqId=transId, origStr=protVar.origStr, docId=protVar.docId,
                        vcfPos=vcfPos, vcfRef=vcfRef, vcfAlt=vcfAlt)
                    rnaVars.append(rnaVar)
            else:
                cdnaEnd = cdnaStart + len(origDnaSeq)
                rnaVar = VariantDescription(mutType=protVar.mutType, seqType="dna", start=cdnaStart, end=cdnaEnd, \
                    origSeq=origDnaSeq.upper(), mutSeq=None, seqId=transId, origStr=protVar.origStr, docId=protVar.docId)
                rnaVars.append(rnaVar)

        return rnaVars

    def dnaAtCodingPos(self, refseqId, start, end, expectAa):
        """ 
        get nucleotide at CODING position in refseqId, check against expected aa
        also return positions on cdna
        """
        logger.debug("Paranoia check: making sure that codons from %d-%d in %s correspond to %s" %
            (start, end, refseqId, expectAa))
        cdsStart = self.geneData.getCdsStart(str(refseqId))
        nuclStart = cdsStart + (3 * start)
        nuclEnd = nuclStart + 3 * (end - start)
        cdnaSeq = self.geneData.getSeq(refseqId)
        nuclSeq = cdnaSeq[nuclStart:nuclEnd]
        foundAa = translate(nuclSeq)
        logger.debug("dnaAtCodingPos: cdsStart=%d, start=%d, end=%d, nuclStart=%d, nuclEnd=%d, nuclSeq=%s, foundAa=%s" % (cdsStart, start, end, nuclStart, nuclEnd, nuclSeq, foundAa))
        # logger.debug("CDS start is %d, nucl pos is %d, codon is %s" % (cdsStart, nuclStart, nuclSeq))
        if not (expectAa == None or foundAa == expectAa):
            return None, None, None
        return nuclSeq, nuclStart, nuclEnd


    def findValidPsls(self, variant):
        if variant.mutType == "splicing":
            logger.info("getting valid PSLs for splicing variant")
            return self.findValidSplicingPsls(variant.seqId, variant)
        return self.geneData.getRefseqPsls(variant.seqId)

    def mapToGenome(self, rnaVars, bedName):
        " map to genome from refseq, remove duplicate results, return as 12-tuple (=BED) "
        beds = []
        for rnaVar in rnaVars:
            logger.debug("Mapping rnaVar %s:%d-%d (offset %d) to genome" % (rnaVar.seqId, rnaVar.start, rnaVar.end, rnaVar.offset))
            # get psl
            psls = self.findValidPsls(rnaVar)
            cdsStart = self.geneData.getCdsStart(rnaVar.seqId)
            variantCodonStartNumber = int((rnaVar.start - cdsStart) / 3.0) + 1
            if len(psls) == 0:
                logger.warn("No mapping for %s, skipping variant" % rnaVar.seqId)
                continue
            maker = pslMapBed.PslMapBedMaker()
            for psl in psls:
                maker.clear()
                # map rna var through psl
                start = rnaVar.start
                end = rnaVar.end
                maker.mapQuery(psl, start, end)

                bed = maker.getBed(name=bedName)
                if bed == None:
                    logger.debug("found mapping psl but nothing was mapped")
                    continue
                # bed.append(re.sub('\s', '_', rnaVar.origStr))
                bed.append(rnaVar.origStr)
                assert rnaVar.offset == 0
                bedChrom = bed[0]
                bedStart = int(bed[1])
                bedEnd = int(bed[2])

                vcfPos = None
                vcfRef = None
                vcfAlt = None
                strand = bed[5]
                if rnaVar.seqType == "dna" and rnaVar.mutType == "sub":
                    vcfPos = str(bed[2])
                    if rnaVar.vcfPos is not None:
                        maker.clear()
                        maker.mapQuery(psl, rnaVar.vcfPos-1, rnaVar.vcfPos)
                        vcfBed = maker.getBed("vcf")
                        if vcfBed is None:
                            vcfPos = None
                        else:
                            logger.info("vcf bed: %s" % str(vcfBed))
                            logger.info("actual bed: %s" % str(bed))
                            vcfPos = str(int(vcfBed[1]) + 1)
                            # assert abs(int(vcfBed[1]) + 1 - int(bedStart)) <= 5
                    if rnaVar.vcfRef is not None or rnaVar.vcfAlt is not None:
                        vcfRef = str(rnaVar.vcfRef)
                        vcfAlt = str(rnaVar.vcfAlt)
                    else:
                        if rnaVar.origSeq and rnaVar.mutSeq:
                            # for weird reasons, we often get "+-" for the plus strand
                            if strand == "+":
                                vcfRef = rnaVar.origSeq.upper()
                                vcfAlt = rnaVar.mutSeq.upper()
                            elif strand == "-":
                                vcfRef = maxbio.revComp(rnaVar.origSeq).upper()
                                vcfAlt = maxbio.revComp(rnaVar.mutSeq).upper()
                elif rnaVar.seqType == "dna" and rnaVar.mutType == "del":
                    vcfPos = bedStart
                    vcfRef = str(self.hg19[bedChrom][bedStart - 1:bedEnd]).upper()
                    vcfAlt = str(self.hg19[bedChrom][bedStart - 1:bedStart]).upper()
                elif rnaVar.seqType == "dna" and rnaVar.mutType == "ins" and rnaVar.mutSeq:
                    if strand == "+":
                        vcfPos = bedStart
                        vcfRef = str(self.hg19[bedChrom][vcfPos - 1:vcfPos]).upper()
                        vcfAlt = vcfRef + rnaVar.mutSeq.upper()
                    elif strand == "-":
                        vcfPos = bedStart + 1
                        vcfRef = str(self.hg19[bedChrom][vcfPos - 1:vcfPos]).upper()
                        vcfAlt = vcfRef + maxbio.revComp(rnaVar.mutSeq).upper()
                exonIntron = self.findExonIntron(rnaVar.seqId, bedStart, bedEnd)
                bed.extend([str(vcfPos), str(vcfRef), str(vcfAlt),
                            rnaVar.seqId + ":" + str(variantCodonStartNumber) + ":" + exonIntron])
                logger.debug("Got bed: %s" % str(bed))
                beds.append(bed)
        return beds


    def findExonIntron(self, seqId, variantStart, variantEnd):
        rv_list = []
        seqId = seqId.split('.')[0]
        if seqId in self.geneData.intronChrStartEndStrands:
            for i, (chrom, start, end, strand) in enumerate(self.geneData.intronChrStartEndStrands[seqId]):
                if set(range(start, end)).intersection(set(range(variantStart, variantEnd))):
                    rv_list.append("i" + str(i + 1))
        if seqId in self.geneData.exonStartEnds:
            for i, (start, end) in enumerate(self.geneData.exonStartEnds[seqId]):
                if set(range(start, end)).intersection(set(range(variantStart, variantEnd))):
                    rv_list.append("e" + str(i + 1))
        return ",".join(rv_list)

    def isSplicingSeqCorrect(self, seqId, variant):
        if not seqId.startswith('NM_'):
            logger.info("Splicing variants can only handle NM_ sequences, not %s" % (seqId))
            return False, []
        validPsls = self.findValidSplicingPsls(seqId, variant)
        logger.info("Found no valid splicing psls for variant %s" % str(variant))
        cdsStart = self.geneData.getCdsStart(seqId)
        rv = []
        for psl in validPsls:
            maker = pslMapBed.PslMapBedMaker()
            maker.clear()
            variantCodonStartNumber = int(variant.start / 3.0) + 1
            maker.mapQuery(psl, variant.start + cdsStart - 1, variant.end + cdsStart - 1)
            bed = maker.getBed(name=variant.docId)
            if bed == None:
                logger.info("bed is None for psl %s (variant %s)" % (str(psl), variant.origStr))
                continue
            # bed.append(re.sub('\s', '_', variant.origStr))
            bed.append(variant.origStr)
            chrom = bed[0]
            start = int(bed[1])
            end = int(bed[2])
            if psl.getTStrand() == "+":
                start += variant.offset
                end += variant.offset
            elif psl.getTStrand() == "-":
                start -= variant.offset
                end -= variant.offset
            else:
                assert False

            logger.info("splicing variant mapped to %s\t%d\t%d (variant %s)" % (chrom, start, end, variant.origStr))
            genomic_nucleotide = str(self.hg19[chrom][start:end]).upper()
            vcfRef = genomic_nucleotide
            vcfAlt = None
            if psl.getTStrand() == "+":
                vcfAlt = variant.mutSeq.upper()
            elif psl.getTStrand() == "-":
                genomic_nucleotide = maxbio.revComp(genomic_nucleotide).upper()
                vcfAlt = maxbio.revComp(variant.mutSeq).upper()
            else:
                assert False
            logger.info("genomic nucleotide at this position: %s" % (genomic_nucleotide))
            if variant.origSeq.upper() == genomic_nucleotide.upper():
                logger.info("found splicing variant")
                bed[1] = str(start)
                bed[2] = str(end)
                bed[6] = str(start)
                bed[7] = str(end)
                bed.append(str(end))
                bed.append(vcfRef)
                bed.append(vcfAlt)
                exonIntron = self.findExonIntron(seqId, start, end)
                bed.append(seqId + ":" + str(variantCodonStartNumber) + ":" + exonIntron)
                rv.append(bed)
        if rv:
            return True, rv
        return False, rv

    def findValidSplicingPsls(self, seqId, variant):
        validPsls = []
        psls = self.geneData.getRefseqPsls(seqId)
        cdsStart = self.geneData.getCdsStart(seqId)
        if cdsStart is None:
            logger.info("Don't have cdsStart for sequence %s" % seqId)
            return []
        vStart = variant.start + cdsStart
        for psl in psls:
            found = False
            logger.info("psl for seqId %s: %s (variant %s)" % (seqId, str(psl), variant.origStr))
            for pslBlock in psl.blocks:
                logger.info("variant start: %d ; pslBlock start: %d , end: %d" % (vStart, pslBlock.qStart, pslBlock.qEnd))
                if variant.offset > 0:
                    if vStart == pslBlock.qEnd:
                        found = True
                        break
                elif variant.offset < 0:
                    # + 1 seems to be correct here
                    if vStart == pslBlock.qStart + 1:
                        found = True
                        break
            if found:
                validPsls.append(psl)
        logger.info("Number of valid psls for variant %s: %d" % (str(variant), len(validPsls)))
        return validPsls

    def isIVSSubSeqCorrect(self, seqId, variant):
        if not seqId.startswith('NM_'):
            logger.info("Splicing variants can only handle NM_ sequences, not %s" % (seqId))
            return False, []
        seqId_base = seqId.split('.')[0]
        if seqId_base not in self.geneData.intronChrStartEndStrands:
            logger.info("Don't have introns for sequence %s" % seqId_base)
            return False, []
        intronChrStartEndStrands = self.geneData.intronChrStartEndStrands[seqId_base]
        if len(intronChrStartEndStrands) <= variant.ivsNumber:
            logger.info("Sequence %s does not have enough introns" % seqId_base)
            return False, []
        rv = []
        chrom, intronStart, intronEnd, strand = intronChrStartEndStrands[variant.ivsNumber - 1]
        variantCodonStartNumber = self.ivsToVariantCodonStartNumber(seqId, variant.ivsNumber, strand)
        if variantCodonStartNumber is None:
            return False, rv
        if strand == "+":
            if variant.offset >= 0:
                variantStart = intronStart - 1 + variant.offset
                variantEnd = intronStart + variant.offset
            else:
                variantStart = intronEnd + variant.offset  # here offset < 0
                variantEnd = intronEnd + 1 + variant.offset
        elif strand == "-":
            if variant.offset >= 0:
                variantStart = intronEnd - variant.offset
                variantEnd = intronEnd - variant.offset + 1
            else:
                variantStart = intronStart - 1 - variant.offset  # here offset < 0
                variantEnd = intronStart - variant.offset
        else:
            assert False, strand
        logger.info("Setting IVS variant to %s:%d-%d for seq ID %s, orig %s" % (chrom, variantStart, variantEnd, seqId_base, variant.origStr))

        genomic_nucleotide = str(self.hg19[chrom][variantStart:variantEnd]).upper()
        vcfRef = genomic_nucleotide
        vcfAlt = None
        if strand == "+":
            vcfAlt = variant.mutSeq.upper()
        elif strand == "-":
            vcfAlt = maxbio.revComp(variant.mutSeq).upper()
            genomic_nucleotide = maxbio.revComp(genomic_nucleotide).upper()
        else:
            assert False, strand
        logger.info("genomic nucleotide at this position: %s" % (genomic_nucleotide))
        if variant.origSeq.upper() == genomic_nucleotide.upper():
            logger.info("found splicing variant")
            exonIntron = "i" + str(variant.ivsNumber)
            bed = [chrom, str(variantStart), str(variantEnd), None, \
                   str(variantEnd - variantStart), strand, str(variantStart), \
                   str(variantEnd), '0', '1', str(variantEnd - variantStart), '0', \
                   # re.sub('\s', '_', variant.origStr), \
                   variant.origStr,
                   str(variantEnd), vcfRef,
                   vcfAlt, seqId + ":" + str(variantCodonStartNumber) + ":" + exonIntron]
            rv.append(bed)
        if rv:
            return True, rv
        return False, rv

    def ivsToVariantCodonStartNumber(self, seqId, ivsNumber, strand):
        exonStartEnds = self.geneData.exonStartEnds[seqId.split('.')[0]][:ivsNumber]
        s = sum(abs(x[1] - x[0]) for x in exonStartEnds)
        cdsStart = self.geneData.getCdsStart(seqId)
        if cdsStart is None:
            return None
        s -= cdsStart
        return int(s / 3.0) + 1

if __name__ == "__main__":
    logger.basicConfig(level=logger.DEBUG)
    import doctest
    doctest.testmod()
