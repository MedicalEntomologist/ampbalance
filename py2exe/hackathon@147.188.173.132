import pysam

class DtwToSam(object):
    
    def __init__(self, rname = "" , outfile = "-", path_to_template_model_file = "", path_to_complement_model_file = ""):
        self.rname = rname
        self.path_to_template_model_file = path_to_template_model_file
        self.path_to_complement_model_file = path_to_complement_model_file
        self.outfile = outfile
        self.open()

    def open(self):
        header = { 'HD': {'VN': '1.0'},
            'SQ': [{'LN': 0, 'SN': self.rname}] , 
            'CO' : ['template_model:' + self.path_to_template_model_file,
                    'complement_model:' + self.path_to_complement_model_file ] 
            }        
        self.sam_file =  pysam.AlignmentFile(self.outfile, "wh", header=header)
        
    def close(self):
        self.sam_file.close()

    def convert(self, aa):
        self.sam_file.write(aa.aligned_segment)

class DtwAlignment(object):

    def __init__(self, event_align_list,filename = None, path_to_read = "", start_clipping = 0, end_clipping = 0):
        self.event_align_list = event_align_list
        self.path_to_read = path_to_read
        self.filename = filename
        self._cigar_string = ""
        self._cigar = []
        self.start_clipping = start_clipping
        self.end_clipping = end_clipping

    @property 
    def rname(self):
        return self.event_align_list[0].ref_name

    @property 
    def flag(self):
        if self.is_reversed:
            return 16
        else:
            return 0

    @property 
    def pos(self):
        return min([e.ref_pos for e in self.event_align_list]) + 1

    @property 
    def is_reversed(self):
        return self.event_align_list[0].is_reversed

    @property 
    def mapq(self):
        return 0
    
    def _generate_cigar_string(self):
        self._prepend_clipping()
        self._previous_ref_pos = -1
        self._previous_query_pos = -1
        self._current_op = "M"
        self._current_op_int = 0
        self._current_count = 0
        for ea in self.event_align_list:
            if self._ref_pos_has_moved(ea) and self._query_pos_has_moved(ea):
                ## M
                self._iterate_M()
            elif self._ref_pos_has_moved(ea) and not self._query_pos_has_moved(ea):
                ## D
                self._iterate_D()
            elif not self._ref_pos_has_moved(ea) and self._query_pos_has_moved(ea):
                ## I
                self._iterate_I()
            else:
                raise ValueError("I don't think repeated rows should occur")
            self._previous_query_pos = ea.query_pos 
            self._previous_ref_pos = ea.ref_pos
        self._write_current_op()
        self._append_clipping()

    def _ref_pos_has_moved(self, ea):
        return ea.ref_pos != self._previous_ref_pos
        
    def _query_pos_has_moved(self, ea):
        return ea.query_pos != self._previous_query_pos

    def _prepend_clipping(self):
        if self.start_clipping:
            self._cigar_string += "%i%s" % (self.start_clipping, "S")
            self._cigar.append((4 , self.start_clipping))

    def _append_clipping(self):
        if self.end_clipping:
            self._cigar_string += "%i%s" % (self.end_clipping, "S")
            self._cigar.append((4 , self.end_clipping))        

    def _iterate_M(self):
        if self._current_op == "M":
            self._current_count += 1
        else:
            self._write_current_op()
            self._current_op = "M"
            self._current_op_int = 0
            self._current_count = 1

    def _iterate_I(self):
        if self._current_op == "I":
            self._current_count += 1
        else:
            self._write_current_op()
            self._current_op = "I"
            self._current_op_int = 1
            self._current_count = 1

    def _iterate_D(self):
        if self._current_op == "D":
            self._current_count += 1
        else:
            self._write_current_op()
            self._current_op = "D"
            self._current_op_int = 2
            self._current_count = 1    

    def _write_current_op(self):
        self._cigar_string += "%i%s" % (self._current_count, self._current_op)
        self._cigar.append((self._current_op_int, self._current_count))

    @property 
    def cigar_string(self):
        if not self._cigar_string:
            self._generate_cigar_string()
            return self._cigar_string
        else:
            return self._cigar_string

    @property
    def cigar(self):
        if not self._cigar:
            self._generate_cigar_string()
            return self._cigar
        else:
            return self._cigar

    @property 
    def rnext(self):
        return "*"

    @property 
    def pnext(self):
        return 0

    @property 
    def tlen(self):
        return 0

    @property 
    def seq(self):
        return "*"

    @property 
    def qual(self):
        return "*" 

    @property   
    def is_increasing(self):
        return all(x<=y for x, y in zip(self.positions, self.positions[1:]))

    @property 
    def es(self):
        if self.is_increasing:
            return 1
        else:
            return -1

    @property 
    def positions(self):
        return [ae.ref_pos + 1 for ae in self.event_align_list]

    @property 
    def aligned_segment(self):
        s = pysam.AlignedSegment()
        s.query_name = str(self.filename)
        s.query_sequence= self.seq
        s.flag = self.flag
        s.reference_start = self.pos
        s.mapping_quality = self.mapq
        s.cigar = self.cigar
        s.next_reference_start = self.pnext
        s.template_length=self.tlen
        s.query_qualities = pysam.fromQualityString(self.qual)
        s.tags = (("ES", self.es),
          ("RP", self.path_to_read))
        return s        
                                                             
class AlignedEvent(object):

    def __init__(self, aligned_event):
        self.aligned_event = aligned_event
        self.query_pos = int(aligned_event[0])
        self.ref_pos = int(aligned_event[1])
        self.dist = float(aligned_event[2])
        self.ref_name = aligned_event[3]
        self.f_r = aligned_event[4]

    @property 
    def is_reversed(self):
        return self.f_r == "R" 
    