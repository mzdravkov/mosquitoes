import logging

from Bio.SearchIO.BlatIO import BlatPslParser
from Bio.SearchIO import Hit, QueryResult
from Bio.SearchIO.BlatIO import _create_hsp

class SafeBlatParser(BlatPslParser):

    def _parse_qresult(self):
        """
        NOTE: this is the same as the original _parse_qresult method
        but modified to handle exceptions when parsing the HSP objects
        and log warnings instead of raising exceptions.

        Yield QueryResult objects (PRIVATE)."""
        # state values, determines what to do for each line
        state_EOF = 0
        state_QRES_NEW = 1
        state_QRES_SAME = 3
        state_HIT_NEW = 2
        state_HIT_SAME = 4
        # initial dummy values
        qres_state = None
        file_state = None
        cur_qid, cur_hid = None, None
        prev_qid, prev_hid = None, None
        cur, prev = None, None
        hit_list, hsp_list = [], []

        while True:
            # store previous line's parsed values for all lines after the first
            if cur is not None:
                prev = cur
                prev_qid = cur_qid
                prev_hid = cur_hid
            # only parse the result row if it's not EOF
            if self.line:
                try:
                    cur = self._parse_row()
                except ValueError as e:
                    logging.warning("Value error in parsing BLAT result in file %s, row '%s'., (%s)",
                                    self.handle.name, self.line, e)
                    continue
                cur_qid = cur["qname"]
                cur_hid = cur["tname"]
            else:
                file_state = state_EOF
                # mock values, since we have nothing to parse
                cur_qid, cur_hid = None, None

            # get the state of hit and qresult
            if prev_qid != cur_qid:
                qres_state = state_QRES_NEW
            else:
                qres_state = state_QRES_SAME
            # new hits are hits with different ids or hits in a new qresult
            if prev_hid != cur_hid or qres_state == state_QRES_NEW:
                hit_state = state_HIT_NEW
            else:
                hit_state = state_HIT_SAME

            if prev is not None:
                # create fragment and HSP and set their attributes
                try:
                    hsp = _create_hsp(prev_hid, prev_qid, prev)
                except AssertionError as e:
                    logging.warning("Assertion error in parsing BLAT result in file %s, query %s, target %s. (%s)",
                                    self.handle.name, cur_qid, cur_hid, e)
                    continue

                hsp_list.append(hsp)

                if hit_state == state_HIT_NEW:
                    # create Hit and set its attributes
                    hit = Hit(hsp_list)
                    hit.seq_len = prev["tsize"]
                    hit_list.append(hit)
                    hsp_list = []

                # create qresult and yield if we're at a new qresult or at EOF
                if qres_state == state_QRES_NEW or file_state == state_EOF:
                    qresult = QueryResult(id=prev_qid)
                    for hit in hit_list:
                        qresult.absorb(hit)
                    qresult.seq_len = prev["qsize"]
                    yield qresult
                    # if we're at EOF, break
                    if file_state == state_EOF:
                        break
                    hit_list = []

            self.line = self.handle.readline()
    