#! python
from ClassTest import *
from BALL import *

START_TEST('<class name>', '$Id: GlobalTypes_test.py,v 1.1.2.1 2005/11/04 08:07:08 oliver Exp $')

CHECK('Distance')
d = Distance()
RESULT()

CHECK('Handle')
h = Handle()
RESULT

CHECK('Index')
i = Index()
RESULT

CHECK('Size')
s = Size()
RESULT

CHECK('Time')
t = Time()
RESULT

CHECK('HashIndex')
h = HashIndex()
RESULT

CHECK('Position')
h = Position()
RESULT

CHECK('Real')
r = Real()
RESULT

CHECK('DoubleRead')
h = DoubleReal()
RESULT

CHECK('Property')
h = Property()
RESULT

CHECK('ErrorCode')
h = ErrorCode()
RESULT

CHECK('Byte')
h = Byte()
RESULT

CHECK('PointerSizeInt')
h = PointerSizeInt()
RESULT

CHECK('LongSize')
h = LongSize()
RESULT

CHECK('LongIndex')
h = LongIndex()
RESULT

END_TEST()
