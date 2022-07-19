import ast
import re

'''
Takes matrix expressions and makes an independent c++ file representing those expressions
Additionally, creates a makefile which can be used to compile the resulting c++ into an executable file

***USAGE***
MatExpr2Gmp.py <algorithm name> <input file> <cycle length>

See the Compiler class for expression syntax and formatting
'''

class Compiler:
	'''Converts matrix expressions to c++ gmp code
	
	Note: Name the matrix variable in the expressions a single capital letter, such as 'M'
	Matrix multiplication is given by '@'
	Entrywise multiplication is given by '*'
	The convert to diagonal matrix function is given by 'd()'
	Exponentiation shorthand: Variables of the form M# are shorthand for M@M@M...
	Hadamard product shorthand: Variables of the form M_# are shorthand for M*M*M...
	
	Example: M@d(M2)*M_3
	'''	
		
	@staticmethod
	def FromFile(fileName, algName):
		'''Compiles newline-delimited expressions into C++
		Note: ALL non-blank lines are read, except if the first two characters are "//", allowing comments'''
		
		exprList = []
		desc = ''
		with open(fileName, 'r') as f:
			for line in f:
				if line == '\n':
					continue
				if line[0:2] == '//':
					desc += line[2:]
					continue
				line = line.replace('\n', '')
				exprList.append(line)
		if desc and desc[-1] == '\n':
			desc = desc[:-1]
		return Compiler(exprList, algName, desc)
	
	def __init__(self, exprList, algName, algDesc):
		self.basisSize = len(exprList)
		if not exprList:
			print('Error: empty expression list')
		match = re.search('[A-Z]', exprList[0])
		self.varName = exprList[0][match.start()]
		self.numTempMatrices = 2
		self.traceTempUsed = False
		algDesc = algDesc and '/* ' + algDesc + ' */\n' or ''
		self.Code = ('{funcDesc}void {funcName}(const Matrix<mpz_t>& {0}, MPList<mpq_t>& terms){{' + '\n\tMatrix<mpz_t> _T1({0}.Rows, {0}.Columns);' + 'Matrix<mpz_t> _T2({0}.Rows, {0}.Columns);mpz_t term, _t1, _t2;mpz_inits(term, _t1, _t2, NULL);').format(self.varName, len(exprList), funcName=algName, funcDesc=algDesc)
		
		# handle shorthand
		powerPattern = re.compile('%s\\d+' % self.varName)
		powerSet = set()
		hadamardPattern = re.compile('%s_\\d+' % self.varName)
		hadamardSet = set()
		for expr in exprList:
			powerSet |= {match for match in powerPattern.findall(expr)}
			hadamardSet |= {match for match in hadamardPattern.findall(expr)}
		powers = sorted(list(map(lambda x: int(x[1:]), powerSet)))
		if powers:
			prev = self.varName
			self.Code += '\n\t'
			for i in range(2, powers[-1] + 1):
				if i in powers:
					storage = self.varName + str(i)
					self.Code += 'Matrix<mpz_t> {1}({0}.Rows, {0}.Columns);'.format(self.varName, storage)
				else:
					storage = '_T1' if i%2 == 0 else '_T2'
				self.Code += '{0}.RightMultiply({1}, {2});'.format(self.varName, storage, prev)
				prev = storage
		hadamards = sorted(list(map(lambda x: int(x[2:]), hadamardSet)))
		if hadamards:
			prev = self.varName
			self.Code += '\n\t'
			for i in range(2, hadamards[-1] + 1):
				if i in hadamards:
					storage = self.varName + '_' + str(i)
					self.Code += 'Matrix<mpz_t> {1}({0}.Rows, {0}.Columns);'.format(self.varName, storage)
				else:
					storage = '_T1' if i%2 == 0 else '_T2'
				self.Code += '{0}.MultiplyEntrywise({1}, {2});'.format(self.varName, storage, prev)
				prev = storage
		
		# compile expressions
		for i,expr in enumerate(exprList):
			self.availableTempMatrices = [('_T%d' % (n + 1)) for n in reversed(range(self.numTempMatrices))] # pop to reserve; push to release
			self.traceTempUsed = False
			self.nameStack = []
			self.Code += '\n\t// %s\n\t' % expr # adds comment of expression before code
			root = ast.parse(expr, mode='single')
			w = self.walk(root) # walk may add to self.Code (through addTempMatrix), so store returned code temporarily before appending to self.Code
			self.Code += w
			lastMatrix = self.nameStack.pop()
			if self.traceTempUsed: # lastMatrix in this case should be _t1 (a trace scalar)
				self.Code += 'mpq_set_z(terms[%d], _t1);' % (i,)
			else:
				self.Code += '%s.Trace(term);mpq_set_z(terms[%d], term);' % (lastMatrix, i)
		self.Code = self.Code.replace(';', ';\n\t') # indent and space properly
		self.Code += '\n}'
	
	def walk(self, node):
		code = ''
		if type(node) is ast.Name:
			self.nameStack.append(node.id)
		elif type(node) is ast.BinOp:
			code += self.walk(node.left)
			code += self.walk(node.right)
			code += self.evalBinOp(node.op)
		elif type(node) is ast.UnaryOp:
			print('No unary operators are currently supported. Attempted to use: ', node)
			# code += self.walk(node.operand)
			# evalUnaryOp
		elif type(node) is ast.Call:
			name = node.func.id
			code += self.walk(node.args[0])
			if name == 'd':
				code += self.diag()
			elif name == 'tr':
				code += self.trace()
			else:
				print('Unrecognized Call %s' % name)
		else: 
			for c in ast.iter_child_nodes(node):
				code += self.walk(c)
		return code
		
	def addTempMatrix(self): 
		self.numTempMatrices += 1
		self.availableTempMatrices.append('_T%d' % self.numTempMatrices)
		self.Code += 'Matrix<mpz_t> ' + self.availableTempMatrices[-1] + '({0}.Rows, {0}.Columns);'.format(self.varName)
		
	def tryReleaseTempMatrix(self, matName):
		if matName[0:2] == '_T' and matName[2].isdigit():
			self.availableTempMatrices.append(matName)
			
	def diag(self):
		if len(self.availableTempMatrices) == 0:
			self.addTempMatrix()
		temp = self.availableTempMatrices.pop()
		operand = self.nameStack.pop()
		self.nameStack.append(temp)
		self.tryReleaseTempMatrix(operand)
		return '%s.GetDiagonal(%s);' % (operand, temp)

	def trace(self):
		if self.traceTempUsed:
			temp = '_t2'
		else:
			temp = '_t1'
			self.traceTempUsed = True
		operand = self.nameStack.pop()
		self.nameStack.append(temp)
		self.tryReleaseTempMatrix(operand)
		return "%s.Trace(%s);" % (operand, temp)
	
	def evalBinOp(self, binop):
		if type(binop) is ast.Mult and self.nameStack[-1][:2] == '_t':
			# scalar trace multiplication
			r = self.nameStack.pop()
			return 'mpz_mul({0}, {0}, {1});'.format(self.nameStack[-1], r)
		t,l,r = self.binopSetup()
		if type(binop) is ast.Mult:
			# Note: scaling a matrix by a constant is not supported
			return self.multiplyEntrywise(t,l,r)
		elif type(binop) is ast.MatMult:
			return self.rightMultiply(t,l,r)
		else:
			print(binop, ' is not supported!')	
	
	def binopSetup(self):
		if len(self.availableTempMatrices) == 0:
			self.addTempMatrix()
		t = self.availableTempMatrices.pop()
		r = self.nameStack.pop()
		l = self.nameStack.pop()
		self.tryReleaseTempMatrix(r)
		self.tryReleaseTempMatrix(l)
		return t, l, r
	
	def rightMultiply(self, temp, left, right):
		self.nameStack.append(temp)
		return '{}.RightMultiply({}, {});'.format(left, temp, right)
		
	def multiplyEntrywise(self, temp, left, right):
		self.nameStack.append(temp)
		return '{}.MultiplyEntrywise({}, {});'.format(left, temp, right)

import sys
if len(sys.argv) != 5:
	print('Usage error: %s <algorithm name> <input file> <cycle length> <basis graphs file>' % sys.argv[0])
	sys.exit()
c = Compiler.FromFile(sys.argv[2], sys.argv[1])

try:
	f = open(sys.argv[1] + '.cpp', 'r')
	f.close()
	f = open('make_' + sys.argv[1], 'r')
	f.close()
except FileNotFoundError:
	pass
else:
	overwrite = input('%s already exists. Overwrite? y/n ' % sys.argv[1])
	if overwrite != 'y':
		sys.exit()

with open(sys.argv[1] + '.cpp', 'w') as o:
	o.write('''#include <gmp.h>
#include <cstdio>
#include "Matrix.h"
#include "MPList.h"
#include <vector>
#include <fstream>

using namespace std;

void {0}(const Matrix<mpz_t>& M, MPList<mpq_t>& terms);

int main(int argc, char* argv[]){{
	int subGraphSize = {1};
	int subGraphIndex = ; // if the subgraph is a basis graph, what index it is. -1 for none
	int dimension = {2}; // how many connected even subgraphs (matrix expressions) of size 'subGraphSize'
	Matrix<mpq_t> system = Matrix<mpq_t>(dimension, dimension + 1);
	MPList<mpq_t> terms = MPList<mpq_t>(dimension); 
	Matrix<mpz_t> cycle = Matrix<mpz_t>::FromCycle(subGraphSize);

	ifstream inFile("{3}");
	string line;
	for (int i = 0; i < dimension; i++) {{
		do {{
			getline(inFile, line);
		}} while (line.find_first_not_of(" ") == string::npos);
		Matrix<mpz_t> basisGraph = Matrix<mpz_t>::FromLineHollowSymmetric(line);
		Matrix<mpz_t> combined = cycle.Combine(basisGraph, false);
		printf("-----iteration %d-----\\n", i+1);

		{0}(combined, terms);
		terms.Print();
		for (int j = 0; j < dimension; j++){{
			mpq_set(system.Index(i, j), terms[j]);
		}}

		mpq_set_ui(system.Index(i, dimension), i == subGraphIndex ? 2 : 1, 1);
	}}
	inFile.close();

	printf("\\n\\n");
	system.Print();
	printf("\\n\\n");
	vector<int> pivots;
	Matrix<mpq_t>::RowEchelonForm(system, pivots);
	system.BackSubstitution(pivots);
	system.Print();

	return 0;
}}

'''.format(sys.argv[1], sys.argv[3], c.basisSize, sys.argv[4]))
	o.write(c.Code)

with open('make_' + sys.argv[1], 'w') as o:
	o.write(
'''CXXFLAGS = -Wall -std=c++14 -O3
LIBS = -lgmpxx -lgmp

{0}: {0}.o mpImpl.o
	$(CXX) {0}.o mpImpl.o $(CXXFLAGS) $(LIBS) -o $@

{0}.o: {0}.cpp Matrix.h MPList.h
	$(CXX) -c {0}.cpp $(CXXFLAGS)

mpImpl.o: mpImpl.h mpImpl.cpp
	make
'''.format(sys.argv[1]))
