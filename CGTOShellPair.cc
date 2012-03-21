#include "CGTOShellPair.hh"
#include "quillmath.hh"

CGTOShellPair::CGTOShellPair(const CGTOShell& shA, const CGTOShell& shB):
	_shA(shA), _shB(shB),
	_weights(shA.weights().replicate(1, shB.size()).rowwise()
		* shB.weights().transpose()),
	_widthsA(shA.widths().replicate(1, shB.size())),
	_widthsB(shB.widths().transpose().replicate(shA.size(), 1)),
	_widths_sum(_widthsA + _widthsB),
	_widths_red(_widthsA * _widthsB / _widths_sum),
	_gauss_red((-(shA.center()-shB.center()).squaredNorm() * _widths_red).qexp()
		/ _widths_sum),
	_hinv_widths(0.5*_widths_sum.inverse()),
	_lsum(shA.lsum() + shB.lsum())
{
	for (int i = 0; i < 3; ++i)
	{
		_P[i] = (widthsA()*shA.center(i) + widthsB()*shB.center(i))
			/ widthsSum();
	}
}
