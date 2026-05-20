#include "dephem/EphemerisRelease.h"

dph::EphemerisRelease::EphemerisRelease(const std::string& binaryFilePath)
{
	// ������������� ���������� ����������:
	clear();

	// ����������� ���� � �����:
	m_binaryFilePath = binaryFilePath;

	// �������� �����:
	m_binaryFileStream.open(m_binaryFilePath.c_str(), std::ios::binary);

	// ���� ������?
	bool isFileOpen = m_binaryFileStream.is_open();

	if (isFileOpen)
	{
		readAndPackData();

		if (isDataCorrect())
		{
			m_ready = true;
		}
		else
		{
			clear();
		}
	}
	else
	{
		m_binaryFilePath.clear();
	}
}

dph::EphemerisRelease::EphemerisRelease(const EphemerisRelease& other)
{
	if (other.m_ready)
	{
		copyHere(other);

		if (isDataCorrect())
		{
			m_ready = true;
		}
		else
		{
			m_ready = false;

			clear();
		}
	}
}

dph::EphemerisRelease& dph::EphemerisRelease::operator=(const EphemerisRelease& other)
{
	if (other.m_ready)
	{
		clear();
		copyHere(other);

		if (isDataCorrect())
		{
			m_ready = true;
		}
		else
		{
			m_ready = false;

			clear();
		}
	}

	return *this;
}

dph::EphemerisRelease::~EphemerisRelease()
{
	m_binaryFileStream.close();
}

void dph::EphemerisRelease::calculateBody(unsigned calculationResult,
	unsigned targetBody, unsigned centerBody, double JED, double* resultArray) const
{
	// ���������� �������� ����������:
	// -------------------------------
	//	- calculationResult:
	//		1 - �������� �������� ������-�������,
	//		2 - �������� �������� ������� ���������
	//		����������: ��������� �������� �� dph::Calculate.
	//
	//	- targetBody, centerBody:
	//		------------------------------------
	//		������	��������
	//		------------------------------------
	//		1		�������� 
	//		2		������
	//		3		�����
	//		4		����
	//		5		������
	//		6		������
	//		7		����
	//		8		������
	//		9		������
	//		10		����
	//		11		������
	//		12		��������� ��������� �������
	//		13		��������� ������� �����-����
	//		------------------------------------
	//		����������: ��������� �������� �� dph::Body.
	//
	//	- JED:
	//		JED ������ ������������ ����������: [m_startDate : m_endDate].
	//
	//	- resultArray:
	//		�� ������������ ��������� �����, ����� ����������� ������ ������� ��� 
	//		���������� ���������� ����������. �� ������ ���� ������� ����������.


	//������� ������������ ��� ������� ������:
	if (this->m_ready == false)
	{
		return;
	}
	else if (calculationResult > 1)
	{
		return;
	}
	else if (targetBody == 0 || centerBody == 0)
	{
		return;
	}
	else if (targetBody > 13 || centerBody > 13)
	{
		return;
	}
	else if (JED < m_startDate || JED > m_endDate)
	{
		return;
	}
	else if (resultArray == NULL)
	{
		return;
	}

	// ���������� ��������� ���������:
	unsigned componentsCount = calculationResult == Calculate::STATE ? 6 : 3;

	// ����� �������� ���������� � ����������� �� ���������� �������� � ������������ ����:
	if (targetBody == centerBody)
	{
		// ������ #1 : ������� ���� �������� �����������.
		//
		// ����������� �������� ������� ������.

		// ���������� ������� ������:
		std::memset(resultArray, 0, sizeof(double) * componentsCount);
	}
	else if (targetBody == Body::SSBARY || centerBody == Body::SSBARY)
	{
		// ������ #2: ������� ��� ����������� ����� �������� ��������� ��������� �������.
		//
		// ��� ��� ��� ������ calculateBase ��� ��� ���������� ������ ������������ ��������� ��,
		// �� ���������� ��������� ���������� ������� ����. � ������, ���� ������� ����� ��������
		// ��� ��������� ��, �� ������������ "����������" ������ ������� ����.

		// ������ ����, ��� �� �������� ����������� ��:
		unsigned notSSBARY = targetBody == Body::SSBARY ? centerBody : targetBody;

		// ����� ������ ���������� � ����������� �� ����:
		switch (notSSBARY)
		{
		case Body::EARTH: calculateBaseEarth(JED, calculationResult, resultArray);	break;
		case Body::MOON: calculateBaseMoon(JED, calculationResult, resultArray);		break;
		case Body::EMBARY: calculateBaseItem(2, JED, calculationResult, resultArray);	break;
		default: calculateBaseItem(notSSBARY - 1, JED, calculationResult, resultArray);
		}

		// ���� ��������� �� �������� ������� �����, �� ������������ "����������" ������:
		if (targetBody == Body::SSBARY)
		{
			for (unsigned i = 0; i < componentsCount; ++i)
			{
				resultArray[i] = -resultArray[i];
			}
		}
	}
	else if (targetBody * centerBody == 30 && targetBody + centerBody == 13)
	{
		// ������ #3 : ������� � ����������� ������ �������� ����� � ���� (��� ���� � �����).
		//
		// � ���� ������ ���������� �������� �������� ��������� ���� ������������ ����� (������� 
		// ������� #9 (�� ����). � ������, ���� ������� ����� �������� �����, �� ������������
		// "���������� ������".

		// ��������� ������-������� (��� ������� ���������) ���� ������������ �����:
		calculateBaseItem(9, JED, calculationResult, resultArray);

		// ���� ������� ����� �������� �����, �� ������������ "����������" ������.
		if (targetBody == Body::EARTH)
		{
			for (unsigned i = 0; i < componentsCount; ++i)
			{
				resultArray[i] = -resultArray[i];
			}
		}
	}
	else
	{
		// ������ #4 : ��� ��������� ���������� ���.
		//
		// ������� ����������� �������� ������������ ���� ������������ ���������� ��, 
		// ����� - ��������. ����������� �������� ������� ����� �������� ������������ ���� �
		// ��������. 

		// ������ ��� ������������ ����:
		double centerBodyArray[6];

		// ��� ��������:
		for (unsigned i = 0; i <= 1; ++i)
		{
			// ����������� ������� � ������� � ����������� �� ������ ��������.
			// i == 0 : ������ � ����������� �����.
			// i == 1 : ������ � ������� �����.
			unsigned currentBodyIndex = i == 0 ? centerBody : targetBody;
			double* currentArray = i == 0 ? centerBodyArray : resultArray;

			// ����� ������ ���������� � ����������� �� ����:
			switch (currentBodyIndex)
			{
			case Body::EARTH: calculateBaseEarth(JED, calculationResult, currentArray);	break;
			case Body::MOON: calculateBaseMoon(JED, calculationResult, currentArray);	break;
			case Body::EMBARY: calculateBaseItem(2, JED, calculationResult, currentArray);	break;
			default: calculateBaseItem(currentBodyIndex - 1, JED, calculationResult, currentArray);
			}
		}

		// ������� ����� �������� ������������ � �������� ����:
		for (unsigned i = 0; i < componentsCount; ++i)
		{
			resultArray[i] -= centerBodyArray[i];
		}
	}
}

void dph::EphemerisRelease::calculateOther(unsigned calculationResult,
	unsigned otherItem, double JED,
	double* resultArray) const
{
	// ���������� �������� ����������:
	// -------------------------------
	//	- calculationResult:
	//		1 - �������� ������������ ��������,
	//		2 - �������� ������������ �������� � ��� (��) ����������� ������� �������.
	//		����������: ��������� �������� �� dph::Calculate.
	//	
	//	- other Item:
	//		----------------------------------------------------------------		
	//		������	��������
	//		----------------------------------------------------------------
	//		14		������ ������� �� �������� � ���������� (������ IAU 1980)	
	//		15		�������� ������ ������
	//		16		������� �������� ������ ������
	//		17		TT - TDB (� ������ �����).
	//		----------------------------------------------------------------
	//		����������: ��������� �������� �� dph::Other.
	//
	//	- JED:
	//		JED ������ ������������ ����������: [m_startDate : m_endDate].
	//
	//	- resultArray:
	//		�� ������������ ��������� �����, ����� ����������� ������ ������� ��� 
	//		���������� ���������� ����������. �� ������ ���� ������� ����������.

	//������� ������������ ��� ������� ������:
	if (this->m_ready == false)
	{
		return;
	}
	else if (calculationResult > 1)
	{
		return;
	}
	else if (otherItem < 14 || otherItem > 17)
	{
		return;
	}
	else if (JED < m_startDate || JED > m_endDate)
	{
		return;
	}
	else if (resultArray == NULL)
	{
		return;
	}
	else
	{
		calculateBaseItem(otherItem - 3, JED, calculationResult, resultArray);
	}
}


bool dph::EphemerisRelease::isReady() const
{
	return m_ready;
}

double dph::EphemerisRelease::startDate() const
{
	return m_startDate;
}

double dph::EphemerisRelease::endDate() const
{
	return m_endDate;
}

uint32_t dph::EphemerisRelease::releaseIndex() const
{
	return m_releaseIndex;
}

const std::string& dph::EphemerisRelease::releaseLabel() const
{
	return m_releaseLabel;
}

double dph::EphemerisRelease::constant(const std::string& constantName) const
{
	if (m_ready == false)
	{
		return 0.0;
	}
	else if (constantName == "AU")
	{
		return m_au;
	}
	else if (constantName == "EMRAT")
	{
		return m_emrat;
	}
	else if (constantName == "DENUM")
	{
		return m_releaseIndex;
	}
	else
	{
		return m_constants.find(constantName)->second;
	}
}

std::string dph::EphemerisRelease::cutBackSpaces(const char* charArray, size_t arraySize)
{
	for (size_t i = arraySize - 1; i > 0; --i)
	{
		if (charArray[i] == ' ' && charArray[i - 1] != ' ')
		{
			return std::string(charArray, i);
		}
	}

	return std::string(charArray, arraySize);
}

void dph::EphemerisRelease::clear()
{
	m_ready = false;

	m_binaryFilePath.clear();
	m_binaryFileStream.close();

	m_releaseLabel.clear();
	m_releaseIndex = 0;
	m_startDate = 0.0;
	m_endDate = 0.0;
	m_blockTimeSpan = 0.0;
	std::memset(m_keys, 0, sizeof(m_keys));
	m_au = 0.0;
	m_emrat = 0.0;
	std::map<std::string, double>().swap(m_constants);	// SWAP TRICK

	m_blocksCount = 0;
	m_ncoeff = 0;
	m_dimensionFit = 0;
	m_blockSize_bytes = 0;

	std::vector<double>().swap(m_buffer);	// SWAP TRICK
	std::vector<double>(1).swap(m_poly);		// SWAP TRICK
	std::vector<double>(2).swap(m_dpoly);	// SWAP TRICK

	m_poly[0] = 1;
	m_dpoly[0] = 0;
	m_dpoly[1] = 1;
}

void dph::EphemerisRelease::copyHere(const EphemerisRelease& other)
{
	// ������������ �:
	//	- ����������� �����������.
	//	- �������� �����������.

	m_ready = other.m_ready;

	m_binaryFilePath = other.m_binaryFilePath;

	m_binaryFileStream.close();
	m_binaryFileStream.open(other.m_binaryFilePath.c_str(), std::ios::binary);

	m_releaseLabel = other.m_releaseLabel;
	m_releaseIndex = other.m_releaseIndex;
	m_startDate = other.m_startDate;
	m_endDate = other.m_endDate;
	m_blockTimeSpan = other.m_blockTimeSpan;
	std::memcpy(m_keys, other.m_keys, sizeof(m_keys));
	m_au = other.m_au;
	m_emrat = other.m_emrat;
	m_constants = other.m_constants;

	m_blocksCount = other.m_blocksCount;
	m_ncoeff = other.m_ncoeff;
	m_emrat2 = other.m_emrat2;
	m_dimensionFit = other.m_dimensionFit;
	m_blockSize_bytes = other.m_blockSize_bytes;

	m_buffer = other.m_buffer;
	m_poly = other.m_poly;
	m_dpoly = other.m_poly;
}

void dph::EphemerisRelease::readAndPackData()
{
	// ������� ��� ������ ���������� �� �����:
	char	releaseLabel_buffer[RLS_LABELS_COUNT][RLS_LABEL_SIZE];	// �����. ���. � �������.
	char	constantsNames_buffer[CCOUNT_MAX_NEW][CNAME_SIZE];		// ����� ��������.
	double	constantsValues_buffer[CCOUNT_MAX_NEW];					// �������� ��������.

	// ���������� �������� � ����� ��������:
	uint32_t constantsCount;
	// ------------------------------------- ������ ����� ------------------------------------- //

	m_binaryFileStream.seekg(0, std::ios::beg);
	m_binaryFileStream.read((char*)&releaseLabel_buffer, RLS_LABEL_SIZE * RLS_LABELS_COUNT);
	m_binaryFileStream.read((char*)&constantsNames_buffer, CNAME_SIZE * CCOUNT_MAX_OLD);
	m_binaryFileStream.read((char*)&m_startDate, 8);
	m_binaryFileStream.read((char*)&m_endDate, 8);
	m_binaryFileStream.read((char*)&m_blockTimeSpan, 8);
	m_binaryFileStream.read((char*)&constantsCount, 4);
	m_binaryFileStream.read((char*)&m_au, 8);
	m_binaryFileStream.read((char*)&m_emrat, 8);
	m_binaryFileStream.read((char*)&m_keys, (12 * 3) * 4);
	m_binaryFileStream.read((char*)&m_releaseIndex, 4);
	m_binaryFileStream.read((char*)&m_keys[12], (3) * 4);

	// ������ �������������� ��������:
	if (constantsCount > 400)
	{
		// ���������� �������������� ��������:
		size_t extraConstantsCount = constantsCount - CCOUNT_MAX_OLD;

		m_binaryFileStream.read((char*)&constantsNames_buffer[CCOUNT_MAX_OLD],
			extraConstantsCount * CNAME_SIZE);
	}

	// ������ �������������� ������:
	m_binaryFileStream.read((char*)&m_keys[13], (3 * 2) * 4);

	// ������� ncoeff (���������� ������������� � �����):
	m_ncoeff = 2;
	for (int i = 0; i < 15; ++i)
	{
		// ���������� ��������� ��� ���������� ��������:
		int comp = i == 11 ? 2 : i == 14 ? 1 : 3;
		m_ncoeff += comp * m_keys[i][1] * m_keys[i][2];
	}

	// ������� � ����� � ����������� � �� ������:	
	if (constantsCount <= CCOUNT_MAX_NEW)
	{
		m_binaryFileStream.seekg(m_ncoeff * 8, std::ios::beg);
		m_binaryFileStream.read((char*)&constantsValues_buffer, constantsCount * 8);
	}


	// -------------------- �������������� � �������� ��������� ���������� --------------------- // 

	// ������������ ����� ����� ���������� � �������:
	for (size_t i = 0; i < RLS_LABELS_COUNT; ++i)
	{
		m_releaseLabel += cutBackSpaces(releaseLabel_buffer[i], RLS_LABEL_SIZE);
		m_releaseLabel += '\n';
	}
	m_releaseLabel;

	// ���������� ���������� m_constants ������� � ���������� ��������:
	if (constantsCount > 0 && constantsCount <= CCOUNT_MAX_NEW)
	{
		for (uint32_t i = 0; i < constantsCount; ++i)
		{
			std::string constantName = cutBackSpaces(constantsNames_buffer[i], CNAME_SIZE);
			m_constants[constantName] = constantsValues_buffer[i];
		}
	}

	// �������������� ����������:
	additionalCalculations();
}

void dph::EphemerisRelease::additionalCalculations()
{
	// ����������� ���. ������������� ��� ������ � �����������:
	m_emrat2 = 1 / (1 + m_emrat);
	m_dimensionFit = 1 / (43200 * m_blockTimeSpan);

	// ����������� ���������� ������ � ����������:
	m_blocksCount = size_t((m_endDate - m_startDate) / m_blockTimeSpan);

	// ������� ������������� ���������� ��������� � �������:
	size_t maxPolynomsCount = 0;
	for (int i = 0; i < 15; ++i)
	{
		if (m_keys[i][1] > maxPolynomsCount)
		{
			maxPolynomsCount = m_keys[i][1];
		}
	}

	// ����������� ������� ����� � ������:
	m_blockSize_bytes = m_ncoeff * sizeof(double);

	// �������������� ������ � ��������:
	m_buffer.resize(m_ncoeff);
	m_poly.resize(maxPolynomsCount);
	m_dpoly.resize(maxPolynomsCount);
}

bool dph::EphemerisRelease::isDataCorrect() const
{
	// � ������ ������ ����������� ������ �� ���������, �������
	// ����� �������� ��������������� �� ���������� �������� ���������, 
	// ���������� � ������� ��������.	

	if (m_binaryFileStream.is_open() == false)			return false;	// ������ �������� �����.
	if (m_startDate >= m_endDate)						return false;
	if (m_blockTimeSpan == 0)							return false;
	if ((m_endDate - m_startDate) < m_blockTimeSpan)	return false;
	if (m_emrat == 0)									return false;
	if (m_ncoeff == 0)									return false;

	if (check_blocksDates() == false)					return false;

	return true;
}

bool dph::EphemerisRelease::check_blocksDates() const
{
	// ����� ������� ����� � �������������� � �����:
	size_t firstBlockAdress = m_blockSize_bytes * 2;

	// ������� � ������� �����:
	m_binaryFileStream.seekg(firstBlockAdress, std::ios::beg);

	// �������� ����� ������� ����� ������ ���� ������ �������������:
	size_t subBlockOffset = (m_ncoeff - 2) * sizeof(double);

	for (size_t blockIndex = 0; blockIndex < m_blocksCount; ++blockIndex)
	{
		// ������ ��� ������ ������ ���� ������������� �� �������� �����:
		double blockDates[2] = { 0.0, 0.0 };

		// ������:
		m_binaryFileStream.read((char*)&blockDates, sizeof(blockDates));

		// ��������, ������� ������ ����:
		double blockStartDate = m_startDate + blockIndex * m_blockTimeSpan;
		double blockEndDate = blockStartDate + m_blockTimeSpan;

		if (blockDates[0] != blockStartDate || blockDates[1] != blockEndDate)
		{
			return false;
		}

		// ������� � ���������� �����:
		m_binaryFileStream.seekg(subBlockOffset, std::ios::cur);
	}

	return true;
}

void dph::EphemerisRelease::fillBuffer(size_t block_num) const
{
	size_t adress = (2 + block_num) * m_blockSize_bytes;

	m_binaryFileStream.seekg(adress, std::ios::beg);

	m_binaryFileStream.read((char*)&m_buffer[0], (m_ncoeff) * 8);
}

void dph::EphemerisRelease::interpolatePosition(unsigned baseItemIndex, double normalizedTime,
	const double* coeffArray, unsigned componentsCount, double* resultArray) const
{
	// ����������� �������� ���������� ������������� �� ����������:
	uint32_t cpec = m_keys[baseItemIndex][1];

	// ��������������� ���������� ��������� (���������� �� ����):
	m_poly[1] = normalizedTime;

	// ���������� ��������� (���������� �� ����):
	for (uint32_t i = 2; i < cpec; ++i)
	{
		m_poly[i] = 2 * normalizedTime * m_poly[i - 1] - m_poly[i - 2];
	}

	// ��������� ������� ���������� ����������:
	memset(resultArray, 0, sizeof(double) * componentsCount);

	// ���������� ���������:
	for (unsigned i = 0; i < componentsCount; ++i)
	{
		for (uint32_t j = 0; j < cpec; ++j)
		{
			resultArray[i] += m_poly[j] * coeffArray[i * cpec + j];
		}
	}
}

void dph::EphemerisRelease::interpolateState(unsigned baseItemIndex, double normalizedTime,
	const double* coeffArray, unsigned componentsCount, double* resultArray) const
{
	// ����������� �������� ���������� ������������� �� ����������:
	uint32_t cpec = m_keys[baseItemIndex][1];

	// ��������������� ���������� ��������� (���������� �� ����):
	m_poly[1] = normalizedTime;
	m_poly[2] = 2 * normalizedTime * normalizedTime - 1;
	m_dpoly[2] = 4 * normalizedTime;

	// ���������� ��������� (���������� �� ����):
	for (uint32_t i = 3; i < cpec; ++i)
	{
		m_poly[i] = 2 * normalizedTime * m_poly[i - 1] - m_poly[i - 2];
		m_dpoly[i] = 2 * m_poly[i - 1] + 2 * normalizedTime * m_dpoly[i - 1] - m_dpoly[i - 2];
	}

	// ��������� ������� ���������� ����������:
	memset(resultArray, 0, sizeof(double) * componentsCount * 2);

	// ����������� ���������� ��� ���������� �����������:
	double derivative_units = m_keys[baseItemIndex][2] * m_dimensionFit;

	// ���������� ���������:
	for (unsigned i = 0; i < componentsCount; ++i)
	{
		for (uint32_t j = 0; j < cpec; ++j, ++coeffArray)
		{
			resultArray[i] += m_poly[j] * *coeffArray;
			resultArray[i + componentsCount] += m_dpoly[j] * *coeffArray;
		}

		resultArray[i + componentsCount] *= derivative_units;
	}
}

void dph::EphemerisRelease::calculateBaseItem(unsigned baseItemIndex, double JED,
	unsigned calculationResult, double* resultArray) const
{
	// ���������� �������� ���������� ����������:
	//	[1]	baseItemIndex - ������ �������� �������� ������� (�� ����).
	// 
	//						��������� ������� ��������� �������
	//		-------------------------------------------------------------------
	//		������	������������
	//		-------------------------------------------------------------------
	//		0		Mercury
	//		1		Venus
	//		2		Earth-Moon barycenter
	//		3		Mars
	//		4		Jupiter
	//		5		Saturn
	//		6		Uranus
	//		7		Neptune
	//		8		Pluto
	//		9		Moon (geocentric)
	//		10		Sun
	//		11		Earth Nutations in longitude and obliquity (IAU 1980 model)
	//		12		Lunar mantle libration
	//		13		Lunar mantle angular velocity
	//		14		TT-TDB (at geocenter)
	//		-------------------------------------------------------------------
	//	[2] JED - ������ ������� �� ������� ��������� �������� ��������� ��������.
	//	[3] calculationResult - ������ ���������� ���������� (��. dph::Calculate).
	//	[4] resultArray - ��������� �� ������ ��� ���������� ����������.

	// ��������! 
	// � ���� ���������� ������� ����� ���������� "normalizedTime" � "offset" ����� ��������.

	// ����. ����� ������������ ���� ������ � �������:
	double normalizedTime = (JED - m_startDate) / m_blockTimeSpan;

	// ���������� ����� �����, �����. �������� ���� JED (����� ����� �� normalizedTime):
	size_t offset = static_cast<size_t>(normalizedTime);

	// ���������� ������� �������������� ���������� �����.
	// ���� ��������� ���� ��� � ���� �������, �� �� �� ����������� ��������.
	// m_buffer[0] - ���� ������ �����.
	// m_buffer[1] - ���� ��������� �����.
	if (JED < m_buffer[0] || JED >= m_buffer[1])
	{
		// ���� JED ����� ��������� ���������� ���� ��� ����������, �� ����������� ��������� ����.

		fillBuffer(offset - (JED == m_endDate ? 1 : 0));
	}

	if (JED == m_endDate)
	{
		// ���������� ����� �������� (��������� �������):
		offset = m_keys[baseItemIndex][2] - 1;

		// ����. ����� ������������ �������� (� ��������� �� -1 �� 1):
		normalizedTime = 1;
	}
	else
	{
		// ����. ����� ������������ ���� ���������:
		normalizedTime = (normalizedTime - offset) * m_keys[baseItemIndex][2];

		// ���������� ����� �������� (����� ����� �� normalizedTime):
		offset = static_cast<size_t>(normalizedTime);

		// ����. ����� ������������ �������� (� ��������� �� -1 �� 1):
		normalizedTime = 2 * (normalizedTime - offset) - 1;
	}

	// ���������� ��������� ��� ���������� �������� ��������:
	unsigned componentsCount = baseItemIndex == 11 ? 2 : baseItemIndex == 14 ? 1 : 3;

	// ���������� ����� ������� ������������ � �����:
	int coeff_pos = m_keys[baseItemIndex][0] - 1 + componentsCount * offset * m_keys[baseItemIndex][1];

	// ����� ������ ���������� � ����������� �� ��������� ���������� ����������:
	switch (calculationResult)
	{
	case Calculate::POSITION:
		interpolatePosition(baseItemIndex, normalizedTime, &m_buffer[coeff_pos], componentsCount,
			resultArray);
		break;

	case Calculate::STATE:
		interpolateState(baseItemIndex, normalizedTime, &m_buffer[coeff_pos], componentsCount,
			resultArray);
		break;

	default:
		memset(resultArray, 0, componentsCount * sizeof(double));
	}
}

void dph::EphemerisRelease::calculateBaseEarth(double JED, unsigned calculationResult,
	double* resultArray) const
{
	// ��������� ������-������� (��� ������� ���������) ���������� �������� �����-����
	// ������������ ���������� ��������� �������:
	calculateBaseItem(2, JED, calculationResult, resultArray);

	// ��������� ������-������� (��� ������� ���������) ���� ����������� �����:
	double MoonRelativeEarth[6];
	calculateBaseItem(9, JED, calculationResult, MoonRelativeEarth);

	// ���������� ���������:
	unsigned componentsCount = calculationResult == Calculate::POSITION ? 3 : 6;

	// ���������� ��������� ����� ������������ ���������� ��������� �������:
	for (unsigned i = 0; i < componentsCount; ++i)
	{
		resultArray[i] -= MoonRelativeEarth[i] * m_emrat2;
	}
}

void dph::EphemerisRelease::calculateBaseMoon(double JED, unsigned calculationResult,
	double* resultArray) const
{
	// ��������� ������-������� (��� ������� ���������) ���������� �������� �����-����
	// ������������ ���������� ��������� �������:
	calculateBaseItem(2, JED, calculationResult, resultArray);

	// ��������� ������-������� (��� ������� ���������) ���� ����������� �����:
	double MoonRelativeEarth[6];
	calculateBaseItem(9, JED, calculationResult, MoonRelativeEarth);

	// ���������� ���������:
	unsigned componentsCount = calculationResult == Calculate::POSITION ? 3 : 6;

	// ����������� �������������� ���������:
	for (unsigned i = 0; i < componentsCount; ++i)
	{
		resultArray[i] += MoonRelativeEarth[i] * (1 - m_emrat2);
	}
}

