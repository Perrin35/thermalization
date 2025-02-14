OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1343159) q[0];
sx q[0];
rz(-2.4957823) q[0];
sx q[0];
rz(2.7752152) q[0];
rz(-2.3925048) q[1];
sx q[1];
rz(-1.6269416) q[1];
sx q[1];
rz(-2.9989624) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34829545) q[0];
sx q[0];
rz(-1.0401588) q[0];
sx q[0];
rz(0.9282309) q[0];
x q[1];
rz(0.22444522) q[2];
sx q[2];
rz(-1.4347174) q[2];
sx q[2];
rz(-2.4405757) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2097358) q[1];
sx q[1];
rz(-2.0448579) q[1];
sx q[1];
rz(-0.13412383) q[1];
x q[2];
rz(0.47358114) q[3];
sx q[3];
rz(-1.8627121) q[3];
sx q[3];
rz(1.1137258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.35965219) q[2];
sx q[2];
rz(-2.6617229) q[2];
sx q[2];
rz(-1.5894319) q[2];
rz(-1.3945256) q[3];
sx q[3];
rz(-0.37140578) q[3];
sx q[3];
rz(2.4594405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29779103) q[0];
sx q[0];
rz(-0.60972917) q[0];
sx q[0];
rz(0.30526701) q[0];
rz(1.8318532) q[1];
sx q[1];
rz(-0.42116183) q[1];
sx q[1];
rz(-3.0895244) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6604495) q[0];
sx q[0];
rz(-1.423712) q[0];
sx q[0];
rz(1.6818462) q[0];
x q[1];
rz(0.37392731) q[2];
sx q[2];
rz(-1.8513334) q[2];
sx q[2];
rz(0.032023059) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0525103) q[1];
sx q[1];
rz(-1.7351314) q[1];
sx q[1];
rz(1.3013863) q[1];
rz(-pi) q[2];
rz(0.90272119) q[3];
sx q[3];
rz(-1.4666712) q[3];
sx q[3];
rz(1.6834761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.084595844) q[2];
sx q[2];
rz(-1.8085542) q[2];
sx q[2];
rz(1.4196654) q[2];
rz(-2.2325884) q[3];
sx q[3];
rz(-2.3125068) q[3];
sx q[3];
rz(1.7349294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97515714) q[0];
sx q[0];
rz(-1.9571914) q[0];
sx q[0];
rz(-1.7433521) q[0];
rz(-2.1947529) q[1];
sx q[1];
rz(-2.0473174) q[1];
sx q[1];
rz(1.2907226) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2159444) q[0];
sx q[0];
rz(-0.15488347) q[0];
sx q[0];
rz(2.0808704) q[0];
rz(0.78611908) q[2];
sx q[2];
rz(-1.2183471) q[2];
sx q[2];
rz(-2.4967683) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5980893) q[1];
sx q[1];
rz(-0.072797983) q[1];
sx q[1];
rz(1.1612215) q[1];
rz(-pi) q[2];
rz(0.47770279) q[3];
sx q[3];
rz(-1.7509499) q[3];
sx q[3];
rz(-0.86828647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1281841) q[2];
sx q[2];
rz(-2.7161697) q[2];
sx q[2];
rz(-2.181633) q[2];
rz(0.32402447) q[3];
sx q[3];
rz(-2.6257381) q[3];
sx q[3];
rz(-0.44017756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0078761) q[0];
sx q[0];
rz(-0.77819264) q[0];
sx q[0];
rz(0.20009759) q[0];
rz(-2.5307185) q[1];
sx q[1];
rz(-2.4433544) q[1];
sx q[1];
rz(0.79455882) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66629529) q[0];
sx q[0];
rz(-1.7202416) q[0];
sx q[0];
rz(-1.4578865) q[0];
x q[1];
rz(-2.9616293) q[2];
sx q[2];
rz(-0.90562253) q[2];
sx q[2];
rz(2.6106204) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8243874) q[1];
sx q[1];
rz(-2.7911268) q[1];
sx q[1];
rz(2.0262655) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95501542) q[3];
sx q[3];
rz(-2.1977067) q[3];
sx q[3];
rz(1.1051429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2574629) q[2];
sx q[2];
rz(-0.56371671) q[2];
sx q[2];
rz(1.5719315) q[2];
rz(1.7852768) q[3];
sx q[3];
rz(-1.9808199) q[3];
sx q[3];
rz(-0.13264382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2785579) q[0];
sx q[0];
rz(-2.886241) q[0];
sx q[0];
rz(-2.417946) q[0];
rz(-3.0408995) q[1];
sx q[1];
rz(-1.0483402) q[1];
sx q[1];
rz(-3.0752799) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6479939) q[0];
sx q[0];
rz(-2.1382522) q[0];
sx q[0];
rz(0.074851224) q[0];
rz(-pi) q[1];
rz(-0.91513779) q[2];
sx q[2];
rz(-2.0210638) q[2];
sx q[2];
rz(1.428626) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4206754) q[1];
sx q[1];
rz(-2.2128792) q[1];
sx q[1];
rz(-2.9319619) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8961433) q[3];
sx q[3];
rz(-1.950437) q[3];
sx q[3];
rz(-1.0285447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2012607) q[2];
sx q[2];
rz(-1.2137493) q[2];
sx q[2];
rz(0.21855375) q[2];
rz(0.39441937) q[3];
sx q[3];
rz(-2.4557178) q[3];
sx q[3];
rz(-2.7421537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.047693096) q[0];
sx q[0];
rz(-2.221929) q[0];
sx q[0];
rz(-0.0041051824) q[0];
rz(-0.63402367) q[1];
sx q[1];
rz(-1.5019633) q[1];
sx q[1];
rz(-0.9563458) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8656971) q[0];
sx q[0];
rz(-1.529954) q[0];
sx q[0];
rz(2.9581822) q[0];
x q[1];
rz(-2.6317503) q[2];
sx q[2];
rz(-1.1948473) q[2];
sx q[2];
rz(-2.816566) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33935336) q[1];
sx q[1];
rz(-1.5763073) q[1];
sx q[1];
rz(-2.8637985) q[1];
rz(-pi) q[2];
rz(0.34849747) q[3];
sx q[3];
rz(-1.4607042) q[3];
sx q[3];
rz(1.7736721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.52028209) q[2];
sx q[2];
rz(-1.4352398) q[2];
sx q[2];
rz(0.83462805) q[2];
rz(-2.7708715) q[3];
sx q[3];
rz(-0.70086896) q[3];
sx q[3];
rz(2.4247775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9653559) q[0];
sx q[0];
rz(-0.014572425) q[0];
sx q[0];
rz(0.45005774) q[0];
rz(0.4484446) q[1];
sx q[1];
rz(-0.83201718) q[1];
sx q[1];
rz(-0.73743302) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0496488) q[0];
sx q[0];
rz(-1.3861602) q[0];
sx q[0];
rz(-1.139527) q[0];
rz(-pi) q[1];
rz(-0.040410553) q[2];
sx q[2];
rz(-1.4295661) q[2];
sx q[2];
rz(-0.34858957) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9818283) q[1];
sx q[1];
rz(-2.6609586) q[1];
sx q[1];
rz(-0.7483866) q[1];
rz(-pi) q[2];
rz(0.78011765) q[3];
sx q[3];
rz(-1.3778049) q[3];
sx q[3];
rz(-1.5867725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8488778) q[2];
sx q[2];
rz(-2.541031) q[2];
sx q[2];
rz(0.71504492) q[2];
rz(2.6333366) q[3];
sx q[3];
rz(-1.6786989) q[3];
sx q[3];
rz(1.3983294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7364863) q[0];
sx q[0];
rz(-2.5953601) q[0];
sx q[0];
rz(-0.35838321) q[0];
rz(2.7690923) q[1];
sx q[1];
rz(-1.7820216) q[1];
sx q[1];
rz(-2.701766) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51727879) q[0];
sx q[0];
rz(-0.77560878) q[0];
sx q[0];
rz(-1.6297352) q[0];
rz(-2.7372375) q[2];
sx q[2];
rz(-1.1431627) q[2];
sx q[2];
rz(0.15844395) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29319841) q[1];
sx q[1];
rz(-1.8351136) q[1];
sx q[1];
rz(-1.7536204) q[1];
x q[2];
rz(2.98073) q[3];
sx q[3];
rz(-2.8144112) q[3];
sx q[3];
rz(-0.88674816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0516574) q[2];
sx q[2];
rz(-3.0445485) q[2];
sx q[2];
rz(-2.4329321) q[2];
rz(0.25157252) q[3];
sx q[3];
rz(-0.99938649) q[3];
sx q[3];
rz(-3.1075509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8468903) q[0];
sx q[0];
rz(-1.0346233) q[0];
sx q[0];
rz(0.59144545) q[0];
rz(3.1239608) q[1];
sx q[1];
rz(-2.089274) q[1];
sx q[1];
rz(-3.0406521) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.712698) q[0];
sx q[0];
rz(-2.7980045) q[0];
sx q[0];
rz(1.1593007) q[0];
x q[1];
rz(-1.5315311) q[2];
sx q[2];
rz(-0.93446839) q[2];
sx q[2];
rz(-2.1992342) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4101319) q[1];
sx q[1];
rz(-1.4515775) q[1];
sx q[1];
rz(-1.340926) q[1];
x q[2];
rz(0.89102192) q[3];
sx q[3];
rz(-1.4263065) q[3];
sx q[3];
rz(-1.7107034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4883604) q[2];
sx q[2];
rz(-0.72335744) q[2];
sx q[2];
rz(1.4867268) q[2];
rz(2.9871353) q[3];
sx q[3];
rz(-1.1052701) q[3];
sx q[3];
rz(2.7753593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79009295) q[0];
sx q[0];
rz(-0.14227754) q[0];
sx q[0];
rz(1.0737786) q[0];
rz(2.0692661) q[1];
sx q[1];
rz(-0.25389478) q[1];
sx q[1];
rz(-0.059018746) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.260215) q[0];
sx q[0];
rz(-2.9395967) q[0];
sx q[0];
rz(-2.9750573) q[0];
x q[1];
rz(-0.26359908) q[2];
sx q[2];
rz(-2.0661743) q[2];
sx q[2];
rz(1.4857298) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6087338) q[1];
sx q[1];
rz(-0.97670943) q[1];
sx q[1];
rz(2.2284177) q[1];
rz(-pi) q[2];
rz(0.62982161) q[3];
sx q[3];
rz(-2.1849106) q[3];
sx q[3];
rz(0.30219004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3286288) q[2];
sx q[2];
rz(-1.8495411) q[2];
sx q[2];
rz(0.22895075) q[2];
rz(-1.1936584) q[3];
sx q[3];
rz(-0.3242068) q[3];
sx q[3];
rz(-1.4540023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1345632) q[0];
sx q[0];
rz(-2.3352191) q[0];
sx q[0];
rz(2.4051608) q[0];
rz(-1.3855343) q[1];
sx q[1];
rz(-1.3722739) q[1];
sx q[1];
rz(-1.3442232) q[1];
rz(-1.6177542) q[2];
sx q[2];
rz(-0.53971077) q[2];
sx q[2];
rz(-2.5765606) q[2];
rz(0.28382873) q[3];
sx q[3];
rz(-0.63224878) q[3];
sx q[3];
rz(-2.1716102) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
