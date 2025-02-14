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
rz(-1.7614814) q[0];
sx q[0];
rz(-1.5505646) q[0];
sx q[0];
rz(-2.759759) q[0];
rz(2.7497357) q[1];
sx q[1];
rz(-1.5500103) q[1];
sx q[1];
rz(-0.92441192) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43542433) q[0];
sx q[0];
rz(-2.2924137) q[0];
sx q[0];
rz(-2.3090786) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7713791) q[2];
sx q[2];
rz(-1.5892803) q[2];
sx q[2];
rz(0.12272515) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0142483) q[1];
sx q[1];
rz(-2.6963628) q[1];
sx q[1];
rz(-0.30502747) q[1];
rz(-pi) q[2];
rz(-2.1788039) q[3];
sx q[3];
rz(-0.6198403) q[3];
sx q[3];
rz(-0.40871295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2048637) q[2];
sx q[2];
rz(-1.2932581) q[2];
sx q[2];
rz(1.5748242) q[2];
rz(1.9751679) q[3];
sx q[3];
rz(-2.0765442) q[3];
sx q[3];
rz(-2.4488357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4676056) q[0];
sx q[0];
rz(-2.3335712) q[0];
sx q[0];
rz(1.1760733) q[0];
rz(-3.0740956) q[1];
sx q[1];
rz(-2.564023) q[1];
sx q[1];
rz(2.8228021) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27869895) q[0];
sx q[0];
rz(-0.79813939) q[0];
sx q[0];
rz(0.44095914) q[0];
rz(-2.712557) q[2];
sx q[2];
rz(-1.3554159) q[2];
sx q[2];
rz(-1.693415) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.76004878) q[1];
sx q[1];
rz(-2.4439619) q[1];
sx q[1];
rz(-2.182644) q[1];
x q[2];
rz(-0.12572581) q[3];
sx q[3];
rz(-1.2866402) q[3];
sx q[3];
rz(-1.8234258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0626283) q[2];
sx q[2];
rz(-0.88125172) q[2];
sx q[2];
rz(-2.6641565) q[2];
rz(1.3580953) q[3];
sx q[3];
rz(-2.4886459) q[3];
sx q[3];
rz(-0.0099946578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39182144) q[0];
sx q[0];
rz(-1.3359767) q[0];
sx q[0];
rz(-1.1392449) q[0];
rz(-2.7353752) q[1];
sx q[1];
rz(-1.8832877) q[1];
sx q[1];
rz(2.4635945) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0292873) q[0];
sx q[0];
rz(-0.59354085) q[0];
sx q[0];
rz(1.7303921) q[0];
rz(-pi) q[1];
rz(1.0933206) q[2];
sx q[2];
rz(-0.68308631) q[2];
sx q[2];
rz(-0.25743279) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3542676) q[1];
sx q[1];
rz(-0.72871043) q[1];
sx q[1];
rz(0.72737741) q[1];
rz(-2.1265245) q[3];
sx q[3];
rz(-2.56009) q[3];
sx q[3];
rz(-2.8648928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71538007) q[2];
sx q[2];
rz(-1.6343626) q[2];
sx q[2];
rz(-1.0769843) q[2];
rz(-1.8814603) q[3];
sx q[3];
rz(-1.0271122) q[3];
sx q[3];
rz(-0.25585678) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28427163) q[0];
sx q[0];
rz(-1.1742598) q[0];
sx q[0];
rz(0.76048365) q[0];
rz(0.35176945) q[1];
sx q[1];
rz(-0.71134174) q[1];
sx q[1];
rz(-0.15028353) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4486769) q[0];
sx q[0];
rz(-1.5674855) q[0];
sx q[0];
rz(-3.1396554) q[0];
rz(0.17268659) q[2];
sx q[2];
rz(-1.9479247) q[2];
sx q[2];
rz(1.0972301) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1898415) q[1];
sx q[1];
rz(-0.69165666) q[1];
sx q[1];
rz(-2.1053949) q[1];
x q[2];
rz(-2.6500449) q[3];
sx q[3];
rz(-0.50954362) q[3];
sx q[3];
rz(0.59922892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.030423) q[2];
sx q[2];
rz(-0.22746484) q[2];
sx q[2];
rz(-2.5566768) q[2];
rz(3.0758744) q[3];
sx q[3];
rz(-2.0021555) q[3];
sx q[3];
rz(-2.0871302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2718662) q[0];
sx q[0];
rz(-1.9166742) q[0];
sx q[0];
rz(2.4140893) q[0];
rz(1.8456521) q[1];
sx q[1];
rz(-2.0869052) q[1];
sx q[1];
rz(0.71472275) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22435221) q[0];
sx q[0];
rz(-0.2324129) q[0];
sx q[0];
rz(-2.8281141) q[0];
rz(-pi) q[1];
rz(-1.4552015) q[2];
sx q[2];
rz(-1.8697479) q[2];
sx q[2];
rz(1.40708) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3989563) q[1];
sx q[1];
rz(-0.42544895) q[1];
sx q[1];
rz(2.0320973) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39387624) q[3];
sx q[3];
rz(-2.1112313) q[3];
sx q[3];
rz(0.5462164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1456445) q[2];
sx q[2];
rz(-1.6137292) q[2];
sx q[2];
rz(-1.9963416) q[2];
rz(0.20259419) q[3];
sx q[3];
rz(-0.069772094) q[3];
sx q[3];
rz(-0.32874671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2957434) q[0];
sx q[0];
rz(-0.29243094) q[0];
sx q[0];
rz(2.9737245) q[0];
rz(0.56495086) q[1];
sx q[1];
rz(-1.047537) q[1];
sx q[1];
rz(-1.3289183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.560379) q[0];
sx q[0];
rz(-0.74958505) q[0];
sx q[0];
rz(-2.4102274) q[0];
rz(-pi) q[1];
rz(0.14018709) q[2];
sx q[2];
rz(-1.5221688) q[2];
sx q[2];
rz(2.1339773) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5649453) q[1];
sx q[1];
rz(-2.7912346) q[1];
sx q[1];
rz(0.74585657) q[1];
rz(-pi) q[2];
rz(-2.9181913) q[3];
sx q[3];
rz(-1.3376682) q[3];
sx q[3];
rz(-2.4883771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.38587511) q[2];
sx q[2];
rz(-1.2473325) q[2];
sx q[2];
rz(2.7937549) q[2];
rz(0.31473413) q[3];
sx q[3];
rz(-2.0458524) q[3];
sx q[3];
rz(-2.0033964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9921853) q[0];
sx q[0];
rz(-1.5062165) q[0];
sx q[0];
rz(2.1904679) q[0];
rz(-0.29257193) q[1];
sx q[1];
rz(-0.89076275) q[1];
sx q[1];
rz(-2.1975885) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68444618) q[0];
sx q[0];
rz(-1.5553885) q[0];
sx q[0];
rz(1.8788565) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47796504) q[2];
sx q[2];
rz(-2.0066119) q[2];
sx q[2];
rz(2.8675241) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.71783644) q[1];
sx q[1];
rz(-1.5482117) q[1];
sx q[1];
rz(1.3976239) q[1];
rz(2.8947325) q[3];
sx q[3];
rz(-1.180397) q[3];
sx q[3];
rz(0.91871432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6501179) q[2];
sx q[2];
rz(-0.74345165) q[2];
sx q[2];
rz(-1.7000807) q[2];
rz(3.1169543) q[3];
sx q[3];
rz(-1.8162497) q[3];
sx q[3];
rz(2.1520481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80855075) q[0];
sx q[0];
rz(-1.4104183) q[0];
sx q[0];
rz(-0.1380052) q[0];
rz(-2.0637312) q[1];
sx q[1];
rz(-1.4638008) q[1];
sx q[1];
rz(-0.93625751) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1753006) q[0];
sx q[0];
rz(-1.202824) q[0];
sx q[0];
rz(1.8938246) q[0];
rz(-pi) q[1];
x q[1];
rz(2.594844) q[2];
sx q[2];
rz(-0.90675747) q[2];
sx q[2];
rz(2.4647317) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.030979217) q[1];
sx q[1];
rz(-0.18337164) q[1];
sx q[1];
rz(-1.6513837) q[1];
rz(-pi) q[2];
rz(-1.3055757) q[3];
sx q[3];
rz(-2.239003) q[3];
sx q[3];
rz(2.1611155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.998698) q[2];
sx q[2];
rz(-1.8691209) q[2];
sx q[2];
rz(1.8227089) q[2];
rz(0.12254347) q[3];
sx q[3];
rz(-2.4428941) q[3];
sx q[3];
rz(2.7114649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1634624) q[0];
sx q[0];
rz(-0.4438816) q[0];
sx q[0];
rz(-2.7769856) q[0];
rz(1.7568024) q[1];
sx q[1];
rz(-0.93162799) q[1];
sx q[1];
rz(-2.2273831) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4791661) q[0];
sx q[0];
rz(-1.437133) q[0];
sx q[0];
rz(-2.0019309) q[0];
rz(-3.0905452) q[2];
sx q[2];
rz(-2.406359) q[2];
sx q[2];
rz(0.027054199) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.89489854) q[1];
sx q[1];
rz(-1.7082038) q[1];
sx q[1];
rz(-2.4989884) q[1];
x q[2];
rz(-2.6007341) q[3];
sx q[3];
rz(-1.4424994) q[3];
sx q[3];
rz(-3.1311664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.30283516) q[2];
sx q[2];
rz(-0.5842394) q[2];
sx q[2];
rz(-1.6287104) q[2];
rz(-2.5527939) q[3];
sx q[3];
rz(-1.2062807) q[3];
sx q[3];
rz(-0.65620667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6673679) q[0];
sx q[0];
rz(-0.56741699) q[0];
sx q[0];
rz(2.8571416) q[0];
rz(-2.4868763) q[1];
sx q[1];
rz(-1.7660716) q[1];
sx q[1];
rz(-2.7213352) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66061879) q[0];
sx q[0];
rz(-1.7848564) q[0];
sx q[0];
rz(-1.995489) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1070232) q[2];
sx q[2];
rz(-0.83597413) q[2];
sx q[2];
rz(-2.7157264) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6984562) q[1];
sx q[1];
rz(-1.4520169) q[1];
sx q[1];
rz(-0.08929792) q[1];
rz(-pi) q[2];
rz(0.56575601) q[3];
sx q[3];
rz(-1.2671736) q[3];
sx q[3];
rz(0.88241258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1493211) q[2];
sx q[2];
rz(-2.6565266) q[2];
sx q[2];
rz(0.9672271) q[2];
rz(1.017717) q[3];
sx q[3];
rz(-1.2844205) q[3];
sx q[3];
rz(-0.73219901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16348542) q[0];
sx q[0];
rz(-2.0757984) q[0];
sx q[0];
rz(2.1708873) q[0];
rz(0.12374395) q[1];
sx q[1];
rz(-2.1950304) q[1];
sx q[1];
rz(-1.310941) q[1];
rz(0.32044784) q[2];
sx q[2];
rz(-1.4437699) q[2];
sx q[2];
rz(-0.52158471) q[2];
rz(1.5354034) q[3];
sx q[3];
rz(-1.8826859) q[3];
sx q[3];
rz(-1.7425698) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
