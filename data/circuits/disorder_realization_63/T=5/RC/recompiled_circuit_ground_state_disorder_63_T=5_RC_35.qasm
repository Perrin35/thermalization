OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37501332) q[0];
sx q[0];
rz(4.3809173) q[0];
sx q[0];
rz(8.3349174) q[0];
rz(1.0640979) q[1];
sx q[1];
rz(4.394726) q[1];
sx q[1];
rz(10.328007) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2555465) q[0];
sx q[0];
rz(-0.95661344) q[0];
sx q[0];
rz(0.79844676) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1268812) q[2];
sx q[2];
rz(-0.88764578) q[2];
sx q[2];
rz(-2.8840182) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.049202327) q[1];
sx q[1];
rz(-0.70194879) q[1];
sx q[1];
rz(2.5939221) q[1];
rz(-pi) q[2];
rz(-1.0711477) q[3];
sx q[3];
rz(-1.6766492) q[3];
sx q[3];
rz(1.9377886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4452867) q[2];
sx q[2];
rz(-0.13662766) q[2];
sx q[2];
rz(-2.7083) q[2];
rz(0.28996921) q[3];
sx q[3];
rz(-0.86499298) q[3];
sx q[3];
rz(-0.32052952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9669773) q[0];
sx q[0];
rz(-2.2791635) q[0];
sx q[0];
rz(-1.8509266) q[0];
rz(2.4621452) q[1];
sx q[1];
rz(-2.3652855) q[1];
sx q[1];
rz(-0.89250934) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3474059) q[0];
sx q[0];
rz(-1.5958324) q[0];
sx q[0];
rz(1.7480609) q[0];
x q[1];
rz(-0.50267668) q[2];
sx q[2];
rz(-0.71841769) q[2];
sx q[2];
rz(-0.65332149) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.034169056) q[1];
sx q[1];
rz(-0.71798199) q[1];
sx q[1];
rz(-2.7173244) q[1];
rz(-pi) q[2];
rz(0.99231798) q[3];
sx q[3];
rz(-0.43799339) q[3];
sx q[3];
rz(-0.43644825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9925925) q[2];
sx q[2];
rz(-1.0789824) q[2];
sx q[2];
rz(-2.0089669) q[2];
rz(2.4752786) q[3];
sx q[3];
rz(-1.3601114) q[3];
sx q[3];
rz(1.9168436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3742974) q[0];
sx q[0];
rz(-0.39300028) q[0];
sx q[0];
rz(-2.1597916) q[0];
rz(-0.94217316) q[1];
sx q[1];
rz(-1.3092382) q[1];
sx q[1];
rz(0.91032496) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019429723) q[0];
sx q[0];
rz(-1.0896026) q[0];
sx q[0];
rz(2.0519755) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.40159638) q[2];
sx q[2];
rz(-0.29916468) q[2];
sx q[2];
rz(-0.49002346) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.83214789) q[1];
sx q[1];
rz(-0.4570241) q[1];
sx q[1];
rz(0.8330658) q[1];
rz(-pi) q[2];
rz(2.8998834) q[3];
sx q[3];
rz(-1.1430267) q[3];
sx q[3];
rz(-2.857353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.94074774) q[2];
sx q[2];
rz(-2.7767599) q[2];
sx q[2];
rz(2.9546837) q[2];
rz(2.9777891) q[3];
sx q[3];
rz(-2.2713594) q[3];
sx q[3];
rz(-3.0189309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0201037) q[0];
sx q[0];
rz(-1.7262456) q[0];
sx q[0];
rz(2.4439268) q[0];
rz(-0.54667073) q[1];
sx q[1];
rz(-0.29269871) q[1];
sx q[1];
rz(-1.1950511) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.228926) q[0];
sx q[0];
rz(-0.9879092) q[0];
sx q[0];
rz(-0.41622644) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.706188) q[2];
sx q[2];
rz(-1.7601331) q[2];
sx q[2];
rz(-2.9385081) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.89151284) q[1];
sx q[1];
rz(-1.7367474) q[1];
sx q[1];
rz(3.0338698) q[1];
rz(-pi) q[2];
rz(-1.2863808) q[3];
sx q[3];
rz(-1.6413601) q[3];
sx q[3];
rz(1.1190991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4635072) q[2];
sx q[2];
rz(-1.4231851) q[2];
sx q[2];
rz(2.9441693) q[2];
rz(-3.0366963) q[3];
sx q[3];
rz(-1.1095122) q[3];
sx q[3];
rz(3.0535898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5474434) q[0];
sx q[0];
rz(-2.6056885) q[0];
sx q[0];
rz(-2.1323668) q[0];
rz(0.028845305) q[1];
sx q[1];
rz(-1.6639158) q[1];
sx q[1];
rz(-0.65863329) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.665425) q[0];
sx q[0];
rz(-1.9840249) q[0];
sx q[0];
rz(2.7318003) q[0];
rz(-pi) q[1];
rz(-1.1971742) q[2];
sx q[2];
rz(-2.6251617) q[2];
sx q[2];
rz(-0.85631285) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4411021) q[1];
sx q[1];
rz(-1.4785188) q[1];
sx q[1];
rz(2.9725595) q[1];
rz(1.3405209) q[3];
sx q[3];
rz(-2.252823) q[3];
sx q[3];
rz(1.7868228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7625526) q[2];
sx q[2];
rz(-0.54658824) q[2];
sx q[2];
rz(0.23400447) q[2];
rz(1.0903357) q[3];
sx q[3];
rz(-1.5749616) q[3];
sx q[3];
rz(1.0696629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2284018) q[0];
sx q[0];
rz(-2.985432) q[0];
sx q[0];
rz(1.8432023) q[0];
rz(-0.78650728) q[1];
sx q[1];
rz(-2.2857917) q[1];
sx q[1];
rz(-1.7826805) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6494736) q[0];
sx q[0];
rz(-2.2303061) q[0];
sx q[0];
rz(1.6716206) q[0];
x q[1];
rz(-2.8739019) q[2];
sx q[2];
rz(-1.3548193) q[2];
sx q[2];
rz(-2.9962073) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7715367) q[1];
sx q[1];
rz(-0.0020469804) q[1];
sx q[1];
rz(-2.5084346) q[1];
x q[2];
rz(-0.46192138) q[3];
sx q[3];
rz(-1.2826903) q[3];
sx q[3];
rz(1.896331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7274373) q[2];
sx q[2];
rz(-1.4667908) q[2];
sx q[2];
rz(0.1725014) q[2];
rz(-0.48834673) q[3];
sx q[3];
rz(-1.1427053) q[3];
sx q[3];
rz(-2.02777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97335029) q[0];
sx q[0];
rz(-2.2689447) q[0];
sx q[0];
rz(2.8016222) q[0];
rz(2.8151457) q[1];
sx q[1];
rz(-0.68845922) q[1];
sx q[1];
rz(-1.9065769) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3351072) q[0];
sx q[0];
rz(-0.77101427) q[0];
sx q[0];
rz(-2.5946027) q[0];
x q[1];
rz(-1.9307617) q[2];
sx q[2];
rz(-1.1380929) q[2];
sx q[2];
rz(-0.77322799) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.82385027) q[1];
sx q[1];
rz(-0.47951937) q[1];
sx q[1];
rz(0.01156263) q[1];
rz(-pi) q[2];
rz(2.6458287) q[3];
sx q[3];
rz(-2.7003482) q[3];
sx q[3];
rz(3.038681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1116703) q[2];
sx q[2];
rz(-2.3387574) q[2];
sx q[2];
rz(-0.55533448) q[2];
rz(-0.16767821) q[3];
sx q[3];
rz(-1.5372814) q[3];
sx q[3];
rz(-2.663747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70917201) q[0];
sx q[0];
rz(-2.2116311) q[0];
sx q[0];
rz(-1.1093371) q[0];
rz(-2.0093911) q[1];
sx q[1];
rz(-0.39674509) q[1];
sx q[1];
rz(2.1999377) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3170211) q[0];
sx q[0];
rz(-2.5618658) q[0];
sx q[0];
rz(-1.2767904) q[0];
x q[1];
rz(-1.2852816) q[2];
sx q[2];
rz(-1.7196349) q[2];
sx q[2];
rz(2.6251453) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4577427) q[1];
sx q[1];
rz(-0.93194973) q[1];
sx q[1];
rz(0.15775494) q[1];
x q[2];
rz(1.3579943) q[3];
sx q[3];
rz(-1.9486685) q[3];
sx q[3];
rz(-1.6869643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.38810101) q[2];
sx q[2];
rz(-2.1145623) q[2];
sx q[2];
rz(1.6027742) q[2];
rz(-1.7704891) q[3];
sx q[3];
rz(-1.1857827) q[3];
sx q[3];
rz(-0.051430844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0592773) q[0];
sx q[0];
rz(-0.61510724) q[0];
sx q[0];
rz(-2.1737461) q[0];
rz(-1.8702501) q[1];
sx q[1];
rz(-1.2920734) q[1];
sx q[1];
rz(-0.06180067) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7907352) q[0];
sx q[0];
rz(-2.4567607) q[0];
sx q[0];
rz(2.5032296) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4872753) q[2];
sx q[2];
rz(-2.8192602) q[2];
sx q[2];
rz(-0.92609012) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.73229746) q[1];
sx q[1];
rz(-1.4753072) q[1];
sx q[1];
rz(0.68640253) q[1];
x q[2];
rz(0.7155719) q[3];
sx q[3];
rz(-0.4288097) q[3];
sx q[3];
rz(-2.9675067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1900078) q[2];
sx q[2];
rz(-2.9374359) q[2];
sx q[2];
rz(-0.07587138) q[2];
rz(-1.9742981) q[3];
sx q[3];
rz(-1.1018402) q[3];
sx q[3];
rz(0.30437881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9645914) q[0];
sx q[0];
rz(-2.8622506) q[0];
sx q[0];
rz(2.8569073) q[0];
rz(-0.074706569) q[1];
sx q[1];
rz(-1.9120522) q[1];
sx q[1];
rz(1.4178735) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7265948) q[0];
sx q[0];
rz(-2.3624067) q[0];
sx q[0];
rz(2.8654773) q[0];
rz(-0.14947628) q[2];
sx q[2];
rz(-2.1227395) q[2];
sx q[2];
rz(-1.9665444) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1571341) q[1];
sx q[1];
rz(-2.1878831) q[1];
sx q[1];
rz(-2.0871215) q[1];
x q[2];
rz(-1.6154434) q[3];
sx q[3];
rz(-0.64638019) q[3];
sx q[3];
rz(-2.165739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3384) q[2];
sx q[2];
rz(-1.1407547) q[2];
sx q[2];
rz(-1.2791862) q[2];
rz(-2.159436) q[3];
sx q[3];
rz(-1.4582062) q[3];
sx q[3];
rz(-0.88472432) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90047705) q[0];
sx q[0];
rz(-1.8178839) q[0];
sx q[0];
rz(-1.5771014) q[0];
rz(-2.1263532) q[1];
sx q[1];
rz(-1.992234) q[1];
sx q[1];
rz(2.0043859) q[1];
rz(-3.0334658) q[2];
sx q[2];
rz(-1.8320638) q[2];
sx q[2];
rz(-2.0083357) q[2];
rz(2.5231936) q[3];
sx q[3];
rz(-0.84155819) q[3];
sx q[3];
rz(0.92084259) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
