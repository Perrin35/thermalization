OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3671626) q[0];
sx q[0];
rz(-2.2280333) q[0];
sx q[0];
rz(1.7295184) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(5.6875416) q[1];
sx q[1];
rz(7.7653801) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1971561) q[0];
sx q[0];
rz(-1.4154139) q[0];
sx q[0];
rz(0.42874254) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6719195) q[2];
sx q[2];
rz(-0.28684068) q[2];
sx q[2];
rz(-1.6490205) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5414121) q[1];
sx q[1];
rz(-2.0611079) q[1];
sx q[1];
rz(-2.6580826) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1336018) q[3];
sx q[3];
rz(-1.6562914) q[3];
sx q[3];
rz(-0.29719719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98510629) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(-0.86581725) q[2];
rz(-0.95430294) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(-1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1433379) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(-3.1153733) q[0];
rz(1.6014618) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(-0.96347934) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026982633) q[0];
sx q[0];
rz(-2.5275505) q[0];
sx q[0];
rz(1.5747889) q[0];
rz(-0.33032592) q[2];
sx q[2];
rz(-2.1147554) q[2];
sx q[2];
rz(0.75737539) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1184247) q[1];
sx q[1];
rz(-1.2512565) q[1];
sx q[1];
rz(-2.097514) q[1];
rz(-1.2198592) q[3];
sx q[3];
rz(-1.3388472) q[3];
sx q[3];
rz(1.2082781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6271237) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(0.13452402) q[2];
rz(-0.7450122) q[3];
sx q[3];
rz(-2.9146505) q[3];
sx q[3];
rz(-2.1988595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2117675) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(-0.79743687) q[0];
rz(-1.047661) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(-0.55999666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5138429) q[0];
sx q[0];
rz(-1.3686413) q[0];
sx q[0];
rz(-1.1802799) q[0];
rz(-pi) q[1];
rz(2.838344) q[2];
sx q[2];
rz(-1.5911284) q[2];
sx q[2];
rz(-3.0603527) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8078976) q[1];
sx q[1];
rz(-1.2119319) q[1];
sx q[1];
rz(1.3586504) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3476944) q[3];
sx q[3];
rz(-0.63459914) q[3];
sx q[3];
rz(-1.3018228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3893163) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(0.17253549) q[2];
rz(-0.98207384) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(-1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.003222) q[0];
sx q[0];
rz(-1.0976185) q[0];
sx q[0];
rz(2.8570783) q[0];
rz(0.31670397) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(1.8428615) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38628681) q[0];
sx q[0];
rz(-0.55314976) q[0];
sx q[0];
rz(1.132071) q[0];
rz(-1.5721333) q[2];
sx q[2];
rz(-0.7664116) q[2];
sx q[2];
rz(3.1198451) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.58060927) q[1];
sx q[1];
rz(-1.0099851) q[1];
sx q[1];
rz(-0.22052712) q[1];
rz(-0.11233791) q[3];
sx q[3];
rz(-1.897052) q[3];
sx q[3];
rz(-2.5374967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(0.38468012) q[2];
rz(-2.3875333) q[3];
sx q[3];
rz(-2.0850756) q[3];
sx q[3];
rz(-1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5383179) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.7549365) q[0];
rz(2.9105913) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(-2.8447661) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1332069) q[0];
sx q[0];
rz(-2.1413681) q[0];
sx q[0];
rz(1.5242566) q[0];
x q[1];
rz(0.16634059) q[2];
sx q[2];
rz(-2.3918249) q[2];
sx q[2];
rz(-1.2321842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2424803) q[1];
sx q[1];
rz(-0.63226262) q[1];
sx q[1];
rz(2.2805023) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4279168) q[3];
sx q[3];
rz(-1.5034961) q[3];
sx q[3];
rz(1.0825368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.37830535) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(-0.39247593) q[2];
rz(1.9893507) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(-2.8241482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1257989) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(-0.75138599) q[0];
rz(1.8136576) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(-0.60633916) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6316815) q[0];
sx q[0];
rz(-1.5874377) q[0];
sx q[0];
rz(-0.04583866) q[0];
rz(0.89247993) q[2];
sx q[2];
rz(-1.2391029) q[2];
sx q[2];
rz(1.4075556) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.025758) q[1];
sx q[1];
rz(-2.5587213) q[1];
sx q[1];
rz(-1.238766) q[1];
x q[2];
rz(-1.8984406) q[3];
sx q[3];
rz(-2.9122105) q[3];
sx q[3];
rz(2.4250507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5027344) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(-2.0022557) q[2];
rz(1.4849439) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(0.10425723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6095603) q[0];
sx q[0];
rz(-0.74815265) q[0];
sx q[0];
rz(0.50810057) q[0];
rz(1.5787026) q[1];
sx q[1];
rz(-1.0888313) q[1];
sx q[1];
rz(-2.3513444) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.761844) q[0];
sx q[0];
rz(-1.2756691) q[0];
sx q[0];
rz(2.5248812) q[0];
rz(-pi) q[1];
x q[1];
rz(3.087567) q[2];
sx q[2];
rz(-2.9276491) q[2];
sx q[2];
rz(-2.8097048) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0223169) q[1];
sx q[1];
rz(-1.1806618) q[1];
sx q[1];
rz(1.8632061) q[1];
rz(2.6584133) q[3];
sx q[3];
rz(-2.4348767) q[3];
sx q[3];
rz(2.8578575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0885075) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(-1.7283758) q[2];
rz(-1.6648071) q[3];
sx q[3];
rz(-2.1063185) q[3];
sx q[3];
rz(3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7168032) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(-1.0472263) q[0];
rz(-0.60910243) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.3887127) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4318651) q[0];
sx q[0];
rz(-2.1301529) q[0];
sx q[0];
rz(-3.1064242) q[0];
rz(1.9954761) q[2];
sx q[2];
rz(-1.0815902) q[2];
sx q[2];
rz(-2.3236772) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70506239) q[1];
sx q[1];
rz(-0.092518004) q[1];
sx q[1];
rz(3.0180879) q[1];
rz(-pi) q[2];
rz(-1.1561398) q[3];
sx q[3];
rz(-0.64838833) q[3];
sx q[3];
rz(-0.63601953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1887112) q[2];
sx q[2];
rz(-2.7313576) q[2];
sx q[2];
rz(-2.2593373) q[2];
rz(-1.4011718) q[3];
sx q[3];
rz(-1.975235) q[3];
sx q[3];
rz(1.9410979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8700478) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(1.7154988) q[0];
rz(0.081461279) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(-2.5833599) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8919864) q[0];
sx q[0];
rz(-2.5476646) q[0];
sx q[0];
rz(-0.83017613) q[0];
rz(-pi) q[1];
rz(1.6782645) q[2];
sx q[2];
rz(-1.5170013) q[2];
sx q[2];
rz(1.6966284) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.12802943) q[1];
sx q[1];
rz(-0.82847825) q[1];
sx q[1];
rz(-0.086704266) q[1];
rz(-pi) q[2];
rz(-0.34096034) q[3];
sx q[3];
rz(-0.51135671) q[3];
sx q[3];
rz(-0.90358678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2925064) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(2.9294087) q[2];
rz(-0.21197453) q[3];
sx q[3];
rz(-0.68325716) q[3];
sx q[3];
rz(-1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.6367209) q[0];
sx q[0];
rz(-0.8240521) q[0];
sx q[0];
rz(-1.5378392) q[0];
rz(2.3161855) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(2.6182981) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035318035) q[0];
sx q[0];
rz(-0.9629074) q[0];
sx q[0];
rz(-1.6465228) q[0];
rz(-pi) q[1];
rz(-1.5027572) q[2];
sx q[2];
rz(-0.57300742) q[2];
sx q[2];
rz(2.9266561) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8955252) q[1];
sx q[1];
rz(-1.9290036) q[1];
sx q[1];
rz(-1.9073061) q[1];
rz(-pi) q[2];
rz(1.80199) q[3];
sx q[3];
rz(-1.8095067) q[3];
sx q[3];
rz(-2.0774487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.837073) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(0.26930299) q[2];
rz(2.6473911) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(-1.0860898) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8158648) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(-1.6745463) q[1];
sx q[1];
rz(-2.84927) q[1];
sx q[1];
rz(-1.897859) q[1];
rz(1.6124484) q[2];
sx q[2];
rz(-0.26298444) q[2];
sx q[2];
rz(2.0900805) q[2];
rz(-1.8105103) q[3];
sx q[3];
rz(-0.78835434) q[3];
sx q[3];
rz(2.2314856) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
