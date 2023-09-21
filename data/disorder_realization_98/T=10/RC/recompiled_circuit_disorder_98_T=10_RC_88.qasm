OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37031072) q[0];
sx q[0];
rz(-2.1455278) q[0];
sx q[0];
rz(-2.2709742) q[0];
rz(2.119996) q[1];
sx q[1];
rz(-2.8586913) q[1];
sx q[1];
rz(0.14970782) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5538841) q[0];
sx q[0];
rz(-0.93548488) q[0];
sx q[0];
rz(0.61650886) q[0];
x q[1];
rz(2.3157273) q[2];
sx q[2];
rz(-2.4561433) q[2];
sx q[2];
rz(0.65537383) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.40587273) q[1];
sx q[1];
rz(-2.957203) q[1];
sx q[1];
rz(-1.7806446) q[1];
x q[2];
rz(2.9505694) q[3];
sx q[3];
rz(-2.1130307) q[3];
sx q[3];
rz(2.1477826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8831138) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(2.4374938) q[2];
rz(0.95300931) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927521) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(2.5090704) q[0];
rz(2.6951492) q[1];
sx q[1];
rz(-1.4182785) q[1];
sx q[1];
rz(2.4893563) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95603847) q[0];
sx q[0];
rz(-1.5421876) q[0];
sx q[0];
rz(-1.5584598) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2488238) q[2];
sx q[2];
rz(-2.7724578) q[2];
sx q[2];
rz(0.54112753) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3995041) q[1];
sx q[1];
rz(-1.3474476) q[1];
sx q[1];
rz(-0.29466596) q[1];
rz(0.90918031) q[3];
sx q[3];
rz(-2.069807) q[3];
sx q[3];
rz(0.0024851174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5941045) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(-1.1616421) q[2];
rz(1.9836327) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050425477) q[0];
sx q[0];
rz(-2.7624891) q[0];
sx q[0];
rz(0.85025775) q[0];
rz(-2.6440874) q[1];
sx q[1];
rz(-1.9583227) q[1];
sx q[1];
rz(1.7920378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3916546) q[0];
sx q[0];
rz(-2.2313801) q[0];
sx q[0];
rz(-1.2664938) q[0];
x q[1];
rz(-3.1042728) q[2];
sx q[2];
rz(-1.63675) q[2];
sx q[2];
rz(-1.0730336) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.68223665) q[1];
sx q[1];
rz(-2.5767234) q[1];
sx q[1];
rz(0.45046803) q[1];
rz(-pi) q[2];
rz(-1.1447103) q[3];
sx q[3];
rz(-1.5618556) q[3];
sx q[3];
rz(0.63250354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1094018) q[2];
sx q[2];
rz(-1.9058062) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(-2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(-1.3999456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49144739) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(-1.7310463) q[0];
rz(-2.5097805) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(0.036380336) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0050126652) q[0];
sx q[0];
rz(-0.67626017) q[0];
sx q[0];
rz(-1.9146634) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8965917) q[2];
sx q[2];
rz(-1.0137644) q[2];
sx q[2];
rz(-1.4404802) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42602793) q[1];
sx q[1];
rz(-1.4392816) q[1];
sx q[1];
rz(0.87042602) q[1];
x q[2];
rz(-0.39156885) q[3];
sx q[3];
rz(-2.8885926) q[3];
sx q[3];
rz(2.5391425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2146384) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(1.9533763) q[2];
rz(2.4711117) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(-0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4398414) q[0];
sx q[0];
rz(-1.259946) q[0];
sx q[0];
rz(-0.16648509) q[0];
rz(-0.75603756) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(2.9072445) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4317961) q[0];
sx q[0];
rz(-2.095982) q[0];
sx q[0];
rz(-1.964142) q[0];
x q[1];
rz(-1.7469823) q[2];
sx q[2];
rz(-0.50054769) q[2];
sx q[2];
rz(-1.2622152) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.25917945) q[1];
sx q[1];
rz(-0.83173527) q[1];
sx q[1];
rz(-1.773136) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3393199) q[3];
sx q[3];
rz(-2.8095062) q[3];
sx q[3];
rz(-1.0486697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4328737) q[2];
sx q[2];
rz(-1.1756228) q[2];
sx q[2];
rz(0.22053545) q[2];
rz(-2.7045414) q[3];
sx q[3];
rz(-2.1190937) q[3];
sx q[3];
rz(2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7261312) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(2.3244526) q[0];
rz(2.5754886) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(1.1436499) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7276579) q[0];
sx q[0];
rz(-1.9415138) q[0];
sx q[0];
rz(-1.6137705) q[0];
rz(-pi) q[1];
rz(0.32240378) q[2];
sx q[2];
rz(-2.44256) q[2];
sx q[2];
rz(0.93271819) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.027187849) q[1];
sx q[1];
rz(-1.4804375) q[1];
sx q[1];
rz(-1.1272217) q[1];
x q[2];
rz(2.6753622) q[3];
sx q[3];
rz(-0.56785184) q[3];
sx q[3];
rz(-0.0044435244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3879261) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(2.6521519) q[2];
rz(2.9135381) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(-2.7155546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.56931) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(1.8547159) q[0];
rz(-2.4781748) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(1.2333262) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7466) q[0];
sx q[0];
rz(-1.1867503) q[0];
sx q[0];
rz(-0.011944255) q[0];
x q[1];
rz(0.6255409) q[2];
sx q[2];
rz(-2.1893246) q[2];
sx q[2];
rz(-2.2369838) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7066321) q[1];
sx q[1];
rz(-0.11206493) q[1];
sx q[1];
rz(-0.12607615) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9222021) q[3];
sx q[3];
rz(-1.1099585) q[3];
sx q[3];
rz(-1.8691065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.037501637) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(2.1255778) q[2];
rz(-3.0715023) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12524097) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(-1.6960779) q[0];
rz(-2.9267172) q[1];
sx q[1];
rz(-2.3863249) q[1];
sx q[1];
rz(1.258237) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90342605) q[0];
sx q[0];
rz(-1.705372) q[0];
sx q[0];
rz(2.0057136) q[0];
x q[1];
rz(-2.3636742) q[2];
sx q[2];
rz(-2.3292543) q[2];
sx q[2];
rz(1.9922436) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.82842365) q[1];
sx q[1];
rz(-2.1143882) q[1];
sx q[1];
rz(-1.6924752) q[1];
rz(-pi) q[2];
rz(2.1636837) q[3];
sx q[3];
rz(-0.80855723) q[3];
sx q[3];
rz(2.9582994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4593279) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(-2.5788467) q[2];
rz(-2.941926) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0257618) q[0];
sx q[0];
rz(-2.0985726) q[0];
sx q[0];
rz(2.8905706) q[0];
rz(-0.42731467) q[1];
sx q[1];
rz(-1.2298093) q[1];
sx q[1];
rz(3.1138611) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4902892) q[0];
sx q[0];
rz(-2.2130744) q[0];
sx q[0];
rz(0.62255967) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2559782) q[2];
sx q[2];
rz(-0.72394365) q[2];
sx q[2];
rz(0.9418504) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8393644) q[1];
sx q[1];
rz(-1.4711958) q[1];
sx q[1];
rz(-2.539413) q[1];
rz(-pi) q[2];
rz(0.025891993) q[3];
sx q[3];
rz(-1.504244) q[3];
sx q[3];
rz(2.6641012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.015908265) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(1.1431747) q[2];
rz(2.9987191) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(-0.85723248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7509572) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(-0.45387682) q[0];
rz(2.4699396) q[1];
sx q[1];
rz(-1.6891054) q[1];
sx q[1];
rz(-0.25751105) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2430902) q[0];
sx q[0];
rz(-1.1872963) q[0];
sx q[0];
rz(-2.6570508) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54458877) q[2];
sx q[2];
rz(-2.0047744) q[2];
sx q[2];
rz(1.7314272) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.820206) q[1];
sx q[1];
rz(-0.69637978) q[1];
sx q[1];
rz(0.022298261) q[1];
rz(-pi) q[2];
rz(2.8994843) q[3];
sx q[3];
rz(-1.1032411) q[3];
sx q[3];
rz(-2.9374591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8132849) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(-2.5349687) q[2];
rz(-2.666752) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(-2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7941147) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(-0.22656245) q[1];
sx q[1];
rz(-1.7270052) q[1];
sx q[1];
rz(-2.5546767) q[1];
rz(-2.4702934) q[2];
sx q[2];
rz(-0.57325267) q[2];
sx q[2];
rz(-0.43043955) q[2];
rz(2.0541035) q[3];
sx q[3];
rz(-1.3369505) q[3];
sx q[3];
rz(-1.7575775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
