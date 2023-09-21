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
rz(4.1376576) q[0];
sx q[0];
rz(7.1538038) q[0];
rz(-1.0215966) q[1];
sx q[1];
rz(-0.28290132) q[1];
sx q[1];
rz(2.9918848) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5538841) q[0];
sx q[0];
rz(-0.93548488) q[0];
sx q[0];
rz(2.5250838) q[0];
rz(-pi) q[1];
rz(1.0295463) q[2];
sx q[2];
rz(-2.0143348) q[2];
sx q[2];
rz(-2.8461547) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.61923164) q[1];
sx q[1];
rz(-1.3904966) q[1];
sx q[1];
rz(0.038832263) q[1];
rz(0.19102328) q[3];
sx q[3];
rz(-2.1130307) q[3];
sx q[3];
rz(-2.1477826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25847882) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(-2.4374938) q[2];
rz(2.1885833) q[3];
sx q[3];
rz(-2.1694031) q[3];
sx q[3];
rz(1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54884058) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(0.63252226) q[0];
rz(-0.44644341) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(-2.4893563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3632293) q[0];
sx q[0];
rz(-0.031154545) q[0];
sx q[0];
rz(-0.40701436) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2781497) q[2];
sx q[2];
rz(-1.7990944) q[2];
sx q[2];
rz(-1.6739068) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3398509) q[1];
sx q[1];
rz(-2.7738214) q[1];
sx q[1];
rz(-0.66373177) q[1];
rz(2.5370876) q[3];
sx q[3];
rz(-1.0009871) q[3];
sx q[3];
rz(1.2165716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5474881) q[2];
sx q[2];
rz(-2.1274121) q[2];
sx q[2];
rz(1.9799505) q[2];
rz(-1.15796) q[3];
sx q[3];
rz(-1.0693113) q[3];
sx q[3];
rz(1.6903711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.050425477) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(-0.85025775) q[0];
rz(2.6440874) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(1.7920378) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3916546) q[0];
sx q[0];
rz(-0.91021252) q[0];
sx q[0];
rz(-1.2664938) q[0];
rz(-pi) q[1];
rz(2.0850052) q[2];
sx q[2];
rz(-3.0658256) q[2];
sx q[2];
rz(0.55759341) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.459356) q[1];
sx q[1];
rz(-0.56486928) q[1];
sx q[1];
rz(2.6911246) q[1];
rz(-1.5924256) q[3];
sx q[3];
rz(-2.7154185) q[3];
sx q[3];
rz(2.222995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0321908) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(0.90399495) q[2];
rz(-2.8404625) q[3];
sx q[3];
rz(-1.7588994) q[3];
sx q[3];
rz(-1.3999456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49144739) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(-1.7310463) q[0];
rz(2.5097805) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(3.1052123) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0050126652) q[0];
sx q[0];
rz(-2.4653325) q[0];
sx q[0];
rz(-1.9146634) q[0];
rz(-pi) q[1];
rz(2.6668947) q[2];
sx q[2];
rz(-2.5050852) q[2];
sx q[2];
rz(2.0091025) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8424884) q[1];
sx q[1];
rz(-0.71055382) q[1];
sx q[1];
rz(1.7732265) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39156885) q[3];
sx q[3];
rz(-2.8885926) q[3];
sx q[3];
rz(0.60245017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.92695421) q[2];
sx q[2];
rz(-1.4321233) q[2];
sx q[2];
rz(-1.9533763) q[2];
rz(0.67048091) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(-2.5454583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(-0.16648509) q[0];
rz(2.3855551) q[1];
sx q[1];
rz(-0.88587228) q[1];
sx q[1];
rz(0.23434815) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4011824) q[0];
sx q[0];
rz(-2.4966842) q[0];
sx q[0];
rz(-0.58437225) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3946103) q[2];
sx q[2];
rz(-0.50054769) q[2];
sx q[2];
rz(-1.2622152) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25917945) q[1];
sx q[1];
rz(-0.83173527) q[1];
sx q[1];
rz(-1.3684567) q[1];
rz(-1.8022728) q[3];
sx q[3];
rz(-2.8095062) q[3];
sx q[3];
rz(1.0486697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.70871893) q[2];
sx q[2];
rz(-1.9659698) q[2];
sx q[2];
rz(-0.22053545) q[2];
rz(-0.43705127) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41546145) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(-2.3244526) q[0];
rz(-2.5754886) q[1];
sx q[1];
rz(-1.7931466) q[1];
sx q[1];
rz(1.1436499) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2957942) q[0];
sx q[0];
rz(-2.7685071) q[0];
sx q[0];
rz(3.0315115) q[0];
rz(-pi) q[1];
rz(-2.4684858) q[2];
sx q[2];
rz(-1.7761201) q[2];
sx q[2];
rz(2.7538607) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.027187849) q[1];
sx q[1];
rz(-1.4804375) q[1];
sx q[1];
rz(2.0143709) q[1];
x q[2];
rz(2.6753622) q[3];
sx q[3];
rz(-2.5737408) q[3];
sx q[3];
rz(0.0044435244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3879261) q[2];
sx q[2];
rz(-1.5619229) q[2];
sx q[2];
rz(-0.4894408) q[2];
rz(-2.9135381) q[3];
sx q[3];
rz(-1.8830048) q[3];
sx q[3];
rz(-2.7155546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5722826) q[0];
sx q[0];
rz(-2.4991878) q[0];
sx q[0];
rz(1.8547159) q[0];
rz(0.6634179) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(1.2333262) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39499261) q[0];
sx q[0];
rz(-1.1867503) q[0];
sx q[0];
rz(-0.011944255) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6255409) q[2];
sx q[2];
rz(-2.1893246) q[2];
sx q[2];
rz(0.90460888) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8334956) q[1];
sx q[1];
rz(-1.6819681) q[1];
sx q[1];
rz(1.5566467) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2193905) q[3];
sx q[3];
rz(-1.1099585) q[3];
sx q[3];
rz(1.8691065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.104091) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(2.1255778) q[2];
rz(3.0715023) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(-1.0664553) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12524097) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(-1.6960779) q[0];
rz(2.9267172) q[1];
sx q[1];
rz(-2.3863249) q[1];
sx q[1];
rz(1.8833556) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1930775) q[0];
sx q[0];
rz(-2.6876039) q[0];
sx q[0];
rz(1.2598739) q[0];
x q[1];
rz(0.7779185) q[2];
sx q[2];
rz(-2.3292543) q[2];
sx q[2];
rz(-1.1493491) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80553493) q[1];
sx q[1];
rz(-1.4667257) q[1];
sx q[1];
rz(0.54688262) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6120841) q[3];
sx q[3];
rz(-0.92748517) q[3];
sx q[3];
rz(-2.1852126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.68226472) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(-2.5788467) q[2];
rz(-0.19966666) q[3];
sx q[3];
rz(-0.66185799) q[3];
sx q[3];
rz(-2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
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
rz(-0.02773157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6513034) q[0];
sx q[0];
rz(-0.9285183) q[0];
sx q[0];
rz(2.519033) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2559782) q[2];
sx q[2];
rz(-0.72394365) q[2];
sx q[2];
rz(-2.1997423) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8393644) q[1];
sx q[1];
rz(-1.6703969) q[1];
sx q[1];
rz(2.539413) q[1];
x q[2];
rz(1.6373709) q[3];
sx q[3];
rz(-1.5449617) q[3];
sx q[3];
rz(-1.0915826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.015908265) q[2];
sx q[2];
rz(-0.97970825) q[2];
sx q[2];
rz(-1.1431747) q[2];
rz(-2.9987191) q[3];
sx q[3];
rz(-1.5121127) q[3];
sx q[3];
rz(2.2843602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3906355) q[0];
sx q[0];
rz(-1.9366783) q[0];
sx q[0];
rz(2.6877158) q[0];
rz(-2.4699396) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(2.8840816) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2902381) q[0];
sx q[0];
rz(-2.5332753) q[0];
sx q[0];
rz(-0.71382199) q[0];
rz(0.72980482) q[2];
sx q[2];
rz(-2.4591755) q[2];
sx q[2];
rz(-0.44621106) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9092907) q[1];
sx q[1];
rz(-1.5850987) q[1];
sx q[1];
rz(0.69625744) q[1];
rz(-pi) q[2];
x q[2];
rz(2.014124) q[3];
sx q[3];
rz(-0.52237288) q[3];
sx q[3];
rz(0.70538196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3283078) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(-2.5349687) q[2];
rz(-2.666752) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(-0.84038466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7941147) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(-2.9150302) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(-2.4702934) q[2];
sx q[2];
rz(-0.57325267) q[2];
sx q[2];
rz(-0.43043955) q[2];
rz(1.0971309) q[3];
sx q[3];
rz(-2.6087425) q[3];
sx q[3];
rz(-2.9125924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];