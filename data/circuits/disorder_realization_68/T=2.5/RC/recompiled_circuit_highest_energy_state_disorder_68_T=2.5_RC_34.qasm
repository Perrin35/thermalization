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
rz(1.6990868) q[0];
sx q[0];
rz(-1.8449755) q[0];
sx q[0];
rz(-1.698864) q[0];
rz(2.7891085) q[1];
sx q[1];
rz(-2.6152857) q[1];
sx q[1];
rz(-1.7736645) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5621924) q[0];
sx q[0];
rz(-2.693667) q[0];
sx q[0];
rz(2.5125487) q[0];
x q[1];
rz(-1.5330731) q[2];
sx q[2];
rz(-1.2153271) q[2];
sx q[2];
rz(2.0672353) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.96208159) q[1];
sx q[1];
rz(-1.700811) q[1];
sx q[1];
rz(2.6165753) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75482224) q[3];
sx q[3];
rz(-1.1140545) q[3];
sx q[3];
rz(0.43908248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0736531) q[2];
sx q[2];
rz(-1.0509793) q[2];
sx q[2];
rz(-1.937872) q[2];
rz(2.1667571) q[3];
sx q[3];
rz(-0.9674415) q[3];
sx q[3];
rz(1.8956641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75342733) q[0];
sx q[0];
rz(-2.0797256) q[0];
sx q[0];
rz(0.88428307) q[0];
rz(0.70941225) q[1];
sx q[1];
rz(-1.246289) q[1];
sx q[1];
rz(0.90219227) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0220778) q[0];
sx q[0];
rz(-1.7762707) q[0];
sx q[0];
rz(0.46066649) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96376597) q[2];
sx q[2];
rz(-1.4861408) q[2];
sx q[2];
rz(0.36106685) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.25833669) q[1];
sx q[1];
rz(-2.751103) q[1];
sx q[1];
rz(-0.68895938) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0608474) q[3];
sx q[3];
rz(-1.4328379) q[3];
sx q[3];
rz(-1.6899275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.117729) q[2];
sx q[2];
rz(-2.8459097) q[2];
sx q[2];
rz(1.6531061) q[2];
rz(1.214341) q[3];
sx q[3];
rz(-1.6818654) q[3];
sx q[3];
rz(1.8620209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0073256) q[0];
sx q[0];
rz(-1.7420344) q[0];
sx q[0];
rz(-0.16635995) q[0];
rz(-2.4436489) q[1];
sx q[1];
rz(-0.78015399) q[1];
sx q[1];
rz(1.4770329) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0255942) q[0];
sx q[0];
rz(-0.4420155) q[0];
sx q[0];
rz(1.563795) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2899621) q[2];
sx q[2];
rz(-1.8321494) q[2];
sx q[2];
rz(-0.42119831) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0946413) q[1];
sx q[1];
rz(-1.1186386) q[1];
sx q[1];
rz(0.96479123) q[1];
rz(-pi) q[2];
rz(0.03607492) q[3];
sx q[3];
rz(-2.6429308) q[3];
sx q[3];
rz(0.3968249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.6537689) q[2];
sx q[2];
rz(-1.9705801) q[2];
sx q[2];
rz(0.6146532) q[2];
rz(1.4608308) q[3];
sx q[3];
rz(-1.5644282) q[3];
sx q[3];
rz(-0.041725807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6093269) q[0];
sx q[0];
rz(-1.8121239) q[0];
sx q[0];
rz(-1.6229269) q[0];
rz(2.3329349) q[1];
sx q[1];
rz(-1.8452019) q[1];
sx q[1];
rz(-0.048695806) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0322389) q[0];
sx q[0];
rz(-1.236602) q[0];
sx q[0];
rz(-2.4346057) q[0];
rz(2.8559309) q[2];
sx q[2];
rz(-0.39644074) q[2];
sx q[2];
rz(-2.2224521) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.16597) q[1];
sx q[1];
rz(-0.73117729) q[1];
sx q[1];
rz(-2.380409) q[1];
x q[2];
rz(-0.96541578) q[3];
sx q[3];
rz(-1.9002689) q[3];
sx q[3];
rz(0.62240619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0398756) q[2];
sx q[2];
rz(-0.81284916) q[2];
sx q[2];
rz(1.8854878) q[2];
rz(0.060955437) q[3];
sx q[3];
rz(-2.239949) q[3];
sx q[3];
rz(0.92756334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1561279) q[0];
sx q[0];
rz(-1.1277072) q[0];
sx q[0];
rz(0.10966478) q[0];
rz(-2.1127286) q[1];
sx q[1];
rz(-1.8548465) q[1];
sx q[1];
rz(1.5922155) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76413918) q[0];
sx q[0];
rz(-1.8528588) q[0];
sx q[0];
rz(0.46746032) q[0];
rz(-pi) q[1];
rz(-2.0711871) q[2];
sx q[2];
rz(-1.8922046) q[2];
sx q[2];
rz(-1.6381581) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2031496) q[1];
sx q[1];
rz(-1.1091653) q[1];
sx q[1];
rz(1.0814971) q[1];
rz(-pi) q[2];
rz(3.1072633) q[3];
sx q[3];
rz(-1.6084371) q[3];
sx q[3];
rz(0.08480367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.5515543) q[2];
sx q[2];
rz(-0.2912713) q[2];
sx q[2];
rz(-2.6596098) q[2];
rz(2.7216952) q[3];
sx q[3];
rz(-1.0651383) q[3];
sx q[3];
rz(2.4317252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9919306) q[0];
sx q[0];
rz(-0.88352942) q[0];
sx q[0];
rz(-0.71257198) q[0];
rz(-0.86992162) q[1];
sx q[1];
rz(-0.37418071) q[1];
sx q[1];
rz(3.1178927) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11940322) q[0];
sx q[0];
rz(-1.4076122) q[0];
sx q[0];
rz(0.52363445) q[0];
rz(-pi) q[1];
rz(-2.4392088) q[2];
sx q[2];
rz(-1.1105892) q[2];
sx q[2];
rz(0.61379877) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1501173) q[1];
sx q[1];
rz(-1.5693226) q[1];
sx q[1];
rz(1.7365843) q[1];
rz(-pi) q[2];
rz(-1.7810443) q[3];
sx q[3];
rz(-2.0422404) q[3];
sx q[3];
rz(-1.2725079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.80453834) q[2];
sx q[2];
rz(-2.4961175) q[2];
sx q[2];
rz(-2.0013334) q[2];
rz(-1.1537457) q[3];
sx q[3];
rz(-1.2146344) q[3];
sx q[3];
rz(1.8203189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7158647) q[0];
sx q[0];
rz(-1.3874929) q[0];
sx q[0];
rz(0.36780372) q[0];
rz(-1.6416719) q[1];
sx q[1];
rz(-0.92253128) q[1];
sx q[1];
rz(-2.0466764) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.80847) q[0];
sx q[0];
rz(-2.5980691) q[0];
sx q[0];
rz(2.031206) q[0];
rz(-2.5911683) q[2];
sx q[2];
rz(-1.4604521) q[2];
sx q[2];
rz(-0.36729022) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3907539) q[1];
sx q[1];
rz(-1.8065284) q[1];
sx q[1];
rz(1.5524992) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1672375) q[3];
sx q[3];
rz(-1.5512244) q[3];
sx q[3];
rz(2.7579466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0423923) q[2];
sx q[2];
rz(-1.3513869) q[2];
sx q[2];
rz(-3.0858827) q[2];
rz(-0.67692155) q[3];
sx q[3];
rz(-2.5261295) q[3];
sx q[3];
rz(0.87853471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95213503) q[0];
sx q[0];
rz(-0.5624693) q[0];
sx q[0];
rz(0.96543717) q[0];
rz(-1.2646487) q[1];
sx q[1];
rz(-0.74600428) q[1];
sx q[1];
rz(-0.3269349) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.801689) q[0];
sx q[0];
rz(-1.9205695) q[0];
sx q[0];
rz(-3.0236493) q[0];
rz(-2.6956396) q[2];
sx q[2];
rz(-1.5805681) q[2];
sx q[2];
rz(-2.8251512) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4070774) q[1];
sx q[1];
rz(-2.8663969) q[1];
sx q[1];
rz(0.86517548) q[1];
rz(-1.502874) q[3];
sx q[3];
rz(-1.1227896) q[3];
sx q[3];
rz(-0.94326708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.38291976) q[2];
sx q[2];
rz(-0.70432538) q[2];
sx q[2];
rz(2.3455589) q[2];
rz(0.087513611) q[3];
sx q[3];
rz(-0.60860601) q[3];
sx q[3];
rz(-2.2039738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4336808) q[0];
sx q[0];
rz(-1.7681363) q[0];
sx q[0];
rz(0.32980907) q[0];
rz(3.0682796) q[1];
sx q[1];
rz(-2.4344756) q[1];
sx q[1];
rz(-1.0940394) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0819223) q[0];
sx q[0];
rz(-2.0412894) q[0];
sx q[0];
rz(-1.612603) q[0];
rz(-1.7297001) q[2];
sx q[2];
rz(-2.7393118) q[2];
sx q[2];
rz(2.2251468) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0702522) q[1];
sx q[1];
rz(-1.0254745) q[1];
sx q[1];
rz(1.7252183) q[1];
x q[2];
rz(3.1018267) q[3];
sx q[3];
rz(-1.13969) q[3];
sx q[3];
rz(2.414961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0988934) q[2];
sx q[2];
rz(-2.3697479) q[2];
sx q[2];
rz(2.0043376) q[2];
rz(0.62816652) q[3];
sx q[3];
rz(-0.38574949) q[3];
sx q[3];
rz(2.1342733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6657418) q[0];
sx q[0];
rz(-2.6306212) q[0];
sx q[0];
rz(-1.092859) q[0];
rz(1.8765556) q[1];
sx q[1];
rz(-1.3053514) q[1];
sx q[1];
rz(-1.0337894) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40202478) q[0];
sx q[0];
rz(-2.0962068) q[0];
sx q[0];
rz(-2.5891586) q[0];
rz(-pi) q[1];
x q[1];
rz(2.381778) q[2];
sx q[2];
rz(-1.7489479) q[2];
sx q[2];
rz(0.15195981) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6792843) q[1];
sx q[1];
rz(-2.1888715) q[1];
sx q[1];
rz(1.0584153) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3218811) q[3];
sx q[3];
rz(-2.9985709) q[3];
sx q[3];
rz(1.7211357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.25896245) q[2];
sx q[2];
rz(-0.7236824) q[2];
sx q[2];
rz(0.21931973) q[2];
rz(0.80260459) q[3];
sx q[3];
rz(-1.8290627) q[3];
sx q[3];
rz(-2.820211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5103067) q[0];
sx q[0];
rz(-1.318537) q[0];
sx q[0];
rz(1.0986811) q[0];
rz(-0.17768271) q[1];
sx q[1];
rz(-1.0164574) q[1];
sx q[1];
rz(-0.16998092) q[1];
rz(0.86831696) q[2];
sx q[2];
rz(-1.6196031) q[2];
sx q[2];
rz(-2.6283787) q[2];
rz(3.0232676) q[3];
sx q[3];
rz(-1.7601624) q[3];
sx q[3];
rz(2.7047529) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
