OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(0.19302364) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(2.4584682) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0177512) q[0];
sx q[0];
rz(-0.089028247) q[0];
sx q[0];
rz(0.71319367) q[0];
x q[1];
rz(1.921466) q[2];
sx q[2];
rz(-1.7979243) q[2];
sx q[2];
rz(2.7820058) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9477387) q[1];
sx q[1];
rz(-0.96809909) q[1];
sx q[1];
rz(-2.1147453) q[1];
rz(1.0115511) q[3];
sx q[3];
rz(-2.4992001) q[3];
sx q[3];
rz(-2.0410283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0156988) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(-2.0430298) q[2];
rz(1.0788318) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(-1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1612448) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(0.20794491) q[0];
rz(-0.57693276) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(-1.4651441) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10509051) q[0];
sx q[0];
rz(-0.72421342) q[0];
sx q[0];
rz(0.0037395517) q[0];
x q[1];
rz(-2.0287839) q[2];
sx q[2];
rz(-0.53456842) q[2];
sx q[2];
rz(-2.2528258) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0456902) q[1];
sx q[1];
rz(-2.0369139) q[1];
sx q[1];
rz(2.5065266) q[1];
rz(-pi) q[2];
rz(-1.915669) q[3];
sx q[3];
rz(-0.85540918) q[3];
sx q[3];
rz(-1.2777002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1229822) q[2];
sx q[2];
rz(-0.98615042) q[2];
sx q[2];
rz(-3.0318276) q[2];
rz(-0.62260735) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(-1.4573147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-0.85997471) q[0];
sx q[0];
rz(0.51613581) q[0];
rz(-2.5667045) q[1];
sx q[1];
rz(-0.92620414) q[1];
sx q[1];
rz(-0.80054545) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24414177) q[0];
sx q[0];
rz(-1.7139009) q[0];
sx q[0];
rz(-1.1600526) q[0];
rz(-2.3163047) q[2];
sx q[2];
rz(-2.2824259) q[2];
sx q[2];
rz(-2.3724144) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.70558207) q[1];
sx q[1];
rz(-1.5758621) q[1];
sx q[1];
rz(0.41298053) q[1];
rz(-pi) q[2];
rz(-0.3018467) q[3];
sx q[3];
rz(-0.89248025) q[3];
sx q[3];
rz(-2.9523926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4642554) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(1.3595954) q[2];
rz(-2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(-1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1490705) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(0.46491369) q[0];
rz(0.34856302) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(-2.0565313) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5497919) q[0];
sx q[0];
rz(-1.8013957) q[0];
sx q[0];
rz(2.6016298) q[0];
rz(-pi) q[1];
rz(2.1539139) q[2];
sx q[2];
rz(-2.7104125) q[2];
sx q[2];
rz(2.0476066) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19238732) q[1];
sx q[1];
rz(-2.9678223) q[1];
sx q[1];
rz(-2.5237571) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9784847) q[3];
sx q[3];
rz(-2.1424322) q[3];
sx q[3];
rz(-2.98416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1365635) q[2];
sx q[2];
rz(-1.0483142) q[2];
sx q[2];
rz(-0.68112779) q[2];
rz(2.629771) q[3];
sx q[3];
rz(-0.32326439) q[3];
sx q[3];
rz(0.26369035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27914771) q[0];
sx q[0];
rz(-0.537916) q[0];
sx q[0];
rz(-1.7329247) q[0];
rz(0.43235835) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(2.1599105) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9156993) q[0];
sx q[0];
rz(-2.3908273) q[0];
sx q[0];
rz(-2.4322926) q[0];
rz(-pi) q[1];
rz(-1.328674) q[2];
sx q[2];
rz(-1.7244581) q[2];
sx q[2];
rz(-0.15649934) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9202068) q[1];
sx q[1];
rz(-2.033794) q[1];
sx q[1];
rz(-1.4057926) q[1];
rz(0.97335191) q[3];
sx q[3];
rz(-1.5095599) q[3];
sx q[3];
rz(0.51613584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1061873) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(2.126746) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(-1.3283407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4846102) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(-0.37242517) q[0];
rz(-1.8107481) q[1];
sx q[1];
rz(-2.1061888) q[1];
sx q[1];
rz(0.16528027) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0342456) q[0];
sx q[0];
rz(-3.0421071) q[0];
sx q[0];
rz(2.8427567) q[0];
rz(-pi) q[1];
rz(-2.4539102) q[2];
sx q[2];
rz(-1.624794) q[2];
sx q[2];
rz(-0.24535594) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.11152553) q[1];
sx q[1];
rz(-2.1038052) q[1];
sx q[1];
rz(2.0433321) q[1];
x q[2];
rz(0.23541707) q[3];
sx q[3];
rz(-1.132292) q[3];
sx q[3];
rz(0.32270839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2281987) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(1.4432663) q[2];
rz(2.7741487) q[3];
sx q[3];
rz(-1.8668709) q[3];
sx q[3];
rz(2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3570324) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(-0.34926397) q[0];
rz(0.7473942) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(-2.4051037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1897141) q[0];
sx q[0];
rz(-1.5218381) q[0];
sx q[0];
rz(-3.1244856) q[0];
rz(-1.2320802) q[2];
sx q[2];
rz(-1.806353) q[2];
sx q[2];
rz(0.92600694) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.10662096) q[1];
sx q[1];
rz(-1.7502516) q[1];
sx q[1];
rz(-0.76645318) q[1];
rz(-2.442939) q[3];
sx q[3];
rz(-0.62056382) q[3];
sx q[3];
rz(-0.26649775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.13654576) q[2];
sx q[2];
rz(-1.2282635) q[2];
sx q[2];
rz(2.6100256) q[2];
rz(-2.452204) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(-1.846107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625967) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(-0.80728665) q[1];
sx q[1];
rz(-1.1785945) q[1];
sx q[1];
rz(-2.2198026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73496504) q[0];
sx q[0];
rz(-1.4359183) q[0];
sx q[0];
rz(-0.69358967) q[0];
x q[1];
rz(-1.8838896) q[2];
sx q[2];
rz(-2.3790092) q[2];
sx q[2];
rz(1.8956172) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8340048) q[1];
sx q[1];
rz(-1.7760135) q[1];
sx q[1];
rz(-2.8819041) q[1];
x q[2];
rz(1.3619625) q[3];
sx q[3];
rz(-0.71912557) q[3];
sx q[3];
rz(1.5987087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2395997) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(-1.9699338) q[2];
rz(1.7840067) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(1.4484891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8885324) q[0];
sx q[0];
rz(-1.1760412) q[0];
sx q[0];
rz(1.6660447) q[0];
rz(1.5400003) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(-2.4618861) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70521351) q[0];
sx q[0];
rz(-2.0853015) q[0];
sx q[0];
rz(-1.8041457) q[0];
rz(-1.7858511) q[2];
sx q[2];
rz(-2.4348223) q[2];
sx q[2];
rz(-0.099345318) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.745984) q[1];
sx q[1];
rz(-1.2320227) q[1];
sx q[1];
rz(-2.4376274) q[1];
rz(-0.9548095) q[3];
sx q[3];
rz(-2.1534352) q[3];
sx q[3];
rz(-2.5481176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1264964) q[2];
sx q[2];
rz(-1.6588147) q[2];
sx q[2];
rz(-0.91040197) q[2];
rz(2.4662468) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(-2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.05242059) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(-0.7146548) q[0];
rz(0.71406281) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(-1.9649327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0929716) q[0];
sx q[0];
rz(-1.4074568) q[0];
sx q[0];
rz(3.0924348) q[0];
rz(2.2628225) q[2];
sx q[2];
rz(-0.76239097) q[2];
sx q[2];
rz(2.9825485) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1706108) q[1];
sx q[1];
rz(-2.5787528) q[1];
sx q[1];
rz(-1.7112205) q[1];
rz(-pi) q[2];
rz(-2.9858399) q[3];
sx q[3];
rz(-0.51625801) q[3];
sx q[3];
rz(-1.932365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.24361336) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(2.0142377) q[2];
rz(2.7838498) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(1.0424785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7174299) q[0];
sx q[0];
rz(-1.2107727) q[0];
sx q[0];
rz(-2.6834224) q[0];
rz(0.39623109) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(0.73666112) q[2];
sx q[2];
rz(-1.7792637) q[2];
sx q[2];
rz(-1.7532495) q[2];
rz(1.5626004) q[3];
sx q[3];
rz(-1.2200439) q[3];
sx q[3];
rz(-0.72190819) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
