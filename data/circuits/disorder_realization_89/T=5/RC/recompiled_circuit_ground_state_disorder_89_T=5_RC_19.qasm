OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2616413) q[0];
sx q[0];
rz(-0.14544848) q[0];
sx q[0];
rz(0.26785904) q[0];
rz(-1.5675867) q[1];
sx q[1];
rz(-0.1685473) q[1];
sx q[1];
rz(0.57810098) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6689396) q[0];
sx q[0];
rz(-0.90757209) q[0];
sx q[0];
rz(2.7884463) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4431345) q[2];
sx q[2];
rz(-1.6751117) q[2];
sx q[2];
rz(1.8836752) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4207124) q[1];
sx q[1];
rz(-2.0053177) q[1];
sx q[1];
rz(1.1195393) q[1];
rz(1.3760482) q[3];
sx q[3];
rz(-0.71310242) q[3];
sx q[3];
rz(-0.76228599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9428923) q[2];
sx q[2];
rz(-0.41112021) q[2];
sx q[2];
rz(1.4760419) q[2];
rz(2.5623411) q[3];
sx q[3];
rz(-1.9871291) q[3];
sx q[3];
rz(-2.4756685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2838374) q[0];
sx q[0];
rz(-1.0455766) q[0];
sx q[0];
rz(3.0864518) q[0];
rz(-0.91446963) q[1];
sx q[1];
rz(-1.6998467) q[1];
sx q[1];
rz(-1.9257911) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037579868) q[0];
sx q[0];
rz(-1.5079588) q[0];
sx q[0];
rz(1.5893905) q[0];
rz(-pi) q[1];
rz(2.6138805) q[2];
sx q[2];
rz(-1.1368504) q[2];
sx q[2];
rz(-0.82810706) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3275274) q[1];
sx q[1];
rz(-1.9694929) q[1];
sx q[1];
rz(1.5507803) q[1];
x q[2];
rz(2.6917879) q[3];
sx q[3];
rz(-1.6946892) q[3];
sx q[3];
rz(1.9270437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2083644) q[2];
sx q[2];
rz(-2.6231982) q[2];
sx q[2];
rz(-2.5210157) q[2];
rz(-1.4536475) q[3];
sx q[3];
rz(-1.0317289) q[3];
sx q[3];
rz(1.0769963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.8715912) q[0];
sx q[0];
rz(-0.33137614) q[0];
sx q[0];
rz(1.510386) q[0];
rz(-0.67391467) q[1];
sx q[1];
rz(-1.8464512) q[1];
sx q[1];
rz(-1.6253701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6166068) q[0];
sx q[0];
rz(-1.5459037) q[0];
sx q[0];
rz(-2.1014433) q[0];
x q[1];
rz(2.7582232) q[2];
sx q[2];
rz(-0.96770937) q[2];
sx q[2];
rz(-0.71345873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4067799) q[1];
sx q[1];
rz(-1.3441097) q[1];
sx q[1];
rz(-0.43155687) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2760157) q[3];
sx q[3];
rz(-1.0350482) q[3];
sx q[3];
rz(-1.2212041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38500938) q[2];
sx q[2];
rz(-1.5417121) q[2];
sx q[2];
rz(2.4760683) q[2];
rz(-2.3982128) q[3];
sx q[3];
rz(-2.0205108) q[3];
sx q[3];
rz(-0.22751787) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3709563) q[0];
sx q[0];
rz(-1.0524858) q[0];
sx q[0];
rz(1.8413405) q[0];
rz(-3.0216253) q[1];
sx q[1];
rz(-1.080546) q[1];
sx q[1];
rz(-1.0467451) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20489731) q[0];
sx q[0];
rz(-0.012210695) q[0];
sx q[0];
rz(1.8050675) q[0];
rz(-0.8530944) q[2];
sx q[2];
rz(-1.2893558) q[2];
sx q[2];
rz(2.6996343) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1018104) q[1];
sx q[1];
rz(-2.6521195) q[1];
sx q[1];
rz(1.6881949) q[1];
x q[2];
rz(-0.98174121) q[3];
sx q[3];
rz(-1.9509754) q[3];
sx q[3];
rz(-1.5721934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7946502) q[2];
sx q[2];
rz(-1.0618989) q[2];
sx q[2];
rz(2.3700355) q[2];
rz(1.0509342) q[3];
sx q[3];
rz(-2.1080878) q[3];
sx q[3];
rz(1.1933901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0070888) q[0];
sx q[0];
rz(-1.6935885) q[0];
sx q[0];
rz(-2.336308) q[0];
rz(1.443642) q[1];
sx q[1];
rz(-2.3256358) q[1];
sx q[1];
rz(-0.03646341) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59234658) q[0];
sx q[0];
rz(-1.0039465) q[0];
sx q[0];
rz(-0.057307883) q[0];
rz(0.068425325) q[2];
sx q[2];
rz(-1.003661) q[2];
sx q[2];
rz(0.65862818) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.441022) q[1];
sx q[1];
rz(-2.6614958) q[1];
sx q[1];
rz(-1.455485) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.24954911) q[3];
sx q[3];
rz(-0.59246906) q[3];
sx q[3];
rz(-1.0147926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61520758) q[2];
sx q[2];
rz(-1.1247164) q[2];
sx q[2];
rz(2.5649694) q[2];
rz(-1.5455101) q[3];
sx q[3];
rz(-0.85448623) q[3];
sx q[3];
rz(0.57687783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0312408) q[0];
sx q[0];
rz(-1.4875655) q[0];
sx q[0];
rz(-2.5417969) q[0];
rz(-0.80232969) q[1];
sx q[1];
rz(-2.0498514) q[1];
sx q[1];
rz(-0.64782992) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0698034) q[0];
sx q[0];
rz(-2.2043071) q[0];
sx q[0];
rz(-2.5833552) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6036649) q[2];
sx q[2];
rz(-0.21015659) q[2];
sx q[2];
rz(1.5446784) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5831754) q[1];
sx q[1];
rz(-0.98924556) q[1];
sx q[1];
rz(1.2152332) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4418026) q[3];
sx q[3];
rz(-1.0465396) q[3];
sx q[3];
rz(-1.4211224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.99737793) q[2];
sx q[2];
rz(-0.37313676) q[2];
sx q[2];
rz(-1.0972265) q[2];
rz(-1.0002452) q[3];
sx q[3];
rz(-2.2179243) q[3];
sx q[3];
rz(-1.8358811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10709396) q[0];
sx q[0];
rz(-1.1544363) q[0];
sx q[0];
rz(-2.9423998) q[0];
rz(1.2983324) q[1];
sx q[1];
rz(-0.36111626) q[1];
sx q[1];
rz(-2.4995506) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6605153) q[0];
sx q[0];
rz(-0.84781983) q[0];
sx q[0];
rz(-2.8882746) q[0];
rz(-pi) q[1];
x q[1];
rz(2.726905) q[2];
sx q[2];
rz(-2.0609591) q[2];
sx q[2];
rz(1.8275593) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.43993801) q[1];
sx q[1];
rz(-1.3716193) q[1];
sx q[1];
rz(1.6526196) q[1];
x q[2];
rz(-2.0208259) q[3];
sx q[3];
rz(-1.3880669) q[3];
sx q[3];
rz(2.8830626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9127427) q[2];
sx q[2];
rz(-2.0759089) q[2];
sx q[2];
rz(-1.9672811) q[2];
rz(-0.92285389) q[3];
sx q[3];
rz(-0.43761161) q[3];
sx q[3];
rz(-2.9722884) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6628424) q[0];
sx q[0];
rz(-1.6599382) q[0];
sx q[0];
rz(-0.99916512) q[0];
rz(2.303458) q[1];
sx q[1];
rz(-2.2331388) q[1];
sx q[1];
rz(-1.025544) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9860172) q[0];
sx q[0];
rz(-2.3087569) q[0];
sx q[0];
rz(-0.9132333) q[0];
rz(0.77825688) q[2];
sx q[2];
rz(-2.7822251) q[2];
sx q[2];
rz(-3.0148413) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2237812) q[1];
sx q[1];
rz(-1.1780699) q[1];
sx q[1];
rz(2.1070943) q[1];
x q[2];
rz(2.5426264) q[3];
sx q[3];
rz(-2.3271797) q[3];
sx q[3];
rz(-1.502445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.020236882) q[2];
sx q[2];
rz(-1.3573703) q[2];
sx q[2];
rz(-0.21719246) q[2];
rz(1.1413261) q[3];
sx q[3];
rz(-0.3749899) q[3];
sx q[3];
rz(-2.2264437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7104257) q[0];
sx q[0];
rz(-1.2942261) q[0];
sx q[0];
rz(-3.1173832) q[0];
rz(2.0503893) q[1];
sx q[1];
rz(-2.0228491) q[1];
sx q[1];
rz(-0.81407636) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0776802) q[0];
sx q[0];
rz(-1.6565431) q[0];
sx q[0];
rz(2.4299939) q[0];
rz(-0.9466968) q[2];
sx q[2];
rz(-1.1128508) q[2];
sx q[2];
rz(0.68489546) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.947727) q[1];
sx q[1];
rz(-1.5720688) q[1];
sx q[1];
rz(2.3296859) q[1];
rz(1.3423479) q[3];
sx q[3];
rz(-1.486393) q[3];
sx q[3];
rz(0.24385246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7476864) q[2];
sx q[2];
rz(-2.2795491) q[2];
sx q[2];
rz(-1.6115335) q[2];
rz(2.090442) q[3];
sx q[3];
rz(-1.3408778) q[3];
sx q[3];
rz(-3.0401201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.53973389) q[0];
sx q[0];
rz(-2.1506385) q[0];
sx q[0];
rz(-0.47716004) q[0];
rz(1.5513783) q[1];
sx q[1];
rz(-1.4789707) q[1];
sx q[1];
rz(1.307055) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4068749) q[0];
sx q[0];
rz(-2.2860323) q[0];
sx q[0];
rz(0.060013219) q[0];
rz(2.3315694) q[2];
sx q[2];
rz(-0.56213435) q[2];
sx q[2];
rz(0.13371828) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2714398) q[1];
sx q[1];
rz(-2.296827) q[1];
sx q[1];
rz(2.4072263) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96316506) q[3];
sx q[3];
rz(-1.213671) q[3];
sx q[3];
rz(-0.70588338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9755379) q[2];
sx q[2];
rz(-2.1259978) q[2];
sx q[2];
rz(0.75237742) q[2];
rz(0.11263975) q[3];
sx q[3];
rz(-2.967716) q[3];
sx q[3];
rz(1.1627452) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024121506) q[0];
sx q[0];
rz(-1.7565256) q[0];
sx q[0];
rz(1.1275445) q[0];
rz(-0.37359259) q[1];
sx q[1];
rz(-1.5126546) q[1];
sx q[1];
rz(1.0027813) q[1];
rz(1.9890979) q[2];
sx q[2];
rz(-1.9807182) q[2];
sx q[2];
rz(2.9302927) q[2];
rz(-0.74736377) q[3];
sx q[3];
rz(-1.4229645) q[3];
sx q[3];
rz(0.67887282) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
