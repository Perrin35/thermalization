OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.81808972) q[0];
sx q[0];
rz(-0.40677318) q[0];
sx q[0];
rz(0.50954252) q[0];
rz(-0.31515631) q[1];
sx q[1];
rz(-2.9480313) q[1];
sx q[1];
rz(1.0055746) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5339342) q[0];
sx q[0];
rz(-2.9276121) q[0];
sx q[0];
rz(0.35983742) q[0];
x q[1];
rz(2.2174913) q[2];
sx q[2];
rz(-1.8515203) q[2];
sx q[2];
rz(-1.678987) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.79709681) q[1];
sx q[1];
rz(-2.2386947) q[1];
sx q[1];
rz(0.75842661) q[1];
rz(-pi) q[2];
rz(2.3629684) q[3];
sx q[3];
rz(-2.0247685) q[3];
sx q[3];
rz(3.0954645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84844184) q[2];
sx q[2];
rz(-1.2520496) q[2];
sx q[2];
rz(-2.0330009) q[2];
rz(2.0075924) q[3];
sx q[3];
rz(-2.0847376) q[3];
sx q[3];
rz(-0.081324287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9950681) q[0];
sx q[0];
rz(-2.0159371) q[0];
sx q[0];
rz(-2.8296237) q[0];
rz(-0.23090714) q[1];
sx q[1];
rz(-2.0377908) q[1];
sx q[1];
rz(-1.1280967) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9474079) q[0];
sx q[0];
rz(-1.1028521) q[0];
sx q[0];
rz(0.2121966) q[0];
x q[1];
rz(2.6083469) q[2];
sx q[2];
rz(-0.74138481) q[2];
sx q[2];
rz(0.35517755) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9018394) q[1];
sx q[1];
rz(-1.1247207) q[1];
sx q[1];
rz(-1.447999) q[1];
x q[2];
rz(2.3771068) q[3];
sx q[3];
rz(-1.2448881) q[3];
sx q[3];
rz(-0.21316646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4565304) q[2];
sx q[2];
rz(-1.637746) q[2];
sx q[2];
rz(0.064229639) q[2];
rz(1.8188933) q[3];
sx q[3];
rz(-0.6367681) q[3];
sx q[3];
rz(0.88671154) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5927758) q[0];
sx q[0];
rz(-0.23574695) q[0];
sx q[0];
rz(-1.2491666) q[0];
rz(-0.33117548) q[1];
sx q[1];
rz(-1.1391897) q[1];
sx q[1];
rz(1.9452728) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6702658) q[0];
sx q[0];
rz(-1.571079) q[0];
sx q[0];
rz(3.1413881) q[0];
x q[1];
rz(0.24748374) q[2];
sx q[2];
rz(-1.5759517) q[2];
sx q[2];
rz(-0.86901074) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1332449) q[1];
sx q[1];
rz(-0.82723325) q[1];
sx q[1];
rz(0.72669795) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9587956) q[3];
sx q[3];
rz(-0.85577589) q[3];
sx q[3];
rz(-1.9388461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5195878) q[2];
sx q[2];
rz(-0.76498166) q[2];
sx q[2];
rz(0.91274846) q[2];
rz(-1.7031472) q[3];
sx q[3];
rz(-1.5104975) q[3];
sx q[3];
rz(1.8057757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66252935) q[0];
sx q[0];
rz(-2.0344069) q[0];
sx q[0];
rz(-2.4712439) q[0];
rz(0.66728512) q[1];
sx q[1];
rz(-2.1332462) q[1];
sx q[1];
rz(2.537312) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1924698) q[0];
sx q[0];
rz(-1.6141103) q[0];
sx q[0];
rz(-0.73661026) q[0];
rz(-pi) q[1];
rz(0.38197869) q[2];
sx q[2];
rz(-1.4264709) q[2];
sx q[2];
rz(-1.195418) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6091023) q[1];
sx q[1];
rz(-0.33342842) q[1];
sx q[1];
rz(-0.87581265) q[1];
rz(1.3210345) q[3];
sx q[3];
rz(-2.3951924) q[3];
sx q[3];
rz(1.6370893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0303354) q[2];
sx q[2];
rz(-1.5680485) q[2];
sx q[2];
rz(-0.37720171) q[2];
rz(0.12353573) q[3];
sx q[3];
rz(-1.7081407) q[3];
sx q[3];
rz(1.7326573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94443026) q[0];
sx q[0];
rz(-2.5640709) q[0];
sx q[0];
rz(-2.8422728) q[0];
rz(-2.0206644) q[1];
sx q[1];
rz(-1.9363554) q[1];
sx q[1];
rz(2.1790806) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4967921) q[0];
sx q[0];
rz(-1.776665) q[0];
sx q[0];
rz(-1.0394761) q[0];
rz(2.148284) q[2];
sx q[2];
rz(-0.74314144) q[2];
sx q[2];
rz(0.20692839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.968041) q[1];
sx q[1];
rz(-2.3690201) q[1];
sx q[1];
rz(2.3501758) q[1];
rz(-pi) q[2];
rz(-2.4351746) q[3];
sx q[3];
rz(-2.4101541) q[3];
sx q[3];
rz(-2.4789435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0943429) q[2];
sx q[2];
rz(-0.80545682) q[2];
sx q[2];
rz(-0.94364014) q[2];
rz(2.0761679) q[3];
sx q[3];
rz(-1.3506972) q[3];
sx q[3];
rz(-2.0431199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0672673) q[0];
sx q[0];
rz(-1.4370947) q[0];
sx q[0];
rz(-3.1332916) q[0];
rz(0.51757327) q[1];
sx q[1];
rz(-2.4868496) q[1];
sx q[1];
rz(1.5932721) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7784938) q[0];
sx q[0];
rz(-0.35176793) q[0];
sx q[0];
rz(-1.5009053) q[0];
rz(-1.8528884) q[2];
sx q[2];
rz(-2.05915) q[2];
sx q[2];
rz(2.2887678) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.949985) q[1];
sx q[1];
rz(-1.9298501) q[1];
sx q[1];
rz(2.784347) q[1];
rz(-pi) q[2];
rz(1.0128138) q[3];
sx q[3];
rz(-2.4468385) q[3];
sx q[3];
rz(-1.5149124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8325309) q[2];
sx q[2];
rz(-1.6522202) q[2];
sx q[2];
rz(1.8050516) q[2];
rz(3.0858223) q[3];
sx q[3];
rz(-2.1499108) q[3];
sx q[3];
rz(2.706004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5830773) q[0];
sx q[0];
rz(-0.4774839) q[0];
sx q[0];
rz(2.4537295) q[0];
rz(1.4631924) q[1];
sx q[1];
rz(-0.72931591) q[1];
sx q[1];
rz(-2.8942143) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1557013) q[0];
sx q[0];
rz(-1.3257605) q[0];
sx q[0];
rz(1.758453) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65785758) q[2];
sx q[2];
rz(-2.0867072) q[2];
sx q[2];
rz(2.6929457) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.069889594) q[1];
sx q[1];
rz(-2.1462198) q[1];
sx q[1];
rz(2.6399122) q[1];
rz(-pi) q[2];
rz(-3.0159146) q[3];
sx q[3];
rz(-1.8018556) q[3];
sx q[3];
rz(-0.30320689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6972202) q[2];
sx q[2];
rz(-1.8070544) q[2];
sx q[2];
rz(-2.7174301) q[2];
rz(-1.0423202) q[3];
sx q[3];
rz(-2.3764231) q[3];
sx q[3];
rz(-0.55646363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1384077) q[0];
sx q[0];
rz(-0.98588949) q[0];
sx q[0];
rz(0.61088046) q[0];
rz(-0.25018397) q[1];
sx q[1];
rz(-1.7411722) q[1];
sx q[1];
rz(2.6649323) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8308423) q[0];
sx q[0];
rz(-1.0888958) q[0];
sx q[0];
rz(3.0974814) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.119461) q[2];
sx q[2];
rz(-0.20344606) q[2];
sx q[2];
rz(0.97285336) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.51591251) q[1];
sx q[1];
rz(-1.1704186) q[1];
sx q[1];
rz(2.6039586) q[1];
rz(-2.1564756) q[3];
sx q[3];
rz(-1.9352311) q[3];
sx q[3];
rz(-1.6062615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1499947) q[2];
sx q[2];
rz(-1.0341045) q[2];
sx q[2];
rz(-2.4617713) q[2];
rz(-0.46197915) q[3];
sx q[3];
rz(-1.446412) q[3];
sx q[3];
rz(1.6405039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3221472) q[0];
sx q[0];
rz(-1.7663904) q[0];
sx q[0];
rz(-2.9875901) q[0];
rz(2.9601861) q[1];
sx q[1];
rz(-2.0712974) q[1];
sx q[1];
rz(3.0214686) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87163602) q[0];
sx q[0];
rz(-1.1678732) q[0];
sx q[0];
rz(-3.0771034) q[0];
x q[1];
rz(0.93379069) q[2];
sx q[2];
rz(-2.8178038) q[2];
sx q[2];
rz(-1.1188913) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.97717211) q[1];
sx q[1];
rz(-1.394632) q[1];
sx q[1];
rz(-0.099343012) q[1];
rz(0.52773962) q[3];
sx q[3];
rz(-1.4459918) q[3];
sx q[3];
rz(-3.0230056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4274365) q[2];
sx q[2];
rz(-1.3849881) q[2];
sx q[2];
rz(-0.45905054) q[2];
rz(-2.6311724) q[3];
sx q[3];
rz(-0.84158689) q[3];
sx q[3];
rz(-2.4418805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089461483) q[0];
sx q[0];
rz(-0.94877807) q[0];
sx q[0];
rz(-2.7556038) q[0];
rz(1.5486859) q[1];
sx q[1];
rz(-2.6451151) q[1];
sx q[1];
rz(-1.5923502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-1.5992461) q[0];
sx q[0];
rz(-1.7398119) q[0];
x q[1];
rz(2.5086918) q[2];
sx q[2];
rz(-2.0582681) q[2];
sx q[2];
rz(0.097214708) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0299227) q[1];
sx q[1];
rz(-0.1055461) q[1];
sx q[1];
rz(-1.8284069) q[1];
rz(-pi) q[2];
rz(1.4172961) q[3];
sx q[3];
rz(-1.5046538) q[3];
sx q[3];
rz(-1.9754174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3647032) q[2];
sx q[2];
rz(-0.93901712) q[2];
sx q[2];
rz(1.4466064) q[2];
rz(2.1052836) q[3];
sx q[3];
rz(-0.88007897) q[3];
sx q[3];
rz(3.115263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.885289) q[0];
sx q[0];
rz(-0.97923179) q[0];
sx q[0];
rz(-1.6543065) q[0];
rz(-0.22008315) q[1];
sx q[1];
rz(-1.1047803) q[1];
sx q[1];
rz(1.6688375) q[1];
rz(1.8840811) q[2];
sx q[2];
rz(-0.85713119) q[2];
sx q[2];
rz(-1.5766889) q[2];
rz(-0.50095273) q[3];
sx q[3];
rz(-1.8617478) q[3];
sx q[3];
rz(-0.72248722) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
