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
rz(-7.3047819) q[1];
sx q[1];
rz(2.8586913) q[1];
sx q[1];
rz(18.999264) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5877085) q[0];
sx q[0];
rz(-0.93548488) q[0];
sx q[0];
rz(0.61650886) q[0];
rz(-2.6354191) q[2];
sx q[2];
rz(-2.0548327) q[2];
sx q[2];
rz(1.6137705) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.183061) q[1];
sx q[1];
rz(-1.6089988) q[1];
sx q[1];
rz(-1.7512291) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9505694) q[3];
sx q[3];
rz(-1.028562) q[3];
sx q[3];
rz(0.99381002) q[3];
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
rz(-0.70409888) q[2];
rz(-0.95300931) q[3];
sx q[3];
rz(-0.97218958) q[3];
sx q[3];
rz(1.7378418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54884058) q[0];
sx q[0];
rz(-1.623818) q[0];
sx q[0];
rz(2.5090704) q[0];
rz(0.44644341) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(2.4893563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95603847) q[0];
sx q[0];
rz(-1.5421876) q[0];
sx q[0];
rz(1.5831328) q[0];
rz(-0.89276887) q[2];
sx q[2];
rz(-0.36913482) q[2];
sx q[2];
rz(0.54112753) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.10416874) q[1];
sx q[1];
rz(-1.2836604) q[1];
sx q[1];
rz(1.3377405) q[1];
x q[2];
rz(-0.90918031) q[3];
sx q[3];
rz(-2.069807) q[3];
sx q[3];
rz(-0.0024851174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5941045) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(1.1616421) q[2];
rz(-1.9836327) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(-1.4512216) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0911672) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(2.2913349) q[0];
rz(0.49750528) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(1.3495548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3916546) q[0];
sx q[0];
rz(-0.91021252) q[0];
sx q[0];
rz(-1.8750989) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6367958) q[2];
sx q[2];
rz(-1.5335576) q[2];
sx q[2];
rz(0.50022349) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2022484) q[1];
sx q[1];
rz(-1.0679686) q[1];
sx q[1];
rz(-1.8400251) q[1];
rz(1.9968824) q[3];
sx q[3];
rz(-1.5797371) q[3];
sx q[3];
rz(-0.63250354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0321908) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(-2.2375977) q[2];
rz(-0.30113014) q[3];
sx q[3];
rz(-1.3826933) q[3];
sx q[3];
rz(1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6501453) q[0];
sx q[0];
rz(-0.97147816) q[0];
sx q[0];
rz(-1.4105463) q[0];
rz(-0.63181216) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(-0.036380336) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.13658) q[0];
sx q[0];
rz(-0.67626017) q[0];
sx q[0];
rz(-1.9146634) q[0];
rz(-pi) q[1];
rz(-0.47469791) q[2];
sx q[2];
rz(-2.5050852) q[2];
sx q[2];
rz(2.0091025) q[2];
rz(-pi) q[3];
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
rz(1.6691469) q[3];
sx q[3];
rz(-1.8042759) q[3];
sx q[3];
rz(-1.0055621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.92695421) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(-1.9533763) q[2];
rz(2.4711117) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(-0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7017512) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(2.9751076) q[0];
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
rz(-2.7974632) q[0];
sx q[0];
rz(-1.9088233) q[0];
sx q[0];
rz(2.5812838) q[0];
x q[1];
rz(-3.0460065) q[2];
sx q[2];
rz(-1.0786973) q[2];
sx q[2];
rz(-1.0620067) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.078951051) q[3];
sx q[3];
rz(-1.2478932) q[3];
sx q[3];
rz(2.3372646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.70871893) q[2];
sx q[2];
rz(-1.9659698) q[2];
sx q[2];
rz(0.22053545) q[2];
rz(0.43705127) q[3];
sx q[3];
rz(-1.022499) q[3];
sx q[3];
rz(-2.3760858) q[3];
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
rz(-2.7261312) q[0];
sx q[0];
rz(-2.8227865) q[0];
sx q[0];
rz(-2.3244526) q[0];
rz(-0.56610402) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(1.1436499) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9691539) q[0];
sx q[0];
rz(-1.5307431) q[0];
sx q[0];
rz(0.37102951) q[0];
rz(-pi) q[1];
rz(1.8311062) q[2];
sx q[2];
rz(-2.227265) q[2];
sx q[2];
rz(-1.3442163) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3559349) q[1];
sx q[1];
rz(-0.4520843) q[1];
sx q[1];
rz(-1.7788586) q[1];
rz(-pi) q[2];
rz(2.6236344) q[3];
sx q[3];
rz(-1.8149788) q[3];
sx q[3];
rz(1.1740008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3879261) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(0.4894408) q[2];
rz(0.22805452) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(2.7155546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5722826) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(-1.2868767) q[0];
rz(-0.6634179) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(-1.2333262) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7466) q[0];
sx q[0];
rz(-1.1867503) q[0];
sx q[0];
rz(0.011944255) q[0];
rz(-2.5160518) q[2];
sx q[2];
rz(-0.95226804) q[2];
sx q[2];
rz(-0.90460888) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7066321) q[1];
sx q[1];
rz(-3.0295277) q[1];
sx q[1];
rz(3.0155165) q[1];
rz(-pi) q[2];
rz(1.2193905) q[3];
sx q[3];
rz(-1.1099585) q[3];
sx q[3];
rz(-1.8691065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.037501637) q[2];
sx q[2];
rz(-0.23510322) q[2];
sx q[2];
rz(1.0160149) q[2];
rz(-3.0715023) q[3];
sx q[3];
rz(-1.2066119) q[3];
sx q[3];
rz(-2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12524097) q[0];
sx q[0];
rz(-2.4588983) q[0];
sx q[0];
rz(1.4455147) q[0];
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
rz(-0.94851516) q[0];
sx q[0];
rz(-0.45398871) q[0];
sx q[0];
rz(1.2598739) q[0];
rz(-pi) q[1];
rz(-2.3636742) q[2];
sx q[2];
rz(-0.81233835) q[2];
sx q[2];
rz(1.1493491) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.313169) q[1];
sx q[1];
rz(-1.0272044) q[1];
sx q[1];
rz(-1.4491175) q[1];
rz(2.2860252) q[3];
sx q[3];
rz(-1.9868402) q[3];
sx q[3];
rz(-0.95203979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4593279) q[2];
sx q[2];
rz(-0.1846281) q[2];
sx q[2];
rz(0.56274596) q[2];
rz(-0.19966666) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(-1.1220042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(-2.8905706) q[0];
rz(2.714278) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(0.02773157) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77596091) q[0];
sx q[0];
rz(-2.2790475) q[0];
sx q[0];
rz(-2.2328949) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2698195) q[2];
sx q[2];
rz(-1.36424) q[2];
sx q[2];
rz(-2.2733462) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7293538) q[1];
sx q[1];
rz(-2.5322399) q[1];
sx q[1];
rz(2.9669697) q[1];
rz(-pi) q[2];
rz(-3.1157007) q[3];
sx q[3];
rz(-1.504244) q[3];
sx q[3];
rz(-0.47749146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1256844) q[2];
sx q[2];
rz(-0.97970825) q[2];
sx q[2];
rz(-1.998418) q[2];
rz(-0.14287359) q[3];
sx q[3];
rz(-1.5121127) q[3];
sx q[3];
rz(-2.2843602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3906355) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(0.45387682) q[0];
rz(2.4699396) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(0.25751105) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8513545) q[0];
sx q[0];
rz(-2.5332753) q[0];
sx q[0];
rz(2.4277707) q[0];
rz(-pi) q[1];
rz(-1.0742513) q[2];
sx q[2];
rz(-1.0814582) q[2];
sx q[2];
rz(0.41000965) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.232302) q[1];
sx q[1];
rz(-1.5850987) q[1];
sx q[1];
rz(2.4453352) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8994843) q[3];
sx q[3];
rz(-2.0383516) q[3];
sx q[3];
rz(-2.9374591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3283078) q[2];
sx q[2];
rz(-1.4537145) q[2];
sx q[2];
rz(0.60662398) q[2];
rz(-2.666752) q[3];
sx q[3];
rz(-0.94687051) q[3];
sx q[3];
rz(-2.301208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.3474779) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(-2.9150302) q[1];
sx q[1];
rz(-1.4145874) q[1];
sx q[1];
rz(0.58691595) q[1];
rz(0.46802855) q[2];
sx q[2];
rz(-1.2266908) q[2];
sx q[2];
rz(-2.5897349) q[2];
rz(-1.0971309) q[3];
sx q[3];
rz(-0.53285014) q[3];
sx q[3];
rz(0.22900029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
