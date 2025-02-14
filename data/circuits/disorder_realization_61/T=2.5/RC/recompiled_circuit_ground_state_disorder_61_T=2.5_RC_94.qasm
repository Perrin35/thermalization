OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9899848) q[0];
sx q[0];
rz(-2.2014872) q[0];
sx q[0];
rz(-0.54036933) q[0];
rz(-2.6481533) q[1];
sx q[1];
rz(-0.72055888) q[1];
sx q[1];
rz(0.61385733) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6709124) q[0];
sx q[0];
rz(-1.7550091) q[0];
sx q[0];
rz(-0.6349011) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0059284728) q[2];
sx q[2];
rz(-1.2790171) q[2];
sx q[2];
rz(2.8242982) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2062224) q[1];
sx q[1];
rz(-0.79806483) q[1];
sx q[1];
rz(-0.3978637) q[1];
x q[2];
rz(-1.5491539) q[3];
sx q[3];
rz(-1.3236681) q[3];
sx q[3];
rz(2.6789042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2962239) q[2];
sx q[2];
rz(-1.2426528) q[2];
sx q[2];
rz(0.63429147) q[2];
rz(-2.2058709) q[3];
sx q[3];
rz(-2.8959385) q[3];
sx q[3];
rz(0.74265695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3212386) q[0];
sx q[0];
rz(-1.5411493) q[0];
sx q[0];
rz(0.4441922) q[0];
rz(0.76002899) q[1];
sx q[1];
rz(-2.0030231) q[1];
sx q[1];
rz(0.98145032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0652577) q[0];
sx q[0];
rz(-1.0629356) q[0];
sx q[0];
rz(1.9064376) q[0];
x q[1];
rz(2.9600328) q[2];
sx q[2];
rz(-1.4019792) q[2];
sx q[2];
rz(-0.09240514) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.76933512) q[1];
sx q[1];
rz(-0.44722873) q[1];
sx q[1];
rz(-2.9370537) q[1];
x q[2];
rz(-1.3055824) q[3];
sx q[3];
rz(-1.9899559) q[3];
sx q[3];
rz(-0.1923628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0280219) q[2];
sx q[2];
rz(-2.5271723) q[2];
sx q[2];
rz(-1.9363972) q[2];
rz(-2.6049854) q[3];
sx q[3];
rz(-1.2380995) q[3];
sx q[3];
rz(1.7369778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9179012) q[0];
sx q[0];
rz(-1.2811998) q[0];
sx q[0];
rz(-0.83876383) q[0];
rz(-0.63610786) q[1];
sx q[1];
rz(-1.4898224) q[1];
sx q[1];
rz(3.1210693) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9560709) q[0];
sx q[0];
rz(-0.98617109) q[0];
sx q[0];
rz(-2.4654075) q[0];
rz(-2.0383561) q[2];
sx q[2];
rz(-2.1439432) q[2];
sx q[2];
rz(-1.3064885) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.62333362) q[1];
sx q[1];
rz(-1.1578528) q[1];
sx q[1];
rz(-1.5790126) q[1];
rz(-pi) q[2];
rz(-2.0666615) q[3];
sx q[3];
rz(-0.32204667) q[3];
sx q[3];
rz(-2.5888458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8013578) q[2];
sx q[2];
rz(-1.6131718) q[2];
sx q[2];
rz(0.94432962) q[2];
rz(-1.0229735) q[3];
sx q[3];
rz(-1.2228271) q[3];
sx q[3];
rz(2.003722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2706547) q[0];
sx q[0];
rz(-1.9671257) q[0];
sx q[0];
rz(1.9842072) q[0];
rz(1.4986787) q[1];
sx q[1];
rz(-1.6694371) q[1];
sx q[1];
rz(1.1882943) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.37876) q[0];
sx q[0];
rz(-0.68183696) q[0];
sx q[0];
rz(-2.3781611) q[0];
rz(1.6158197) q[2];
sx q[2];
rz(-2.1189711) q[2];
sx q[2];
rz(0.5612824) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3664361) q[1];
sx q[1];
rz(-0.85298733) q[1];
sx q[1];
rz(-0.78182533) q[1];
x q[2];
rz(2.9784936) q[3];
sx q[3];
rz(-1.2360611) q[3];
sx q[3];
rz(1.3770449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0231861) q[2];
sx q[2];
rz(-1.1772757) q[2];
sx q[2];
rz(-2.4510621) q[2];
rz(1.1150507) q[3];
sx q[3];
rz(-0.77459049) q[3];
sx q[3];
rz(1.640865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0951776) q[0];
sx q[0];
rz(-1.4354118) q[0];
sx q[0];
rz(1.0272367) q[0];
rz(-1.3014303) q[1];
sx q[1];
rz(-0.65168989) q[1];
sx q[1];
rz(2.8177736) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0831858) q[0];
sx q[0];
rz(-2.2290285) q[0];
sx q[0];
rz(0.51256521) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1236046) q[2];
sx q[2];
rz(-1.5587285) q[2];
sx q[2];
rz(-0.34876212) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.71621543) q[1];
sx q[1];
rz(-1.0310804) q[1];
sx q[1];
rz(-1.1951237) q[1];
rz(-1.845076) q[3];
sx q[3];
rz(-0.52805942) q[3];
sx q[3];
rz(-1.5825281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7258437) q[2];
sx q[2];
rz(-0.21012935) q[2];
sx q[2];
rz(-2.6591163) q[2];
rz(1.1680565) q[3];
sx q[3];
rz(-1.9246293) q[3];
sx q[3];
rz(2.4647958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93550682) q[0];
sx q[0];
rz(-2.2914903) q[0];
sx q[0];
rz(-2.2555943) q[0];
rz(0.63367263) q[1];
sx q[1];
rz(-0.88992563) q[1];
sx q[1];
rz(0.69127965) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.752744) q[0];
sx q[0];
rz(-0.92557478) q[0];
sx q[0];
rz(-1.6165125) q[0];
rz(1.5070314) q[2];
sx q[2];
rz(-0.87997961) q[2];
sx q[2];
rz(-1.7036167) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4659652) q[1];
sx q[1];
rz(-1.6576997) q[1];
sx q[1];
rz(2.3244832) q[1];
x q[2];
rz(-2.3522934) q[3];
sx q[3];
rz(-1.3914445) q[3];
sx q[3];
rz(-1.2937677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1401691) q[2];
sx q[2];
rz(-2.3052577) q[2];
sx q[2];
rz(0.26958618) q[2];
rz(-2.1932898) q[3];
sx q[3];
rz(-1.6092665) q[3];
sx q[3];
rz(1.3986826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.36169323) q[0];
sx q[0];
rz(-0.43710709) q[0];
sx q[0];
rz(-1.0070739) q[0];
rz(-0.40564793) q[1];
sx q[1];
rz(-0.59527731) q[1];
sx q[1];
rz(2.5880623) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8419475) q[0];
sx q[0];
rz(-1.8252234) q[0];
sx q[0];
rz(-0.05924745) q[0];
rz(-pi) q[1];
rz(2.5895025) q[2];
sx q[2];
rz(-2.0497397) q[2];
sx q[2];
rz(2.8450784) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0515159) q[1];
sx q[1];
rz(-1.1184268) q[1];
sx q[1];
rz(0.69454792) q[1];
rz(-pi) q[2];
x q[2];
rz(0.062640142) q[3];
sx q[3];
rz(-2.7609112) q[3];
sx q[3];
rz(-2.9306987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3125399) q[2];
sx q[2];
rz(-1.092814) q[2];
sx q[2];
rz(1.9276169) q[2];
rz(-2.35516) q[3];
sx q[3];
rz(-1.7460456) q[3];
sx q[3];
rz(0.023155183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19602747) q[0];
sx q[0];
rz(-0.29226154) q[0];
sx q[0];
rz(1.8096402) q[0];
rz(0.61344433) q[1];
sx q[1];
rz(-2.1211801) q[1];
sx q[1];
rz(0.44949284) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1163131) q[0];
sx q[0];
rz(-2.6465102) q[0];
sx q[0];
rz(-1.3950924) q[0];
rz(1.8138936) q[2];
sx q[2];
rz(-0.34533325) q[2];
sx q[2];
rz(-1.9497046) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.66287884) q[1];
sx q[1];
rz(-0.7906853) q[1];
sx q[1];
rz(1.3776758) q[1];
rz(-1.6291796) q[3];
sx q[3];
rz(-1.8280142) q[3];
sx q[3];
rz(-1.2793878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.31676644) q[2];
sx q[2];
rz(-0.68244857) q[2];
sx q[2];
rz(2.5332434) q[2];
rz(-1.1634722) q[3];
sx q[3];
rz(-1.3694265) q[3];
sx q[3];
rz(2.5782862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1083531) q[0];
sx q[0];
rz(-0.9110564) q[0];
sx q[0];
rz(-0.60229993) q[0];
rz(-0.96744084) q[1];
sx q[1];
rz(-0.93092218) q[1];
sx q[1];
rz(1.1036576) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3924343) q[0];
sx q[0];
rz(-1.9604248) q[0];
sx q[0];
rz(0.97831877) q[0];
rz(-2.5050312) q[2];
sx q[2];
rz(-1.5590073) q[2];
sx q[2];
rz(-1.7280886) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1638318) q[1];
sx q[1];
rz(-1.1196616) q[1];
sx q[1];
rz(-3.1402863) q[1];
x q[2];
rz(-1.4673442) q[3];
sx q[3];
rz(-0.58052968) q[3];
sx q[3];
rz(-0.79642297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.935282) q[2];
sx q[2];
rz(-1.3056359) q[2];
sx q[2];
rz(-0.3375816) q[2];
rz(-0.28389367) q[3];
sx q[3];
rz(-2.4604535) q[3];
sx q[3];
rz(-2.9547227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.5481446) q[0];
sx q[0];
rz(-1.633506) q[0];
sx q[0];
rz(-3.0737851) q[0];
rz(-1.9920805) q[1];
sx q[1];
rz(-1.5510473) q[1];
sx q[1];
rz(-2.6670719) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15856811) q[0];
sx q[0];
rz(-1.8095542) q[0];
sx q[0];
rz(-0.013834133) q[0];
rz(-pi) q[1];
rz(1.7065918) q[2];
sx q[2];
rz(-1.3320005) q[2];
sx q[2];
rz(0.39993024) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1637516) q[1];
sx q[1];
rz(-3.0545898) q[1];
sx q[1];
rz(-1.6496501) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69863221) q[3];
sx q[3];
rz(-2.1967874) q[3];
sx q[3];
rz(-2.0076942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6600251) q[2];
sx q[2];
rz(-2.2109172) q[2];
sx q[2];
rz(2.068326) q[2];
rz(2.977071) q[3];
sx q[3];
rz(-1.8142895) q[3];
sx q[3];
rz(2.8209414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58707033) q[0];
sx q[0];
rz(-1.5759435) q[0];
sx q[0];
rz(-1.6389621) q[0];
rz(-1.5837689) q[1];
sx q[1];
rz(-2.0638034) q[1];
sx q[1];
rz(2.6149909) q[1];
rz(-0.64939349) q[2];
sx q[2];
rz(-0.55772256) q[2];
sx q[2];
rz(-1.100308) q[2];
rz(0.14231331) q[3];
sx q[3];
rz(-1.6518136) q[3];
sx q[3];
rz(1.470737) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
