OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8060057) q[0];
sx q[0];
rz(-0.94526362) q[0];
sx q[0];
rz(2.6160016) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(-0.90484172) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68569505) q[0];
sx q[0];
rz(-2.064961) q[0];
sx q[0];
rz(2.6463638) q[0];
rz(-pi) q[1];
rz(0.22613871) q[2];
sx q[2];
rz(-1.3790501) q[2];
sx q[2];
rz(2.6297671) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9368254) q[1];
sx q[1];
rz(-1.966404) q[1];
sx q[1];
rz(-2.8863465) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0337898) q[3];
sx q[3];
rz(-2.783943) q[3];
sx q[3];
rz(1.3057749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9378172) q[2];
sx q[2];
rz(-1.4062466) q[2];
sx q[2];
rz(-3.0453483) q[2];
rz(-2.105666) q[3];
sx q[3];
rz(-0.38714287) q[3];
sx q[3];
rz(2.9878785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.7006943) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(-2.3764215) q[0];
rz(-1.8493429) q[1];
sx q[1];
rz(-2.6563829) q[1];
sx q[1];
rz(-0.66295019) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6843296) q[0];
sx q[0];
rz(-1.508679) q[0];
sx q[0];
rz(1.5044042) q[0];
x q[1];
rz(1.0639973) q[2];
sx q[2];
rz(-2.4231744) q[2];
sx q[2];
rz(1.4425299) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.502683) q[1];
sx q[1];
rz(-2.0109634) q[1];
sx q[1];
rz(-1.6619976) q[1];
x q[2];
rz(0.29942056) q[3];
sx q[3];
rz(-0.50588183) q[3];
sx q[3];
rz(-2.8477856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.092676) q[2];
sx q[2];
rz(-2.0930223) q[2];
sx q[2];
rz(-2.8125787) q[2];
rz(-0.66550231) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(-1.3177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.5927785) q[0];
sx q[0];
rz(-0.22286649) q[0];
sx q[0];
rz(0.22234017) q[0];
rz(2.1242583) q[1];
sx q[1];
rz(-2.4203114) q[1];
sx q[1];
rz(-0.51868784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9496574) q[0];
sx q[0];
rz(-1.9266832) q[0];
sx q[0];
rz(-2.3200672) q[0];
rz(-pi) q[1];
rz(-0.33883314) q[2];
sx q[2];
rz(-0.28140861) q[2];
sx q[2];
rz(-2.0237405) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.31836244) q[1];
sx q[1];
rz(-1.3981817) q[1];
sx q[1];
rz(0.77962064) q[1];
x q[2];
rz(-3.1268901) q[3];
sx q[3];
rz(-3.0869752) q[3];
sx q[3];
rz(-2.2811449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2805933) q[2];
sx q[2];
rz(-0.36182797) q[2];
sx q[2];
rz(-0.068543531) q[2];
rz(-2.5391501) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(-3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.65072) q[0];
sx q[0];
rz(-0.77712178) q[0];
sx q[0];
rz(2.9673476) q[0];
rz(0.53025591) q[1];
sx q[1];
rz(-1.6590051) q[1];
sx q[1];
rz(0.51309103) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.355277) q[0];
sx q[0];
rz(-2.516541) q[0];
sx q[0];
rz(-3.0062208) q[0];
rz(-2.1999173) q[2];
sx q[2];
rz(-0.86949124) q[2];
sx q[2];
rz(-0.47052449) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0780371) q[1];
sx q[1];
rz(-0.47680285) q[1];
sx q[1];
rz(-0.30492353) q[1];
x q[2];
rz(1.7771878) q[3];
sx q[3];
rz(-1.5406113) q[3];
sx q[3];
rz(0.15011945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47485581) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(-0.22988698) q[2];
rz(-0.41904467) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(2.6823147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4693562) q[0];
sx q[0];
rz(-2.3903963) q[0];
sx q[0];
rz(2.4647734) q[0];
rz(0.49304402) q[1];
sx q[1];
rz(-0.9489916) q[1];
sx q[1];
rz(-0.61606032) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84232932) q[0];
sx q[0];
rz(-0.063532524) q[0];
sx q[0];
rz(-2.1701943) q[0];
rz(2.8370503) q[2];
sx q[2];
rz(-0.38123044) q[2];
sx q[2];
rz(-2.9074557) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7580326) q[1];
sx q[1];
rz(-0.57764232) q[1];
sx q[1];
rz(-0.010096117) q[1];
rz(0.89415278) q[3];
sx q[3];
rz(-1.812462) q[3];
sx q[3];
rz(-1.579293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.73413509) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(0.25536728) q[2];
rz(1.6051965) q[3];
sx q[3];
rz(-2.1891749) q[3];
sx q[3];
rz(0.77409625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(2.7191294) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(2.6690924) q[0];
rz(2.6155112) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(2.2568259) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038267604) q[0];
sx q[0];
rz(-1.190435) q[0];
sx q[0];
rz(-0.079770712) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0753291) q[2];
sx q[2];
rz(-2.931086) q[2];
sx q[2];
rz(-1.2103684) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5342418) q[1];
sx q[1];
rz(-0.51041767) q[1];
sx q[1];
rz(-2.1900602) q[1];
rz(-1.1690508) q[3];
sx q[3];
rz(-0.570795) q[3];
sx q[3];
rz(1.4626383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.74449599) q[2];
sx q[2];
rz(-2.3366191) q[2];
sx q[2];
rz(2.8302144) q[2];
rz(-1.3686251) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(-2.6122724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41806528) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(-3.080522) q[0];
rz(-3.1014077) q[1];
sx q[1];
rz(-1.9804852) q[1];
sx q[1];
rz(0.73289245) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.830025) q[0];
sx q[0];
rz(-2.4550779) q[0];
sx q[0];
rz(-0.65450432) q[0];
rz(0.86874666) q[2];
sx q[2];
rz(-1.947543) q[2];
sx q[2];
rz(1.6422611) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1317132) q[1];
sx q[1];
rz(-0.94505802) q[1];
sx q[1];
rz(3.1335319) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2039127) q[3];
sx q[3];
rz(-3.0266579) q[3];
sx q[3];
rz(0.68655187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0682893) q[2];
sx q[2];
rz(-1.1969593) q[2];
sx q[2];
rz(0.51458365) q[2];
rz(1.2375281) q[3];
sx q[3];
rz(-2.5585744) q[3];
sx q[3];
rz(2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.085389) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(-2.432166) q[0];
rz(-1.6363232) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(0.27871305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73662169) q[0];
sx q[0];
rz(-1.706089) q[0];
sx q[0];
rz(-0.21185974) q[0];
rz(-pi) q[1];
rz(2.7630373) q[2];
sx q[2];
rz(-0.59213973) q[2];
sx q[2];
rz(0.90422599) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4019926) q[1];
sx q[1];
rz(-2.7584689) q[1];
sx q[1];
rz(-2.1869704) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7940984) q[3];
sx q[3];
rz(-1.4947065) q[3];
sx q[3];
rz(3.1267816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8401106) q[2];
sx q[2];
rz(-1.8436517) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(-2.2864443) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(-0.40518951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1466325) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(0.96889281) q[0];
rz(-2.6682207) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(1.999058) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0363077) q[0];
sx q[0];
rz(-1.4799911) q[0];
sx q[0];
rz(-1.275195) q[0];
x q[1];
rz(-1.5683453) q[2];
sx q[2];
rz(-0.8674538) q[2];
sx q[2];
rz(3.0424812) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0513623) q[1];
sx q[1];
rz(-1.4453332) q[1];
sx q[1];
rz(-0.54162234) q[1];
rz(-pi) q[2];
rz(-2.9386018) q[3];
sx q[3];
rz(-0.91068017) q[3];
sx q[3];
rz(0.70860329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2450976) q[2];
sx q[2];
rz(-1.921804) q[2];
sx q[2];
rz(-2.1208105) q[2];
rz(0.3237237) q[3];
sx q[3];
rz(-2.3886069) q[3];
sx q[3];
rz(-3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1616515) q[0];
sx q[0];
rz(-3.1136944) q[0];
sx q[0];
rz(0.7014057) q[0];
rz(0.91570634) q[1];
sx q[1];
rz(-1.0083102) q[1];
sx q[1];
rz(1.9030301) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68574821) q[0];
sx q[0];
rz(-2.0861174) q[0];
sx q[0];
rz(-0.2483764) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3886016) q[2];
sx q[2];
rz(-1.6082616) q[2];
sx q[2];
rz(2.3245036) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.30675948) q[1];
sx q[1];
rz(-0.66654897) q[1];
sx q[1];
rz(-1.761761) q[1];
x q[2];
rz(3.0276887) q[3];
sx q[3];
rz(-2.1602727) q[3];
sx q[3];
rz(-1.1978428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4548268) q[2];
sx q[2];
rz(-1.0964311) q[2];
sx q[2];
rz(0.043838538) q[2];
rz(-1.2010126) q[3];
sx q[3];
rz(-2.4062556) q[3];
sx q[3];
rz(-1.003585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6702406) q[0];
sx q[0];
rz(-0.72605194) q[0];
sx q[0];
rz(-1.3656021) q[0];
rz(0.025370601) q[1];
sx q[1];
rz(-1.8352958) q[1];
sx q[1];
rz(-1.8713554) q[1];
rz(-2.8364137) q[2];
sx q[2];
rz(-1.9532433) q[2];
sx q[2];
rz(0.80079186) q[2];
rz(1.3958037) q[3];
sx q[3];
rz(-2.3621515) q[3];
sx q[3];
rz(3.0099517) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];