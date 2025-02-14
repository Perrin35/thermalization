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
rz(0.42143917) q[0];
sx q[0];
rz(5.0285664) q[0];
sx q[0];
rz(9.9528735) q[0];
rz(2.7948607) q[1];
sx q[1];
rz(-1.3277227) q[1];
sx q[1];
rz(0.2833856) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91862667) q[0];
sx q[0];
rz(-0.98969995) q[0];
sx q[0];
rz(1.6870935) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6706477) q[2];
sx q[2];
rz(-2.2631553) q[2];
sx q[2];
rz(-0.67519855) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8979075) q[1];
sx q[1];
rz(-2.7449611) q[1];
sx q[1];
rz(2.4325788) q[1];
x q[2];
rz(1.3259726) q[3];
sx q[3];
rz(-2.7283165) q[3];
sx q[3];
rz(0.63140376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1957207) q[2];
sx q[2];
rz(-1.6678145) q[2];
sx q[2];
rz(-2.8742068) q[2];
rz(-0.17313677) q[3];
sx q[3];
rz(-1.1370398) q[3];
sx q[3];
rz(-0.6828298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6807569) q[0];
sx q[0];
rz(-2.2956678) q[0];
sx q[0];
rz(-0.2980921) q[0];
rz(2.8744892) q[1];
sx q[1];
rz(-0.88228455) q[1];
sx q[1];
rz(-2.2411236) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1698477) q[0];
sx q[0];
rz(-1.3296223) q[0];
sx q[0];
rz(1.3154047) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7212333) q[2];
sx q[2];
rz(-1.4263881) q[2];
sx q[2];
rz(-0.23061801) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3648405) q[1];
sx q[1];
rz(-0.89438182) q[1];
sx q[1];
rz(-1.8218356) q[1];
x q[2];
rz(-0.44690981) q[3];
sx q[3];
rz(-0.70820184) q[3];
sx q[3];
rz(0.59659905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.71821) q[2];
sx q[2];
rz(-0.98576468) q[2];
sx q[2];
rz(-0.59744376) q[2];
rz(-1.0159703) q[3];
sx q[3];
rz(-1.6705325) q[3];
sx q[3];
rz(-1.5684675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51282561) q[0];
sx q[0];
rz(-0.91575423) q[0];
sx q[0];
rz(0.065091982) q[0];
rz(1.1029296) q[1];
sx q[1];
rz(-1.256559) q[1];
sx q[1];
rz(0.38806134) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2723615) q[0];
sx q[0];
rz(-0.80565208) q[0];
sx q[0];
rz(2.8153573) q[0];
x q[1];
rz(0.36696649) q[2];
sx q[2];
rz(-1.6410368) q[2];
sx q[2];
rz(2.5052469) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6522576) q[1];
sx q[1];
rz(-1.9804238) q[1];
sx q[1];
rz(2.9698305) q[1];
rz(-1.8592836) q[3];
sx q[3];
rz(-1.7955209) q[3];
sx q[3];
rz(2.2109179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1661561) q[2];
sx q[2];
rz(-1.0338975) q[2];
sx q[2];
rz(-0.088509716) q[2];
rz(-2.4837808) q[3];
sx q[3];
rz(-2.7031873) q[3];
sx q[3];
rz(1.8678166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6863962) q[0];
sx q[0];
rz(-0.96298591) q[0];
sx q[0];
rz(-2.2121867) q[0];
rz(1.9724253) q[1];
sx q[1];
rz(-2.2678352) q[1];
sx q[1];
rz(-0.1263617) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7469661) q[0];
sx q[0];
rz(-1.639059) q[0];
sx q[0];
rz(1.9771876) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54985993) q[2];
sx q[2];
rz(-2.7603995) q[2];
sx q[2];
rz(1.3709244) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11495464) q[1];
sx q[1];
rz(-1.8787787) q[1];
sx q[1];
rz(-0.46197173) q[1];
rz(-pi) q[2];
rz(-2.5764719) q[3];
sx q[3];
rz(-0.45101803) q[3];
sx q[3];
rz(2.7489611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8440642) q[2];
sx q[2];
rz(-1.0928096) q[2];
sx q[2];
rz(-1.7139942) q[2];
rz(1.1045688) q[3];
sx q[3];
rz(-1.5449056) q[3];
sx q[3];
rz(-2.4692718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70668689) q[0];
sx q[0];
rz(-0.66449419) q[0];
sx q[0];
rz(-1.8988761) q[0];
rz(-0.13936123) q[1];
sx q[1];
rz(-0.31529537) q[1];
sx q[1];
rz(-0.79162663) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7032486) q[0];
sx q[0];
rz(-1.7243963) q[0];
sx q[0];
rz(-3.1036882) q[0];
x q[1];
rz(-0.56812079) q[2];
sx q[2];
rz(-1.4213398) q[2];
sx q[2];
rz(-2.7338701) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5593223) q[1];
sx q[1];
rz(-0.30836654) q[1];
sx q[1];
rz(1.4919623) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3982806) q[3];
sx q[3];
rz(-2.5316187) q[3];
sx q[3];
rz(-2.0016746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7255154) q[2];
sx q[2];
rz(-1.5308341) q[2];
sx q[2];
rz(1.7804954) q[2];
rz(2.1040037) q[3];
sx q[3];
rz(-1.4525843) q[3];
sx q[3];
rz(0.008755412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0233651) q[0];
sx q[0];
rz(-2.8814377) q[0];
sx q[0];
rz(-0.17876974) q[0];
rz(-0.44890064) q[1];
sx q[1];
rz(-1.9763549) q[1];
sx q[1];
rz(1.9662205) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7473874) q[0];
sx q[0];
rz(-0.76524599) q[0];
sx q[0];
rz(-3.0470362) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3162546) q[2];
sx q[2];
rz(-0.78331982) q[2];
sx q[2];
rz(2.6412727) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0302802) q[1];
sx q[1];
rz(-2.2759109) q[1];
sx q[1];
rz(0.45781237) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39534335) q[3];
sx q[3];
rz(-0.55231491) q[3];
sx q[3];
rz(2.0022587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7281404) q[2];
sx q[2];
rz(-0.38273013) q[2];
sx q[2];
rz(0.21391301) q[2];
rz(-1.7119857) q[3];
sx q[3];
rz(-1.6702646) q[3];
sx q[3];
rz(-1.8179998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0837309) q[0];
sx q[0];
rz(-1.8566751) q[0];
sx q[0];
rz(0.65823746) q[0];
rz(1.7103051) q[1];
sx q[1];
rz(-0.88600102) q[1];
sx q[1];
rz(-1.5586982) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76097902) q[0];
sx q[0];
rz(-1.3833933) q[0];
sx q[0];
rz(2.6949469) q[0];
rz(-pi) q[1];
rz(1.8020242) q[2];
sx q[2];
rz(-1.2716846) q[2];
sx q[2];
rz(2.0264152) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4728707) q[1];
sx q[1];
rz(-1.8897561) q[1];
sx q[1];
rz(1.215056) q[1];
rz(-1.962108) q[3];
sx q[3];
rz(-2.1150622) q[3];
sx q[3];
rz(-2.5044092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.3805286) q[2];
sx q[2];
rz(-1.5498127) q[2];
sx q[2];
rz(-1.2065678) q[2];
rz(-1.0208463) q[3];
sx q[3];
rz(-2.6267093) q[3];
sx q[3];
rz(-0.90768874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.9970488) q[0];
sx q[0];
rz(-0.21268614) q[0];
sx q[0];
rz(-0.3120684) q[0];
rz(-0.32876217) q[1];
sx q[1];
rz(-1.6203251) q[1];
sx q[1];
rz(-0.75334466) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67344224) q[0];
sx q[0];
rz(-0.99981013) q[0];
sx q[0];
rz(-2.2500413) q[0];
rz(3.1198959) q[2];
sx q[2];
rz(-1.0602131) q[2];
sx q[2];
rz(1.7006602) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4792013) q[1];
sx q[1];
rz(-1.1117142) q[1];
sx q[1];
rz(0.69028141) q[1];
rz(1.0114134) q[3];
sx q[3];
rz(-2.0205343) q[3];
sx q[3];
rz(1.442329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3030777) q[2];
sx q[2];
rz(-1.7531771) q[2];
sx q[2];
rz(0.77084368) q[2];
rz(-2.5924957) q[3];
sx q[3];
rz(-1.2884459) q[3];
sx q[3];
rz(-0.87636605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2751806) q[0];
sx q[0];
rz(-2.488945) q[0];
sx q[0];
rz(-1.8021679) q[0];
rz(-1.132698) q[1];
sx q[1];
rz(-0.39342543) q[1];
sx q[1];
rz(-3.0963617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1605555) q[0];
sx q[0];
rz(-2.0038031) q[0];
sx q[0];
rz(-2.9287522) q[0];
rz(-pi) q[1];
rz(1.8884097) q[2];
sx q[2];
rz(-0.3738974) q[2];
sx q[2];
rz(-1.31391) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.40881904) q[1];
sx q[1];
rz(-1.7013927) q[1];
sx q[1];
rz(-0.25326107) q[1];
rz(-pi) q[2];
rz(1.0186152) q[3];
sx q[3];
rz(-1.3613627) q[3];
sx q[3];
rz(-2.6554573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.905978) q[2];
sx q[2];
rz(-0.9367885) q[2];
sx q[2];
rz(2.3255685) q[2];
rz(-1.7488545) q[3];
sx q[3];
rz(-1.1002898) q[3];
sx q[3];
rz(1.5131153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29869646) q[0];
sx q[0];
rz(-0.5492292) q[0];
sx q[0];
rz(1.5604875) q[0];
rz(2.9807978) q[1];
sx q[1];
rz(-1.1409047) q[1];
sx q[1];
rz(-0.92170942) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35730001) q[0];
sx q[0];
rz(-0.18942115) q[0];
sx q[0];
rz(1.6209391) q[0];
rz(0.99803136) q[2];
sx q[2];
rz(-2.0819527) q[2];
sx q[2];
rz(-0.17705189) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3474011) q[1];
sx q[1];
rz(-0.32116613) q[1];
sx q[1];
rz(-1.0845409) q[1];
rz(-pi) q[2];
rz(1.7245737) q[3];
sx q[3];
rz(-1.3594846) q[3];
sx q[3];
rz(2.7696424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86772743) q[2];
sx q[2];
rz(-0.92550698) q[2];
sx q[2];
rz(1.072139) q[2];
rz(-0.22465651) q[3];
sx q[3];
rz(-1.4916689) q[3];
sx q[3];
rz(2.2255161) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12162019) q[0];
sx q[0];
rz(-2.2157123) q[0];
sx q[0];
rz(-2.8257688) q[0];
rz(-3.0674122) q[1];
sx q[1];
rz(-1.0646432) q[1];
sx q[1];
rz(-0.021312996) q[1];
rz(-1.5431719) q[2];
sx q[2];
rz(-0.78747126) q[2];
sx q[2];
rz(-3.1046545) q[2];
rz(1.9613135) q[3];
sx q[3];
rz(-2.2979679) q[3];
sx q[3];
rz(-1.2818492) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
