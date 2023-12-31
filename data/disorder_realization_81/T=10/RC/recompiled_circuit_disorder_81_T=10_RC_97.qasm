OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73683357) q[0];
sx q[0];
rz(-1.3614549) q[0];
sx q[0];
rz(1.7629495) q[0];
rz(-0.8575851) q[1];
sx q[1];
rz(-1.4839988) q[1];
sx q[1];
rz(0.4508957) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2338976) q[0];
sx q[0];
rz(-1.4820218) q[0];
sx q[0];
rz(-1.5855473) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11159201) q[2];
sx q[2];
rz(-1.0942232) q[2];
sx q[2];
rz(3.1389719) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8494107) q[1];
sx q[1];
rz(-1.0725478) q[1];
sx q[1];
rz(-1.317418) q[1];
rz(-pi) q[2];
rz(1.4536742) q[3];
sx q[3];
rz(-2.4544567) q[3];
sx q[3];
rz(2.5205034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9840055) q[2];
sx q[2];
rz(-1.6820587) q[2];
sx q[2];
rz(-2.297304) q[2];
rz(-2.700581) q[3];
sx q[3];
rz(-0.35566548) q[3];
sx q[3];
rz(-0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5490897) q[0];
sx q[0];
rz(-1.2298158) q[0];
sx q[0];
rz(-2.8785008) q[0];
rz(-2.198055) q[1];
sx q[1];
rz(-0.5967921) q[1];
sx q[1];
rz(1.1862322) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3620136) q[0];
sx q[0];
rz(-1.626726) q[0];
sx q[0];
rz(-1.5895784) q[0];
x q[1];
rz(-2.8430014) q[2];
sx q[2];
rz(-2.9332187) q[2];
sx q[2];
rz(-1.6076455) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8801404) q[1];
sx q[1];
rz(-1.6958691) q[1];
sx q[1];
rz(-2.2373881) q[1];
rz(1.3783185) q[3];
sx q[3];
rz(-0.88497439) q[3];
sx q[3];
rz(1.6845077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1295604) q[2];
sx q[2];
rz(-1.0027145) q[2];
sx q[2];
rz(-1.1594695) q[2];
rz(2.7705079) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(0.31093591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3804669) q[0];
sx q[0];
rz(-1.1300056) q[0];
sx q[0];
rz(2.3348715) q[0];
rz(0.21356788) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(-0.82021964) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7263111) q[0];
sx q[0];
rz(-1.4584686) q[0];
sx q[0];
rz(-2.2583654) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9633031) q[2];
sx q[2];
rz(-2.6143392) q[2];
sx q[2];
rz(0.96780992) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0318451) q[1];
sx q[1];
rz(-2.577563) q[1];
sx q[1];
rz(-0.86947039) q[1];
rz(-pi) q[2];
rz(-1.1060171) q[3];
sx q[3];
rz(-1.4651863) q[3];
sx q[3];
rz(-1.5231903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8308668) q[2];
sx q[2];
rz(-1.5006289) q[2];
sx q[2];
rz(-0.93079981) q[2];
rz(-2.9860949) q[3];
sx q[3];
rz(-1.5036539) q[3];
sx q[3];
rz(-0.29155198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.61313066) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(2.2303175) q[0];
rz(-2.7032734) q[1];
sx q[1];
rz(-1.8194018) q[1];
sx q[1];
rz(1.8211676) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5116918) q[0];
sx q[0];
rz(-1.1612079) q[0];
sx q[0];
rz(-3.1303309) q[0];
rz(2.690372) q[2];
sx q[2];
rz(-1.7570474) q[2];
sx q[2];
rz(0.95552432) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3418386) q[1];
sx q[1];
rz(-2.2463887) q[1];
sx q[1];
rz(-1.9104596) q[1];
x q[2];
rz(0.96418013) q[3];
sx q[3];
rz(-0.88084953) q[3];
sx q[3];
rz(-2.6026158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0358255) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(-2.7992115) q[2];
rz(-0.17677447) q[3];
sx q[3];
rz(-0.43313679) q[3];
sx q[3];
rz(-1.140973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3115561) q[0];
sx q[0];
rz(-0.72768584) q[0];
sx q[0];
rz(2.2763021) q[0];
rz(-1.226549) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(-1.3006166) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0375992) q[0];
sx q[0];
rz(-0.35498699) q[0];
sx q[0];
rz(0.07261891) q[0];
rz(2.956203) q[2];
sx q[2];
rz(-2.7692147) q[2];
sx q[2];
rz(-2.7951954) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.865766) q[1];
sx q[1];
rz(-1.8991125) q[1];
sx q[1];
rz(-2.5715716) q[1];
x q[2];
rz(-1.424765) q[3];
sx q[3];
rz(-1.8304123) q[3];
sx q[3];
rz(-2.8433593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1317923) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(0.47719964) q[2];
rz(0.19208433) q[3];
sx q[3];
rz(-1.6936857) q[3];
sx q[3];
rz(-0.93311667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.3451097) q[0];
sx q[0];
rz(-2.5273297) q[0];
sx q[0];
rz(3.1298424) q[0];
rz(-0.55039644) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(1.5884429) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9827305) q[0];
sx q[0];
rz(-1.9059062) q[0];
sx q[0];
rz(-1.1984675) q[0];
x q[1];
rz(-1.5137709) q[2];
sx q[2];
rz(-0.36643039) q[2];
sx q[2];
rz(2.1232405) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1482684) q[1];
sx q[1];
rz(-0.71422186) q[1];
sx q[1];
rz(-0.12970129) q[1];
rz(-0.074514975) q[3];
sx q[3];
rz(-0.86935589) q[3];
sx q[3];
rz(-0.45267347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6340296) q[2];
sx q[2];
rz(-0.66036779) q[2];
sx q[2];
rz(-1.2825512) q[2];
rz(1.3698618) q[3];
sx q[3];
rz(-1.7462574) q[3];
sx q[3];
rz(-1.1184568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58105528) q[0];
sx q[0];
rz(-2.9736309) q[0];
sx q[0];
rz(-2.4643331) q[0];
rz(0.15180763) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(-0.97704926) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1099694) q[0];
sx q[0];
rz(-0.97951802) q[0];
sx q[0];
rz(-1.4800319) q[0];
rz(-3.0095519) q[2];
sx q[2];
rz(-1.5706976) q[2];
sx q[2];
rz(-1.4591109) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3358826) q[1];
sx q[1];
rz(-1.4937595) q[1];
sx q[1];
rz(1.2309993) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3903923) q[3];
sx q[3];
rz(-2.6874472) q[3];
sx q[3];
rz(-0.62869149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3892422) q[2];
sx q[2];
rz(-0.82169473) q[2];
sx q[2];
rz(-1.0127257) q[2];
rz(-1.1879454) q[3];
sx q[3];
rz(-1.0725189) q[3];
sx q[3];
rz(0.48721203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.0174719) q[0];
sx q[0];
rz(-0.033360632) q[0];
sx q[0];
rz(-2.4429328) q[0];
rz(-2.0195122) q[1];
sx q[1];
rz(-0.84609234) q[1];
sx q[1];
rz(-1.8922071) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.218924) q[0];
sx q[0];
rz(-1.7814753) q[0];
sx q[0];
rz(-1.6181437) q[0];
rz(-pi) q[1];
rz(-2.6685733) q[2];
sx q[2];
rz(-2.0447391) q[2];
sx q[2];
rz(-0.75616403) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3932712) q[1];
sx q[1];
rz(-1.5652579) q[1];
sx q[1];
rz(-3.1194411) q[1];
rz(2.2035355) q[3];
sx q[3];
rz(-1.3360099) q[3];
sx q[3];
rz(-2.7126922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32968783) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(1.1784941) q[2];
rz(1.4568436) q[3];
sx q[3];
rz(-1.0624351) q[3];
sx q[3];
rz(-2.7594574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4998528) q[0];
sx q[0];
rz(-1.7620182) q[0];
sx q[0];
rz(1.8485803) q[0];
rz(-1.7199843) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(0.59757772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7062089) q[0];
sx q[0];
rz(-1.4357114) q[0];
sx q[0];
rz(3.1027017) q[0];
rz(-1.045268) q[2];
sx q[2];
rz(-2.0647486) q[2];
sx q[2];
rz(-1.2509105) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4328879) q[1];
sx q[1];
rz(-1.27379) q[1];
sx q[1];
rz(1.9087285) q[1];
rz(-pi) q[2];
rz(-0.58595539) q[3];
sx q[3];
rz(-1.6349941) q[3];
sx q[3];
rz(-2.5776598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.22275816) q[2];
sx q[2];
rz(-1.6901878) q[2];
sx q[2];
rz(1.9082327) q[2];
rz(2.2402066) q[3];
sx q[3];
rz(-3.021535) q[3];
sx q[3];
rz(1.4982769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0937061) q[0];
sx q[0];
rz(-0.77195764) q[0];
sx q[0];
rz(-3.1179324) q[0];
rz(-0.95611447) q[1];
sx q[1];
rz(-1.8319943) q[1];
sx q[1];
rz(-2.4694209) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2490847) q[0];
sx q[0];
rz(-1.5605643) q[0];
sx q[0];
rz(3.1363048) q[0];
rz(-0.51151885) q[2];
sx q[2];
rz(-1.5626972) q[2];
sx q[2];
rz(-0.17009232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.086239554) q[1];
sx q[1];
rz(-1.9085911) q[1];
sx q[1];
rz(2.4243381) q[1];
x q[2];
rz(-0.43283312) q[3];
sx q[3];
rz(-1.4308617) q[3];
sx q[3];
rz(0.39633745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51222926) q[2];
sx q[2];
rz(-1.9146634) q[2];
sx q[2];
rz(2.771634) q[2];
rz(1.6379179) q[3];
sx q[3];
rz(-0.88589293) q[3];
sx q[3];
rz(-1.9406208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5794012) q[0];
sx q[0];
rz(-0.36407064) q[0];
sx q[0];
rz(-1.9343485) q[0];
rz(-2.4178986) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(-1.9696708) q[2];
sx q[2];
rz(-1.4234067) q[2];
sx q[2];
rz(1.1248551) q[2];
rz(-1.1643812) q[3];
sx q[3];
rz(-1.8462528) q[3];
sx q[3];
rz(0.13312199) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
