OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.70513201) q[0];
sx q[0];
rz(-2.5897265) q[0];
sx q[0];
rz(-0.021835672) q[0];
rz(-0.39437374) q[1];
sx q[1];
rz(-1.6819277) q[1];
sx q[1];
rz(0.2149166) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.874274) q[0];
sx q[0];
rz(-1.8288757) q[0];
sx q[0];
rz(1.2866856) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8432437) q[2];
sx q[2];
rz(-0.48765182) q[2];
sx q[2];
rz(1.2717441) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.51337459) q[1];
sx q[1];
rz(-0.77925032) q[1];
sx q[1];
rz(0.90374225) q[1];
x q[2];
rz(-2.05902) q[3];
sx q[3];
rz(-2.2256652) q[3];
sx q[3];
rz(1.0893351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(-2.5773876) q[2];
rz(1.7764067) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(-1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0974225) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(-0.92798293) q[0];
rz(1.9762951) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(-0.89675084) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1363163) q[0];
sx q[0];
rz(-1.3747842) q[0];
sx q[0];
rz(-0.85640237) q[0];
rz(2.9727544) q[2];
sx q[2];
rz(-1.642792) q[2];
sx q[2];
rz(-2.8746586) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7795418) q[1];
sx q[1];
rz(-0.78501399) q[1];
sx q[1];
rz(1.7464459) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41665839) q[3];
sx q[3];
rz(-2.2778802) q[3];
sx q[3];
rz(-1.7435031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.26560489) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(-2.1014452) q[2];
rz(-1.4552207) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74137694) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(1.0282015) q[0];
rz(-1.0785412) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(0.43513402) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2745278) q[0];
sx q[0];
rz(-1.2194249) q[0];
sx q[0];
rz(1.4562796) q[0];
rz(-pi) q[1];
rz(-2.8912796) q[2];
sx q[2];
rz(-1.4280983) q[2];
sx q[2];
rz(2.755969) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.16776925) q[1];
sx q[1];
rz(-2.6288189) q[1];
sx q[1];
rz(1.2330526) q[1];
rz(-2.7828091) q[3];
sx q[3];
rz(-2.3453418) q[3];
sx q[3];
rz(-2.1735454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(2.8386774) q[2];
rz(1.3251925) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19642297) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(0.91745013) q[0];
rz(-0.67287412) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(0.26487574) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6661975) q[0];
sx q[0];
rz(-0.65984939) q[0];
sx q[0];
rz(-3.092479) q[0];
rz(-pi) q[1];
rz(2.3189544) q[2];
sx q[2];
rz(-2.7441141) q[2];
sx q[2];
rz(0.17695225) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3600196) q[1];
sx q[1];
rz(-2.2847166) q[1];
sx q[1];
rz(-0.54846958) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75986741) q[3];
sx q[3];
rz(-2.2025975) q[3];
sx q[3];
rz(2.2894273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.36310568) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(1.5650361) q[2];
rz(1.0270843) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(-2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.251579) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(-0.47873163) q[0];
rz(2.1084673) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(-0.95265257) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.032741) q[0];
sx q[0];
rz(-1.6881588) q[0];
sx q[0];
rz(-0.54876859) q[0];
rz(-pi) q[1];
rz(-0.41575899) q[2];
sx q[2];
rz(-1.8766878) q[2];
sx q[2];
rz(-1.8005467) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.488727) q[1];
sx q[1];
rz(-1.9619202) q[1];
sx q[1];
rz(2.8847242) q[1];
x q[2];
rz(-1.5289375) q[3];
sx q[3];
rz(-2.0484945) q[3];
sx q[3];
rz(2.4822513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4218563) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(2.6110113) q[2];
rz(-1.4060219) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(0.83166844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56753165) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(-1.5166327) q[0];
rz(1.3051055) q[1];
sx q[1];
rz(-1.3508947) q[1];
sx q[1];
rz(-0.17257246) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1369551) q[0];
sx q[0];
rz(-1.9872268) q[0];
sx q[0];
rz(3.1266771) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2595915) q[2];
sx q[2];
rz(-1.1278369) q[2];
sx q[2];
rz(-1.7829347) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.29218201) q[1];
sx q[1];
rz(-1.9841571) q[1];
sx q[1];
rz(-0.11489111) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63038007) q[3];
sx q[3];
rz(-1.3833589) q[3];
sx q[3];
rz(2.558625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3036348) q[2];
sx q[2];
rz(-1.6875608) q[2];
sx q[2];
rz(-1.1266358) q[2];
rz(-0.78222328) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(-1.8036028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28850266) q[0];
sx q[0];
rz(-0.30650109) q[0];
sx q[0];
rz(0.66147584) q[0];
rz(0.96039564) q[1];
sx q[1];
rz(-1.7405225) q[1];
sx q[1];
rz(0.75659928) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7717168) q[0];
sx q[0];
rz(-1.6635832) q[0];
sx q[0];
rz(1.7232399) q[0];
rz(-pi) q[1];
rz(2.7301894) q[2];
sx q[2];
rz(-1.7559933) q[2];
sx q[2];
rz(-1.7298557) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7927671) q[1];
sx q[1];
rz(-1.7276689) q[1];
sx q[1];
rz(-1.8030333) q[1];
x q[2];
rz(-1.6301304) q[3];
sx q[3];
rz(-1.3795128) q[3];
sx q[3];
rz(0.64185601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1371655) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(0.19443092) q[2];
rz(-0.91313177) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(-3.074926) q[0];
rz(0.32456675) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(-2.1527122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38935223) q[0];
sx q[0];
rz(-1.9703431) q[0];
sx q[0];
rz(0.35895343) q[0];
rz(-pi) q[1];
rz(-1.5867932) q[2];
sx q[2];
rz(-1.4774067) q[2];
sx q[2];
rz(-2.8430251) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5576396) q[1];
sx q[1];
rz(-1.156731) q[1];
sx q[1];
rz(-0.072903452) q[1];
rz(-pi) q[2];
rz(-2.3905121) q[3];
sx q[3];
rz(-2.4820648) q[3];
sx q[3];
rz(0.86576033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4618335) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(-0.40763339) q[2];
rz(-0.76861012) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(-0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75893629) q[0];
sx q[0];
rz(-1.8396682) q[0];
sx q[0];
rz(2.5323903) q[0];
rz(-0.095104782) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(-0.87337714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0926889) q[0];
sx q[0];
rz(-0.19177076) q[0];
sx q[0];
rz(1.2094686) q[0];
x q[1];
rz(-1.4344425) q[2];
sx q[2];
rz(-1.7169723) q[2];
sx q[2];
rz(1.2863976) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1904618) q[1];
sx q[1];
rz(-1.8776263) q[1];
sx q[1];
rz(2.5224586) q[1];
rz(0.43975131) q[3];
sx q[3];
rz(-1.0381178) q[3];
sx q[3];
rz(0.38910481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1197027) q[2];
sx q[2];
rz(-0.78803524) q[2];
sx q[2];
rz(2.4592887) q[2];
rz(-0.37426379) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69797126) q[0];
sx q[0];
rz(-1.3598096) q[0];
sx q[0];
rz(2.1886254) q[0];
rz(0.8264181) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(-1.7451161) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6122702) q[0];
sx q[0];
rz(-1.4850052) q[0];
sx q[0];
rz(2.8180608) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7637485) q[2];
sx q[2];
rz(-2.4477738) q[2];
sx q[2];
rz(-0.56570429) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5212621) q[1];
sx q[1];
rz(-1.2863103) q[1];
sx q[1];
rz(1.8950589) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8793328) q[3];
sx q[3];
rz(-1.6256623) q[3];
sx q[3];
rz(-0.39615397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46618) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(-0.15979016) q[2];
rz(-2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(-0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0970584) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(-0.13327577) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(-1.3143905) q[2];
sx q[2];
rz(-2.495043) q[2];
sx q[2];
rz(-1.3174353) q[2];
rz(-1.521048) q[3];
sx q[3];
rz(-1.0377025) q[3];
sx q[3];
rz(2.514537) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];