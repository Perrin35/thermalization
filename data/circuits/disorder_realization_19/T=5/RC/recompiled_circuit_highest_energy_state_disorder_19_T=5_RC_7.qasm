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
rz(0.49993604) q[0];
sx q[0];
rz(1.9498107) q[0];
sx q[0];
rz(10.796588) q[0];
rz(1.6504047) q[1];
sx q[1];
rz(-0.86725441) q[1];
sx q[1];
rz(0.39630085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21090487) q[0];
sx q[0];
rz(-2.9172724) q[0];
sx q[0];
rz(2.7773662) q[0];
x q[1];
rz(-1.195186) q[2];
sx q[2];
rz(-1.7306788) q[2];
sx q[2];
rz(-2.218046) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5057136) q[1];
sx q[1];
rz(-1.6649655) q[1];
sx q[1];
rz(-1.0416743) q[1];
rz(2.7462148) q[3];
sx q[3];
rz(-1.3002743) q[3];
sx q[3];
rz(0.45365712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1402011) q[2];
sx q[2];
rz(-2.3346257) q[2];
sx q[2];
rz(1.3414471) q[2];
rz(-2.0724824) q[3];
sx q[3];
rz(-0.1461229) q[3];
sx q[3];
rz(1.3816381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5694201) q[0];
sx q[0];
rz(-0.91745806) q[0];
sx q[0];
rz(-0.65895748) q[0];
rz(-3.068889) q[1];
sx q[1];
rz(-2.0854988) q[1];
sx q[1];
rz(1.2603849) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1457659) q[0];
sx q[0];
rz(-0.46016177) q[0];
sx q[0];
rz(1.3868679) q[0];
rz(-pi) q[1];
rz(2.6953893) q[2];
sx q[2];
rz(-1.1006736) q[2];
sx q[2];
rz(-2.0044495) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.61140984) q[1];
sx q[1];
rz(-1.5059222) q[1];
sx q[1];
rz(2.8758509) q[1];
x q[2];
rz(0.96003344) q[3];
sx q[3];
rz(-1.472578) q[3];
sx q[3];
rz(-1.4032422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8947123) q[2];
sx q[2];
rz(-0.74328819) q[2];
sx q[2];
rz(1.7903719) q[2];
rz(2.5249935) q[3];
sx q[3];
rz(-1.0298046) q[3];
sx q[3];
rz(-0.29115796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5093444) q[0];
sx q[0];
rz(-2.6743439) q[0];
sx q[0];
rz(-0.21009357) q[0];
rz(0.083077438) q[1];
sx q[1];
rz(-1.9030842) q[1];
sx q[1];
rz(-3.074379) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34850304) q[0];
sx q[0];
rz(-1.5572892) q[0];
sx q[0];
rz(1.5587139) q[0];
x q[1];
rz(-2.1234197) q[2];
sx q[2];
rz(-2.6913096) q[2];
sx q[2];
rz(-2.1076815) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.018118) q[1];
sx q[1];
rz(-0.90318455) q[1];
sx q[1];
rz(-0.43161824) q[1];
rz(0.70560734) q[3];
sx q[3];
rz(-1.5021245) q[3];
sx q[3];
rz(-2.9182383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3136966) q[2];
sx q[2];
rz(-1.6500429) q[2];
sx q[2];
rz(-1.776604) q[2];
rz(2.7374173) q[3];
sx q[3];
rz(-1.8242691) q[3];
sx q[3];
rz(-1.9425758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2760524) q[0];
sx q[0];
rz(-1.7409538) q[0];
sx q[0];
rz(2.6633967) q[0];
rz(0.95282355) q[1];
sx q[1];
rz(-0.15004221) q[1];
sx q[1];
rz(2.731954) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3296369) q[0];
sx q[0];
rz(-2.4853737) q[0];
sx q[0];
rz(0.043270525) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6686001) q[2];
sx q[2];
rz(-1.1162283) q[2];
sx q[2];
rz(1.821988) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9497021) q[1];
sx q[1];
rz(-0.97964761) q[1];
sx q[1];
rz(-0.77456919) q[1];
rz(-pi) q[2];
rz(1.8919935) q[3];
sx q[3];
rz(-1.7068591) q[3];
sx q[3];
rz(2.4709159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5459368) q[2];
sx q[2];
rz(-1.7031534) q[2];
sx q[2];
rz(1.1379918) q[2];
rz(-0.99500895) q[3];
sx q[3];
rz(-2.5840839) q[3];
sx q[3];
rz(3.0549684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74395162) q[0];
sx q[0];
rz(-1.8331563) q[0];
sx q[0];
rz(-2.8635136) q[0];
rz(-1.8771578) q[1];
sx q[1];
rz(-1.8587298) q[1];
sx q[1];
rz(1.5637195) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20533097) q[0];
sx q[0];
rz(-2.7561128) q[0];
sx q[0];
rz(1.9140713) q[0];
rz(-pi) q[1];
rz(0.78863849) q[2];
sx q[2];
rz(-1.83162) q[2];
sx q[2];
rz(2.3344085) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.75018308) q[1];
sx q[1];
rz(-0.58779189) q[1];
sx q[1];
rz(2.9255735) q[1];
x q[2];
rz(2.0966846) q[3];
sx q[3];
rz(-2.9581986) q[3];
sx q[3];
rz(2.6332651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9331253) q[2];
sx q[2];
rz(-1.6105904) q[2];
sx q[2];
rz(2.7739286) q[2];
rz(2.0465046) q[3];
sx q[3];
rz(-1.0235267) q[3];
sx q[3];
rz(0.17525214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4472189) q[0];
sx q[0];
rz(-2.8419438) q[0];
sx q[0];
rz(-1.8023941) q[0];
rz(-0.12204349) q[1];
sx q[1];
rz(-1.3587147) q[1];
sx q[1];
rz(0.36062127) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1674408) q[0];
sx q[0];
rz(-1.9527148) q[0];
sx q[0];
rz(-2.4171305) q[0];
rz(-pi) q[1];
rz(-1.1507658) q[2];
sx q[2];
rz(-1.5345062) q[2];
sx q[2];
rz(2.0535713) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5721467) q[1];
sx q[1];
rz(-1.1317557) q[1];
sx q[1];
rz(-2.8296058) q[1];
rz(3.120947) q[3];
sx q[3];
rz(-0.25293487) q[3];
sx q[3];
rz(-1.3731352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3660761) q[2];
sx q[2];
rz(-1.6785494) q[2];
sx q[2];
rz(-2.0029081) q[2];
rz(2.8876997) q[3];
sx q[3];
rz(-0.64520276) q[3];
sx q[3];
rz(-2.3937288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426386) q[0];
sx q[0];
rz(-0.88976088) q[0];
sx q[0];
rz(0.97578543) q[0];
rz(1.1294533) q[1];
sx q[1];
rz(-0.44810805) q[1];
sx q[1];
rz(1.7879558) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8944446) q[0];
sx q[0];
rz(-1.8631051) q[0];
sx q[0];
rz(2.4185989) q[0];
x q[1];
rz(-0.079965429) q[2];
sx q[2];
rz(-1.1478394) q[2];
sx q[2];
rz(-2.1048196) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6401939) q[1];
sx q[1];
rz(-2.7796449) q[1];
sx q[1];
rz(0.37935235) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72809345) q[3];
sx q[3];
rz(-0.59641664) q[3];
sx q[3];
rz(2.8567258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5740616) q[2];
sx q[2];
rz(-2.3109544) q[2];
sx q[2];
rz(-2.2243824) q[2];
rz(1.5926825) q[3];
sx q[3];
rz(-2.8062688) q[3];
sx q[3];
rz(-0.77939051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2061283) q[0];
sx q[0];
rz(-1.6213106) q[0];
sx q[0];
rz(3.1167378) q[0];
rz(-1.8553597) q[1];
sx q[1];
rz(-1.7846466) q[1];
sx q[1];
rz(1.6110427) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8899621) q[0];
sx q[0];
rz(-2.0109091) q[0];
sx q[0];
rz(1.1391231) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3127682) q[2];
sx q[2];
rz(-2.240643) q[2];
sx q[2];
rz(-0.19425288) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3973786) q[1];
sx q[1];
rz(-1.2894003) q[1];
sx q[1];
rz(-1.534627) q[1];
rz(-pi) q[2];
rz(-3.0026818) q[3];
sx q[3];
rz(-1.362097) q[3];
sx q[3];
rz(3.1155966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1207235) q[2];
sx q[2];
rz(-2.1376231) q[2];
sx q[2];
rz(1.8570159) q[2];
rz(-2.8203216) q[3];
sx q[3];
rz(-0.64906859) q[3];
sx q[3];
rz(-2.9200714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016841737) q[0];
sx q[0];
rz(-1.9945972) q[0];
sx q[0];
rz(-2.2220213) q[0];
rz(2.549767) q[1];
sx q[1];
rz(-1.0706736) q[1];
sx q[1];
rz(-2.8048973) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29558173) q[0];
sx q[0];
rz(-2.8407556) q[0];
sx q[0];
rz(1.1016125) q[0];
rz(2.8317802) q[2];
sx q[2];
rz(-2.7618119) q[2];
sx q[2];
rz(-1.7526232) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.10178653) q[1];
sx q[1];
rz(-2.617604) q[1];
sx q[1];
rz(0.26079674) q[1];
x q[2];
rz(-1.0180264) q[3];
sx q[3];
rz(-2.8447731) q[3];
sx q[3];
rz(-0.33949131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0261953) q[2];
sx q[2];
rz(-0.43033174) q[2];
sx q[2];
rz(-0.78286147) q[2];
rz(3.03249) q[3];
sx q[3];
rz(-2.1249873) q[3];
sx q[3];
rz(-1.2469863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03610177) q[0];
sx q[0];
rz(-1.4694659) q[0];
sx q[0];
rz(-0.92392695) q[0];
rz(2.6149514) q[1];
sx q[1];
rz(-1.8662607) q[1];
sx q[1];
rz(2.142876) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0786506) q[0];
sx q[0];
rz(-1.6237769) q[0];
sx q[0];
rz(2.9880301) q[0];
rz(1.8417712) q[2];
sx q[2];
rz(-0.30108115) q[2];
sx q[2];
rz(-1.4512514) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0793887) q[1];
sx q[1];
rz(-0.98323373) q[1];
sx q[1];
rz(-0.04157898) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4630354) q[3];
sx q[3];
rz(-1.863409) q[3];
sx q[3];
rz(1.7998526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1363137) q[2];
sx q[2];
rz(-0.88999358) q[2];
sx q[2];
rz(-0.33373731) q[2];
rz(-2.2809095) q[3];
sx q[3];
rz(-1.656683) q[3];
sx q[3];
rz(3.1349643) q[3];
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
rz(0.46398791) q[0];
sx q[0];
rz(-1.9515568) q[0];
sx q[0];
rz(-2.4757181) q[0];
rz(-2.5869276) q[1];
sx q[1];
rz(-1.8634836) q[1];
sx q[1];
rz(-2.9992933) q[1];
rz(2.7016976) q[2];
sx q[2];
rz(-1.3439726) q[2];
sx q[2];
rz(-1.7769565) q[2];
rz(0.71698112) q[3];
sx q[3];
rz(-2.3752799) q[3];
sx q[3];
rz(0.3175288) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
