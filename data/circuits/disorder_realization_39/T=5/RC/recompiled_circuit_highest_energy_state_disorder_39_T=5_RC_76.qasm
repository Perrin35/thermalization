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
rz(0.14035913) q[0];
sx q[0];
rz(4.601534) q[0];
sx q[0];
rz(10.277716) q[0];
rz(1.9536904) q[1];
sx q[1];
rz(3.3878769) q[1];
sx q[1];
rz(9.090957) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.649619) q[0];
sx q[0];
rz(-0.57081693) q[0];
sx q[0];
rz(1.094857) q[0];
x q[1];
rz(-0.63392459) q[2];
sx q[2];
rz(-1.156154) q[2];
sx q[2];
rz(-2.7716605) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.546963) q[1];
sx q[1];
rz(-1.3876042) q[1];
sx q[1];
rz(-1.4180844) q[1];
rz(0.6972972) q[3];
sx q[3];
rz(-3.0099031) q[3];
sx q[3];
rz(-3.0856709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9217801) q[2];
sx q[2];
rz(-0.87679902) q[2];
sx q[2];
rz(-0.63966695) q[2];
rz(-2.2031247) q[3];
sx q[3];
rz(-2.7191021) q[3];
sx q[3];
rz(-1.870702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4939782) q[0];
sx q[0];
rz(-3.0632126) q[0];
sx q[0];
rz(-1.6139503) q[0];
rz(2.899462) q[1];
sx q[1];
rz(-2.1277728) q[1];
sx q[1];
rz(0.72227824) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5764389) q[0];
sx q[0];
rz(-1.1761888) q[0];
sx q[0];
rz(0.37532401) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3039178) q[2];
sx q[2];
rz(-0.79395771) q[2];
sx q[2];
rz(-2.4030952) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7178044) q[1];
sx q[1];
rz(-0.77149888) q[1];
sx q[1];
rz(-2.3256734) q[1];
rz(-1.8037968) q[3];
sx q[3];
rz(-0.63574857) q[3];
sx q[3];
rz(-2.014117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0904514) q[2];
sx q[2];
rz(-0.34153667) q[2];
sx q[2];
rz(-3.0549468) q[2];
rz(2.5610949) q[3];
sx q[3];
rz(-1.2248421) q[3];
sx q[3];
rz(2.1897924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3314303) q[0];
sx q[0];
rz(-2.0482752) q[0];
sx q[0];
rz(1.6828368) q[0];
rz(-3.0259865) q[1];
sx q[1];
rz(-1.1980201) q[1];
sx q[1];
rz(0.16673949) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0519823) q[0];
sx q[0];
rz(-1.3638743) q[0];
sx q[0];
rz(0.99159411) q[0];
rz(-pi) q[1];
rz(0.71311976) q[2];
sx q[2];
rz(-2.6067039) q[2];
sx q[2];
rz(-2.5120017) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.48651153) q[1];
sx q[1];
rz(-1.1919199) q[1];
sx q[1];
rz(-0.98030555) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9882713) q[3];
sx q[3];
rz(-0.75107952) q[3];
sx q[3];
rz(0.73634232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.705767) q[2];
sx q[2];
rz(-2.9220118) q[2];
sx q[2];
rz(1.1337918) q[2];
rz(0.052915834) q[3];
sx q[3];
rz(-1.8688801) q[3];
sx q[3];
rz(2.4165966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.13919203) q[0];
sx q[0];
rz(-2.5097558) q[0];
sx q[0];
rz(-2.8644526) q[0];
rz(-0.78284043) q[1];
sx q[1];
rz(-1.6354086) q[1];
sx q[1];
rz(-0.35071075) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6445054) q[0];
sx q[0];
rz(-0.86928029) q[0];
sx q[0];
rz(0.10414609) q[0];
rz(-3.0530351) q[2];
sx q[2];
rz(-0.32677256) q[2];
sx q[2];
rz(-1.8984924) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.86736401) q[1];
sx q[1];
rz(-1.1038934) q[1];
sx q[1];
rz(2.2850249) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1923741) q[3];
sx q[3];
rz(-0.11648341) q[3];
sx q[3];
rz(1.0939712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0413401) q[2];
sx q[2];
rz(-1.2790054) q[2];
sx q[2];
rz(-1.9913541) q[2];
rz(-0.1746812) q[3];
sx q[3];
rz(-0.29398578) q[3];
sx q[3];
rz(-0.090350769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1133076) q[0];
sx q[0];
rz(-2.568013) q[0];
sx q[0];
rz(-1.9453402) q[0];
rz(-0.47736564) q[1];
sx q[1];
rz(-0.82672516) q[1];
sx q[1];
rz(0.54944077) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7350205) q[0];
sx q[0];
rz(-1.0112678) q[0];
sx q[0];
rz(2.4491549) q[0];
rz(-pi) q[1];
rz(-0.35713335) q[2];
sx q[2];
rz(-2.5741842) q[2];
sx q[2];
rz(-1.0536989) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0116033) q[1];
sx q[1];
rz(-1.450256) q[1];
sx q[1];
rz(-1.1611564) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0277689) q[3];
sx q[3];
rz(-1.1871871) q[3];
sx q[3];
rz(-2.3758604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4598733) q[2];
sx q[2];
rz(-2.7698066) q[2];
sx q[2];
rz(2.8968774) q[2];
rz(-1.5185482) q[3];
sx q[3];
rz(-2.0270429) q[3];
sx q[3];
rz(-3.0486619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2731584) q[0];
sx q[0];
rz(-2.1491829) q[0];
sx q[0];
rz(0.42700818) q[0];
rz(1.2097516) q[1];
sx q[1];
rz(-1.6170343) q[1];
sx q[1];
rz(-0.0078113656) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2526379) q[0];
sx q[0];
rz(-1.3547621) q[0];
sx q[0];
rz(0.81058575) q[0];
rz(-pi) q[1];
rz(1.5862238) q[2];
sx q[2];
rz(-0.54718218) q[2];
sx q[2];
rz(1.5146966) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5558034) q[1];
sx q[1];
rz(-0.92098707) q[1];
sx q[1];
rz(1.5073677) q[1];
rz(-2.2676007) q[3];
sx q[3];
rz(-1.7808508) q[3];
sx q[3];
rz(1.7396613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4012287) q[2];
sx q[2];
rz(-0.87330356) q[2];
sx q[2];
rz(1.1391696) q[2];
rz(0.27979699) q[3];
sx q[3];
rz(-1.8179025) q[3];
sx q[3];
rz(-1.4626224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6497659) q[0];
sx q[0];
rz(-0.38182807) q[0];
sx q[0];
rz(-2.5873798) q[0];
rz(0.062049374) q[1];
sx q[1];
rz(-0.65826145) q[1];
sx q[1];
rz(2.2668692) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5420008) q[0];
sx q[0];
rz(-2.3440147) q[0];
sx q[0];
rz(2.9717658) q[0];
rz(-pi) q[1];
rz(3.1360097) q[2];
sx q[2];
rz(-1.5783596) q[2];
sx q[2];
rz(-0.39952229) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8043912) q[1];
sx q[1];
rz(-0.60742765) q[1];
sx q[1];
rz(-2.0359705) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0482385) q[3];
sx q[3];
rz(-2.4191609) q[3];
sx q[3];
rz(2.5337766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41500652) q[2];
sx q[2];
rz(-1.0174624) q[2];
sx q[2];
rz(-2.2028108) q[2];
rz(2.4472661) q[3];
sx q[3];
rz(-2.4735579) q[3];
sx q[3];
rz(0.15650775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
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
rz(2.0701377) q[0];
sx q[0];
rz(-2.5818765) q[0];
sx q[0];
rz(0.30858421) q[0];
rz(-1.4831108) q[1];
sx q[1];
rz(-1.7959271) q[1];
sx q[1];
rz(-2.9187091) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73327979) q[0];
sx q[0];
rz(-2.2594497) q[0];
sx q[0];
rz(-2.0613502) q[0];
rz(-pi) q[1];
rz(0.22391386) q[2];
sx q[2];
rz(-1.4707047) q[2];
sx q[2];
rz(1.8318286) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.4207238) q[1];
sx q[1];
rz(-1.6096186) q[1];
sx q[1];
rz(1.3255408) q[1];
rz(-2.127366) q[3];
sx q[3];
rz(-2.1750919) q[3];
sx q[3];
rz(2.6984071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8181204) q[2];
sx q[2];
rz(-2.1151395) q[2];
sx q[2];
rz(0.66413122) q[2];
rz(-1.1431665) q[3];
sx q[3];
rz(-1.4799456) q[3];
sx q[3];
rz(1.0955048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7185862) q[0];
sx q[0];
rz(-1.9976595) q[0];
sx q[0];
rz(-3.123172) q[0];
rz(-2.9711235) q[1];
sx q[1];
rz(-1.2455995) q[1];
sx q[1];
rz(2.9439994) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7391846) q[0];
sx q[0];
rz(-1.5614334) q[0];
sx q[0];
rz(-1.739517) q[0];
rz(-pi) q[1];
rz(1.9371447) q[2];
sx q[2];
rz(-1.6197259) q[2];
sx q[2];
rz(-1.4885224) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8325628) q[1];
sx q[1];
rz(-1.4352928) q[1];
sx q[1];
rz(1.871528) q[1];
rz(2.7457947) q[3];
sx q[3];
rz(-1.5185818) q[3];
sx q[3];
rz(-2.8767916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1175055) q[2];
sx q[2];
rz(-1.0988289) q[2];
sx q[2];
rz(2.7433024) q[2];
rz(-2.6309218) q[3];
sx q[3];
rz(-1.6528249) q[3];
sx q[3];
rz(0.85810703) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2354105) q[0];
sx q[0];
rz(-2.372083) q[0];
sx q[0];
rz(-0.19943516) q[0];
rz(2.8681352) q[1];
sx q[1];
rz(-2.7077935) q[1];
sx q[1];
rz(2.1308897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0293297) q[0];
sx q[0];
rz(-0.71319095) q[0];
sx q[0];
rz(2.7351456) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23538537) q[2];
sx q[2];
rz(-1.590325) q[2];
sx q[2];
rz(0.69695401) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.005474394) q[1];
sx q[1];
rz(-1.6228383) q[1];
sx q[1];
rz(-2.4050557) q[1];
rz(-pi) q[2];
rz(-1.2730678) q[3];
sx q[3];
rz(-1.862251) q[3];
sx q[3];
rz(3.1133661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.79591862) q[2];
sx q[2];
rz(-2.5359539) q[2];
sx q[2];
rz(-2.3460713) q[2];
rz(-2.7703721) q[3];
sx q[3];
rz(-1.755654) q[3];
sx q[3];
rz(-0.25446874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026074792) q[0];
sx q[0];
rz(-1.4160897) q[0];
sx q[0];
rz(-1.8427451) q[0];
rz(-0.55116354) q[1];
sx q[1];
rz(-1.9438585) q[1];
sx q[1];
rz(1.9895947) q[1];
rz(1.3462832) q[2];
sx q[2];
rz(-0.85805744) q[2];
sx q[2];
rz(-0.35408264) q[2];
rz(2.1195246) q[3];
sx q[3];
rz(-2.7014501) q[3];
sx q[3];
rz(-0.81934495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
