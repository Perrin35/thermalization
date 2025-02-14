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
rz(0.87822479) q[0];
sx q[0];
rz(-2.4189334) q[0];
sx q[0];
rz(1.0455796) q[0];
rz(-3.1388406) q[1];
sx q[1];
rz(-1.6834919) q[1];
sx q[1];
rz(-2.4207065) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2963639) q[0];
sx q[0];
rz(-0.90334821) q[0];
sx q[0];
rz(3.0240701) q[0];
rz(-pi) q[1];
rz(3.0903417) q[2];
sx q[2];
rz(-1.379671) q[2];
sx q[2];
rz(-2.4911936) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.73511368) q[1];
sx q[1];
rz(-1.8476474) q[1];
sx q[1];
rz(0.39445725) q[1];
x q[2];
rz(-3.0009427) q[3];
sx q[3];
rz(-0.74797927) q[3];
sx q[3];
rz(1.5281531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8568933) q[2];
sx q[2];
rz(-1.928494) q[2];
sx q[2];
rz(-2.6846679) q[2];
rz(-2.4871155) q[3];
sx q[3];
rz(-2.579253) q[3];
sx q[3];
rz(-2.7175609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.65861312) q[0];
sx q[0];
rz(-2.787866) q[0];
sx q[0];
rz(-2.5335627) q[0];
rz(1.1990625) q[1];
sx q[1];
rz(-2.6741195) q[1];
sx q[1];
rz(-0.79493585) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011493135) q[0];
sx q[0];
rz(-1.8602295) q[0];
sx q[0];
rz(-0.35057803) q[0];
rz(-pi) q[1];
rz(-1.5420004) q[2];
sx q[2];
rz(-1.8378647) q[2];
sx q[2];
rz(1.1065337) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39917377) q[1];
sx q[1];
rz(-1.292242) q[1];
sx q[1];
rz(0.919085) q[1];
rz(-pi) q[2];
rz(-2.7057418) q[3];
sx q[3];
rz(-2.4727949) q[3];
sx q[3];
rz(-0.28025249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6752211) q[2];
sx q[2];
rz(-1.0543062) q[2];
sx q[2];
rz(1.6611453) q[2];
rz(0.36888567) q[3];
sx q[3];
rz(-1.8921013) q[3];
sx q[3];
rz(2.4108385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9166301) q[0];
sx q[0];
rz(-2.1376305) q[0];
sx q[0];
rz(0.68600255) q[0];
rz(1.8983768) q[1];
sx q[1];
rz(-2.9153283) q[1];
sx q[1];
rz(2.539198) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7334977) q[0];
sx q[0];
rz(-1.0801093) q[0];
sx q[0];
rz(-2.6564834) q[0];
rz(-pi) q[1];
rz(-3.1092293) q[2];
sx q[2];
rz(-1.6113069) q[2];
sx q[2];
rz(1.4537741) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6210052) q[1];
sx q[1];
rz(-0.57656724) q[1];
sx q[1];
rz(-1.3919049) q[1];
rz(-pi) q[2];
rz(2.8152594) q[3];
sx q[3];
rz(-2.0179393) q[3];
sx q[3];
rz(0.23995072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7556222) q[2];
sx q[2];
rz(-1.8381511) q[2];
sx q[2];
rz(-1.8899567) q[2];
rz(-2.5899467) q[3];
sx q[3];
rz(-2.9988204) q[3];
sx q[3];
rz(-2.2940206) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30152339) q[0];
sx q[0];
rz(-1.2395549) q[0];
sx q[0];
rz(0.10313343) q[0];
rz(0.81941191) q[1];
sx q[1];
rz(-2.2401919) q[1];
sx q[1];
rz(-2.588396) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8545733) q[0];
sx q[0];
rz(-0.4724882) q[0];
sx q[0];
rz(2.2437139) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0290506) q[2];
sx q[2];
rz(-2.5248387) q[2];
sx q[2];
rz(2.665131) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8242632) q[1];
sx q[1];
rz(-1.1019442) q[1];
sx q[1];
rz(0.15293185) q[1];
x q[2];
rz(1.3534773) q[3];
sx q[3];
rz(-0.44063755) q[3];
sx q[3];
rz(-2.8348451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72054046) q[2];
sx q[2];
rz(-2.4002176) q[2];
sx q[2];
rz(-0.53483024) q[2];
rz(-0.77830166) q[3];
sx q[3];
rz(-2.4141267) q[3];
sx q[3];
rz(-2.9261869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(2.4541723) q[0];
sx q[0];
rz(-1.9947616) q[0];
sx q[0];
rz(0.73690328) q[0];
rz(-0.058628254) q[1];
sx q[1];
rz(-0.77719378) q[1];
sx q[1];
rz(-0.58707213) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0995057) q[0];
sx q[0];
rz(-1.5387968) q[0];
sx q[0];
rz(-0.93670292) q[0];
rz(-pi) q[1];
rz(-0.063063085) q[2];
sx q[2];
rz(-1.3418256) q[2];
sx q[2];
rz(-1.582043) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78898326) q[1];
sx q[1];
rz(-2.8375344) q[1];
sx q[1];
rz(1.6795781) q[1];
rz(-pi) q[2];
rz(-2.8361017) q[3];
sx q[3];
rz(-1.2661714) q[3];
sx q[3];
rz(-0.50034278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29207841) q[2];
sx q[2];
rz(-0.92486113) q[2];
sx q[2];
rz(-2.4936131) q[2];
rz(-0.61602151) q[3];
sx q[3];
rz(-1.5024622) q[3];
sx q[3];
rz(-0.047671635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7306526) q[0];
sx q[0];
rz(-0.76496449) q[0];
sx q[0];
rz(2.1867645) q[0];
rz(-3.094589) q[1];
sx q[1];
rz(-0.67263293) q[1];
sx q[1];
rz(1.0161374) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2178041) q[0];
sx q[0];
rz(-2.1013484) q[0];
sx q[0];
rz(2.1235882) q[0];
rz(0.59784933) q[2];
sx q[2];
rz(-2.3232004) q[2];
sx q[2];
rz(-2.1482836) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.816427) q[1];
sx q[1];
rz(-1.0652958) q[1];
sx q[1];
rz(-1.5542404) q[1];
rz(-0.079214752) q[3];
sx q[3];
rz(-2.5159249) q[3];
sx q[3];
rz(0.46712671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.59812349) q[2];
sx q[2];
rz(-1.8757966) q[2];
sx q[2];
rz(-0.050696105) q[2];
rz(-3.0905881) q[3];
sx q[3];
rz(-1.8404605) q[3];
sx q[3];
rz(0.72371975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26492587) q[0];
sx q[0];
rz(-2.6267138) q[0];
sx q[0];
rz(1.4798973) q[0];
rz(-1.121608) q[1];
sx q[1];
rz(-1.2622967) q[1];
sx q[1];
rz(-0.46953896) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4001842) q[0];
sx q[0];
rz(-0.11708507) q[0];
sx q[0];
rz(1.1598806) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15596822) q[2];
sx q[2];
rz(-1.9591589) q[2];
sx q[2];
rz(-2.591382) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.6958298) q[1];
sx q[1];
rz(-1.1053021) q[1];
sx q[1];
rz(0.19449921) q[1];
rz(-pi) q[2];
rz(-2.2121053) q[3];
sx q[3];
rz(-0.90228122) q[3];
sx q[3];
rz(0.53169441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7349714) q[2];
sx q[2];
rz(-1.5512286) q[2];
sx q[2];
rz(2.0242019) q[2];
rz(-2.0495074) q[3];
sx q[3];
rz(-1.4039682) q[3];
sx q[3];
rz(-2.7058069) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0206873) q[0];
sx q[0];
rz(-1.1130604) q[0];
sx q[0];
rz(0.06509617) q[0];
rz(2.8096325) q[1];
sx q[1];
rz(-1.8354225) q[1];
sx q[1];
rz(1.8665705) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1972361) q[0];
sx q[0];
rz(-1.580339) q[0];
sx q[0];
rz(-1.8002323) q[0];
rz(-3.0049999) q[2];
sx q[2];
rz(-2.0099735) q[2];
sx q[2];
rz(-2.3802649) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.79679497) q[1];
sx q[1];
rz(-1.5196504) q[1];
sx q[1];
rz(2.4586492) q[1];
rz(-pi) q[2];
rz(-0.80294369) q[3];
sx q[3];
rz(-1.6979473) q[3];
sx q[3];
rz(-1.2252327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8518105) q[2];
sx q[2];
rz(-1.6921356) q[2];
sx q[2];
rz(-0.97274485) q[2];
rz(0.70074493) q[3];
sx q[3];
rz(-2.2795129) q[3];
sx q[3];
rz(-0.8249445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.7026611) q[0];
sx q[0];
rz(-0.40533608) q[0];
sx q[0];
rz(2.7324556) q[0];
rz(1.9955955) q[1];
sx q[1];
rz(-1.0640249) q[1];
sx q[1];
rz(1.8786059) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0008903) q[0];
sx q[0];
rz(-1.5626199) q[0];
sx q[0];
rz(-3.1414638) q[0];
x q[1];
rz(1.9648329) q[2];
sx q[2];
rz(-1.3140643) q[2];
sx q[2];
rz(0.22025157) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4597171) q[1];
sx q[1];
rz(-0.72725216) q[1];
sx q[1];
rz(1.9161112) q[1];
x q[2];
rz(1.658012) q[3];
sx q[3];
rz(-2.1391807) q[3];
sx q[3];
rz(-1.5414814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9001749) q[2];
sx q[2];
rz(-1.2215542) q[2];
sx q[2];
rz(0.22132604) q[2];
rz(-2.6545702) q[3];
sx q[3];
rz(-2.2723787) q[3];
sx q[3];
rz(-1.2041436) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0814334) q[0];
sx q[0];
rz(-2.8964323) q[0];
sx q[0];
rz(-1.2231476) q[0];
rz(-3.0606015) q[1];
sx q[1];
rz(-1.6209737) q[1];
sx q[1];
rz(-1.5222668) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.318383) q[0];
sx q[0];
rz(-1.2391187) q[0];
sx q[0];
rz(1.3964425) q[0];
x q[1];
rz(2.1606101) q[2];
sx q[2];
rz(-2.207617) q[2];
sx q[2];
rz(2.3825519) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.63357089) q[1];
sx q[1];
rz(-2.7765818) q[1];
sx q[1];
rz(0.99442039) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6468148) q[3];
sx q[3];
rz(-1.4162546) q[3];
sx q[3];
rz(3.0734594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12485518) q[2];
sx q[2];
rz(-1.3913466) q[2];
sx q[2];
rz(-0.69046956) q[2];
rz(-0.27103439) q[3];
sx q[3];
rz(-2.6764968) q[3];
sx q[3];
rz(-1.5439699) q[3];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4616213) q[0];
sx q[0];
rz(-0.93628708) q[0];
sx q[0];
rz(-0.17967889) q[0];
rz(0.24196729) q[1];
sx q[1];
rz(-2.2503743) q[1];
sx q[1];
rz(1.888884) q[1];
rz(-2.7015637) q[2];
sx q[2];
rz(-2.6806705) q[2];
sx q[2];
rz(-0.086967998) q[2];
rz(-0.88884066) q[3];
sx q[3];
rz(-1.325229) q[3];
sx q[3];
rz(-1.6390507) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
