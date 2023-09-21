OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.93958062) q[0];
sx q[0];
rz(-0.35020819) q[0];
sx q[0];
rz(2.7749618) q[0];
rz(-2.2740013) q[1];
sx q[1];
rz(-0.64414135) q[1];
sx q[1];
rz(1.6860513) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35858425) q[0];
sx q[0];
rz(-1.9134247) q[0];
sx q[0];
rz(-1.1230099) q[0];
x q[1];
rz(-0.26379649) q[2];
sx q[2];
rz(-2.2552239) q[2];
sx q[2];
rz(0.85927187) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.75016025) q[1];
sx q[1];
rz(-2.5176628) q[1];
sx q[1];
rz(-2.0425914) q[1];
rz(1.5815758) q[3];
sx q[3];
rz(-1.618715) q[3];
sx q[3];
rz(-0.77424327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.036844) q[2];
sx q[2];
rz(-1.3354744) q[2];
sx q[2];
rz(2.74995) q[2];
rz(0.13970217) q[3];
sx q[3];
rz(-2.468686) q[3];
sx q[3];
rz(3.1203549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4166819) q[0];
sx q[0];
rz(-1.4044489) q[0];
sx q[0];
rz(-1.8649944) q[0];
rz(2.2712767) q[1];
sx q[1];
rz(-1.566193) q[1];
sx q[1];
rz(1.2044027) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6635839) q[0];
sx q[0];
rz(-1.1482129) q[0];
sx q[0];
rz(2.708486) q[0];
rz(1.3538829) q[2];
sx q[2];
rz(-2.4241944) q[2];
sx q[2];
rz(0.72620981) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7960297) q[1];
sx q[1];
rz(-1.2064762) q[1];
sx q[1];
rz(-1.6771392) q[1];
x q[2];
rz(-2.8562206) q[3];
sx q[3];
rz(-0.91802363) q[3];
sx q[3];
rz(0.6245581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5045972) q[2];
sx q[2];
rz(-2.6837139) q[2];
sx q[2];
rz(1.2191999) q[2];
rz(-0.35456625) q[3];
sx q[3];
rz(-2.6597326) q[3];
sx q[3];
rz(-1.5725117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59042674) q[0];
sx q[0];
rz(-1.5274436) q[0];
sx q[0];
rz(0.61808008) q[0];
rz(0.016050054) q[1];
sx q[1];
rz(-0.87688223) q[1];
sx q[1];
rz(-1.9504257) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58142692) q[0];
sx q[0];
rz(-1.418857) q[0];
sx q[0];
rz(1.8877958) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70706681) q[2];
sx q[2];
rz(-2.7756049) q[2];
sx q[2];
rz(2.6235839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6776059) q[1];
sx q[1];
rz(-1.364344) q[1];
sx q[1];
rz(2.9264005) q[1];
rz(-pi) q[2];
rz(0.54179811) q[3];
sx q[3];
rz(-1.2216611) q[3];
sx q[3];
rz(0.19402129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.96737343) q[2];
sx q[2];
rz(-2.1790395) q[2];
sx q[2];
rz(-1.013914) q[2];
rz(1.3570471) q[3];
sx q[3];
rz(-0.8299399) q[3];
sx q[3];
rz(-1.0323662) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8106666) q[0];
sx q[0];
rz(-1.7174915) q[0];
sx q[0];
rz(-2.9372835) q[0];
rz(1.3775685) q[1];
sx q[1];
rz(-1.7267449) q[1];
sx q[1];
rz(2.2185982) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8078634) q[0];
sx q[0];
rz(-2.6834052) q[0];
sx q[0];
rz(2.3966167) q[0];
rz(-3.0976474) q[2];
sx q[2];
rz(-1.5825795) q[2];
sx q[2];
rz(0.95566434) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.35034414) q[1];
sx q[1];
rz(-1.6748168) q[1];
sx q[1];
rz(-1.9700325) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.022425671) q[3];
sx q[3];
rz(-0.86655819) q[3];
sx q[3];
rz(1.1426329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4541645) q[2];
sx q[2];
rz(-1.5565846) q[2];
sx q[2];
rz(-0.50951177) q[2];
rz(-2.4984958) q[3];
sx q[3];
rz(-2.0431079) q[3];
sx q[3];
rz(-2.174214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61280695) q[0];
sx q[0];
rz(-1.2061773) q[0];
sx q[0];
rz(-2.2221785) q[0];
rz(2.1557504) q[1];
sx q[1];
rz(-1.7835833) q[1];
sx q[1];
rz(1.3607508) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0589941) q[0];
sx q[0];
rz(-1.3816125) q[0];
sx q[0];
rz(2.9181705) q[0];
rz(0.76061337) q[2];
sx q[2];
rz(-2.0509655) q[2];
sx q[2];
rz(1.7825356) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67140019) q[1];
sx q[1];
rz(-2.6365286) q[1];
sx q[1];
rz(-2.7027674) q[1];
rz(3.0642062) q[3];
sx q[3];
rz(-2.1688528) q[3];
sx q[3];
rz(0.070904562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3551066) q[2];
sx q[2];
rz(-1.3856709) q[2];
sx q[2];
rz(-3.0267267) q[2];
rz(0.45587513) q[3];
sx q[3];
rz(-1.3331648) q[3];
sx q[3];
rz(1.280064) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961261) q[0];
sx q[0];
rz(-1.8969314) q[0];
sx q[0];
rz(-1.8967569) q[0];
rz(-0.70760977) q[1];
sx q[1];
rz(-1.6376303) q[1];
sx q[1];
rz(0.27522603) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9305785) q[0];
sx q[0];
rz(-2.0654581) q[0];
sx q[0];
rz(-2.7633694) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49595828) q[2];
sx q[2];
rz(-2.1514116) q[2];
sx q[2];
rz(2.1883287) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.77995283) q[1];
sx q[1];
rz(-1.0859023) q[1];
sx q[1];
rz(1.1760902) q[1];
rz(-pi) q[2];
rz(0.082645881) q[3];
sx q[3];
rz(-2.189889) q[3];
sx q[3];
rz(-1.9649399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.39367166) q[2];
sx q[2];
rz(-1.5966406) q[2];
sx q[2];
rz(2.6829524) q[2];
rz(0.7115055) q[3];
sx q[3];
rz(-2.4578874) q[3];
sx q[3];
rz(2.5512364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8608619) q[0];
sx q[0];
rz(-1.9545398) q[0];
sx q[0];
rz(1.3635427) q[0];
rz(-1.4200312) q[1];
sx q[1];
rz(-2.6224711) q[1];
sx q[1];
rz(-2.4218959) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62455432) q[0];
sx q[0];
rz(-1.8554243) q[0];
sx q[0];
rz(-2.4863003) q[0];
rz(-pi) q[1];
rz(0.72549707) q[2];
sx q[2];
rz(-2.2273387) q[2];
sx q[2];
rz(3.0692435) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.81855782) q[1];
sx q[1];
rz(-0.54669387) q[1];
sx q[1];
rz(3.0562889) q[1];
x q[2];
rz(2.8119874) q[3];
sx q[3];
rz(-2.4334987) q[3];
sx q[3];
rz(-3.0081188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.65934962) q[2];
sx q[2];
rz(-1.6094001) q[2];
sx q[2];
rz(-0.68816319) q[2];
rz(3.0996389) q[3];
sx q[3];
rz(-2.0418906) q[3];
sx q[3];
rz(2.8614614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2280837) q[0];
sx q[0];
rz(-1.9788454) q[0];
sx q[0];
rz(2.9550609) q[0];
rz(-2.5371011) q[1];
sx q[1];
rz(-2.1291321) q[1];
sx q[1];
rz(1.4950745) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.625531) q[0];
sx q[0];
rz(-0.4058668) q[0];
sx q[0];
rz(0.56674515) q[0];
x q[1];
rz(2.9096793) q[2];
sx q[2];
rz(-2.4319318) q[2];
sx q[2];
rz(2.2556925) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1868134) q[1];
sx q[1];
rz(-2.8615132) q[1];
sx q[1];
rz(-2.7571452) q[1];
rz(-0.80963482) q[3];
sx q[3];
rz(-1.0757043) q[3];
sx q[3];
rz(1.1822869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8994393) q[2];
sx q[2];
rz(-0.95726761) q[2];
sx q[2];
rz(-0.014766679) q[2];
rz(1.8642558) q[3];
sx q[3];
rz(-1.3093964) q[3];
sx q[3];
rz(1.4130672) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9799141) q[0];
sx q[0];
rz(-2.4665687) q[0];
sx q[0];
rz(2.263608) q[0];
rz(0.19628482) q[1];
sx q[1];
rz(-1.9613962) q[1];
sx q[1];
rz(1.3605798) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1454865) q[0];
sx q[0];
rz(-1.6084451) q[0];
sx q[0];
rz(2.1426175) q[0];
rz(-pi) q[1];
rz(-3.0948823) q[2];
sx q[2];
rz(-2.197406) q[2];
sx q[2];
rz(-1.4057297) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0226375) q[1];
sx q[1];
rz(-1.2345018) q[1];
sx q[1];
rz(0.35204661) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1514879) q[3];
sx q[3];
rz(-1.7125704) q[3];
sx q[3];
rz(1.7820953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.517841) q[2];
sx q[2];
rz(-1.5072284) q[2];
sx q[2];
rz(-0.8927792) q[2];
rz(-0.039285224) q[3];
sx q[3];
rz(-1.4792484) q[3];
sx q[3];
rz(2.3915496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2980625) q[0];
sx q[0];
rz(-0.003304464) q[0];
sx q[0];
rz(-0.046534006) q[0];
rz(0.90905601) q[1];
sx q[1];
rz(-2.2687056) q[1];
sx q[1];
rz(-2.4216901) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2324632) q[0];
sx q[0];
rz(-1.4913519) q[0];
sx q[0];
rz(1.3264873) q[0];
rz(1.5762395) q[2];
sx q[2];
rz(-1.2619702) q[2];
sx q[2];
rz(2.960161) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.1137005) q[1];
sx q[1];
rz(-2.8942331) q[1];
sx q[1];
rz(1.5075831) q[1];
rz(-0.7474483) q[3];
sx q[3];
rz(-0.8851074) q[3];
sx q[3];
rz(1.7789343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4204734) q[2];
sx q[2];
rz(-1.7420008) q[2];
sx q[2];
rz(3.1151248) q[2];
rz(-1.3039533) q[3];
sx q[3];
rz(-0.81726685) q[3];
sx q[3];
rz(1.2560237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.7883041) q[0];
sx q[0];
rz(-1.371654) q[0];
sx q[0];
rz(-1.3375682) q[0];
rz(-0.2164671) q[1];
sx q[1];
rz(-1.6995866) q[1];
sx q[1];
rz(1.235984) q[1];
rz(1.886006) q[2];
sx q[2];
rz(-2.253058) q[2];
sx q[2];
rz(-0.81894973) q[2];
rz(0.98417102) q[3];
sx q[3];
rz(-1.5979206) q[3];
sx q[3];
rz(0.4978705) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];