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
rz(0.69819063) q[0];
sx q[0];
rz(-0.31261045) q[0];
sx q[0];
rz(2.1460549) q[0];
rz(-2.0350463) q[1];
sx q[1];
rz(-2.3732329) q[1];
sx q[1];
rz(0.33222693) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26259188) q[0];
sx q[0];
rz(-2.342412) q[0];
sx q[0];
rz(0.57882092) q[0];
rz(-pi) q[1];
x q[1];
rz(2.034445) q[2];
sx q[2];
rz(-1.4476748) q[2];
sx q[2];
rz(2.5277918) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2549887) q[1];
sx q[1];
rz(-2.0333148) q[1];
sx q[1];
rz(-0.47420331) q[1];
x q[2];
rz(0.51789051) q[3];
sx q[3];
rz(-2.0841408) q[3];
sx q[3];
rz(2.1012517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.60653162) q[2];
sx q[2];
rz(-1.0142832) q[2];
sx q[2];
rz(3.0895341) q[2];
rz(-2.0590797) q[3];
sx q[3];
rz(-2.1942997) q[3];
sx q[3];
rz(-1.989971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.624619) q[0];
sx q[0];
rz(-0.032051429) q[0];
sx q[0];
rz(-0.33729851) q[0];
rz(0.15268046) q[1];
sx q[1];
rz(-2.3744507) q[1];
sx q[1];
rz(-1.5229567) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3169976) q[0];
sx q[0];
rz(-2.499541) q[0];
sx q[0];
rz(2.8606728) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27045336) q[2];
sx q[2];
rz(-2.4861397) q[2];
sx q[2];
rz(1.7846817) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.26451787) q[1];
sx q[1];
rz(-1.3691069) q[1];
sx q[1];
rz(-1.1796477) q[1];
rz(0.60373276) q[3];
sx q[3];
rz(-1.4430178) q[3];
sx q[3];
rz(-1.4462245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1107669) q[2];
sx q[2];
rz(-0.36588565) q[2];
sx q[2];
rz(-2.7552674) q[2];
rz(-1.857916) q[3];
sx q[3];
rz(-1.9648896) q[3];
sx q[3];
rz(1.2181237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074742643) q[0];
sx q[0];
rz(-2.061494) q[0];
sx q[0];
rz(-1.0190438) q[0];
rz(-2.3661803) q[1];
sx q[1];
rz(-1.2762028) q[1];
sx q[1];
rz(0.48616854) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2371444) q[0];
sx q[0];
rz(-1.6021361) q[0];
sx q[0];
rz(-1.5678568) q[0];
rz(2.3535291) q[2];
sx q[2];
rz(-1.9275086) q[2];
sx q[2];
rz(-2.8887088) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0224494) q[1];
sx q[1];
rz(-2.0538985) q[1];
sx q[1];
rz(1.3234322) q[1];
rz(-pi) q[2];
rz(-1.8649549) q[3];
sx q[3];
rz(-1.2614864) q[3];
sx q[3];
rz(0.33634392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7458618) q[2];
sx q[2];
rz(-1.1660601) q[2];
sx q[2];
rz(-0.80500785) q[2];
rz(-0.65195596) q[3];
sx q[3];
rz(-2.0396353) q[3];
sx q[3];
rz(-2.2074047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66251278) q[0];
sx q[0];
rz(-1.3211687) q[0];
sx q[0];
rz(1.3385734) q[0];
rz(-1.5486807) q[1];
sx q[1];
rz(-0.70392307) q[1];
sx q[1];
rz(-0.42565638) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073624728) q[0];
sx q[0];
rz(-1.5626161) q[0];
sx q[0];
rz(0.091853022) q[0];
rz(-pi) q[1];
rz(-1.6792308) q[2];
sx q[2];
rz(-2.8320356) q[2];
sx q[2];
rz(2.5029687) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7158) q[1];
sx q[1];
rz(-1.9746026) q[1];
sx q[1];
rz(-2.3753662) q[1];
x q[2];
rz(0.32069499) q[3];
sx q[3];
rz(-2.3983722) q[3];
sx q[3];
rz(-0.37494613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.272133) q[2];
sx q[2];
rz(-1.774186) q[2];
sx q[2];
rz(0.41932219) q[2];
rz(-1.9648633) q[3];
sx q[3];
rz(-0.71507016) q[3];
sx q[3];
rz(0.56593219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45593539) q[0];
sx q[0];
rz(-2.3929907) q[0];
sx q[0];
rz(-1.5022044) q[0];
rz(1.4337076) q[1];
sx q[1];
rz(-1.8610443) q[1];
sx q[1];
rz(-0.40275231) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1486738) q[0];
sx q[0];
rz(-2.4687827) q[0];
sx q[0];
rz(-1.1514173) q[0];
x q[1];
rz(-3.0482015) q[2];
sx q[2];
rz(-0.56677188) q[2];
sx q[2];
rz(3.0393578) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.50583831) q[1];
sx q[1];
rz(-0.16316667) q[1];
sx q[1];
rz(2.5058305) q[1];
x q[2];
rz(0.66688315) q[3];
sx q[3];
rz(-1.6630238) q[3];
sx q[3];
rz(0.86574829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6413642) q[2];
sx q[2];
rz(-1.0470752) q[2];
sx q[2];
rz(-2.8889528) q[2];
rz(0.80444515) q[3];
sx q[3];
rz(-0.52144709) q[3];
sx q[3];
rz(2.6430602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8135524) q[0];
sx q[0];
rz(-3.033162) q[0];
sx q[0];
rz(0.88388467) q[0];
rz(0.49993316) q[1];
sx q[1];
rz(-1.4686613) q[1];
sx q[1];
rz(-0.51441851) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.894569) q[0];
sx q[0];
rz(-1.2366364) q[0];
sx q[0];
rz(2.1804125) q[0];
x q[1];
rz(-2.7398501) q[2];
sx q[2];
rz(-0.56890024) q[2];
sx q[2];
rz(0.57050975) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0757671) q[1];
sx q[1];
rz(-1.3525241) q[1];
sx q[1];
rz(1.9820007) q[1];
x q[2];
rz(2.3522369) q[3];
sx q[3];
rz(-1.5522) q[3];
sx q[3];
rz(-2.9531053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1540692) q[2];
sx q[2];
rz(-1.6776626) q[2];
sx q[2];
rz(-0.12518159) q[2];
rz(0.82031885) q[3];
sx q[3];
rz(-1.9671974) q[3];
sx q[3];
rz(-0.34631795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5888551) q[0];
sx q[0];
rz(-0.6468361) q[0];
sx q[0];
rz(-3.0297025) q[0];
rz(-1.1397859) q[1];
sx q[1];
rz(-2.9264989) q[1];
sx q[1];
rz(-1.2824167) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70045602) q[0];
sx q[0];
rz(-2.4917291) q[0];
sx q[0];
rz(1.1101686) q[0];
rz(-pi) q[1];
rz(1.4256023) q[2];
sx q[2];
rz(-1.34158) q[2];
sx q[2];
rz(-2.9289233) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0541898) q[1];
sx q[1];
rz(-2.3034181) q[1];
sx q[1];
rz(-2.0778632) q[1];
x q[2];
rz(2.4573648) q[3];
sx q[3];
rz(-1.6155532) q[3];
sx q[3];
rz(-1.0563204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.68975) q[2];
sx q[2];
rz(-0.30561438) q[2];
sx q[2];
rz(1.3172733) q[2];
rz(3.1210323) q[3];
sx q[3];
rz(-0.33187425) q[3];
sx q[3];
rz(-2.1338972) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9356215) q[0];
sx q[0];
rz(-0.99140778) q[0];
sx q[0];
rz(0.29888612) q[0];
rz(-1.0609974) q[1];
sx q[1];
rz(-1.7540878) q[1];
sx q[1];
rz(1.908173) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73288146) q[0];
sx q[0];
rz(-1.2378843) q[0];
sx q[0];
rz(-3.0391418) q[0];
rz(1.685254) q[2];
sx q[2];
rz(-1.6117764) q[2];
sx q[2];
rz(-1.51454) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.778925) q[1];
sx q[1];
rz(-0.92955076) q[1];
sx q[1];
rz(-2.769313) q[1];
rz(-pi) q[2];
rz(-2.4806907) q[3];
sx q[3];
rz(-1.9993625) q[3];
sx q[3];
rz(2.1235583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.79248205) q[2];
sx q[2];
rz(-2.0439549) q[2];
sx q[2];
rz(-0.20723542) q[2];
rz(-0.66050291) q[3];
sx q[3];
rz(-2.1410172) q[3];
sx q[3];
rz(-1.8628666) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5527363) q[0];
sx q[0];
rz(-1.4016466) q[0];
sx q[0];
rz(1.2218342) q[0];
rz(-0.51512042) q[1];
sx q[1];
rz(-3.0312067) q[1];
sx q[1];
rz(2.5882904) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3344655) q[0];
sx q[0];
rz(-0.57274023) q[0];
sx q[0];
rz(-1.3329026) q[0];
rz(-pi) q[1];
rz(1.4261999) q[2];
sx q[2];
rz(-0.42486496) q[2];
sx q[2];
rz(-2.1230842) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8043704) q[1];
sx q[1];
rz(-0.17717277) q[1];
sx q[1];
rz(-1.5371662) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5946383) q[3];
sx q[3];
rz(-0.90180221) q[3];
sx q[3];
rz(3.075595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.33636609) q[2];
sx q[2];
rz(-1.1104501) q[2];
sx q[2];
rz(-2.6160348) q[2];
rz(-1.6622274) q[3];
sx q[3];
rz(-2.1879304) q[3];
sx q[3];
rz(0.76013887) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0568327) q[0];
sx q[0];
rz(-2.3469717) q[0];
sx q[0];
rz(0.42924616) q[0];
rz(-1.5224573) q[1];
sx q[1];
rz(-1.2702962) q[1];
sx q[1];
rz(2.2116908) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4711205) q[0];
sx q[0];
rz(-1.7848624) q[0];
sx q[0];
rz(-1.2250958) q[0];
x q[1];
rz(0.11008628) q[2];
sx q[2];
rz(-1.7890245) q[2];
sx q[2];
rz(2.9492117) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.20969756) q[1];
sx q[1];
rz(-1.731809) q[1];
sx q[1];
rz(-1.8568532) q[1];
rz(3.0449379) q[3];
sx q[3];
rz(-0.55516637) q[3];
sx q[3];
rz(0.47278178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4055736) q[2];
sx q[2];
rz(-0.47601998) q[2];
sx q[2];
rz(2.948577) q[2];
rz(0.080282601) q[3];
sx q[3];
rz(-2.2716227) q[3];
sx q[3];
rz(-1.75753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3748462) q[0];
sx q[0];
rz(-1.9369047) q[0];
sx q[0];
rz(1.9932224) q[0];
rz(2.171352) q[1];
sx q[1];
rz(-1.2087676) q[1];
sx q[1];
rz(-1.2420775) q[1];
rz(-1.8850897) q[2];
sx q[2];
rz(-2.3705924) q[2];
sx q[2];
rz(2.3938897) q[2];
rz(3.0474256) q[3];
sx q[3];
rz(-1.5945024) q[3];
sx q[3];
rz(1.9645346) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
