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
rz(1.053831) q[0];
sx q[0];
rz(3.6917917) q[0];
sx q[0];
rz(10.79296) q[0];
rz(-0.47222459) q[1];
sx q[1];
rz(-3.0250186) q[1];
sx q[1];
rz(-0.20224686) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8919173) q[0];
sx q[0];
rz(-1.37171) q[0];
sx q[0];
rz(2.9249935) q[0];
rz(-pi) q[1];
rz(-0.80533452) q[2];
sx q[2];
rz(-2.4676358) q[2];
sx q[2];
rz(-1.0588624) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3766915) q[1];
sx q[1];
rz(-2.4623532) q[1];
sx q[1];
rz(2.0308073) q[1];
rz(-2.0368882) q[3];
sx q[3];
rz(-1.5370318) q[3];
sx q[3];
rz(-1.0068144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54174417) q[2];
sx q[2];
rz(-0.84459633) q[2];
sx q[2];
rz(-2.2005626) q[2];
rz(-2.8372676) q[3];
sx q[3];
rz(-0.86186886) q[3];
sx q[3];
rz(2.8743675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32352725) q[0];
sx q[0];
rz(-2.7637988) q[0];
sx q[0];
rz(0.7793119) q[0];
rz(1.7794973) q[1];
sx q[1];
rz(-2.7423488) q[1];
sx q[1];
rz(-1.13387) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2529566) q[0];
sx q[0];
rz(-2.0042814) q[0];
sx q[0];
rz(-3.033328) q[0];
x q[1];
rz(-0.46102779) q[2];
sx q[2];
rz(-0.13153464) q[2];
sx q[2];
rz(1.5054838) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.459455) q[1];
sx q[1];
rz(-1.3054779) q[1];
sx q[1];
rz(-3.0299761) q[1];
rz(-pi) q[2];
rz(0.00098382424) q[3];
sx q[3];
rz(-1.0615225) q[3];
sx q[3];
rz(0.60722095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4506932) q[2];
sx q[2];
rz(-1.0319812) q[2];
sx q[2];
rz(2.9449985) q[2];
rz(-2.172566) q[3];
sx q[3];
rz(-1.6220379) q[3];
sx q[3];
rz(-1.903418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75152385) q[0];
sx q[0];
rz(-0.5539493) q[0];
sx q[0];
rz(-0.73177904) q[0];
rz(0.36830184) q[1];
sx q[1];
rz(-2.383547) q[1];
sx q[1];
rz(0.36645737) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55933773) q[0];
sx q[0];
rz(-0.32024239) q[0];
sx q[0];
rz(0.59424627) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.75080504) q[2];
sx q[2];
rz(-2.8158203) q[2];
sx q[2];
rz(0.70394403) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71354564) q[1];
sx q[1];
rz(-0.69880345) q[1];
sx q[1];
rz(2.2627463) q[1];
rz(-2.5440574) q[3];
sx q[3];
rz(-0.96216494) q[3];
sx q[3];
rz(0.017692117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9590108) q[2];
sx q[2];
rz(-2.3704447) q[2];
sx q[2];
rz(-0.89243531) q[2];
rz(-2.6871032) q[3];
sx q[3];
rz(-2.317704) q[3];
sx q[3];
rz(0.23196001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7682122) q[0];
sx q[0];
rz(-0.1425655) q[0];
sx q[0];
rz(2.7349803) q[0];
rz(-0.82798249) q[1];
sx q[1];
rz(-1.7384638) q[1];
sx q[1];
rz(0.83555317) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31553651) q[0];
sx q[0];
rz(-0.18504194) q[0];
sx q[0];
rz(-1.3547784) q[0];
rz(-pi) q[1];
rz(1.2242975) q[2];
sx q[2];
rz(-1.76148) q[2];
sx q[2];
rz(-2.4257223) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9367076) q[1];
sx q[1];
rz(-0.46006535) q[1];
sx q[1];
rz(-2.8060444) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1855679) q[3];
sx q[3];
rz(-1.6611757) q[3];
sx q[3];
rz(-0.65672311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7416209) q[2];
sx q[2];
rz(-2.8394832) q[2];
sx q[2];
rz(-0.92140222) q[2];
rz(-1.7737927) q[3];
sx q[3];
rz(-2.1801345) q[3];
sx q[3];
rz(0.14149806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88857404) q[0];
sx q[0];
rz(-1.0971917) q[0];
sx q[0];
rz(-0.15884037) q[0];
rz(-1.6244434) q[1];
sx q[1];
rz(-0.30888638) q[1];
sx q[1];
rz(1.3295757) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81013524) q[0];
sx q[0];
rz(-1.5525426) q[0];
sx q[0];
rz(1.5601331) q[0];
rz(-pi) q[1];
rz(1.9407523) q[2];
sx q[2];
rz(-0.59365497) q[2];
sx q[2];
rz(-1.4482738) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.5097039) q[1];
sx q[1];
rz(-1.77138) q[1];
sx q[1];
rz(1.1840118) q[1];
rz(-pi) q[2];
rz(1.8908126) q[3];
sx q[3];
rz(-2.5502) q[3];
sx q[3];
rz(-1.0726014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.57881957) q[2];
sx q[2];
rz(-0.69391888) q[2];
sx q[2];
rz(-2.5895183) q[2];
rz(-0.79365802) q[3];
sx q[3];
rz(-2.5439883) q[3];
sx q[3];
rz(2.1551267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3119222) q[0];
sx q[0];
rz(-2.2770918) q[0];
sx q[0];
rz(-0.70575869) q[0];
rz(-0.13748473) q[1];
sx q[1];
rz(-2.4973713) q[1];
sx q[1];
rz(2.4286043) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1437837) q[0];
sx q[0];
rz(-1.9558055) q[0];
sx q[0];
rz(1.817817) q[0];
rz(-pi) q[1];
rz(0.51432864) q[2];
sx q[2];
rz(-1.1900969) q[2];
sx q[2];
rz(-1.2027539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9116152) q[1];
sx q[1];
rz(-1.4739081) q[1];
sx q[1];
rz(-1.1712892) q[1];
rz(-pi) q[2];
rz(-1.8637795) q[3];
sx q[3];
rz(-1.6821386) q[3];
sx q[3];
rz(0.35929832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2116427) q[2];
sx q[2];
rz(-0.86923081) q[2];
sx q[2];
rz(-2.8575274) q[2];
rz(-3.0197213) q[3];
sx q[3];
rz(-0.28212306) q[3];
sx q[3];
rz(-1.6596644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26509869) q[0];
sx q[0];
rz(-0.6186741) q[0];
sx q[0];
rz(2.4342243) q[0];
rz(-0.44037285) q[1];
sx q[1];
rz(-0.71763867) q[1];
sx q[1];
rz(-1.7792938) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4934568) q[0];
sx q[0];
rz(-1.3246944) q[0];
sx q[0];
rz(1.2854693) q[0];
rz(0.42476219) q[2];
sx q[2];
rz(-0.889689) q[2];
sx q[2];
rz(-0.29007888) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.728546) q[1];
sx q[1];
rz(-0.73592454) q[1];
sx q[1];
rz(-0.82637431) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23162095) q[3];
sx q[3];
rz(-2.6745788) q[3];
sx q[3];
rz(1.2488332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.97544396) q[2];
sx q[2];
rz(-2.7232309) q[2];
sx q[2];
rz(-2.5725906) q[2];
rz(-2.1042018) q[3];
sx q[3];
rz(-0.85809696) q[3];
sx q[3];
rz(0.4666127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5910832) q[0];
sx q[0];
rz(-2.789848) q[0];
sx q[0];
rz(1.9367223) q[0];
rz(0.51271802) q[1];
sx q[1];
rz(-2.3725489) q[1];
sx q[1];
rz(0.18096322) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9009243) q[0];
sx q[0];
rz(-0.5597194) q[0];
sx q[0];
rz(-2.3842035) q[0];
rz(0.36278533) q[2];
sx q[2];
rz(-0.37172592) q[2];
sx q[2];
rz(0.065601018) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9258283) q[1];
sx q[1];
rz(-1.3604394) q[1];
sx q[1];
rz(-2.1390361) q[1];
x q[2];
rz(-1.4970786) q[3];
sx q[3];
rz(-1.4604919) q[3];
sx q[3];
rz(-3.0347787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.60244954) q[2];
sx q[2];
rz(-0.64843833) q[2];
sx q[2];
rz(-0.096972018) q[2];
rz(-0.55475956) q[3];
sx q[3];
rz(-1.3223038) q[3];
sx q[3];
rz(-0.35840148) q[3];
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
rz(2.6612514) q[0];
sx q[0];
rz(-2.3619409) q[0];
sx q[0];
rz(-0.10777792) q[0];
rz(2.5906471) q[1];
sx q[1];
rz(-2.7007553) q[1];
sx q[1];
rz(-1.617618) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11883987) q[0];
sx q[0];
rz(-1.700239) q[0];
sx q[0];
rz(-1.3368827) q[0];
rz(-pi) q[1];
rz(0.6438821) q[2];
sx q[2];
rz(-1.1538299) q[2];
sx q[2];
rz(-1.8089) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2366166) q[1];
sx q[1];
rz(-2.0324685) q[1];
sx q[1];
rz(-3.0932337) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49160853) q[3];
sx q[3];
rz(-1.6171241) q[3];
sx q[3];
rz(-2.2674048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.23065755) q[2];
sx q[2];
rz(-0.50369889) q[2];
sx q[2];
rz(-1.1304193) q[2];
rz(2.9039827) q[3];
sx q[3];
rz(-2.7311324) q[3];
sx q[3];
rz(2.3835278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0216574) q[0];
sx q[0];
rz(-0.1788685) q[0];
sx q[0];
rz(-2.2005431) q[0];
rz(2.7811116) q[1];
sx q[1];
rz(-1.7345813) q[1];
sx q[1];
rz(-1.2233268) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1038641) q[0];
sx q[0];
rz(-2.1489804) q[0];
sx q[0];
rz(1.8072771) q[0];
rz(-pi) q[1];
rz(2.1603196) q[2];
sx q[2];
rz(-0.86893493) q[2];
sx q[2];
rz(1.5379932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.66050038) q[1];
sx q[1];
rz(-1.6410488) q[1];
sx q[1];
rz(0.41892799) q[1];
rz(-0.19630614) q[3];
sx q[3];
rz(-2.098791) q[3];
sx q[3];
rz(1.7332026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2337522) q[2];
sx q[2];
rz(-1.3774104) q[2];
sx q[2];
rz(-3.0506548) q[2];
rz(-0.20445538) q[3];
sx q[3];
rz(-2.3369868) q[3];
sx q[3];
rz(1.0819775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7627056) q[0];
sx q[0];
rz(-1.611041) q[0];
sx q[0];
rz(1.8086717) q[0];
rz(1.4270205) q[1];
sx q[1];
rz(-0.49976977) q[1];
sx q[1];
rz(1.5907092) q[1];
rz(2.0841523) q[2];
sx q[2];
rz(-2.3427137) q[2];
sx q[2];
rz(-2.8165934) q[2];
rz(-1.6283725) q[3];
sx q[3];
rz(-1.3575469) q[3];
sx q[3];
rz(-1.224106) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
