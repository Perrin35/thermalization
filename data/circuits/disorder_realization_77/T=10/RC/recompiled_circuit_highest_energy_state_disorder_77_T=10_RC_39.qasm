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
rz(-0.51802975) q[0];
sx q[0];
rz(-2.0002444) q[0];
sx q[0];
rz(-0.040123392) q[0];
rz(0.24815458) q[1];
sx q[1];
rz(5.2207898) q[1];
sx q[1];
rz(9.3407486) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9548708) q[0];
sx q[0];
rz(-1.499641) q[0];
sx q[0];
rz(1.6297718) q[0];
x q[1];
rz(2.7366287) q[2];
sx q[2];
rz(-0.42879196) q[2];
sx q[2];
rz(2.2909209) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5539816) q[1];
sx q[1];
rz(-1.1003255) q[1];
sx q[1];
rz(-1.379983) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5658251) q[3];
sx q[3];
rz(-1.218349) q[3];
sx q[3];
rz(-1.5653597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6930801) q[2];
sx q[2];
rz(-1.7757519) q[2];
sx q[2];
rz(2.8903294) q[2];
rz(1.9692028) q[3];
sx q[3];
rz(-2.0386212) q[3];
sx q[3];
rz(-3.0917061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4548816) q[0];
sx q[0];
rz(-0.48191163) q[0];
sx q[0];
rz(1.7676109) q[0];
rz(2.8107367) q[1];
sx q[1];
rz(-1.4127981) q[1];
sx q[1];
rz(2.8673598) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86048638) q[0];
sx q[0];
rz(-1.8869068) q[0];
sx q[0];
rz(-1.4962026) q[0];
rz(-pi) q[1];
rz(1.577967) q[2];
sx q[2];
rz(-2.104055) q[2];
sx q[2];
rz(2.8523142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9668558) q[1];
sx q[1];
rz(-1.2464332) q[1];
sx q[1];
rz(1.3423389) q[1];
x q[2];
rz(3.0586309) q[3];
sx q[3];
rz(-1.9881691) q[3];
sx q[3];
rz(2.7132332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8809044) q[2];
sx q[2];
rz(-1.8131249) q[2];
sx q[2];
rz(-2.0534959) q[2];
rz(1.043383) q[3];
sx q[3];
rz(-0.16845307) q[3];
sx q[3];
rz(-2.8823631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.008721) q[0];
sx q[0];
rz(-2.6368124) q[0];
sx q[0];
rz(-1.1016499) q[0];
rz(0.77611008) q[1];
sx q[1];
rz(-1.2932237) q[1];
sx q[1];
rz(-2.4993842) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6519847) q[0];
sx q[0];
rz(-2.5074158) q[0];
sx q[0];
rz(2.339667) q[0];
rz(1.3234336) q[2];
sx q[2];
rz(-0.86630922) q[2];
sx q[2];
rz(-3.041628) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.57802618) q[1];
sx q[1];
rz(-1.1456923) q[1];
sx q[1];
rz(-3.0277287) q[1];
rz(2.1884584) q[3];
sx q[3];
rz(-2.1336305) q[3];
sx q[3];
rz(-1.9632393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8316101) q[2];
sx q[2];
rz(-1.9881366) q[2];
sx q[2];
rz(-0.93008053) q[2];
rz(2.3275404) q[3];
sx q[3];
rz(-1.1907153) q[3];
sx q[3];
rz(1.2306151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26152465) q[0];
sx q[0];
rz(-1.1898758) q[0];
sx q[0];
rz(0.50450182) q[0];
rz(-0.91744676) q[1];
sx q[1];
rz(-2.4023299) q[1];
sx q[1];
rz(2.9333072) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84346164) q[0];
sx q[0];
rz(-1.7195452) q[0];
sx q[0];
rz(2.896033) q[0];
rz(0.71983003) q[2];
sx q[2];
rz(-1.3583169) q[2];
sx q[2];
rz(1.302254) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5985344) q[1];
sx q[1];
rz(-2.0428791) q[1];
sx q[1];
rz(-0.92453875) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3127906) q[3];
sx q[3];
rz(-1.0258706) q[3];
sx q[3];
rz(1.2598002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8248799) q[2];
sx q[2];
rz(-2.0471639) q[2];
sx q[2];
rz(1.971395) q[2];
rz(-0.96585387) q[3];
sx q[3];
rz(-1.0713157) q[3];
sx q[3];
rz(2.9962208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.589094) q[0];
sx q[0];
rz(-0.67411244) q[0];
sx q[0];
rz(0.66459769) q[0];
rz(0.26847863) q[1];
sx q[1];
rz(-1.45739) q[1];
sx q[1];
rz(-2.1427515) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12832175) q[0];
sx q[0];
rz(-1.7230464) q[0];
sx q[0];
rz(2.1571742) q[0];
rz(-pi) q[1];
rz(-1.1469969) q[2];
sx q[2];
rz(-1.4487954) q[2];
sx q[2];
rz(1.854055) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84643572) q[1];
sx q[1];
rz(-1.4017516) q[1];
sx q[1];
rz(-3.0205932) q[1];
rz(-2.2625974) q[3];
sx q[3];
rz(-1.5546518) q[3];
sx q[3];
rz(3.0825577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9773679) q[2];
sx q[2];
rz(-0.94027844) q[2];
sx q[2];
rz(-3.0244136) q[2];
rz(-3.0818648) q[3];
sx q[3];
rz(-1.9220587) q[3];
sx q[3];
rz(-0.81502771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3672459) q[0];
sx q[0];
rz(-2.7278439) q[0];
sx q[0];
rz(-1.2600979) q[0];
rz(-0.11481181) q[1];
sx q[1];
rz(-1.3453307) q[1];
sx q[1];
rz(-0.095349163) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0336825) q[0];
sx q[0];
rz(-2.9386407) q[0];
sx q[0];
rz(0.072921948) q[0];
x q[1];
rz(2.0748069) q[2];
sx q[2];
rz(-1.4246449) q[2];
sx q[2];
rz(-2.604217) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5231875) q[1];
sx q[1];
rz(-1.9908991) q[1];
sx q[1];
rz(-1.0805431) q[1];
rz(-pi) q[2];
rz(-2.0982025) q[3];
sx q[3];
rz(-1.5800484) q[3];
sx q[3];
rz(-1.9394138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0872588) q[2];
sx q[2];
rz(-0.91780535) q[2];
sx q[2];
rz(-2.4799662) q[2];
rz(2.3719487) q[3];
sx q[3];
rz(-1.4504284) q[3];
sx q[3];
rz(-0.86696398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0041644) q[0];
sx q[0];
rz(-0.34808174) q[0];
sx q[0];
rz(2.3837756) q[0];
rz(-0.025253145) q[1];
sx q[1];
rz(-0.39552894) q[1];
sx q[1];
rz(-1.1753561) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3536153) q[0];
sx q[0];
rz(-0.60871658) q[0];
sx q[0];
rz(2.489021) q[0];
rz(-0.48904959) q[2];
sx q[2];
rz(-2.1903775) q[2];
sx q[2];
rz(-1.3330158) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0569799) q[1];
sx q[1];
rz(-2.2723722) q[1];
sx q[1];
rz(0.52053156) q[1];
rz(-1.5937012) q[3];
sx q[3];
rz(-1.9996907) q[3];
sx q[3];
rz(-2.1351086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.56388277) q[2];
sx q[2];
rz(-1.6174199) q[2];
sx q[2];
rz(1.2152524) q[2];
rz(0.75872129) q[3];
sx q[3];
rz(-1.4164378) q[3];
sx q[3];
rz(-1.5569713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6587081) q[0];
sx q[0];
rz(-1.2889129) q[0];
sx q[0];
rz(2.7523852) q[0];
rz(3.0534741) q[1];
sx q[1];
rz(-1.7033109) q[1];
sx q[1];
rz(1.5315936) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9154598) q[0];
sx q[0];
rz(-0.32308871) q[0];
sx q[0];
rz(2.3996572) q[0];
rz(-pi) q[1];
x q[1];
rz(0.66122239) q[2];
sx q[2];
rz(-1.2031789) q[2];
sx q[2];
rz(-1.6501381) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4952505) q[1];
sx q[1];
rz(-1.6227437) q[1];
sx q[1];
rz(2.6990128) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1346899) q[3];
sx q[3];
rz(-2.3824771) q[3];
sx q[3];
rz(1.5576943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4072998) q[2];
sx q[2];
rz(-2.5278842) q[2];
sx q[2];
rz(1.3235922) q[2];
rz(1.7981516) q[3];
sx q[3];
rz(-0.7898134) q[3];
sx q[3];
rz(2.8048803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4765428) q[0];
sx q[0];
rz(-0.90710586) q[0];
sx q[0];
rz(-0.46515775) q[0];
rz(3.0910885) q[1];
sx q[1];
rz(-0.89021325) q[1];
sx q[1];
rz(-2.6506298) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57632724) q[0];
sx q[0];
rz(-2.6748228) q[0];
sx q[0];
rz(0.76961036) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7247164) q[2];
sx q[2];
rz(-1.2668005) q[2];
sx q[2];
rz(-0.15862578) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95066324) q[1];
sx q[1];
rz(-1.0196389) q[1];
sx q[1];
rz(-2.3930156) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8601362) q[3];
sx q[3];
rz(-1.6773684) q[3];
sx q[3];
rz(-0.64079277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4831873) q[2];
sx q[2];
rz(-2.5453973) q[2];
sx q[2];
rz(2.2779706) q[2];
rz(2.9577799) q[3];
sx q[3];
rz(-0.96537662) q[3];
sx q[3];
rz(-2.1553154) q[3];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42613906) q[0];
sx q[0];
rz(-1.5197536) q[0];
sx q[0];
rz(2.5157628) q[0];
rz(-2.4848056) q[1];
sx q[1];
rz(-0.96581179) q[1];
sx q[1];
rz(-1.2084557) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6656553) q[0];
sx q[0];
rz(-1.6724419) q[0];
sx q[0];
rz(1.7045506) q[0];
rz(-pi) q[1];
rz(2.7371128) q[2];
sx q[2];
rz(-2.4411429) q[2];
sx q[2];
rz(-0.03534783) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0565666) q[1];
sx q[1];
rz(-1.657227) q[1];
sx q[1];
rz(1.5734476) q[1];
rz(-pi) q[2];
rz(-2.5700611) q[3];
sx q[3];
rz(-1.4213955) q[3];
sx q[3];
rz(2.7461441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5586231) q[2];
sx q[2];
rz(-1.308459) q[2];
sx q[2];
rz(-2.6226131) q[2];
rz(2.6988622) q[3];
sx q[3];
rz(-1.3230007) q[3];
sx q[3];
rz(-1.7980661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2124355) q[0];
sx q[0];
rz(-2.1514308) q[0];
sx q[0];
rz(1.9707752) q[0];
rz(-0.15405542) q[1];
sx q[1];
rz(-0.93692056) q[1];
sx q[1];
rz(-0.9137203) q[1];
rz(-3.0689756) q[2];
sx q[2];
rz(-2.6238863) q[2];
sx q[2];
rz(-3.1095502) q[2];
rz(-2.7967706) q[3];
sx q[3];
rz(-0.92690868) q[3];
sx q[3];
rz(2.3186701) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
