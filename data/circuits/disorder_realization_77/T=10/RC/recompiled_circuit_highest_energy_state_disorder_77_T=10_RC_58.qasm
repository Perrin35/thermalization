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
rz(-1.0623955) q[1];
sx q[1];
rz(-0.084029347) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18672189) q[0];
sx q[0];
rz(-1.6419517) q[0];
sx q[0];
rz(1.5118208) q[0];
rz(-pi) q[1];
rz(-2.7366287) q[2];
sx q[2];
rz(-2.7128007) q[2];
sx q[2];
rz(-0.85067174) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1511615) q[1];
sx q[1];
rz(-0.50499454) q[1];
sx q[1];
rz(-2.7846365) q[1];
rz(-2.5474882) q[3];
sx q[3];
rz(-0.66451529) q[3];
sx q[3];
rz(-2.6474109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6930801) q[2];
sx q[2];
rz(-1.7757519) q[2];
sx q[2];
rz(-2.8903294) q[2];
rz(1.9692028) q[3];
sx q[3];
rz(-2.0386212) q[3];
sx q[3];
rz(0.049886543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6867111) q[0];
sx q[0];
rz(-0.48191163) q[0];
sx q[0];
rz(1.3739817) q[0];
rz(-0.33085597) q[1];
sx q[1];
rz(-1.7287946) q[1];
sx q[1];
rz(-2.8673598) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4080547) q[0];
sx q[0];
rz(-1.499905) q[0];
sx q[0];
rz(-2.8246585) q[0];
rz(-pi) q[1];
x q[1];
rz(1.577967) q[2];
sx q[2];
rz(-1.0375377) q[2];
sx q[2];
rz(-2.8523142) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.80503217) q[1];
sx q[1];
rz(-2.7471886) q[1];
sx q[1];
rz(0.59275643) q[1];
rz(1.1521454) q[3];
sx q[3];
rz(-1.4949706) q[3];
sx q[3];
rz(1.1087429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8809044) q[2];
sx q[2];
rz(-1.8131249) q[2];
sx q[2];
rz(-2.0534959) q[2];
rz(-2.0982096) q[3];
sx q[3];
rz(-0.16845307) q[3];
sx q[3];
rz(-2.8823631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1328717) q[0];
sx q[0];
rz(-0.50478029) q[0];
sx q[0];
rz(-2.0399427) q[0];
rz(-0.77611008) q[1];
sx q[1];
rz(-1.848369) q[1];
sx q[1];
rz(-2.4993842) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74321824) q[0];
sx q[0];
rz(-1.1461598) q[0];
sx q[0];
rz(2.0570799) q[0];
rz(-pi) q[1];
rz(1.8181591) q[2];
sx q[2];
rz(-0.86630922) q[2];
sx q[2];
rz(-0.099964634) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8485198) q[1];
sx q[1];
rz(-2.7024033) q[1];
sx q[1];
rz(-1.8166914) q[1];
rz(0.65861838) q[3];
sx q[3];
rz(-1.0589387) q[3];
sx q[3];
rz(3.1114674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8316101) q[2];
sx q[2];
rz(-1.9881366) q[2];
sx q[2];
rz(-2.2115121) q[2];
rz(2.3275404) q[3];
sx q[3];
rz(-1.9508773) q[3];
sx q[3];
rz(-1.2306151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.26152465) q[0];
sx q[0];
rz(-1.9517169) q[0];
sx q[0];
rz(2.6370908) q[0];
rz(-2.2241459) q[1];
sx q[1];
rz(-0.73926273) q[1];
sx q[1];
rz(-0.20828542) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69021001) q[0];
sx q[0];
rz(-1.3280032) q[0];
sx q[0];
rz(-1.4175178) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8253352) q[2];
sx q[2];
rz(-0.74511792) q[2];
sx q[2];
rz(3.1090914) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.69729561) q[1];
sx q[1];
rz(-1.0047067) q[1];
sx q[1];
rz(2.5725911) q[1];
x q[2];
rz(-0.39842968) q[3];
sx q[3];
rz(-0.59729415) q[3];
sx q[3];
rz(0.78890991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8248799) q[2];
sx q[2];
rz(-1.0944288) q[2];
sx q[2];
rz(-1.971395) q[2];
rz(-0.96585387) q[3];
sx q[3];
rz(-1.0713157) q[3];
sx q[3];
rz(-0.14537183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55249864) q[0];
sx q[0];
rz(-2.4674802) q[0];
sx q[0];
rz(-2.476995) q[0];
rz(2.873114) q[1];
sx q[1];
rz(-1.45739) q[1];
sx q[1];
rz(2.1427515) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3420606) q[0];
sx q[0];
rz(-2.149509) q[0];
sx q[0];
rz(0.18216746) q[0];
rz(-pi) q[1];
rz(-1.9945958) q[2];
sx q[2];
rz(-1.6927973) q[2];
sx q[2];
rz(1.854055) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4376862) q[1];
sx q[1];
rz(-1.4515299) q[1];
sx q[1];
rz(1.7410623) q[1];
rz(-pi) q[2];
x q[2];
rz(0.020962997) q[3];
sx q[3];
rz(-0.87910324) q[3];
sx q[3];
rz(1.5251336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9773679) q[2];
sx q[2];
rz(-0.94027844) q[2];
sx q[2];
rz(-0.11717907) q[2];
rz(0.059727877) q[3];
sx q[3];
rz(-1.2195339) q[3];
sx q[3];
rz(-2.3265649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77434671) q[0];
sx q[0];
rz(-2.7278439) q[0];
sx q[0];
rz(1.8814948) q[0];
rz(-0.11481181) q[1];
sx q[1];
rz(-1.7962619) q[1];
sx q[1];
rz(0.095349163) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60854429) q[0];
sx q[0];
rz(-1.5854821) q[0];
sx q[0];
rz(2.9391654) q[0];
x q[1];
rz(2.9750455) q[2];
sx q[2];
rz(-2.0689365) q[2];
sx q[2];
rz(-2.0280251) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5231875) q[1];
sx q[1];
rz(-1.1506936) q[1];
sx q[1];
rz(-2.0610496) q[1];
x q[2];
rz(-1.5891777) q[3];
sx q[3];
rz(-0.52747969) q[3];
sx q[3];
rz(0.35273409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0872588) q[2];
sx q[2];
rz(-2.2237873) q[2];
sx q[2];
rz(-2.4799662) q[2];
rz(0.76964393) q[3];
sx q[3];
rz(-1.4504284) q[3];
sx q[3];
rz(-2.2746287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1374283) q[0];
sx q[0];
rz(-0.34808174) q[0];
sx q[0];
rz(-0.75781703) q[0];
rz(-0.025253145) q[1];
sx q[1];
rz(-0.39552894) q[1];
sx q[1];
rz(-1.1753561) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7879774) q[0];
sx q[0];
rz(-0.60871658) q[0];
sx q[0];
rz(2.489021) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1532159) q[2];
sx q[2];
rz(-0.76887956) q[2];
sx q[2];
rz(1.066756) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.35835207) q[1];
sx q[1];
rz(-0.84642998) q[1];
sx q[1];
rz(1.0388166) q[1];
rz(-pi) q[2];
rz(0.42899363) q[3];
sx q[3];
rz(-1.5499664) q[3];
sx q[3];
rz(0.57383895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5777099) q[2];
sx q[2];
rz(-1.6174199) q[2];
sx q[2];
rz(-1.2152524) q[2];
rz(2.3828714) q[3];
sx q[3];
rz(-1.7251549) q[3];
sx q[3];
rz(-1.5569713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6587081) q[0];
sx q[0];
rz(-1.2889129) q[0];
sx q[0];
rz(-2.7523852) q[0];
rz(-3.0534741) q[1];
sx q[1];
rz(-1.7033109) q[1];
sx q[1];
rz(1.6099991) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9154598) q[0];
sx q[0];
rz(-2.8185039) q[0];
sx q[0];
rz(0.74193546) q[0];
rz(-pi) q[1];
rz(-0.56014748) q[2];
sx q[2];
rz(-2.3986926) q[2];
sx q[2];
rz(-2.7882238) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1080146) q[1];
sx q[1];
rz(-0.44541767) q[1];
sx q[1];
rz(0.12081318) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0069028) q[3];
sx q[3];
rz(-0.75911555) q[3];
sx q[3];
rz(-1.5838983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7342928) q[2];
sx q[2];
rz(-2.5278842) q[2];
sx q[2];
rz(1.3235922) q[2];
rz(1.7981516) q[3];
sx q[3];
rz(-2.3517793) q[3];
sx q[3];
rz(0.33671236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4765428) q[0];
sx q[0];
rz(-0.90710586) q[0];
sx q[0];
rz(-2.6764349) q[0];
rz(-0.05050412) q[1];
sx q[1];
rz(-2.2513794) q[1];
sx q[1];
rz(-0.49096289) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5652654) q[0];
sx q[0];
rz(-0.46676985) q[0];
sx q[0];
rz(-0.76961036) q[0];
rz(-pi) q[1];
rz(-1.4168763) q[2];
sx q[2];
rz(-1.2668005) q[2];
sx q[2];
rz(0.15862578) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9742187) q[1];
sx q[1];
rz(-0.95229665) q[1];
sx q[1];
rz(2.2688686) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36763962) q[3];
sx q[3];
rz(-0.30045569) q[3];
sx q[3];
rz(-2.5641172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4831873) q[2];
sx q[2];
rz(-0.5961954) q[2];
sx q[2];
rz(-0.8636221) q[2];
rz(-0.18381271) q[3];
sx q[3];
rz(-0.96537662) q[3];
sx q[3];
rz(0.98627728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.42613906) q[0];
sx q[0];
rz(-1.621839) q[0];
sx q[0];
rz(0.62582985) q[0];
rz(-0.65678701) q[1];
sx q[1];
rz(-2.1757809) q[1];
sx q[1];
rz(1.9331369) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6656553) q[0];
sx q[0];
rz(-1.6724419) q[0];
sx q[0];
rz(1.4370421) q[0];
x q[1];
rz(0.65932806) q[2];
sx q[2];
rz(-1.314333) q[2];
sx q[2];
rz(-1.9224482) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0565666) q[1];
sx q[1];
rz(-1.4843656) q[1];
sx q[1];
rz(1.5734476) q[1];
rz(-pi) q[2];
rz(-1.7478862) q[3];
sx q[3];
rz(-1.0064126) q[3];
sx q[3];
rz(1.0799112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5829696) q[2];
sx q[2];
rz(-1.308459) q[2];
sx q[2];
rz(0.51897955) q[2];
rz(-0.44273043) q[3];
sx q[3];
rz(-1.3230007) q[3];
sx q[3];
rz(1.3435266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.9291572) q[0];
sx q[0];
rz(-0.99016187) q[0];
sx q[0];
rz(-1.1708175) q[0];
rz(0.15405542) q[1];
sx q[1];
rz(-2.2046721) q[1];
sx q[1];
rz(2.2278723) q[1];
rz(-0.51657233) q[2];
sx q[2];
rz(-1.5348829) q[2];
sx q[2];
rz(1.665967) q[2];
rz(0.89755015) q[3];
sx q[3];
rz(-1.297045) q[3];
sx q[3];
rz(0.53551077) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
