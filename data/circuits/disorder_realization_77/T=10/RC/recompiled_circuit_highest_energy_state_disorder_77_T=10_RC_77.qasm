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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18672189) q[0];
sx q[0];
rz(-1.6419517) q[0];
sx q[0];
rz(1.5118208) q[0];
x q[1];
rz(0.3977836) q[2];
sx q[2];
rz(-1.7353463) q[2];
sx q[2];
rz(-2.0497422) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2457468) q[1];
sx q[1];
rz(-1.7406643) q[1];
sx q[1];
rz(0.47791153) q[1];
x q[2];
rz(-2.5658251) q[3];
sx q[3];
rz(-1.9232436) q[3];
sx q[3];
rz(1.5653597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6930801) q[2];
sx q[2];
rz(-1.3658407) q[2];
sx q[2];
rz(2.8903294) q[2];
rz(-1.9692028) q[3];
sx q[3];
rz(-1.1029714) q[3];
sx q[3];
rz(0.049886543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4548816) q[0];
sx q[0];
rz(-2.659681) q[0];
sx q[0];
rz(1.7676109) q[0];
rz(0.33085597) q[1];
sx q[1];
rz(-1.7287946) q[1];
sx q[1];
rz(2.8673598) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73353798) q[0];
sx q[0];
rz(-1.499905) q[0];
sx q[0];
rz(-2.8246585) q[0];
rz(-0.53326989) q[2];
sx q[2];
rz(-1.5769713) q[2];
sx q[2];
rz(1.285163) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3365605) q[1];
sx q[1];
rz(-0.39440409) q[1];
sx q[1];
rz(-2.5488362) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9894473) q[3];
sx q[3];
rz(-1.6466221) q[3];
sx q[3];
rz(-2.0328498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8809044) q[2];
sx q[2];
rz(-1.3284677) q[2];
sx q[2];
rz(2.0534959) q[2];
rz(-2.0982096) q[3];
sx q[3];
rz(-2.9731396) q[3];
sx q[3];
rz(2.8823631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1328717) q[0];
sx q[0];
rz(-2.6368124) q[0];
sx q[0];
rz(1.1016499) q[0];
rz(0.77611008) q[1];
sx q[1];
rz(-1.848369) q[1];
sx q[1];
rz(-0.64220846) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6519847) q[0];
sx q[0];
rz(-0.63417681) q[0];
sx q[0];
rz(0.80192566) q[0];
rz(-pi) q[1];
x q[1];
rz(0.28046728) q[2];
sx q[2];
rz(-0.73958042) q[2];
sx q[2];
rz(0.47175874) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2930729) q[1];
sx q[1];
rz(-2.7024033) q[1];
sx q[1];
rz(1.8166914) q[1];
rz(-pi) q[2];
rz(0.74263708) q[3];
sx q[3];
rz(-2.3315695) q[3];
sx q[3];
rz(-2.1049316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3099826) q[2];
sx q[2];
rz(-1.9881366) q[2];
sx q[2];
rz(-2.2115121) q[2];
rz(0.81405226) q[3];
sx q[3];
rz(-1.9508773) q[3];
sx q[3];
rz(1.2306151) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26152465) q[0];
sx q[0];
rz(-1.9517169) q[0];
sx q[0];
rz(2.6370908) q[0];
rz(2.2241459) q[1];
sx q[1];
rz(-2.4023299) q[1];
sx q[1];
rz(2.9333072) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8802279) q[0];
sx q[0];
rz(-0.28631908) q[0];
sx q[0];
rz(2.5891735) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2913877) q[2];
sx q[2];
rz(-0.87051755) q[2];
sx q[2];
rz(2.6902187) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.69729561) q[1];
sx q[1];
rz(-1.0047067) q[1];
sx q[1];
rz(2.5725911) q[1];
rz(-pi) q[2];
rz(-1.3127906) q[3];
sx q[3];
rz(-1.0258706) q[3];
sx q[3];
rz(1.2598002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8248799) q[2];
sx q[2];
rz(-2.0471639) q[2];
sx q[2];
rz(-1.971395) q[2];
rz(-0.96585387) q[3];
sx q[3];
rz(-2.070277) q[3];
sx q[3];
rz(0.14537183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.589094) q[0];
sx q[0];
rz(-0.67411244) q[0];
sx q[0];
rz(-0.66459769) q[0];
rz(2.873114) q[1];
sx q[1];
rz(-1.6842027) q[1];
sx q[1];
rz(-2.1427515) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3420606) q[0];
sx q[0];
rz(-0.99208368) q[0];
sx q[0];
rz(-0.18216746) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8605609) q[2];
sx q[2];
rz(-2.701607) q[2];
sx q[2];
rz(-0.54674613) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4376862) q[1];
sx q[1];
rz(-1.4515299) q[1];
sx q[1];
rz(-1.7410623) q[1];
rz(-pi) q[2];
rz(-3.1206297) q[3];
sx q[3];
rz(-0.87910324) q[3];
sx q[3];
rz(-1.616459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9773679) q[2];
sx q[2];
rz(-2.2013142) q[2];
sx q[2];
rz(0.11717907) q[2];
rz(-3.0818648) q[3];
sx q[3];
rz(-1.9220587) q[3];
sx q[3];
rz(-0.81502771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77434671) q[0];
sx q[0];
rz(-0.4137488) q[0];
sx q[0];
rz(1.8814948) q[0];
rz(0.11481181) q[1];
sx q[1];
rz(-1.7962619) q[1];
sx q[1];
rz(-0.095349163) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60854429) q[0];
sx q[0];
rz(-1.5854821) q[0];
sx q[0];
rz(0.20242723) q[0];
x q[1];
rz(2.0748069) q[2];
sx q[2];
rz(-1.4246449) q[2];
sx q[2];
rz(-2.604217) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8417175) q[1];
sx q[1];
rz(-2.5073194) q[1];
sx q[1];
rz(2.3298765) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.010706832) q[3];
sx q[3];
rz(-1.043415) q[3];
sx q[3];
rz(2.7675865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.054333869) q[2];
sx q[2];
rz(-2.2237873) q[2];
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
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1374283) q[0];
sx q[0];
rz(-2.7935109) q[0];
sx q[0];
rz(0.75781703) q[0];
rz(0.025253145) q[1];
sx q[1];
rz(-0.39552894) q[1];
sx q[1];
rz(1.1753561) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7987418) q[0];
sx q[0];
rz(-1.9254058) q[0];
sx q[0];
rz(-2.6358428) q[0];
rz(-pi) q[1];
rz(-0.48904959) q[2];
sx q[2];
rz(-0.95121517) q[2];
sx q[2];
rz(-1.8085768) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3010294) q[1];
sx q[1];
rz(-1.1811273) q[1];
sx q[1];
rz(2.3430166) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5478915) q[3];
sx q[3];
rz(-1.1419019) q[3];
sx q[3];
rz(-2.1351086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.56388277) q[2];
sx q[2];
rz(-1.5241728) q[2];
sx q[2];
rz(1.9263402) q[2];
rz(0.75872129) q[3];
sx q[3];
rz(-1.4164378) q[3];
sx q[3];
rz(1.5846213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4828846) q[0];
sx q[0];
rz(-1.8526798) q[0];
sx q[0];
rz(-0.38920745) q[0];
rz(-3.0534741) q[1];
sx q[1];
rz(-1.4382818) q[1];
sx q[1];
rz(-1.6099991) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22613285) q[0];
sx q[0];
rz(-2.8185039) q[0];
sx q[0];
rz(-2.3996572) q[0];
rz(-pi) q[1];
rz(0.56014748) q[2];
sx q[2];
rz(-0.74290007) q[2];
sx q[2];
rz(-2.7882238) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.64634215) q[1];
sx q[1];
rz(-1.6227437) q[1];
sx q[1];
rz(-2.6990128) q[1];
rz(-1.1346899) q[3];
sx q[3];
rz(-0.75911555) q[3];
sx q[3];
rz(1.5576943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4072998) q[2];
sx q[2];
rz(-0.6137085) q[2];
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
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4765428) q[0];
sx q[0];
rz(-2.2344868) q[0];
sx q[0];
rz(-2.6764349) q[0];
rz(3.0910885) q[1];
sx q[1];
rz(-2.2513794) q[1];
sx q[1];
rz(2.6506298) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5652654) q[0];
sx q[0];
rz(-2.6748228) q[0];
sx q[0];
rz(-2.3719823) q[0];
x q[1];
rz(0.45456205) q[2];
sx q[2];
rz(-0.33966065) q[2];
sx q[2];
rz(0.31955921) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1909294) q[1];
sx q[1];
rz(-1.0196389) q[1];
sx q[1];
rz(-0.74857708) q[1];
rz(-pi) q[2];
rz(0.36763962) q[3];
sx q[3];
rz(-0.30045569) q[3];
sx q[3];
rz(-0.57747546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4831873) q[2];
sx q[2];
rz(-0.5961954) q[2];
sx q[2];
rz(-2.2779706) q[2];
rz(-2.9577799) q[3];
sx q[3];
rz(-2.176216) q[3];
sx q[3];
rz(-2.1553154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42613906) q[0];
sx q[0];
rz(-1.5197536) q[0];
sx q[0];
rz(0.62582985) q[0];
rz(0.65678701) q[1];
sx q[1];
rz(-0.96581179) q[1];
sx q[1];
rz(1.9331369) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0330808) q[0];
sx q[0];
rz(-1.4377366) q[0];
sx q[0];
rz(-0.10255524) q[0];
x q[1];
rz(0.65932806) q[2];
sx q[2];
rz(-1.8272597) q[2];
sx q[2];
rz(-1.2191445) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0543218) q[1];
sx q[1];
rz(-3.0551214) q[1];
sx q[1];
rz(-0.030589624) q[1];
x q[2];
rz(1.3937065) q[3];
sx q[3];
rz(-1.0064126) q[3];
sx q[3];
rz(1.0799112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5829696) q[2];
sx q[2];
rz(-1.308459) q[2];
sx q[2];
rz(0.51897955) q[2];
rz(-0.44273043) q[3];
sx q[3];
rz(-1.3230007) q[3];
sx q[3];
rz(-1.7980661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9291572) q[0];
sx q[0];
rz(-2.1514308) q[0];
sx q[0];
rz(1.9707752) q[0];
rz(-2.9875372) q[1];
sx q[1];
rz(-2.2046721) q[1];
sx q[1];
rz(2.2278723) q[1];
rz(-1.6120934) q[2];
sx q[2];
rz(-2.0870024) q[2];
sx q[2];
rz(0.11556297) q[2];
rz(-2.2440425) q[3];
sx q[3];
rz(-1.297045) q[3];
sx q[3];
rz(0.53551077) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
