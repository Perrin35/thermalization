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
rz(2.0410886) q[0];
sx q[0];
rz(2.4325844) q[0];
sx q[0];
rz(11.772052) q[0];
rz(-2.2494443) q[1];
sx q[1];
rz(-1.4798857) q[1];
sx q[1];
rz(-2.7977112) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00044879221) q[0];
sx q[0];
rz(-1.4836932) q[0];
sx q[0];
rz(-0.25229519) q[0];
rz(-0.8240118) q[2];
sx q[2];
rz(-0.67979797) q[2];
sx q[2];
rz(2.4275017) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0176489) q[1];
sx q[1];
rz(-2.9363465) q[1];
sx q[1];
rz(-0.68428792) q[1];
rz(-0.12012336) q[3];
sx q[3];
rz(-1.7645482) q[3];
sx q[3];
rz(-1.4258949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6227365) q[2];
sx q[2];
rz(-1.3750261) q[2];
sx q[2];
rz(-2.1905017) q[2];
rz(-0.91894346) q[3];
sx q[3];
rz(-0.94075847) q[3];
sx q[3];
rz(2.9520891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6491991) q[0];
sx q[0];
rz(-0.7890107) q[0];
sx q[0];
rz(1.5845818) q[0];
rz(2.4015265) q[1];
sx q[1];
rz(-2.2941755) q[1];
sx q[1];
rz(1.3620296) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1213367) q[0];
sx q[0];
rz(-2.0411939) q[0];
sx q[0];
rz(-1.008995) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9744423) q[2];
sx q[2];
rz(-1.6975308) q[2];
sx q[2];
rz(2.6069146) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6891514) q[1];
sx q[1];
rz(-1.2106306) q[1];
sx q[1];
rz(-2.4634201) q[1];
rz(-pi) q[2];
x q[2];
rz(0.091097761) q[3];
sx q[3];
rz(-2.280844) q[3];
sx q[3];
rz(1.2595081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7035199) q[2];
sx q[2];
rz(-1.9813931) q[2];
sx q[2];
rz(-2.4959219) q[2];
rz(-0.034363184) q[3];
sx q[3];
rz(-0.86302257) q[3];
sx q[3];
rz(3.0784975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7647758) q[0];
sx q[0];
rz(-0.21012981) q[0];
sx q[0];
rz(2.3576376) q[0];
rz(-2.8249373) q[1];
sx q[1];
rz(-1.6461992) q[1];
sx q[1];
rz(2.1276316) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72958845) q[0];
sx q[0];
rz(-0.17590287) q[0];
sx q[0];
rz(0.094502016) q[0];
x q[1];
rz(-1.1523109) q[2];
sx q[2];
rz(-2.6259702) q[2];
sx q[2];
rz(1.9299233) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5875776) q[1];
sx q[1];
rz(-1.3237751) q[1];
sx q[1];
rz(2.8783911) q[1];
rz(-0.52041913) q[3];
sx q[3];
rz(-1.0380967) q[3];
sx q[3];
rz(-1.0315438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.97518146) q[2];
sx q[2];
rz(-1.760249) q[2];
sx q[2];
rz(0.4057917) q[2];
rz(-0.25519145) q[3];
sx q[3];
rz(-1.8813335) q[3];
sx q[3];
rz(0.7757799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-2.880421) q[0];
sx q[0];
rz(-0.71735993) q[0];
sx q[0];
rz(-0.98947155) q[0];
rz(-0.62670341) q[1];
sx q[1];
rz(-0.53755212) q[1];
sx q[1];
rz(0.42868838) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1470448) q[0];
sx q[0];
rz(-2.328985) q[0];
sx q[0];
rz(-1.7269686) q[0];
rz(-pi) q[1];
rz(-1.4312885) q[2];
sx q[2];
rz(-1.2409455) q[2];
sx q[2];
rz(-2.863186) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3544412) q[1];
sx q[1];
rz(-1.0686744) q[1];
sx q[1];
rz(0.67014613) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4267202) q[3];
sx q[3];
rz(-1.360641) q[3];
sx q[3];
rz(1.6870267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.89877146) q[2];
sx q[2];
rz(-3.0148744) q[2];
sx q[2];
rz(0.80647331) q[2];
rz(1.1558695) q[3];
sx q[3];
rz(-0.75771302) q[3];
sx q[3];
rz(2.9466467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72215801) q[0];
sx q[0];
rz(-2.2787155) q[0];
sx q[0];
rz(-0.3062329) q[0];
rz(-0.58079863) q[1];
sx q[1];
rz(-0.88169801) q[1];
sx q[1];
rz(-2.4961848) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78272143) q[0];
sx q[0];
rz(-2.5497265) q[0];
sx q[0];
rz(-2.84821) q[0];
x q[1];
rz(-3.0642411) q[2];
sx q[2];
rz(-2.0708551) q[2];
sx q[2];
rz(-1.7986705) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.66845208) q[1];
sx q[1];
rz(-0.73669556) q[1];
sx q[1];
rz(2.2948242) q[1];
x q[2];
rz(1.3401806) q[3];
sx q[3];
rz(-1.5531887) q[3];
sx q[3];
rz(1.0215774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.2885651) q[2];
sx q[2];
rz(-1.7662851) q[2];
sx q[2];
rz(-2.1993401) q[2];
rz(-2.6906768) q[3];
sx q[3];
rz(-0.98603606) q[3];
sx q[3];
rz(0.26421419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0838202) q[0];
sx q[0];
rz(-1.8696996) q[0];
sx q[0];
rz(0.94495946) q[0];
rz(-0.68929976) q[1];
sx q[1];
rz(-2.0336475) q[1];
sx q[1];
rz(2.0384929) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8980628) q[0];
sx q[0];
rz(-1.8193428) q[0];
sx q[0];
rz(2.4175279) q[0];
rz(1.9951773) q[2];
sx q[2];
rz(-1.7327899) q[2];
sx q[2];
rz(-1.4206518) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0755495) q[1];
sx q[1];
rz(-2.4166921) q[1];
sx q[1];
rz(-2.3182475) q[1];
rz(0.50838105) q[3];
sx q[3];
rz(-1.0219921) q[3];
sx q[3];
rz(-0.8027975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5276706) q[2];
sx q[2];
rz(-0.066378243) q[2];
sx q[2];
rz(1.684368) q[2];
rz(-0.97113222) q[3];
sx q[3];
rz(-1.5338273) q[3];
sx q[3];
rz(-3.0253313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3986452) q[0];
sx q[0];
rz(-1.8859517) q[0];
sx q[0];
rz(-0.1980814) q[0];
rz(-0.48175851) q[1];
sx q[1];
rz(-1.2629197) q[1];
sx q[1];
rz(-1.6162965) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85563816) q[0];
sx q[0];
rz(-1.5795991) q[0];
sx q[0];
rz(1.5670526) q[0];
rz(-pi) q[1];
rz(3.1379596) q[2];
sx q[2];
rz(-1.1143408) q[2];
sx q[2];
rz(-2.9287287) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9349667) q[1];
sx q[1];
rz(-0.68892067) q[1];
sx q[1];
rz(0.72241131) q[1];
rz(-pi) q[2];
rz(-2.0644998) q[3];
sx q[3];
rz(-1.6564661) q[3];
sx q[3];
rz(2.6926797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5109479) q[2];
sx q[2];
rz(-1.7787361) q[2];
sx q[2];
rz(-2.2044568) q[2];
rz(-2.7142094) q[3];
sx q[3];
rz(-2.1338227) q[3];
sx q[3];
rz(-1.5510604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9464924) q[0];
sx q[0];
rz(-1.6747403) q[0];
sx q[0];
rz(1.7505769) q[0];
rz(-1.6777978) q[1];
sx q[1];
rz(-0.59278667) q[1];
sx q[1];
rz(-0.52267271) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14870384) q[0];
sx q[0];
rz(-0.33451053) q[0];
sx q[0];
rz(-1.3687031) q[0];
x q[1];
rz(0.14334579) q[2];
sx q[2];
rz(-2.0269577) q[2];
sx q[2];
rz(-0.83269115) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1795038) q[1];
sx q[1];
rz(-0.25280646) q[1];
sx q[1];
rz(-2.4180555) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0955047) q[3];
sx q[3];
rz(-1.3541835) q[3];
sx q[3];
rz(-3.0334186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3576144) q[2];
sx q[2];
rz(-1.7366333) q[2];
sx q[2];
rz(0.092255175) q[2];
rz(2.9914894) q[3];
sx q[3];
rz(-1.152252) q[3];
sx q[3];
rz(-2.2805555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0263696) q[0];
sx q[0];
rz(-2.8404591) q[0];
sx q[0];
rz(1.2979771) q[0];
rz(0.84833974) q[1];
sx q[1];
rz(-1.4804761) q[1];
sx q[1];
rz(-0.27378219) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0299579) q[0];
sx q[0];
rz(-0.23914251) q[0];
sx q[0];
rz(2.9947623) q[0];
x q[1];
rz(1.8452497) q[2];
sx q[2];
rz(-2.2355766) q[2];
sx q[2];
rz(2.279778) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50934511) q[1];
sx q[1];
rz(-1.6654105) q[1];
sx q[1];
rz(2.3791887) q[1];
x q[2];
rz(1.3630834) q[3];
sx q[3];
rz(-1.9915607) q[3];
sx q[3];
rz(1.8150235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0691284) q[2];
sx q[2];
rz(-1.0651361) q[2];
sx q[2];
rz(-2.1799942) q[2];
rz(-1.0570863) q[3];
sx q[3];
rz(-1.0870533) q[3];
sx q[3];
rz(0.45007733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7708275) q[0];
sx q[0];
rz(-1.2130883) q[0];
sx q[0];
rz(0.83592498) q[0];
rz(-1.4113374) q[1];
sx q[1];
rz(-2.5524499) q[1];
sx q[1];
rz(-0.020847598) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1040346) q[0];
sx q[0];
rz(-1.7171894) q[0];
sx q[0];
rz(0.030265042) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7078988) q[2];
sx q[2];
rz(-1.1011656) q[2];
sx q[2];
rz(3.0933703) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.97562645) q[1];
sx q[1];
rz(-1.8667815) q[1];
sx q[1];
rz(1.8459666) q[1];
x q[2];
rz(-1.3278058) q[3];
sx q[3];
rz(-1.3919481) q[3];
sx q[3];
rz(-2.230165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6622322) q[2];
sx q[2];
rz(-1.6785097) q[2];
sx q[2];
rz(1.642891) q[2];
rz(-2.6920476) q[3];
sx q[3];
rz(-0.42732626) q[3];
sx q[3];
rz(2.8699365) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0679006) q[0];
sx q[0];
rz(-1.3561159) q[0];
sx q[0];
rz(2.7664716) q[0];
rz(-0.66315229) q[1];
sx q[1];
rz(-1.2620435) q[1];
sx q[1];
rz(-0.97558998) q[1];
rz(-1.763487) q[2];
sx q[2];
rz(-1.4507254) q[2];
sx q[2];
rz(-2.7117827) q[2];
rz(-0.4573902) q[3];
sx q[3];
rz(-1.8531696) q[3];
sx q[3];
rz(2.8715092) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
