OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.37031072) q[0];
sx q[0];
rz(-2.1455278) q[0];
sx q[0];
rz(-2.2709742) q[0];
rz(2.119996) q[1];
sx q[1];
rz(-2.8586913) q[1];
sx q[1];
rz(0.14970782) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4275442) q[0];
sx q[0];
rz(-2.2872426) q[0];
sx q[0];
rz(-0.9057522) q[0];
x q[1];
rz(-2.3157273) q[2];
sx q[2];
rz(-0.68544938) q[2];
sx q[2];
rz(0.65537383) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.522361) q[1];
sx q[1];
rz(-1.751096) q[1];
sx q[1];
rz(-3.1027604) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9505694) q[3];
sx q[3];
rz(-2.1130307) q[3];
sx q[3];
rz(2.1477826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8831138) q[2];
sx q[2];
rz(-0.09720619) q[2];
sx q[2];
rz(2.4374938) q[2];
rz(-2.1885833) q[3];
sx q[3];
rz(-2.1694031) q[3];
sx q[3];
rz(-1.4037508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54884058) q[0];
sx q[0];
rz(-1.5177746) q[0];
sx q[0];
rz(-2.5090704) q[0];
rz(0.44644341) q[1];
sx q[1];
rz(-1.7233142) q[1];
sx q[1];
rz(2.4893563) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95603847) q[0];
sx q[0];
rz(-1.599405) q[0];
sx q[0];
rz(-1.5584598) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9035283) q[2];
sx q[2];
rz(-1.2859584) q[2];
sx q[2];
rz(-0.1711947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.10416874) q[1];
sx q[1];
rz(-1.2836604) q[1];
sx q[1];
rz(-1.3377405) q[1];
rz(-pi) q[2];
rz(-0.84516256) q[3];
sx q[3];
rz(-2.3361428) q[3];
sx q[3];
rz(1.0172539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5941045) q[2];
sx q[2];
rz(-1.0141806) q[2];
sx q[2];
rz(-1.1616421) q[2];
rz(-1.9836327) q[3];
sx q[3];
rz(-2.0722814) q[3];
sx q[3];
rz(-1.4512216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050425477) q[0];
sx q[0];
rz(-0.37910351) q[0];
sx q[0];
rz(0.85025775) q[0];
rz(-2.6440874) q[1];
sx q[1];
rz(-1.18327) q[1];
sx q[1];
rz(-1.7920378) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2230167) q[0];
sx q[0];
rz(-2.4239459) q[0];
sx q[0];
rz(-2.7735604) q[0];
rz(-pi) q[1];
rz(-1.5047968) q[2];
sx q[2];
rz(-1.5335576) q[2];
sx q[2];
rz(0.50022349) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2022484) q[1];
sx q[1];
rz(-1.0679686) q[1];
sx q[1];
rz(1.8400251) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9968824) q[3];
sx q[3];
rz(-1.5618556) q[3];
sx q[3];
rz(2.5090891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1094018) q[2];
sx q[2];
rz(-1.2357864) q[2];
sx q[2];
rz(-0.90399495) q[2];
rz(-0.30113014) q[3];
sx q[3];
rz(-1.3826933) q[3];
sx q[3];
rz(1.7416471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
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
rz(2.6501453) q[0];
sx q[0];
rz(-2.1701145) q[0];
sx q[0];
rz(1.4105463) q[0];
rz(2.5097805) q[1];
sx q[1];
rz(-1.8099064) q[1];
sx q[1];
rz(3.1052123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43544337) q[0];
sx q[0];
rz(-2.2009146) q[0];
sx q[0];
rz(0.26421996) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6668947) q[2];
sx q[2];
rz(-2.5050852) q[2];
sx q[2];
rz(-2.0091025) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42602793) q[1];
sx q[1];
rz(-1.4392816) q[1];
sx q[1];
rz(-0.87042602) q[1];
rz(-pi) q[2];
rz(-0.23457228) q[3];
sx q[3];
rz(-1.6664701) q[3];
sx q[3];
rz(-2.5535339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.92695421) q[2];
sx q[2];
rz(-1.7094694) q[2];
sx q[2];
rz(1.1882163) q[2];
rz(-0.67048091) q[3];
sx q[3];
rz(-1.9159578) q[3];
sx q[3];
rz(-0.59613434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4398414) q[0];
sx q[0];
rz(-1.8816467) q[0];
sx q[0];
rz(0.16648509) q[0];
rz(-2.3855551) q[1];
sx q[1];
rz(-2.2557204) q[1];
sx q[1];
rz(-2.9072445) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4317961) q[0];
sx q[0];
rz(-2.095982) q[0];
sx q[0];
rz(-1.964142) q[0];
rz(1.0767897) q[2];
sx q[2];
rz(-1.4865808) q[2];
sx q[2];
rz(-2.6780724) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.25917945) q[1];
sx q[1];
rz(-2.3098574) q[1];
sx q[1];
rz(1.773136) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.078951051) q[3];
sx q[3];
rz(-1.2478932) q[3];
sx q[3];
rz(0.80432804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4328737) q[2];
sx q[2];
rz(-1.9659698) q[2];
sx q[2];
rz(-2.9210572) q[2];
rz(2.7045414) q[3];
sx q[3];
rz(-2.1190937) q[3];
sx q[3];
rz(-2.3760858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41546145) q[0];
sx q[0];
rz(-0.3188062) q[0];
sx q[0];
rz(0.81714001) q[0];
rz(-0.56610402) q[1];
sx q[1];
rz(-1.348446) q[1];
sx q[1];
rz(1.1436499) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9691539) q[0];
sx q[0];
rz(-1.5307431) q[0];
sx q[0];
rz(0.37102951) q[0];
x q[1];
rz(2.8191889) q[2];
sx q[2];
rz(-0.69903261) q[2];
sx q[2];
rz(0.93271819) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3559349) q[1];
sx q[1];
rz(-2.6895084) q[1];
sx q[1];
rz(1.7788586) q[1];
rz(-0.51795824) q[3];
sx q[3];
rz(-1.8149788) q[3];
sx q[3];
rz(1.1740008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3879261) q[2];
sx q[2];
rz(-1.5796698) q[2];
sx q[2];
rz(-0.4894408) q[2];
rz(-2.9135381) q[3];
sx q[3];
rz(-1.2585879) q[3];
sx q[3];
rz(-0.42603809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5722826) q[0];
sx q[0];
rz(-0.64240488) q[0];
sx q[0];
rz(-1.2868767) q[0];
rz(-2.4781748) q[1];
sx q[1];
rz(-1.5723012) q[1];
sx q[1];
rz(-1.9082665) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39499261) q[0];
sx q[0];
rz(-1.1867503) q[0];
sx q[0];
rz(-0.011944255) q[0];
rz(-pi) q[1];
rz(2.5160518) q[2];
sx q[2];
rz(-2.1893246) q[2];
sx q[2];
rz(-0.90460888) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.26112939) q[1];
sx q[1];
rz(-1.556734) q[1];
sx q[1];
rz(3.0304099) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.535378) q[3];
sx q[3];
rz(-2.5698235) q[3];
sx q[3];
rz(-1.9619463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.104091) q[2];
sx q[2];
rz(-2.9064894) q[2];
sx q[2];
rz(2.1255778) q[2];
rz(3.0715023) q[3];
sx q[3];
rz(-1.9349808) q[3];
sx q[3];
rz(-2.0751374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12524097) q[0];
sx q[0];
rz(-0.68269435) q[0];
sx q[0];
rz(1.6960779) q[0];
rz(0.21487543) q[1];
sx q[1];
rz(-2.3863249) q[1];
sx q[1];
rz(1.258237) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1930775) q[0];
sx q[0];
rz(-2.6876039) q[0];
sx q[0];
rz(-1.8817188) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7779185) q[2];
sx q[2];
rz(-2.3292543) q[2];
sx q[2];
rz(-1.9922436) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.59626034) q[1];
sx q[1];
rz(-0.55570554) q[1];
sx q[1];
rz(-2.943379) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6120841) q[3];
sx q[3];
rz(-0.92748517) q[3];
sx q[3];
rz(0.95638004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4593279) q[2];
sx q[2];
rz(-2.9569646) q[2];
sx q[2];
rz(-2.5788467) q[2];
rz(2.941926) q[3];
sx q[3];
rz(-2.4797347) q[3];
sx q[3];
rz(2.0195885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11583081) q[0];
sx q[0];
rz(-1.04302) q[0];
sx q[0];
rz(-0.2510221) q[0];
rz(-0.42731467) q[1];
sx q[1];
rz(-1.9117833) q[1];
sx q[1];
rz(-3.1138611) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8160307) q[0];
sx q[0];
rz(-1.0849909) q[0];
sx q[0];
rz(-2.3150139) q[0];
rz(-pi) q[1];
rz(-0.8717732) q[2];
sx q[2];
rz(-1.36424) q[2];
sx q[2];
rz(2.2733462) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8393644) q[1];
sx q[1];
rz(-1.4711958) q[1];
sx q[1];
rz(0.60217963) q[1];
rz(1.5042217) q[3];
sx q[3];
rz(-1.596631) q[3];
sx q[3];
rz(-1.0915826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1256844) q[2];
sx q[2];
rz(-2.1618844) q[2];
sx q[2];
rz(1.998418) q[2];
rz(0.14287359) q[3];
sx q[3];
rz(-1.6294799) q[3];
sx q[3];
rz(-2.2843602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7509572) q[0];
sx q[0];
rz(-1.2049144) q[0];
sx q[0];
rz(0.45387682) q[0];
rz(-0.67165309) q[1];
sx q[1];
rz(-1.4524873) q[1];
sx q[1];
rz(-2.8840816) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2902381) q[0];
sx q[0];
rz(-0.60831735) q[0];
sx q[0];
rz(-2.4277707) q[0];
rz(2.5970039) q[2];
sx q[2];
rz(-2.0047744) q[2];
sx q[2];
rz(-1.7314272) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.232302) q[1];
sx q[1];
rz(-1.5564939) q[1];
sx q[1];
rz(-0.69625744) q[1];
rz(-2.0503644) q[3];
sx q[3];
rz(-1.7864831) q[3];
sx q[3];
rz(1.8857764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8132849) q[2];
sx q[2];
rz(-1.6878781) q[2];
sx q[2];
rz(2.5349687) q[2];
rz(-2.666752) q[3];
sx q[3];
rz(-2.1947221) q[3];
sx q[3];
rz(-0.84038466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3474779) q[0];
sx q[0];
rz(-1.5119727) q[0];
sx q[0];
rz(2.150362) q[0];
rz(2.9150302) q[1];
sx q[1];
rz(-1.7270052) q[1];
sx q[1];
rz(-2.5546767) q[1];
rz(1.1889585) q[2];
sx q[2];
rz(-1.1321862) q[2];
sx q[2];
rz(1.95375) q[2];
rz(1.0971309) q[3];
sx q[3];
rz(-2.6087425) q[3];
sx q[3];
rz(-2.9125924) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
