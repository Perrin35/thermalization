OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.99958217) q[0];
sx q[0];
rz(-1.002123) q[0];
sx q[0];
rz(2.2440417) q[0];
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(-1.0645359) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25047725) q[0];
sx q[0];
rz(-2.5076137) q[0];
sx q[0];
rz(-1.4527713) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35383309) q[2];
sx q[2];
rz(-0.96828038) q[2];
sx q[2];
rz(1.2474071) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0734288) q[1];
sx q[1];
rz(-0.66232077) q[1];
sx q[1];
rz(-2.8627002) q[1];
rz(-1.3863871) q[3];
sx q[3];
rz(-0.82114906) q[3];
sx q[3];
rz(0.26339312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.66427461) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(2.409639) q[2];
rz(-0.96015635) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(-1.7094973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17523781) q[0];
sx q[0];
rz(-2.2580999) q[0];
sx q[0];
rz(-0.91645855) q[0];
rz(-2.6610999) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(-0.8786456) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7322313) q[0];
sx q[0];
rz(-1.1890113) q[0];
sx q[0];
rz(2.2344927) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36388134) q[2];
sx q[2];
rz(-2.4907787) q[2];
sx q[2];
rz(-1.770307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2616927) q[1];
sx q[1];
rz(-0.42527929) q[1];
sx q[1];
rz(0.17875032) q[1];
x q[2];
rz(-1.3602123) q[3];
sx q[3];
rz(-2.8254291) q[3];
sx q[3];
rz(-2.8663243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8615222) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(-2.4943165) q[2];
rz(0.17368008) q[3];
sx q[3];
rz(-1.1207542) q[3];
sx q[3];
rz(2.98996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52755255) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(2.8955984) q[0];
rz(1.7315158) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(-2.0203967) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66378731) q[0];
sx q[0];
rz(-2.2664547) q[0];
sx q[0];
rz(2.0545161) q[0];
x q[1];
rz(-1.1779184) q[2];
sx q[2];
rz(-1.6298721) q[2];
sx q[2];
rz(-1.3438366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.4157297) q[1];
sx q[1];
rz(-1.9007069) q[1];
sx q[1];
rz(1.8833453) q[1];
rz(-pi) q[2];
rz(-2.9259053) q[3];
sx q[3];
rz(-0.50369278) q[3];
sx q[3];
rz(1.4195201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7401509) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(-2.1901954) q[2];
rz(2.4915063) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(0.29561177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51405108) q[0];
sx q[0];
rz(-2.5294332) q[0];
sx q[0];
rz(-2.3535368) q[0];
rz(-1.0568985) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(-3.025211) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6560293) q[0];
sx q[0];
rz(-2.2040327) q[0];
sx q[0];
rz(-1.5390736) q[0];
rz(0.94159796) q[2];
sx q[2];
rz(-0.59855748) q[2];
sx q[2];
rz(1.7184005) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1981922) q[1];
sx q[1];
rz(-2.0640089) q[1];
sx q[1];
rz(1.1073768) q[1];
rz(-pi) q[2];
rz(-1.5192401) q[3];
sx q[3];
rz(-2.0600852) q[3];
sx q[3];
rz(-1.6980905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6115761) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(-2.8895565) q[2];
rz(-2.7633372) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(-0.3616412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5761121) q[0];
sx q[0];
rz(-1.7385372) q[0];
sx q[0];
rz(0.53260032) q[0];
rz(-1.4415007) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(1.1486357) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2411597) q[0];
sx q[0];
rz(-2.1923089) q[0];
sx q[0];
rz(0.09089367) q[0];
x q[1];
rz(-0.89678905) q[2];
sx q[2];
rz(-1.2090346) q[2];
sx q[2];
rz(2.0828431) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83541218) q[1];
sx q[1];
rz(-0.40459834) q[1];
sx q[1];
rz(-2.6447891) q[1];
x q[2];
rz(-2.7025181) q[3];
sx q[3];
rz(-2.3665161) q[3];
sx q[3];
rz(0.72417688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0107515) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(-1.032069) q[2];
rz(-2.4268835) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(-3.1177974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59153581) q[0];
sx q[0];
rz(-1.2055826) q[0];
sx q[0];
rz(-2.7456039) q[0];
rz(1.4453325) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(2.9352303) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2287184) q[0];
sx q[0];
rz(-1.4836856) q[0];
sx q[0];
rz(-0.74761439) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4762819) q[2];
sx q[2];
rz(-2.1449001) q[2];
sx q[2];
rz(-2.2523508) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5620835) q[1];
sx q[1];
rz(-1.4663896) q[1];
sx q[1];
rz(-0.097192055) q[1];
x q[2];
rz(2.9456375) q[3];
sx q[3];
rz(-1.2489508) q[3];
sx q[3];
rz(-0.67924196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.53987327) q[2];
sx q[2];
rz(-1.0605992) q[2];
sx q[2];
rz(0.27077857) q[2];
rz(0.74832908) q[3];
sx q[3];
rz(-1.7691408) q[3];
sx q[3];
rz(-1.07871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2449743) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(0.57624972) q[0];
rz(-2.9684864) q[1];
sx q[1];
rz(-0.74179596) q[1];
sx q[1];
rz(1.9304088) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.422238) q[0];
sx q[0];
rz(-2.1033759) q[0];
sx q[0];
rz(-2.5371735) q[0];
rz(0.48970512) q[2];
sx q[2];
rz(-1.4486794) q[2];
sx q[2];
rz(-2.7011938) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3019575) q[1];
sx q[1];
rz(-0.17050276) q[1];
sx q[1];
rz(1.9819928) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.022642) q[3];
sx q[3];
rz(-0.65330905) q[3];
sx q[3];
rz(-2.2484231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8048191) q[2];
sx q[2];
rz(-2.4833198) q[2];
sx q[2];
rz(-3.1170735) q[2];
rz(-2.426614) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(1.582675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0054758469) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(-2.1210282) q[0];
rz(0.15469805) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(-1.221009) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85212612) q[0];
sx q[0];
rz(-2.3131436) q[0];
sx q[0];
rz(-2.7702987) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4079354) q[2];
sx q[2];
rz(-1.7454299) q[2];
sx q[2];
rz(2.3101431) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.358566) q[1];
sx q[1];
rz(-1.4604124) q[1];
sx q[1];
rz(-3.0049938) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9302759) q[3];
sx q[3];
rz(-0.7822789) q[3];
sx q[3];
rz(2.4855011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8528379) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(2.4364046) q[2];
rz(-2.5382036) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(-1.6220629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.110638) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(0.13701339) q[0];
rz(0.6048454) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(0.12577122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3868956) q[0];
sx q[0];
rz(-1.3820952) q[0];
sx q[0];
rz(-2.9625091) q[0];
x q[1];
rz(2.8081886) q[2];
sx q[2];
rz(-1.7575022) q[2];
sx q[2];
rz(-0.29078996) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8217433) q[1];
sx q[1];
rz(-2.2935767) q[1];
sx q[1];
rz(-0.38139947) q[1];
rz(0.59154193) q[3];
sx q[3];
rz(-2.7898443) q[3];
sx q[3];
rz(-1.8464309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.13828364) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(2.4323145) q[2];
rz(-0.62018958) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(-2.9848849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9897292) q[0];
sx q[0];
rz(-0.79280889) q[0];
sx q[0];
rz(-1.2003157) q[0];
rz(2.8740846) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(2.0013924) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8126497) q[0];
sx q[0];
rz(-2.86033) q[0];
sx q[0];
rz(-1.2055231) q[0];
rz(-2.033038) q[2];
sx q[2];
rz(-1.8046364) q[2];
sx q[2];
rz(0.086694593) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2503567) q[1];
sx q[1];
rz(-1.48238) q[1];
sx q[1];
rz(2.9524515) q[1];
rz(-0.31302281) q[3];
sx q[3];
rz(-2.3232197) q[3];
sx q[3];
rz(-2.5767874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.37529477) q[2];
sx q[2];
rz(-1.4582429) q[2];
sx q[2];
rz(-0.89861384) q[2];
rz(3.1344154) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(-1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89467775) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(-2.6208411) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(2.649879) q[2];
sx q[2];
rz(-1.1469054) q[2];
sx q[2];
rz(2.7804874) q[2];
rz(-2.5557774) q[3];
sx q[3];
rz(-1.2231493) q[3];
sx q[3];
rz(2.3412658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
