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
rz(2.763971) q[0];
sx q[0];
rz(-0.42855898) q[0];
sx q[0];
rz(-2.9052486) q[0];
rz(-1.6726681) q[1];
sx q[1];
rz(-1.5863215) q[1];
sx q[1];
rz(2.9798887) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16888389) q[0];
sx q[0];
rz(-2.5492269) q[0];
sx q[0];
rz(1.7667207) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9457372) q[2];
sx q[2];
rz(-1.5136949) q[2];
sx q[2];
rz(-1.7498992) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0031944) q[1];
sx q[1];
rz(-1.5439543) q[1];
sx q[1];
rz(-3.1346727) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3507648) q[3];
sx q[3];
rz(-1.9373978) q[3];
sx q[3];
rz(2.0442968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8990495) q[2];
sx q[2];
rz(-2.264302) q[2];
sx q[2];
rz(-0.38929942) q[2];
rz(0.23046514) q[3];
sx q[3];
rz(-3.1233628) q[3];
sx q[3];
rz(-2.5267595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56735754) q[0];
sx q[0];
rz(-2.1901972) q[0];
sx q[0];
rz(-1.6472598) q[0];
rz(-1.5556473) q[1];
sx q[1];
rz(-2.9289398) q[1];
sx q[1];
rz(-2.0071323) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8788505) q[0];
sx q[0];
rz(-1.8244184) q[0];
sx q[0];
rz(2.4154759) q[0];
rz(-pi) q[1];
rz(-0.73573839) q[2];
sx q[2];
rz(-1.2366857) q[2];
sx q[2];
rz(0.10057893) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8061773) q[1];
sx q[1];
rz(-1.9796438) q[1];
sx q[1];
rz(-0.93542904) q[1];
x q[2];
rz(2.3207619) q[3];
sx q[3];
rz(-1.9811673) q[3];
sx q[3];
rz(-0.73843187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.050934164) q[2];
sx q[2];
rz(-0.91973534) q[2];
sx q[2];
rz(1.8402137) q[2];
rz(-1.0743514) q[3];
sx q[3];
rz(-2.8315872) q[3];
sx q[3];
rz(-1.5955135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12304561) q[0];
sx q[0];
rz(-2.8392241) q[0];
sx q[0];
rz(-0.61755919) q[0];
rz(2.0410208) q[1];
sx q[1];
rz(-3.1220084) q[1];
sx q[1];
rz(-0.40357959) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36369187) q[0];
sx q[0];
rz(-3.0263623) q[0];
sx q[0];
rz(1.4794502) q[0];
rz(-pi) q[1];
rz(-0.17346548) q[2];
sx q[2];
rz(-2.6092165) q[2];
sx q[2];
rz(-2.8796822) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3364279) q[1];
sx q[1];
rz(-1.4138828) q[1];
sx q[1];
rz(3.0057231) q[1];
rz(-pi) q[2];
rz(-2.5837901) q[3];
sx q[3];
rz(-1.3645759) q[3];
sx q[3];
rz(2.5286412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6214211) q[2];
sx q[2];
rz(-1.6941864) q[2];
sx q[2];
rz(-2.3604895) q[2];
rz(-2.805294) q[3];
sx q[3];
rz(-1.4399485) q[3];
sx q[3];
rz(-2.4380016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4562562) q[0];
sx q[0];
rz(-2.6254613) q[0];
sx q[0];
rz(1.2180895) q[0];
rz(-1.3457899) q[1];
sx q[1];
rz(-3.1333874) q[1];
sx q[1];
rz(1.2340612) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4505554) q[0];
sx q[0];
rz(-0.19007401) q[0];
sx q[0];
rz(-3.0897753) q[0];
rz(2.7379964) q[2];
sx q[2];
rz(-2.4880313) q[2];
sx q[2];
rz(-0.60774481) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7750596) q[1];
sx q[1];
rz(-1.9294191) q[1];
sx q[1];
rz(2.8721362) q[1];
rz(2.5472121) q[3];
sx q[3];
rz(-3.1295589) q[3];
sx q[3];
rz(-1.8569225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9464843) q[2];
sx q[2];
rz(-0.4230963) q[2];
sx q[2];
rz(-1.1484324) q[2];
rz(0.75442433) q[3];
sx q[3];
rz(-1.1634588) q[3];
sx q[3];
rz(-2.8783126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36963439) q[0];
sx q[0];
rz(-0.088005528) q[0];
sx q[0];
rz(-0.60246402) q[0];
rz(0.98214904) q[1];
sx q[1];
rz(-3.1352477) q[1];
sx q[1];
rz(-2.6314661) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23778039) q[0];
sx q[0];
rz(-2.8717715) q[0];
sx q[0];
rz(0.80679195) q[0];
x q[1];
rz(2.192239) q[2];
sx q[2];
rz(-0.1489677) q[2];
sx q[2];
rz(-2.0003275) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.229631) q[1];
sx q[1];
rz(-1.0385286) q[1];
sx q[1];
rz(2.0144573) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3272133) q[3];
sx q[3];
rz(-1.571221) q[3];
sx q[3];
rz(-2.9546776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0977352) q[2];
sx q[2];
rz(-2.0243702) q[2];
sx q[2];
rz(-1.0350234) q[2];
rz(1.467661) q[3];
sx q[3];
rz(-0.36164713) q[3];
sx q[3];
rz(0.58290946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8408836) q[0];
sx q[0];
rz(-2.1568334) q[0];
sx q[0];
rz(-2.050052) q[0];
rz(-0.40065271) q[1];
sx q[1];
rz(-0.0023829208) q[1];
sx q[1];
rz(-2.5750652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3728722) q[0];
sx q[0];
rz(-0.75516381) q[0];
sx q[0];
rz(-1.7674957) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4479273) q[2];
sx q[2];
rz(-0.98204192) q[2];
sx q[2];
rz(-2.6679469) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2926082) q[1];
sx q[1];
rz(-2.1273534) q[1];
sx q[1];
rz(1.3054331) q[1];
rz(-pi) q[2];
rz(0.57265307) q[3];
sx q[3];
rz(-1.6653456) q[3];
sx q[3];
rz(1.8714286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1953676) q[2];
sx q[2];
rz(-2.3905498) q[2];
sx q[2];
rz(1.5283778) q[2];
rz(1.001312) q[3];
sx q[3];
rz(-0.87037218) q[3];
sx q[3];
rz(-0.24054578) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2840851) q[0];
sx q[0];
rz(-0.79428285) q[0];
sx q[0];
rz(1.1450144) q[0];
rz(1.5223632) q[1];
sx q[1];
rz(-3.122819) q[1];
sx q[1];
rz(-1.1633263) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6451241) q[0];
sx q[0];
rz(-1.860052) q[0];
sx q[0];
rz(3.0933063) q[0];
rz(2.2448213) q[2];
sx q[2];
rz(-1.4538897) q[2];
sx q[2];
rz(-0.74806556) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.39966291) q[1];
sx q[1];
rz(-2.1737639) q[1];
sx q[1];
rz(-0.16027995) q[1];
x q[2];
rz(-1.6895503) q[3];
sx q[3];
rz(-2.2014479) q[3];
sx q[3];
rz(-1.5753559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6254977) q[2];
sx q[2];
rz(-0.98573804) q[2];
sx q[2];
rz(2.7682448) q[2];
rz(-0.18800023) q[3];
sx q[3];
rz(-1.5865654) q[3];
sx q[3];
rz(1.1628304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2242551) q[0];
sx q[0];
rz(-2.6118216) q[0];
sx q[0];
rz(2.8944471) q[0];
rz(2.9180134) q[1];
sx q[1];
rz(-3.1378742) q[1];
sx q[1];
rz(1.5162969) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9947544) q[0];
sx q[0];
rz(-0.13299599) q[0];
sx q[0];
rz(2.07703) q[0];
rz(-0.76659492) q[2];
sx q[2];
rz(-1.8803036) q[2];
sx q[2];
rz(0.37878894) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.030236472) q[1];
sx q[1];
rz(-1.9647536) q[1];
sx q[1];
rz(-1.5535105) q[1];
rz(-pi) q[2];
rz(1.6187632) q[3];
sx q[3];
rz(-0.9771416) q[3];
sx q[3];
rz(-1.1831212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.197) q[2];
sx q[2];
rz(-2.4304515) q[2];
sx q[2];
rz(1.5468583) q[2];
rz(-0.77557766) q[3];
sx q[3];
rz(-1.8577134) q[3];
sx q[3];
rz(-1.6790793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.813886) q[0];
sx q[0];
rz(-2.1729108) q[0];
sx q[0];
rz(1.1073329) q[0];
rz(-2.2164717) q[1];
sx q[1];
rz(-0.0019625891) q[1];
sx q[1];
rz(0.75606871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30413142) q[0];
sx q[0];
rz(-2.5072949) q[0];
sx q[0];
rz(-2.1055806) q[0];
rz(-1.5599361) q[2];
sx q[2];
rz(-0.53087044) q[2];
sx q[2];
rz(0.52299352) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9097164) q[1];
sx q[1];
rz(-1.8617587) q[1];
sx q[1];
rz(-0.58334406) q[1];
x q[2];
rz(-2.5697034) q[3];
sx q[3];
rz(-1.4684621) q[3];
sx q[3];
rz(0.97692797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5209311) q[2];
sx q[2];
rz(-0.88520092) q[2];
sx q[2];
rz(-1.0010285) q[2];
rz(1.3641317) q[3];
sx q[3];
rz(-2.2083211) q[3];
sx q[3];
rz(0.53396839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-3.119568) q[0];
sx q[0];
rz(-1.3725932) q[0];
sx q[0];
rz(-0.46510988) q[0];
rz(1.8067092) q[1];
sx q[1];
rz(-2.7692134) q[1];
sx q[1];
rz(1.5660657) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6380999) q[0];
sx q[0];
rz(-1.3412807) q[0];
sx q[0];
rz(0.044661836) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31120531) q[2];
sx q[2];
rz(-2.0796144) q[2];
sx q[2];
rz(-2.6151163) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5425092) q[1];
sx q[1];
rz(-1.5712106) q[1];
sx q[1];
rz(-3.1387657) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6080721) q[3];
sx q[3];
rz(-1.9498943) q[3];
sx q[3];
rz(-2.9994734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.89432013) q[2];
sx q[2];
rz(-0.044581052) q[2];
sx q[2];
rz(1.0856005) q[2];
rz(1.9191437) q[3];
sx q[3];
rz(-2.6294851) q[3];
sx q[3];
rz(1.2781757) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6923675) q[0];
sx q[0];
rz(-1.7074371) q[0];
sx q[0];
rz(-1.3771124) q[0];
rz(1.5832681) q[1];
sx q[1];
rz(-2.227034) q[1];
sx q[1];
rz(-2.9169678) q[1];
rz(3.0947826) q[2];
sx q[2];
rz(-3.0407314) q[2];
sx q[2];
rz(-2.8032816) q[2];
rz(1.1261945) q[3];
sx q[3];
rz(-1.3570519) q[3];
sx q[3];
rz(-2.9828664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
