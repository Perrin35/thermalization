OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3616537) q[0];
sx q[0];
rz(-3.0335463) q[0];
sx q[0];
rz(-1.894423) q[0];
rz(-0.78020686) q[1];
sx q[1];
rz(8.297774) q[1];
sx q[1];
rz(10.170593) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4853349) q[0];
sx q[0];
rz(-1.7575329) q[0];
sx q[0];
rz(1.223808) q[0];
rz(-3.132944) q[2];
sx q[2];
rz(-1.4100572) q[2];
sx q[2];
rz(2.6313105) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.76872003) q[1];
sx q[1];
rz(-0.24133397) q[1];
sx q[1];
rz(2.0745501) q[1];
x q[2];
rz(-0.28663825) q[3];
sx q[3];
rz(-0.53813808) q[3];
sx q[3];
rz(-1.9236444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1625533) q[2];
sx q[2];
rz(-2.5472842) q[2];
sx q[2];
rz(-1.7621367) q[2];
rz(2.4123794) q[3];
sx q[3];
rz(-2.3766434) q[3];
sx q[3];
rz(-0.929207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5586229) q[0];
sx q[0];
rz(-2.0780777) q[0];
sx q[0];
rz(-2.9003918) q[0];
rz(2.8855715) q[1];
sx q[1];
rz(-2.119901) q[1];
sx q[1];
rz(2.3604438) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1597693) q[0];
sx q[0];
rz(-2.0836813) q[0];
sx q[0];
rz(2.485093) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0623661) q[2];
sx q[2];
rz(-0.56622711) q[2];
sx q[2];
rz(1.5528284) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.95805146) q[1];
sx q[1];
rz(-2.1528917) q[1];
sx q[1];
rz(0.84611012) q[1];
x q[2];
rz(3.097104) q[3];
sx q[3];
rz(-1.4334363) q[3];
sx q[3];
rz(0.024976211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8476734) q[2];
sx q[2];
rz(-1.5896229) q[2];
sx q[2];
rz(-1.5839362) q[2];
rz(-3.1368351) q[3];
sx q[3];
rz(-2.5465951) q[3];
sx q[3];
rz(0.34576542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7541589) q[0];
sx q[0];
rz(-1.4492946) q[0];
sx q[0];
rz(0.19101983) q[0];
rz(2.7905131) q[1];
sx q[1];
rz(-2.7103238) q[1];
sx q[1];
rz(-1.4494337) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3397749) q[0];
sx q[0];
rz(-1.9097985) q[0];
sx q[0];
rz(1.5168241) q[0];
x q[1];
rz(-1.1507785) q[2];
sx q[2];
rz(-1.1598829) q[2];
sx q[2];
rz(1.0622298) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8061132) q[1];
sx q[1];
rz(-1.8361143) q[1];
sx q[1];
rz(2.669304) q[1];
rz(2.6552917) q[3];
sx q[3];
rz(-2.4913906) q[3];
sx q[3];
rz(0.40806684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.2991952) q[2];
sx q[2];
rz(-1.4292052) q[2];
sx q[2];
rz(3.0136285) q[2];
rz(-1.1658824) q[3];
sx q[3];
rz(-2.2539625) q[3];
sx q[3];
rz(0.33795801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2853476) q[0];
sx q[0];
rz(-1.0839533) q[0];
sx q[0];
rz(1.6612843) q[0];
rz(0.081143204) q[1];
sx q[1];
rz(-1.0973884) q[1];
sx q[1];
rz(0.88044423) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9362088) q[0];
sx q[0];
rz(-3.1160389) q[0];
sx q[0];
rz(-0.25746246) q[0];
rz(-pi) q[1];
x q[1];
rz(1.501785) q[2];
sx q[2];
rz(-0.61997783) q[2];
sx q[2];
rz(-0.91199707) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.26759155) q[1];
sx q[1];
rz(-2.1028215) q[1];
sx q[1];
rz(1.6997972) q[1];
x q[2];
rz(2.5838636) q[3];
sx q[3];
rz(-1.5909877) q[3];
sx q[3];
rz(-0.23947434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.60260281) q[2];
sx q[2];
rz(-2.4675641) q[2];
sx q[2];
rz(-1.4269786) q[2];
rz(1.1566628) q[3];
sx q[3];
rz(-1.7159228) q[3];
sx q[3];
rz(2.9837933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2222774) q[0];
sx q[0];
rz(-2.3265525) q[0];
sx q[0];
rz(-0.94006938) q[0];
rz(-2.6433511) q[1];
sx q[1];
rz(-2.2254641) q[1];
sx q[1];
rz(-1.1246276) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7263382) q[0];
sx q[0];
rz(-0.7043411) q[0];
sx q[0];
rz(-1.6319153) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1612879) q[2];
sx q[2];
rz(-1.7589658) q[2];
sx q[2];
rz(1.3649324) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.985253) q[1];
sx q[1];
rz(-0.96846995) q[1];
sx q[1];
rz(-2.1039318) q[1];
rz(-pi) q[2];
rz(-3.1203007) q[3];
sx q[3];
rz(-0.7856889) q[3];
sx q[3];
rz(-1.1710492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70790946) q[2];
sx q[2];
rz(-0.83432546) q[2];
sx q[2];
rz(-1.679812) q[2];
rz(-2.8953569) q[3];
sx q[3];
rz(-0.88053954) q[3];
sx q[3];
rz(-2.06854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6980625) q[0];
sx q[0];
rz(-1.227523) q[0];
sx q[0];
rz(-0.45502934) q[0];
rz(-0.21806923) q[1];
sx q[1];
rz(-2.2360305) q[1];
sx q[1];
rz(2.6630482) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6821176) q[0];
sx q[0];
rz(-0.90930258) q[0];
sx q[0];
rz(-1.4465998) q[0];
rz(-pi) q[1];
rz(-0.63963525) q[2];
sx q[2];
rz(-1.860485) q[2];
sx q[2];
rz(-1.3418496) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.26520935) q[1];
sx q[1];
rz(-0.57683101) q[1];
sx q[1];
rz(-0.93872197) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68799893) q[3];
sx q[3];
rz(-1.08537) q[3];
sx q[3];
rz(-1.345705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1856508) q[2];
sx q[2];
rz(-1.4943244) q[2];
sx q[2];
rz(0.68598023) q[2];
rz(-1.8020804) q[3];
sx q[3];
rz(-1.1687665) q[3];
sx q[3];
rz(1.7495988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5865536) q[0];
sx q[0];
rz(-1.3864484) q[0];
sx q[0];
rz(0.47898022) q[0];
rz(-2.2827177) q[1];
sx q[1];
rz(-0.54194599) q[1];
sx q[1];
rz(-0.10471334) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59706748) q[0];
sx q[0];
rz(-1.1865718) q[0];
sx q[0];
rz(-1.3659507) q[0];
rz(-2.5129086) q[2];
sx q[2];
rz(-2.5631944) q[2];
sx q[2];
rz(0.70369988) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0942451) q[1];
sx q[1];
rz(-0.99818238) q[1];
sx q[1];
rz(-2.8837573) q[1];
rz(-1.0927622) q[3];
sx q[3];
rz(-1.6217188) q[3];
sx q[3];
rz(-1.6855406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0571478) q[2];
sx q[2];
rz(-1.2556602) q[2];
sx q[2];
rz(-0.88195938) q[2];
rz(3.0706578) q[3];
sx q[3];
rz(-1.2627914) q[3];
sx q[3];
rz(-0.96496636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1607745) q[0];
sx q[0];
rz(-0.47176281) q[0];
sx q[0];
rz(-1.2534575) q[0];
rz(2.4313633) q[1];
sx q[1];
rz(-1.1980779) q[1];
sx q[1];
rz(-1.089878) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.010243) q[0];
sx q[0];
rz(-0.88307086) q[0];
sx q[0];
rz(-2.1136081) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2267672) q[2];
sx q[2];
rz(-1.1448749) q[2];
sx q[2];
rz(-2.3461208) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3115499) q[1];
sx q[1];
rz(-1.122992) q[1];
sx q[1];
rz(2.2954659) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23183159) q[3];
sx q[3];
rz(-0.58887312) q[3];
sx q[3];
rz(-0.26118054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39500427) q[2];
sx q[2];
rz(-0.85902625) q[2];
sx q[2];
rz(2.11917) q[2];
rz(-2.5452781) q[3];
sx q[3];
rz(-2.3583581) q[3];
sx q[3];
rz(-0.35308009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8280867) q[0];
sx q[0];
rz(-2.6441898) q[0];
sx q[0];
rz(1.5561546) q[0];
rz(1.4900788) q[1];
sx q[1];
rz(-0.45526344) q[1];
sx q[1];
rz(-2.1447287) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6108902) q[0];
sx q[0];
rz(-0.35323745) q[0];
sx q[0];
rz(-2.9363786) q[0];
x q[1];
rz(1.7559986) q[2];
sx q[2];
rz(-1.5723937) q[2];
sx q[2];
rz(-1.8271556) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.91040033) q[1];
sx q[1];
rz(-0.93437372) q[1];
sx q[1];
rz(0.6271805) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6410651) q[3];
sx q[3];
rz(-2.1033035) q[3];
sx q[3];
rz(0.73126572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8263714) q[2];
sx q[2];
rz(-0.73306495) q[2];
sx q[2];
rz(-0.23615393) q[2];
rz(-0.30596966) q[3];
sx q[3];
rz(-1.1924084) q[3];
sx q[3];
rz(1.4845622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.9295101) q[0];
sx q[0];
rz(-0.65123737) q[0];
sx q[0];
rz(-0.65810743) q[0];
rz(1.8923538) q[1];
sx q[1];
rz(-0.58964261) q[1];
sx q[1];
rz(-2.6499937) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.126287) q[0];
sx q[0];
rz(-1.1231866) q[0];
sx q[0];
rz(-2.5372895) q[0];
rz(-pi) q[1];
rz(-1.7936321) q[2];
sx q[2];
rz(-2.1249944) q[2];
sx q[2];
rz(-1.8775307) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2053615) q[1];
sx q[1];
rz(-0.27589162) q[1];
sx q[1];
rz(3.107168) q[1];
rz(-pi) q[2];
rz(3.0594143) q[3];
sx q[3];
rz(-2.058147) q[3];
sx q[3];
rz(0.10457071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7813985) q[2];
sx q[2];
rz(-2.7935544) q[2];
sx q[2];
rz(-1.6602824) q[2];
rz(-0.71634746) q[3];
sx q[3];
rz(-0.64645386) q[3];
sx q[3];
rz(-2.9334478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.79138712) q[0];
sx q[0];
rz(-1.5443784) q[0];
sx q[0];
rz(-2.7271893) q[0];
rz(2.1898337) q[1];
sx q[1];
rz(-1.2928243) q[1];
sx q[1];
rz(-1.181319) q[1];
rz(-2.4454115) q[2];
sx q[2];
rz(-2.0700818) q[2];
sx q[2];
rz(-2.5250057) q[2];
rz(-1.7813206) q[3];
sx q[3];
rz(-2.4518436) q[3];
sx q[3];
rz(2.8393815) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
