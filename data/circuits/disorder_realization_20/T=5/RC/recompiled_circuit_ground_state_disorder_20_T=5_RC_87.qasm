OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.60431689) q[0];
sx q[0];
rz(-2.9000227) q[0];
sx q[0];
rz(-0.33302745) q[0];
rz(-1.1355407) q[1];
sx q[1];
rz(3.9685213) q[1];
sx q[1];
rz(10.068738) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7111904) q[0];
sx q[0];
rz(-1.2997753) q[0];
sx q[0];
rz(0.94078101) q[0];
x q[1];
rz(-1.7276554) q[2];
sx q[2];
rz(-1.867395) q[2];
sx q[2];
rz(3.0576884) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.17929303) q[1];
sx q[1];
rz(-0.4518309) q[1];
sx q[1];
rz(0.56698124) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1584362) q[3];
sx q[3];
rz(-1.8915911) q[3];
sx q[3];
rz(-3.0123346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6551299) q[2];
sx q[2];
rz(-0.4643521) q[2];
sx q[2];
rz(2.4008524) q[2];
rz(1.5517392) q[3];
sx q[3];
rz(-2.430075) q[3];
sx q[3];
rz(-0.5927425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7084259) q[0];
sx q[0];
rz(-1.1335224) q[0];
sx q[0];
rz(0.11696996) q[0];
rz(-0.50432694) q[1];
sx q[1];
rz(-1.802899) q[1];
sx q[1];
rz(0.28409827) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4000716) q[0];
sx q[0];
rz(-1.0062381) q[0];
sx q[0];
rz(-3.0429521) q[0];
rz(-pi) q[1];
rz(-2.3440222) q[2];
sx q[2];
rz(-1.1402545) q[2];
sx q[2];
rz(2.0883462) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.57840216) q[1];
sx q[1];
rz(-0.12688533) q[1];
sx q[1];
rz(1.2940426) q[1];
x q[2];
rz(0.25155622) q[3];
sx q[3];
rz(-1.8402035) q[3];
sx q[3];
rz(-0.081307383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32610193) q[2];
sx q[2];
rz(-2.2559866) q[2];
sx q[2];
rz(2.3808114) q[2];
rz(-0.86959362) q[3];
sx q[3];
rz(-1.875501) q[3];
sx q[3];
rz(-0.45309711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78075439) q[0];
sx q[0];
rz(-2.0212845) q[0];
sx q[0];
rz(-0.25217062) q[0];
rz(-1.699327) q[1];
sx q[1];
rz(-2.1863054) q[1];
sx q[1];
rz(-0.35626492) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29554825) q[0];
sx q[0];
rz(-1.5112119) q[0];
sx q[0];
rz(0.0061959717) q[0];
rz(0.61932694) q[2];
sx q[2];
rz(-0.64420036) q[2];
sx q[2];
rz(-1.2540115) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.60881847) q[1];
sx q[1];
rz(-1.428481) q[1];
sx q[1];
rz(-0.93138333) q[1];
x q[2];
rz(1.8271133) q[3];
sx q[3];
rz(-1.5689092) q[3];
sx q[3];
rz(-2.6607115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1189271) q[2];
sx q[2];
rz(-1.7513195) q[2];
sx q[2];
rz(1.6290132) q[2];
rz(3.1318943) q[3];
sx q[3];
rz(-1.9111948) q[3];
sx q[3];
rz(0.62140083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.4951204) q[0];
sx q[0];
rz(-2.5992114) q[0];
sx q[0];
rz(-2.8908308) q[0];
rz(-1.3465025) q[1];
sx q[1];
rz(-0.52918068) q[1];
sx q[1];
rz(-1.4432602) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051606962) q[0];
sx q[0];
rz(-2.3653226) q[0];
sx q[0];
rz(-1.2597741) q[0];
rz(-pi) q[1];
rz(-2.2002831) q[2];
sx q[2];
rz(-2.513859) q[2];
sx q[2];
rz(2.8986487) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.52407284) q[1];
sx q[1];
rz(-1.0403324) q[1];
sx q[1];
rz(-0.73442119) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.79899995) q[3];
sx q[3];
rz(-0.36741396) q[3];
sx q[3];
rz(1.1297117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.5556339) q[2];
sx q[2];
rz(-1.1052174) q[2];
sx q[2];
rz(2.8982758) q[2];
rz(-2.5569052) q[3];
sx q[3];
rz(-0.57716113) q[3];
sx q[3];
rz(0.028133597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1481767) q[0];
sx q[0];
rz(-0.36574829) q[0];
sx q[0];
rz(-2.962033) q[0];
rz(-1.2252294) q[1];
sx q[1];
rz(-1.7370677) q[1];
sx q[1];
rz(-0.45825759) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69773795) q[0];
sx q[0];
rz(-2.8354044) q[0];
sx q[0];
rz(-2.8809887) q[0];
rz(2.8446537) q[2];
sx q[2];
rz(-2.3910948) q[2];
sx q[2];
rz(-0.82210449) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.51650713) q[1];
sx q[1];
rz(-1.436215) q[1];
sx q[1];
rz(3.0662886) q[1];
rz(-pi) q[2];
rz(2.2716801) q[3];
sx q[3];
rz(-1.8509838) q[3];
sx q[3];
rz(-2.3126569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7563584) q[2];
sx q[2];
rz(-0.42432722) q[2];
sx q[2];
rz(2.2198086) q[2];
rz(-1.2515757) q[3];
sx q[3];
rz(-1.4082963) q[3];
sx q[3];
rz(-3.0799227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0475622) q[0];
sx q[0];
rz(-0.91218364) q[0];
sx q[0];
rz(0.14866522) q[0];
rz(-2.8246763) q[1];
sx q[1];
rz(-0.47873679) q[1];
sx q[1];
rz(-2.5247578) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6233404) q[0];
sx q[0];
rz(-0.099661544) q[0];
sx q[0];
rz(-0.90604337) q[0];
rz(-2.0111175) q[2];
sx q[2];
rz(-1.8991578) q[2];
sx q[2];
rz(-0.87000123) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9636391) q[1];
sx q[1];
rz(-0.63040367) q[1];
sx q[1];
rz(1.1373345) q[1];
x q[2];
rz(1.0471763) q[3];
sx q[3];
rz(-1.8111808) q[3];
sx q[3];
rz(-2.2464858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.8646912) q[2];
sx q[2];
rz(-1.428705) q[2];
sx q[2];
rz(-3.1332664) q[2];
rz(-2.5152123) q[3];
sx q[3];
rz(-2.4255987) q[3];
sx q[3];
rz(0.00055073784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3061227) q[0];
sx q[0];
rz(-1.9897113) q[0];
sx q[0];
rz(-2.6595111) q[0];
rz(0.029190633) q[1];
sx q[1];
rz(-1.2960641) q[1];
sx q[1];
rz(0.76404244) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4493443) q[0];
sx q[0];
rz(-1.5034165) q[0];
sx q[0];
rz(-1.9420366) q[0];
x q[1];
rz(-2.2305626) q[2];
sx q[2];
rz(-1.5507959) q[2];
sx q[2];
rz(-2.3338712) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.58043874) q[1];
sx q[1];
rz(-2.517546) q[1];
sx q[1];
rz(1.8638205) q[1];
rz(-1.5519616) q[3];
sx q[3];
rz(-0.40098396) q[3];
sx q[3];
rz(0.31664059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.46334106) q[2];
sx q[2];
rz(-1.8222787) q[2];
sx q[2];
rz(2.6574262) q[2];
rz(2.4593027) q[3];
sx q[3];
rz(-0.65410084) q[3];
sx q[3];
rz(-3.0068523) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.084455647) q[0];
sx q[0];
rz(-0.77780044) q[0];
sx q[0];
rz(1.2917668) q[0];
rz(0.46503398) q[1];
sx q[1];
rz(-0.51883042) q[1];
sx q[1];
rz(-2.8344287) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7503459) q[0];
sx q[0];
rz(-2.0758817) q[0];
sx q[0];
rz(0.48663346) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2276344) q[2];
sx q[2];
rz(-1.982053) q[2];
sx q[2];
rz(0.62818254) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7502082) q[1];
sx q[1];
rz(-2.5473875) q[1];
sx q[1];
rz(1.2886402) q[1];
rz(-1.5504595) q[3];
sx q[3];
rz(-2.1957948) q[3];
sx q[3];
rz(-0.53903264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.18053599) q[2];
sx q[2];
rz(-2.7013216) q[2];
sx q[2];
rz(0.72009909) q[2];
rz(2.2911206) q[3];
sx q[3];
rz(-1.1668147) q[3];
sx q[3];
rz(-1.7299962) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20486031) q[0];
sx q[0];
rz(-0.16848773) q[0];
sx q[0];
rz(-2.4334461) q[0];
rz(-1.6234966) q[1];
sx q[1];
rz(-2.203439) q[1];
sx q[1];
rz(-0.35619563) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9954279) q[0];
sx q[0];
rz(-1.08279) q[0];
sx q[0];
rz(-2.0515217) q[0];
rz(0.041650305) q[2];
sx q[2];
rz(-2.5232968) q[2];
sx q[2];
rz(0.30660812) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5070008) q[1];
sx q[1];
rz(-2.2193877) q[1];
sx q[1];
rz(-1.6084987) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1061927) q[3];
sx q[3];
rz(-1.2038017) q[3];
sx q[3];
rz(1.7758241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0492964) q[2];
sx q[2];
rz(-0.80731097) q[2];
sx q[2];
rz(0.88579196) q[2];
rz(-0.93585912) q[3];
sx q[3];
rz(-1.1805308) q[3];
sx q[3];
rz(-0.22127557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.43779272) q[0];
sx q[0];
rz(-3.0942823) q[0];
sx q[0];
rz(1.4495151) q[0];
rz(1.0225147) q[1];
sx q[1];
rz(-2.6578564) q[1];
sx q[1];
rz(-1.6922916) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81240772) q[0];
sx q[0];
rz(-1.158876) q[0];
sx q[0];
rz(0.054188577) q[0];
x q[1];
rz(-1.9841393) q[2];
sx q[2];
rz(-1.9660232) q[2];
sx q[2];
rz(2.4918567) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.30561799) q[1];
sx q[1];
rz(-2.545649) q[1];
sx q[1];
rz(-1.9553595) q[1];
rz(-2.9935097) q[3];
sx q[3];
rz(-1.8335153) q[3];
sx q[3];
rz(0.99276357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.25541043) q[2];
sx q[2];
rz(-1.1899199) q[2];
sx q[2];
rz(2.4667242) q[2];
rz(0.97149649) q[3];
sx q[3];
rz(-2.4392305) q[3];
sx q[3];
rz(-1.6500047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7991199) q[0];
sx q[0];
rz(-1.5958888) q[0];
sx q[0];
rz(2.2842443) q[0];
rz(2.9091861) q[1];
sx q[1];
rz(-1.0127761) q[1];
sx q[1];
rz(1.3565328) q[1];
rz(-2.7276298) q[2];
sx q[2];
rz(-1.0775177) q[2];
sx q[2];
rz(-2.1106363) q[2];
rz(0.84216778) q[3];
sx q[3];
rz(-1.6292059) q[3];
sx q[3];
rz(-0.90730351) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
