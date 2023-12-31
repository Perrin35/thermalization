OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10652868) q[0];
sx q[0];
rz(-1.0892692) q[0];
sx q[0];
rz(0.16103345) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(-2.6452126) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.046145) q[0];
sx q[0];
rz(-0.60337043) q[0];
sx q[0];
rz(-1.3674111) q[0];
rz(-pi) q[1];
rz(2.6167294) q[2];
sx q[2];
rz(-2.7537991) q[2];
sx q[2];
rz(0.85688574) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1094184) q[1];
sx q[1];
rz(-2.2757581) q[1];
sx q[1];
rz(1.6900307) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6707889) q[3];
sx q[3];
rz(-1.9766207) q[3];
sx q[3];
rz(1.0909181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6266142) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(-2.853945) q[2];
rz(1.7488165) q[3];
sx q[3];
rz(-2.1578622) q[3];
sx q[3];
rz(2.0387409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87067938) q[0];
sx q[0];
rz(-0.7551071) q[0];
sx q[0];
rz(-2.5640008) q[0];
rz(-1.5354935) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(1.577852) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6710885) q[0];
sx q[0];
rz(-1.7894063) q[0];
sx q[0];
rz(1.455362) q[0];
rz(-pi) q[1];
rz(0.30303843) q[2];
sx q[2];
rz(-1.3088471) q[2];
sx q[2];
rz(2.9532202) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4430122) q[1];
sx q[1];
rz(-2.6786782) q[1];
sx q[1];
rz(1.7673311) q[1];
rz(-pi) q[2];
rz(2.9868449) q[3];
sx q[3];
rz(-1.4156716) q[3];
sx q[3];
rz(0.89150067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0343798) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(2.0593026) q[2];
rz(1.1632495) q[3];
sx q[3];
rz(-1.2000368) q[3];
sx q[3];
rz(2.5015586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0619693) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(0.18297718) q[0];
rz(-3.0984763) q[1];
sx q[1];
rz(-2.1886107) q[1];
sx q[1];
rz(1.144369) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17830081) q[0];
sx q[0];
rz(-0.68094567) q[0];
sx q[0];
rz(-0.91239022) q[0];
x q[1];
rz(1.794533) q[2];
sx q[2];
rz(-1.9800595) q[2];
sx q[2];
rz(2.4350016) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5971165) q[1];
sx q[1];
rz(-0.58864486) q[1];
sx q[1];
rz(-2.4878923) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8539682) q[3];
sx q[3];
rz(-1.5678741) q[3];
sx q[3];
rz(-2.6401816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.019505067) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(-2.2369177) q[2];
rz(-2.4217862) q[3];
sx q[3];
rz(-2.7147839) q[3];
sx q[3];
rz(-2.3864746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4633789) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(2.4257207) q[0];
rz(2.8158358) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(-0.99266565) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216114) q[0];
sx q[0];
rz(-2.5489759) q[0];
sx q[0];
rz(2.6556334) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86187141) q[2];
sx q[2];
rz(-2.3700691) q[2];
sx q[2];
rz(2.9574403) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6550733) q[1];
sx q[1];
rz(-1.7336978) q[1];
sx q[1];
rz(-1.4763114) q[1];
rz(-pi) q[2];
rz(-2.643814) q[3];
sx q[3];
rz(-2.0911502) q[3];
sx q[3];
rz(-1.9073515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0830393) q[2];
sx q[2];
rz(-0.71648592) q[2];
sx q[2];
rz(-1.4286208) q[2];
rz(2.2955017) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5457526) q[0];
sx q[0];
rz(-1.7359474) q[0];
sx q[0];
rz(-0.76675057) q[0];
rz(1.2402361) q[1];
sx q[1];
rz(-2.2940472) q[1];
sx q[1];
rz(-0.3516745) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85580101) q[0];
sx q[0];
rz(-1.8579146) q[0];
sx q[0];
rz(0.65335269) q[0];
rz(-pi) q[1];
rz(0.50260966) q[2];
sx q[2];
rz(-0.99684274) q[2];
sx q[2];
rz(0.9290907) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1406447) q[1];
sx q[1];
rz(-1.5854156) q[1];
sx q[1];
rz(-0.39477243) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0293803) q[3];
sx q[3];
rz(-0.57313985) q[3];
sx q[3];
rz(1.0669607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9436283) q[2];
sx q[2];
rz(-1.6683656) q[2];
sx q[2];
rz(2.6082805) q[2];
rz(0.42896459) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(1.126948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0241942) q[0];
sx q[0];
rz(-1.0929996) q[0];
sx q[0];
rz(-0.54164106) q[0];
rz(2.533124) q[1];
sx q[1];
rz(-2.922373) q[1];
sx q[1];
rz(-2.0419962) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8894419) q[0];
sx q[0];
rz(-2.5398206) q[0];
sx q[0];
rz(-2.2579231) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1018283) q[2];
sx q[2];
rz(-1.7341988) q[2];
sx q[2];
rz(2.642717) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9601599) q[1];
sx q[1];
rz(-0.90248855) q[1];
sx q[1];
rz(-1.798435) q[1];
rz(2.3533456) q[3];
sx q[3];
rz(-1.3282093) q[3];
sx q[3];
rz(-1.5188252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.6301443) q[2];
sx q[2];
rz(-1.5254598) q[2];
sx q[2];
rz(2.8586094) q[2];
rz(-2.4354368) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(-0.63505665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054984897) q[0];
sx q[0];
rz(-2.1540756) q[0];
sx q[0];
rz(-1.6116066) q[0];
rz(2.4967172) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(2.9842916) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0387715) q[0];
sx q[0];
rz(-0.73909315) q[0];
sx q[0];
rz(3.0022013) q[0];
rz(-pi) q[1];
rz(0.83432014) q[2];
sx q[2];
rz(-1.4093219) q[2];
sx q[2];
rz(-1.8404567) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.339817) q[1];
sx q[1];
rz(-2.1784557) q[1];
sx q[1];
rz(-2.9031309) q[1];
rz(-1.8623452) q[3];
sx q[3];
rz(-2.1778717) q[3];
sx q[3];
rz(2.8214422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1823696) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(-0.83089337) q[2];
rz(-3.0977141) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(3.1307401) q[0];
rz(-0.27663484) q[1];
sx q[1];
rz(-0.77181548) q[1];
sx q[1];
rz(-2.8483134) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0115259) q[0];
sx q[0];
rz(-1.5767326) q[0];
sx q[0];
rz(2.5991711) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3515527) q[2];
sx q[2];
rz(-2.4862137) q[2];
sx q[2];
rz(-0.13742451) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9353232) q[1];
sx q[1];
rz(-2.2895797) q[1];
sx q[1];
rz(-2.1410336) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6777612) q[3];
sx q[3];
rz(-1.9411191) q[3];
sx q[3];
rz(2.2390389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.04348065) q[2];
sx q[2];
rz(-3.0467693) q[2];
sx q[2];
rz(3.0528255) q[2];
rz(0.25012112) q[3];
sx q[3];
rz(-2.0177896) q[3];
sx q[3];
rz(0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1373238) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(2.0776757) q[0];
rz(2.6990199) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(-1.7907422) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3063072) q[0];
sx q[0];
rz(-1.4915691) q[0];
sx q[0];
rz(1.6552734) q[0];
rz(-pi) q[1];
rz(-1.4883792) q[2];
sx q[2];
rz(-0.91225183) q[2];
sx q[2];
rz(1.8694307) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.46637725) q[1];
sx q[1];
rz(-1.5870759) q[1];
sx q[1];
rz(-2.0051763) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1397347) q[3];
sx q[3];
rz(-2.8981879) q[3];
sx q[3];
rz(-0.36153015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.110048) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(2.6943977) q[2];
rz(1.3859008) q[3];
sx q[3];
rz(-0.30062506) q[3];
sx q[3];
rz(0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8266325) q[0];
sx q[0];
rz(-1.5216014) q[0];
sx q[0];
rz(0.78053027) q[0];
rz(2.7138117) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(-2.9719877) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34925941) q[0];
sx q[0];
rz(-1.0997286) q[0];
sx q[0];
rz(0.40339289) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5461966) q[2];
sx q[2];
rz(-1.149017) q[2];
sx q[2];
rz(-1.8116902) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5055713) q[1];
sx q[1];
rz(-0.95658703) q[1];
sx q[1];
rz(3.0843656) q[1];
rz(-pi) q[2];
rz(1.8923558) q[3];
sx q[3];
rz(-1.7064629) q[3];
sx q[3];
rz(-0.78007573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2910989) q[2];
sx q[2];
rz(-1.9591745) q[2];
sx q[2];
rz(-0.69236857) q[2];
rz(2.678357) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(2.0326116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8511843) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(0.51207536) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(-0.26906536) q[2];
sx q[2];
rz(-1.27956) q[2];
sx q[2];
rz(0.051824311) q[2];
rz(-3.0242596) q[3];
sx q[3];
rz(-1.3640079) q[3];
sx q[3];
rz(-0.97466536) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
