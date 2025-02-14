OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.3761223) q[0];
sx q[0];
rz(2.0630615) q[0];
sx q[0];
rz(10.800092) q[0];
rz(-2.7219661) q[1];
sx q[1];
rz(-1.3803177) q[1];
sx q[1];
rz(-0.79545155) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8574816) q[0];
sx q[0];
rz(-1.6499551) q[0];
sx q[0];
rz(2.7678431) q[0];
rz(2.8654557) q[2];
sx q[2];
rz(-2.727446) q[2];
sx q[2];
rz(2.2091231) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4013942) q[1];
sx q[1];
rz(-2.2371082) q[1];
sx q[1];
rz(1.3191651) q[1];
rz(-pi) q[2];
rz(3.0555435) q[3];
sx q[3];
rz(-1.3121735) q[3];
sx q[3];
rz(1.2899959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2835283) q[2];
sx q[2];
rz(-2.8499446) q[2];
sx q[2];
rz(0.38823271) q[2];
rz(-2.9325716) q[3];
sx q[3];
rz(-1.4744604) q[3];
sx q[3];
rz(2.2504263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.46880576) q[0];
sx q[0];
rz(-0.46872941) q[0];
sx q[0];
rz(2.3720473) q[0];
rz(2.1107213) q[1];
sx q[1];
rz(-0.85846916) q[1];
sx q[1];
rz(1.4268202) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1311817) q[0];
sx q[0];
rz(-2.7658434) q[0];
sx q[0];
rz(0.6300169) q[0];
x q[1];
rz(1.6278817) q[2];
sx q[2];
rz(-0.89623129) q[2];
sx q[2];
rz(-0.13007535) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.12995686) q[1];
sx q[1];
rz(-2.4677688) q[1];
sx q[1];
rz(2.7259299) q[1];
rz(-pi) q[2];
rz(1.5094425) q[3];
sx q[3];
rz(-2.4262145) q[3];
sx q[3];
rz(0.74980199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0017121) q[2];
sx q[2];
rz(-2.6931245) q[2];
sx q[2];
rz(1.2343538) q[2];
rz(2.6185696) q[3];
sx q[3];
rz(-2.0165899) q[3];
sx q[3];
rz(0.6559059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4457542) q[0];
sx q[0];
rz(-0.0025175968) q[0];
sx q[0];
rz(-2.5939831) q[0];
rz(2.2721263) q[1];
sx q[1];
rz(-0.97620669) q[1];
sx q[1];
rz(-1.4617823) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48058352) q[0];
sx q[0];
rz(-0.55381539) q[0];
sx q[0];
rz(-0.69032918) q[0];
x q[1];
rz(2.3794054) q[2];
sx q[2];
rz(-1.7584287) q[2];
sx q[2];
rz(-2.4248276) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3981769) q[1];
sx q[1];
rz(-1.2733439) q[1];
sx q[1];
rz(-1.1156082) q[1];
rz(0.84949228) q[3];
sx q[3];
rz(-0.85334001) q[3];
sx q[3];
rz(-0.8161374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4168641) q[2];
sx q[2];
rz(-2.8432196) q[2];
sx q[2];
rz(2.5073063) q[2];
rz(1.0861081) q[3];
sx q[3];
rz(-0.48091286) q[3];
sx q[3];
rz(-2.0079131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.5241765) q[0];
sx q[0];
rz(-0.15707459) q[0];
sx q[0];
rz(-1.6198535) q[0];
rz(-2.1583083) q[1];
sx q[1];
rz(-2.0359928) q[1];
sx q[1];
rz(-3.0599248) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3043609) q[0];
sx q[0];
rz(-2.1262693) q[0];
sx q[0];
rz(-0.13120478) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1278243) q[2];
sx q[2];
rz(-1.9038611) q[2];
sx q[2];
rz(-0.027160732) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5809243) q[1];
sx q[1];
rz(-2.4663975) q[1];
sx q[1];
rz(-0.84431727) q[1];
rz(-pi) q[2];
rz(-2.0606629) q[3];
sx q[3];
rz(-1.6182163) q[3];
sx q[3];
rz(2.1159548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4982768) q[2];
sx q[2];
rz(-1.4772819) q[2];
sx q[2];
rz(0.91900438) q[2];
rz(-0.32215858) q[3];
sx q[3];
rz(-1.6291658) q[3];
sx q[3];
rz(2.7951516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0568986) q[0];
sx q[0];
rz(-1.1321122) q[0];
sx q[0];
rz(-1.3380916) q[0];
rz(-2.7390506) q[1];
sx q[1];
rz(-1.2502547) q[1];
sx q[1];
rz(1.5571627) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65619233) q[0];
sx q[0];
rz(-1.3043405) q[0];
sx q[0];
rz(1.043826) q[0];
rz(-0.7942368) q[2];
sx q[2];
rz(-1.4068687) q[2];
sx q[2];
rz(-1.2212703) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8338523) q[1];
sx q[1];
rz(-0.55451143) q[1];
sx q[1];
rz(0.28077287) q[1];
x q[2];
rz(1.1013056) q[3];
sx q[3];
rz(-1.3302505) q[3];
sx q[3];
rz(-2.8194373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1160139) q[2];
sx q[2];
rz(-1.2776926) q[2];
sx q[2];
rz(-2.4063827) q[2];
rz(-2.0871346) q[3];
sx q[3];
rz(-2.1078883) q[3];
sx q[3];
rz(1.2247156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9758107) q[0];
sx q[0];
rz(-2.7175792) q[0];
sx q[0];
rz(1.8884416) q[0];
rz(-1.5755298) q[1];
sx q[1];
rz(-2.0063446) q[1];
sx q[1];
rz(-2.6083686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0237817) q[0];
sx q[0];
rz(-0.83988777) q[0];
sx q[0];
rz(0.31857849) q[0];
rz(-pi) q[1];
rz(2.8715233) q[2];
sx q[2];
rz(-0.52723072) q[2];
sx q[2];
rz(-0.6747077) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8044843) q[1];
sx q[1];
rz(-1.3829559) q[1];
sx q[1];
rz(-2.9153182) q[1];
rz(-0.56508692) q[3];
sx q[3];
rz(-2.0425519) q[3];
sx q[3];
rz(-0.12607546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.48625654) q[2];
sx q[2];
rz(-0.37898263) q[2];
sx q[2];
rz(0.8302702) q[2];
rz(1.8298979) q[3];
sx q[3];
rz(-2.0528841) q[3];
sx q[3];
rz(-2.8970498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.526392) q[0];
sx q[0];
rz(-1.6561693) q[0];
sx q[0];
rz(0.8335337) q[0];
rz(-0.99710456) q[1];
sx q[1];
rz(-1.763696) q[1];
sx q[1];
rz(-1.5514143) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1520321) q[0];
sx q[0];
rz(-0.44267926) q[0];
sx q[0];
rz(-2.9665806) q[0];
rz(2.1623108) q[2];
sx q[2];
rz(-1.2752609) q[2];
sx q[2];
rz(-2.599567) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.60234705) q[1];
sx q[1];
rz(-1.6825469) q[1];
sx q[1];
rz(1.69126) q[1];
rz(2.6069712) q[3];
sx q[3];
rz(-1.7869014) q[3];
sx q[3];
rz(0.099319746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80018631) q[2];
sx q[2];
rz(-2.1087346) q[2];
sx q[2];
rz(1.7436854) q[2];
rz(2.7684815) q[3];
sx q[3];
rz(-2.2171376) q[3];
sx q[3];
rz(2.1227409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.43847346) q[0];
sx q[0];
rz(-1.6699474) q[0];
sx q[0];
rz(0.31914172) q[0];
rz(-1.9769662) q[1];
sx q[1];
rz(-1.9705557) q[1];
sx q[1];
rz(-3.107531) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3122404) q[0];
sx q[0];
rz(-2.3213648) q[0];
sx q[0];
rz(-1.149631) q[0];
rz(-2.8197779) q[2];
sx q[2];
rz(-2.3452367) q[2];
sx q[2];
rz(-2.1200695) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0314591) q[1];
sx q[1];
rz(-1.5601399) q[1];
sx q[1];
rz(1.5489122) q[1];
rz(-pi) q[2];
rz(0.56896992) q[3];
sx q[3];
rz(-0.60980699) q[3];
sx q[3];
rz(-2.7066674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.065319149) q[2];
sx q[2];
rz(-1.7460145) q[2];
sx q[2];
rz(-3.0714387) q[2];
rz(-1.3089199) q[3];
sx q[3];
rz(-0.85246837) q[3];
sx q[3];
rz(-3.0064228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0604621) q[0];
sx q[0];
rz(-1.9996996) q[0];
sx q[0];
rz(-0.4749701) q[0];
rz(-1.9604856) q[1];
sx q[1];
rz(-2.2468552) q[1];
sx q[1];
rz(2.2244577) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48045298) q[0];
sx q[0];
rz(-2.5256143) q[0];
sx q[0];
rz(2.7333214) q[0];
rz(-pi) q[1];
rz(-1.195329) q[2];
sx q[2];
rz(-0.72403833) q[2];
sx q[2];
rz(-0.53508102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.18439427) q[1];
sx q[1];
rz(-2.9177114) q[1];
sx q[1];
rz(1.9489524) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5417622) q[3];
sx q[3];
rz(-1.2496061) q[3];
sx q[3];
rz(2.7293492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.26943031) q[2];
sx q[2];
rz(-1.3698801) q[2];
sx q[2];
rz(1.5838464) q[2];
rz(0.15548429) q[3];
sx q[3];
rz(-1.0289861) q[3];
sx q[3];
rz(-0.99229971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6601722) q[0];
sx q[0];
rz(-1.476113) q[0];
sx q[0];
rz(2.7760264) q[0];
rz(1.9407326) q[1];
sx q[1];
rz(-1.6114085) q[1];
sx q[1];
rz(-1.0467122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7378254) q[0];
sx q[0];
rz(-2.1008472) q[0];
sx q[0];
rz(-2.5395495) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44136845) q[2];
sx q[2];
rz(-2.0725996) q[2];
sx q[2];
rz(2.8428915) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9871708) q[1];
sx q[1];
rz(-2.1959496) q[1];
sx q[1];
rz(2.465017) q[1];
rz(-2.3975054) q[3];
sx q[3];
rz(-0.73089337) q[3];
sx q[3];
rz(2.6958974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4552348) q[2];
sx q[2];
rz(-2.0265667) q[2];
sx q[2];
rz(0.70351797) q[2];
rz(-2.5840058) q[3];
sx q[3];
rz(-2.360354) q[3];
sx q[3];
rz(2.668146) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8483509) q[0];
sx q[0];
rz(-2.1108755) q[0];
sx q[0];
rz(-2.1847771) q[0];
rz(-1.035607) q[1];
sx q[1];
rz(-2.5685812) q[1];
sx q[1];
rz(-2.0249637) q[1];
rz(2.4734617) q[2];
sx q[2];
rz(-1.1429759) q[2];
sx q[2];
rz(-0.10980724) q[2];
rz(-2.9241548) q[3];
sx q[3];
rz(-1.5296546) q[3];
sx q[3];
rz(-1.1543867) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
