OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3611276) q[0];
sx q[0];
rz(-0.68929231) q[0];
sx q[0];
rz(-0.33049345) q[0];
rz(-2.8117872) q[1];
sx q[1];
rz(-2.2916315) q[1];
sx q[1];
rz(-0.70911521) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.37492) q[0];
sx q[0];
rz(-1.8639038) q[0];
sx q[0];
rz(2.2105182) q[0];
rz(1.8843295) q[2];
sx q[2];
rz(-1.6013813) q[2];
sx q[2];
rz(-0.41536301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4830556) q[1];
sx q[1];
rz(-2.4194948) q[1];
sx q[1];
rz(0.87941054) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0508918) q[3];
sx q[3];
rz(-1.9016148) q[3];
sx q[3];
rz(-1.0370805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4101397) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(-1.5343792) q[2];
rz(2.2065227) q[3];
sx q[3];
rz(-1.8362703) q[3];
sx q[3];
rz(-2.3527761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4193029) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(2.5193135) q[0];
rz(-0.17624804) q[1];
sx q[1];
rz(-1.9269678) q[1];
sx q[1];
rz(-0.91631779) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4692357) q[0];
sx q[0];
rz(-1.0807481) q[0];
sx q[0];
rz(0.83067466) q[0];
x q[1];
rz(1.7158521) q[2];
sx q[2];
rz(-0.74160355) q[2];
sx q[2];
rz(1.8650101) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1044836) q[1];
sx q[1];
rz(-0.47905211) q[1];
sx q[1];
rz(2.3399809) q[1];
rz(-pi) q[2];
rz(0.022577062) q[3];
sx q[3];
rz(-1.1728247) q[3];
sx q[3];
rz(0.28385362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9006485) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(-2.3201578) q[2];
rz(3.1243096) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(-2.3582874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.18773742) q[0];
sx q[0];
rz(-1.5783577) q[0];
sx q[0];
rz(2.1333372) q[0];
rz(3.1058274) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(-0.52454138) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13428155) q[0];
sx q[0];
rz(-1.5410005) q[0];
sx q[0];
rz(-1.7344463) q[0];
x q[1];
rz(2.8183297) q[2];
sx q[2];
rz(-0.70463902) q[2];
sx q[2];
rz(-1.7636253) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.14568612) q[1];
sx q[1];
rz(-1.6325163) q[1];
sx q[1];
rz(-1.593959) q[1];
rz(-1.2479765) q[3];
sx q[3];
rz(-0.4399235) q[3];
sx q[3];
rz(2.3390714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3699469) q[2];
sx q[2];
rz(-1.6235855) q[2];
sx q[2];
rz(-1.6463722) q[2];
rz(1.3211936) q[3];
sx q[3];
rz(-1.7318055) q[3];
sx q[3];
rz(-2.7041919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8111073) q[0];
sx q[0];
rz(-0.91902584) q[0];
sx q[0];
rz(3.0526429) q[0];
rz(-0.51070172) q[1];
sx q[1];
rz(-0.82534868) q[1];
sx q[1];
rz(0.68960062) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8423858) q[0];
sx q[0];
rz(-0.33873522) q[0];
sx q[0];
rz(2.6329106) q[0];
x q[1];
rz(1.3454516) q[2];
sx q[2];
rz(-1.0190522) q[2];
sx q[2];
rz(-1.4622886) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.00011132414) q[1];
sx q[1];
rz(-2.3173216) q[1];
sx q[1];
rz(2.6898756) q[1];
rz(-2.6570286) q[3];
sx q[3];
rz(-1.441941) q[3];
sx q[3];
rz(-1.7948922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4758063) q[2];
sx q[2];
rz(-0.67725724) q[2];
sx q[2];
rz(-0.62292567) q[2];
rz(-1.1359435) q[3];
sx q[3];
rz(-2.9562852) q[3];
sx q[3];
rz(-0.46666551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2899807) q[0];
sx q[0];
rz(-1.8988134) q[0];
sx q[0];
rz(0.583453) q[0];
rz(-1.1460229) q[1];
sx q[1];
rz(-1.4910411) q[1];
sx q[1];
rz(1.4978283) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9153584) q[0];
sx q[0];
rz(-2.3225473) q[0];
sx q[0];
rz(-0.35924964) q[0];
x q[1];
rz(1.4719109) q[2];
sx q[2];
rz(-0.77557287) q[2];
sx q[2];
rz(2.0688187) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0651107) q[1];
sx q[1];
rz(-1.5981734) q[1];
sx q[1];
rz(-2.4379424) q[1];
x q[2];
rz(-2.8605117) q[3];
sx q[3];
rz(-2.6087458) q[3];
sx q[3];
rz(-1.7588774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4862165) q[2];
sx q[2];
rz(-1.5465753) q[2];
sx q[2];
rz(1.5636469) q[2];
rz(-2.2359713) q[3];
sx q[3];
rz(-0.27035299) q[3];
sx q[3];
rz(0.63846987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49801302) q[0];
sx q[0];
rz(-1.6828515) q[0];
sx q[0];
rz(-3.0946099) q[0];
rz(0.14818305) q[1];
sx q[1];
rz(-2.2427185) q[1];
sx q[1];
rz(1.7061589) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5481789) q[0];
sx q[0];
rz(-2.3313064) q[0];
sx q[0];
rz(-0.20472783) q[0];
rz(1.4749182) q[2];
sx q[2];
rz(-0.55170689) q[2];
sx q[2];
rz(2.4193537) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6378577) q[1];
sx q[1];
rz(-1.6956455) q[1];
sx q[1];
rz(1.1605074) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4581497) q[3];
sx q[3];
rz(-2.9023691) q[3];
sx q[3];
rz(-0.84157543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0753714) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(0.071468778) q[2];
rz(1.6890769) q[3];
sx q[3];
rz(-1.444933) q[3];
sx q[3];
rz(-2.215109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8554095) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(-0.57762161) q[0];
rz(1.8619934) q[1];
sx q[1];
rz(-1.3459233) q[1];
sx q[1];
rz(-2.1320027) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44430915) q[0];
sx q[0];
rz(-2.5464006) q[0];
sx q[0];
rz(0.62023456) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7570653) q[2];
sx q[2];
rz(-1.3587917) q[2];
sx q[2];
rz(-2.7146313) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.91612591) q[1];
sx q[1];
rz(-1.9088609) q[1];
sx q[1];
rz(-0.52789968) q[1];
x q[2];
rz(-2.1647251) q[3];
sx q[3];
rz(-2.2403324) q[3];
sx q[3];
rz(0.23564786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0916831) q[2];
sx q[2];
rz(-2.78806) q[2];
sx q[2];
rz(1.9419149) q[2];
rz(-2.6692634) q[3];
sx q[3];
rz(-1.4940558) q[3];
sx q[3];
rz(1.002731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6427479) q[0];
sx q[0];
rz(-1.5002748) q[0];
sx q[0];
rz(0.2391267) q[0];
rz(-2.3587976) q[1];
sx q[1];
rz(-2.6233964) q[1];
sx q[1];
rz(0.4447287) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1040161) q[0];
sx q[0];
rz(-2.4009279) q[0];
sx q[0];
rz(-2.2292577) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4773981) q[2];
sx q[2];
rz(-2.1014629) q[2];
sx q[2];
rz(-1.2150089) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4152894) q[1];
sx q[1];
rz(-0.88978926) q[1];
sx q[1];
rz(-2.9773832) q[1];
rz(-pi) q[2];
rz(-1.4374251) q[3];
sx q[3];
rz(-1.1629472) q[3];
sx q[3];
rz(0.22920242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42177054) q[2];
sx q[2];
rz(-2.3193216) q[2];
sx q[2];
rz(0.81531173) q[2];
rz(-2.1067965) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(1.8243272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3141044) q[0];
sx q[0];
rz(-2.1056471) q[0];
sx q[0];
rz(-0.28717336) q[0];
rz(0.18889591) q[1];
sx q[1];
rz(-0.72383988) q[1];
sx q[1];
rz(-2.8093991) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0848541) q[0];
sx q[0];
rz(-2.2915654) q[0];
sx q[0];
rz(0.42215729) q[0];
x q[1];
rz(-3.095093) q[2];
sx q[2];
rz(-2.2691155) q[2];
sx q[2];
rz(2.7547714) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.48391446) q[1];
sx q[1];
rz(-2.6473443) q[1];
sx q[1];
rz(2.5792522) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6094749) q[3];
sx q[3];
rz(-1.9241153) q[3];
sx q[3];
rz(-1.6649099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2087848) q[2];
sx q[2];
rz(-1.00495) q[2];
sx q[2];
rz(0.83958158) q[2];
rz(1.2906637) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(2.4850142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.055450913) q[0];
sx q[0];
rz(-0.76776004) q[0];
sx q[0];
rz(-1.0593876) q[0];
rz(-2.1620031) q[1];
sx q[1];
rz(-1.9444119) q[1];
sx q[1];
rz(-2.8870781) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5352288) q[0];
sx q[0];
rz(-1.4974721) q[0];
sx q[0];
rz(1.3396157) q[0];
x q[1];
rz(1.326667) q[2];
sx q[2];
rz(-2.0460528) q[2];
sx q[2];
rz(0.82820669) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.75603112) q[1];
sx q[1];
rz(-0.33547151) q[1];
sx q[1];
rz(-2.8940593) q[1];
rz(-0.43838318) q[3];
sx q[3];
rz(-1.66072) q[3];
sx q[3];
rz(0.59059483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0467726) q[2];
sx q[2];
rz(-2.5929055) q[2];
sx q[2];
rz(1.0478896) q[2];
rz(-1.7808328) q[3];
sx q[3];
rz(-0.69793099) q[3];
sx q[3];
rz(0.60539436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.027325252) q[0];
sx q[0];
rz(-2.0226759) q[0];
sx q[0];
rz(-0.080060536) q[0];
rz(0.36021532) q[1];
sx q[1];
rz(-1.6811014) q[1];
sx q[1];
rz(-0.95492687) q[1];
rz(0.24047492) q[2];
sx q[2];
rz(-2.0839811) q[2];
sx q[2];
rz(2.6032084) q[2];
rz(0.9568692) q[3];
sx q[3];
rz(-1.3812243) q[3];
sx q[3];
rz(0.45637043) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
