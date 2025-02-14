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
rz(0.90647107) q[0];
sx q[0];
rz(-1.4426761) q[0];
sx q[0];
rz(0.28191167) q[0];
rz(-2.6126722) q[1];
sx q[1];
rz(-1.6493874) q[1];
sx q[1];
rz(-1.5780916) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96254565) q[0];
sx q[0];
rz(-1.7508649) q[0];
sx q[0];
rz(-1.1992886) q[0];
rz(-pi) q[1];
rz(0.018997832) q[2];
sx q[2];
rz(-2.8667031) q[2];
sx q[2];
rz(1.887448) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7719745) q[1];
sx q[1];
rz(-1.6965515) q[1];
sx q[1];
rz(-0.63196504) q[1];
x q[2];
rz(0.6880349) q[3];
sx q[3];
rz(-1.3322988) q[3];
sx q[3];
rz(1.9515338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6866744) q[2];
sx q[2];
rz(-2.8439971) q[2];
sx q[2];
rz(0.61304027) q[2];
rz(-0.4736627) q[3];
sx q[3];
rz(-1.1981755) q[3];
sx q[3];
rz(1.4325498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7365725) q[0];
sx q[0];
rz(-0.18017811) q[0];
sx q[0];
rz(0.75102425) q[0];
rz(-0.48149064) q[1];
sx q[1];
rz(-2.057169) q[1];
sx q[1];
rz(-2.1717333) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4879726) q[0];
sx q[0];
rz(-2.3832085) q[0];
sx q[0];
rz(0.55413306) q[0];
rz(1.5300691) q[2];
sx q[2];
rz(-1.3551522) q[2];
sx q[2];
rz(-3.0037896) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4837326) q[1];
sx q[1];
rz(-2.2071903) q[1];
sx q[1];
rz(-2.8191393) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0131575) q[3];
sx q[3];
rz(-2.1214888) q[3];
sx q[3];
rz(-2.9912419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.91886175) q[2];
sx q[2];
rz(-1.6657882) q[2];
sx q[2];
rz(-0.53885031) q[2];
rz(-3.1239964) q[3];
sx q[3];
rz(-2.9774057) q[3];
sx q[3];
rz(-0.9084475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4410412) q[0];
sx q[0];
rz(-0.38527641) q[0];
sx q[0];
rz(-3.0803296) q[0];
rz(1.6368658) q[1];
sx q[1];
rz(-0.69925362) q[1];
sx q[1];
rz(1.9452852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17249566) q[0];
sx q[0];
rz(-2.3845256) q[0];
sx q[0];
rz(1.9226546) q[0];
x q[1];
rz(-1.5012791) q[2];
sx q[2];
rz(-1.8325153) q[2];
sx q[2];
rz(1.6659322) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4704901) q[1];
sx q[1];
rz(-1.2467978) q[1];
sx q[1];
rz(-0.98172202) q[1];
rz(-pi) q[2];
rz(2.263271) q[3];
sx q[3];
rz(-1.3008144) q[3];
sx q[3];
rz(0.53966554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.68941826) q[2];
sx q[2];
rz(-1.3070561) q[2];
sx q[2];
rz(-2.7090731) q[2];
rz(-1.0906667) q[3];
sx q[3];
rz(-0.64041036) q[3];
sx q[3];
rz(1.0961016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5525621) q[0];
sx q[0];
rz(-0.93695372) q[0];
sx q[0];
rz(-0.20137782) q[0];
rz(0.51678139) q[1];
sx q[1];
rz(-0.38025451) q[1];
sx q[1];
rz(-2.1929599) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3654914) q[0];
sx q[0];
rz(-2.1946215) q[0];
sx q[0];
rz(0.43419773) q[0];
x q[1];
rz(-0.29251137) q[2];
sx q[2];
rz(-0.88380948) q[2];
sx q[2];
rz(2.3623737) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.779114) q[1];
sx q[1];
rz(-2.6370905) q[1];
sx q[1];
rz(-0.027536784) q[1];
rz(-pi) q[2];
rz(1.473677) q[3];
sx q[3];
rz(-2.0896974) q[3];
sx q[3];
rz(-0.023954602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2739233) q[2];
sx q[2];
rz(-2.7322768) q[2];
sx q[2];
rz(-0.55740994) q[2];
rz(0.080816001) q[3];
sx q[3];
rz(-1.9010474) q[3];
sx q[3];
rz(-1.2062937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
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
rz(-1.4163365) q[0];
sx q[0];
rz(-1.4755604) q[0];
sx q[0];
rz(-2.1112554) q[0];
rz(0.15580767) q[1];
sx q[1];
rz(-2.0741597) q[1];
sx q[1];
rz(2.9373998) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1847398) q[0];
sx q[0];
rz(-2.889688) q[0];
sx q[0];
rz(0.27916081) q[0];
rz(-2.5956018) q[2];
sx q[2];
rz(-2.8936849) q[2];
sx q[2];
rz(-1.966983) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9928007) q[1];
sx q[1];
rz(-0.53314236) q[1];
sx q[1];
rz(2.7270182) q[1];
rz(-pi) q[2];
rz(-2.89661) q[3];
sx q[3];
rz(-2.3534273) q[3];
sx q[3];
rz(0.80313166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.92945176) q[2];
sx q[2];
rz(-1.6484304) q[2];
sx q[2];
rz(-0.41352752) q[2];
rz(2.2996969) q[3];
sx q[3];
rz(-1.3827518) q[3];
sx q[3];
rz(-0.97257417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0464762) q[0];
sx q[0];
rz(-2.0501417) q[0];
sx q[0];
rz(-3.1360151) q[0];
rz(1.3439517) q[1];
sx q[1];
rz(-0.19897142) q[1];
sx q[1];
rz(0.70095789) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98788059) q[0];
sx q[0];
rz(-1.0544262) q[0];
sx q[0];
rz(2.6285281) q[0];
x q[1];
rz(2.1326601) q[2];
sx q[2];
rz(-0.72266662) q[2];
sx q[2];
rz(-2.780513) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.13853874) q[1];
sx q[1];
rz(-1.8459311) q[1];
sx q[1];
rz(-1.6134) q[1];
x q[2];
rz(0.57311975) q[3];
sx q[3];
rz(-2.0175211) q[3];
sx q[3];
rz(-2.5271811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.31412101) q[2];
sx q[2];
rz(-2.8077329) q[2];
sx q[2];
rz(1.1673048) q[2];
rz(-0.88268924) q[3];
sx q[3];
rz(-0.76789951) q[3];
sx q[3];
rz(-2.816443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2348787) q[0];
sx q[0];
rz(-1.119708) q[0];
sx q[0];
rz(-0.085302189) q[0];
rz(2.3872497) q[1];
sx q[1];
rz(-2.6383548) q[1];
sx q[1];
rz(-2.1139961) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7639973) q[0];
sx q[0];
rz(-0.21276424) q[0];
sx q[0];
rz(-0.048325267) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59664388) q[2];
sx q[2];
rz(-1.575257) q[2];
sx q[2];
rz(2.4031134) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3049406) q[1];
sx q[1];
rz(-1.2479559) q[1];
sx q[1];
rz(3.0129827) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2619234) q[3];
sx q[3];
rz(-2.4604049) q[3];
sx q[3];
rz(-0.78024125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1520285) q[2];
sx q[2];
rz(-2.1542408) q[2];
sx q[2];
rz(-2.7222471) q[2];
rz(-0.63465214) q[3];
sx q[3];
rz(-2.3366163) q[3];
sx q[3];
rz(2.0161207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80970508) q[0];
sx q[0];
rz(-0.87600791) q[0];
sx q[0];
rz(-2.4494655) q[0];
rz(-0.57811111) q[1];
sx q[1];
rz(-2.7222241) q[1];
sx q[1];
rz(1.3828166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7958) q[0];
sx q[0];
rz(-1.5032856) q[0];
sx q[0];
rz(3.1221005) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4695669) q[2];
sx q[2];
rz(-1.0254854) q[2];
sx q[2];
rz(1.7500851) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4647411) q[1];
sx q[1];
rz(-0.60982043) q[1];
sx q[1];
rz(-1.1552385) q[1];
rz(-pi) q[2];
rz(-2.2949831) q[3];
sx q[3];
rz(-1.9384406) q[3];
sx q[3];
rz(-1.5049619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3708923) q[2];
sx q[2];
rz(-0.53052491) q[2];
sx q[2];
rz(1.4250866) q[2];
rz(-2.6254081) q[3];
sx q[3];
rz(-2.1631212) q[3];
sx q[3];
rz(-1.2396575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4952963) q[0];
sx q[0];
rz(-1.4343867) q[0];
sx q[0];
rz(-0.37149757) q[0];
rz(-0.83483541) q[1];
sx q[1];
rz(-0.8131665) q[1];
sx q[1];
rz(0.017597839) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8397578) q[0];
sx q[0];
rz(-2.5389414) q[0];
sx q[0];
rz(1.0687123) q[0];
x q[1];
rz(0.73569466) q[2];
sx q[2];
rz(-1.7580877) q[2];
sx q[2];
rz(3.1295619) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0549893) q[1];
sx q[1];
rz(-1.2003683) q[1];
sx q[1];
rz(-2.5309674) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5334355) q[3];
sx q[3];
rz(-1.8377234) q[3];
sx q[3];
rz(-1.9723231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1060433) q[2];
sx q[2];
rz(-0.86842662) q[2];
sx q[2];
rz(3.0412728) q[2];
rz(-2.912168) q[3];
sx q[3];
rz(-2.094163) q[3];
sx q[3];
rz(-0.44982287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68596524) q[0];
sx q[0];
rz(-2.8534511) q[0];
sx q[0];
rz(-0.57383865) q[0];
rz(1.0539184) q[1];
sx q[1];
rz(-2.0951447) q[1];
sx q[1];
rz(-2.4651249) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1944502) q[0];
sx q[0];
rz(-1.8136171) q[0];
sx q[0];
rz(0.28698289) q[0];
rz(1.0208548) q[2];
sx q[2];
rz(-0.46560198) q[2];
sx q[2];
rz(-2.4684722) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6868848) q[1];
sx q[1];
rz(-2.096972) q[1];
sx q[1];
rz(1.673102) q[1];
rz(-pi) q[2];
rz(-1.2419534) q[3];
sx q[3];
rz(-0.55677215) q[3];
sx q[3];
rz(1.1278576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68710589) q[2];
sx q[2];
rz(-2.3837619) q[2];
sx q[2];
rz(1.6472316) q[2];
rz(2.2729661) q[3];
sx q[3];
rz(-0.70610154) q[3];
sx q[3];
rz(0.47006616) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070521991) q[0];
sx q[0];
rz(-1.1916397) q[0];
sx q[0];
rz(1.2217039) q[0];
rz(1.3337878) q[1];
sx q[1];
rz(-1.318327) q[1];
sx q[1];
rz(-0.64073906) q[1];
rz(-2.376198) q[2];
sx q[2];
rz(-1.5049878) q[2];
sx q[2];
rz(2.4743248) q[2];
rz(-3.0972539) q[3];
sx q[3];
rz(-1.7462742) q[3];
sx q[3];
rz(-2.3901226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
