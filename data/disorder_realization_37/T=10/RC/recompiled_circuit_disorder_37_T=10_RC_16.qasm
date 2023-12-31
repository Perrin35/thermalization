OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10387575) q[0];
sx q[0];
rz(1.2020943) q[0];
sx q[0];
rz(10.572875) q[0];
rz(-1.8885053) q[1];
sx q[1];
rz(-0.94068599) q[1];
sx q[1];
rz(1.747945) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0855334) q[0];
sx q[0];
rz(-0.30963184) q[0];
sx q[0];
rz(0.0066030533) q[0];
x q[1];
rz(2.355028) q[2];
sx q[2];
rz(-2.2872891) q[2];
sx q[2];
rz(0.63209817) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4204971) q[1];
sx q[1];
rz(-1.3008529) q[1];
sx q[1];
rz(2.6707778) q[1];
rz(-pi) q[2];
rz(2.1288793) q[3];
sx q[3];
rz(-1.3289684) q[3];
sx q[3];
rz(1.1954607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.48774886) q[2];
sx q[2];
rz(-1.8493435) q[2];
sx q[2];
rz(3.0207108) q[2];
rz(-2.9623048) q[3];
sx q[3];
rz(-0.59569734) q[3];
sx q[3];
rz(-0.15549913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-3.0497465) q[0];
sx q[0];
rz(-2.3738528) q[0];
sx q[0];
rz(3.0088186) q[0];
rz(-1.4615387) q[1];
sx q[1];
rz(-1.5613873) q[1];
sx q[1];
rz(0.24138385) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5854908) q[0];
sx q[0];
rz(-1.9560768) q[0];
sx q[0];
rz(-0.38498621) q[0];
rz(2.6163289) q[2];
sx q[2];
rz(-1.1366476) q[2];
sx q[2];
rz(-2.1936552) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.0028210359) q[1];
sx q[1];
rz(-1.661794) q[1];
sx q[1];
rz(2.2343193) q[1];
x q[2];
rz(-1.7216191) q[3];
sx q[3];
rz(-1.5449459) q[3];
sx q[3];
rz(-2.6005656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1277348) q[2];
sx q[2];
rz(-1.3188136) q[2];
sx q[2];
rz(-2.0347118) q[2];
rz(-1.7539304) q[3];
sx q[3];
rz(-0.51968402) q[3];
sx q[3];
rz(-1.9096411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-2.7805507) q[0];
sx q[0];
rz(-1.1013958) q[0];
sx q[0];
rz(0.74044359) q[0];
rz(2.6904147) q[1];
sx q[1];
rz(-1.2606882) q[1];
sx q[1];
rz(-1.0528475) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2967865) q[0];
sx q[0];
rz(-2.1150981) q[0];
sx q[0];
rz(-1.143572) q[0];
rz(-pi) q[1];
rz(-2.244433) q[2];
sx q[2];
rz(-1.3348012) q[2];
sx q[2];
rz(1.2981851) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.14245089) q[1];
sx q[1];
rz(-1.4186727) q[1];
sx q[1];
rz(-1.692148) q[1];
x q[2];
rz(-1.5925519) q[3];
sx q[3];
rz(-0.43020159) q[3];
sx q[3];
rz(-2.3019626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25049245) q[2];
sx q[2];
rz(-1.6001469) q[2];
sx q[2];
rz(-0.54692522) q[2];
rz(0.28918239) q[3];
sx q[3];
rz(-0.5368036) q[3];
sx q[3];
rz(3.1291936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96994394) q[0];
sx q[0];
rz(-2.8778853) q[0];
sx q[0];
rz(-1.7893715) q[0];
rz(2.9317454) q[1];
sx q[1];
rz(-0.70979697) q[1];
sx q[1];
rz(0.23637493) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073233152) q[0];
sx q[0];
rz(-1.3993565) q[0];
sx q[0];
rz(-1.839848) q[0];
rz(-pi) q[1];
rz(-0.3503372) q[2];
sx q[2];
rz(-2.0603893) q[2];
sx q[2];
rz(1.9621153) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31503662) q[1];
sx q[1];
rz(-2.9737925) q[1];
sx q[1];
rz(1.5312974) q[1];
rz(-pi) q[2];
x q[2];
rz(1.719791) q[3];
sx q[3];
rz(-2.0553203) q[3];
sx q[3];
rz(-1.6001584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9005047) q[2];
sx q[2];
rz(-0.97854486) q[2];
sx q[2];
rz(-0.55348712) q[2];
rz(0.91529804) q[3];
sx q[3];
rz(-1.883029) q[3];
sx q[3];
rz(-2.4966911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47167641) q[0];
sx q[0];
rz(-1.3277418) q[0];
sx q[0];
rz(2.4374403) q[0];
rz(2.0856805) q[1];
sx q[1];
rz(-2.6864955) q[1];
sx q[1];
rz(0.59590894) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0674954) q[0];
sx q[0];
rz(-1.6620381) q[0];
sx q[0];
rz(1.3542961) q[0];
rz(1.709183) q[2];
sx q[2];
rz(-0.90183898) q[2];
sx q[2];
rz(1.578518) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0101498) q[1];
sx q[1];
rz(-1.3859268) q[1];
sx q[1];
rz(-2.7308982) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44645198) q[3];
sx q[3];
rz(-0.40735301) q[3];
sx q[3];
rz(0.27094597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6490877) q[2];
sx q[2];
rz(-1.1143782) q[2];
sx q[2];
rz(-2.5904783) q[2];
rz(2.9344432) q[3];
sx q[3];
rz(-1.8271577) q[3];
sx q[3];
rz(-1.6023887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.9838487) q[0];
sx q[0];
rz(-2.284323) q[0];
sx q[0];
rz(2.1898848) q[0];
rz(2.6668008) q[1];
sx q[1];
rz(-1.9665078) q[1];
sx q[1];
rz(0.29528433) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2518894) q[0];
sx q[0];
rz(-1.9648874) q[0];
sx q[0];
rz(0.15051145) q[0];
x q[1];
rz(-1.4163383) q[2];
sx q[2];
rz(-1.8784349) q[2];
sx q[2];
rz(0.49738202) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30027896) q[1];
sx q[1];
rz(-1.0611371) q[1];
sx q[1];
rz(2.120244) q[1];
x q[2];
rz(1.0601677) q[3];
sx q[3];
rz(-1.5139297) q[3];
sx q[3];
rz(2.1631654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7513912) q[2];
sx q[2];
rz(-1.7753121) q[2];
sx q[2];
rz(-0.93377101) q[2];
rz(-1.6823403) q[3];
sx q[3];
rz(-1.3542342) q[3];
sx q[3];
rz(1.4565844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.075832531) q[0];
sx q[0];
rz(-1.9969143) q[0];
sx q[0];
rz(0.24205762) q[0];
rz(2.4767955) q[1];
sx q[1];
rz(-1.2882065) q[1];
sx q[1];
rz(-2.738293) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2780669) q[0];
sx q[0];
rz(-2.6289872) q[0];
sx q[0];
rz(-0.82786064) q[0];
rz(1.2861916) q[2];
sx q[2];
rz(-0.91214123) q[2];
sx q[2];
rz(1.2496442) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.62337263) q[1];
sx q[1];
rz(-1.5689335) q[1];
sx q[1];
rz(-1.3606447) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1163252) q[3];
sx q[3];
rz(-1.7606887) q[3];
sx q[3];
rz(-2.5437298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91281259) q[2];
sx q[2];
rz(-2.0704724) q[2];
sx q[2];
rz(0.43441233) q[2];
rz(-2.1408391) q[3];
sx q[3];
rz(-0.36608168) q[3];
sx q[3];
rz(2.6312857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0077165724) q[0];
sx q[0];
rz(-0.0061329734) q[0];
sx q[0];
rz(-2.6469321) q[0];
rz(-1.6330632) q[1];
sx q[1];
rz(-1.6260908) q[1];
sx q[1];
rz(-0.60044926) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6143215) q[0];
sx q[0];
rz(-1.1724171) q[0];
sx q[0];
rz(0.50171731) q[0];
x q[1];
rz(-0.026293228) q[2];
sx q[2];
rz(-2.0956989) q[2];
sx q[2];
rz(2.3412447) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.820103) q[1];
sx q[1];
rz(-0.21807018) q[1];
sx q[1];
rz(2.5597508) q[1];
rz(-pi) q[2];
rz(-0.59372254) q[3];
sx q[3];
rz(-0.36168081) q[3];
sx q[3];
rz(-2.7152674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1373458) q[2];
sx q[2];
rz(-2.1477951) q[2];
sx q[2];
rz(0.11631575) q[2];
rz(2.7159193) q[3];
sx q[3];
rz(-0.95723546) q[3];
sx q[3];
rz(-11*pi/12) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15329926) q[0];
sx q[0];
rz(-2.9635552) q[0];
sx q[0];
rz(1.4784038) q[0];
rz(2.2019745) q[1];
sx q[1];
rz(-1.3213108) q[1];
sx q[1];
rz(0.41752648) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4672887) q[0];
sx q[0];
rz(-2.1018565) q[0];
sx q[0];
rz(-0.051961016) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3067901) q[2];
sx q[2];
rz(-1.4679113) q[2];
sx q[2];
rz(-1.7057982) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.49839834) q[1];
sx q[1];
rz(-2.8865951) q[1];
sx q[1];
rz(-2.4229089) q[1];
x q[2];
rz(1.9741251) q[3];
sx q[3];
rz(-2.6780193) q[3];
sx q[3];
rz(1.4232672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6614762) q[2];
sx q[2];
rz(-0.63642234) q[2];
sx q[2];
rz(2.647906) q[2];
rz(-2.4168329) q[3];
sx q[3];
rz(-1.9829491) q[3];
sx q[3];
rz(-0.31989583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17393728) q[0];
sx q[0];
rz(-2.4854361) q[0];
sx q[0];
rz(2.4560112) q[0];
rz(-0.29742345) q[1];
sx q[1];
rz(-0.23700266) q[1];
sx q[1];
rz(-1.1313653) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5269276) q[0];
sx q[0];
rz(-1.0173431) q[0];
sx q[0];
rz(1.9795896) q[0];
x q[1];
rz(2.4360043) q[2];
sx q[2];
rz(-1.0552647) q[2];
sx q[2];
rz(0.74598344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3494306) q[1];
sx q[1];
rz(-1.7948479) q[1];
sx q[1];
rz(0.38778023) q[1];
x q[2];
rz(2.0496619) q[3];
sx q[3];
rz(-2.1225956) q[3];
sx q[3];
rz(2.5294876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.82548213) q[2];
sx q[2];
rz(-1.9448514) q[2];
sx q[2];
rz(-0.70739174) q[2];
rz(-0.80983821) q[3];
sx q[3];
rz(-2.4707068) q[3];
sx q[3];
rz(-2.9530318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512882) q[0];
sx q[0];
rz(-1.123883) q[0];
sx q[0];
rz(-0.71258769) q[0];
rz(0.21223016) q[1];
sx q[1];
rz(-1.4490912) q[1];
sx q[1];
rz(2.6279411) q[1];
rz(2.1420508) q[2];
sx q[2];
rz(-2.0959601) q[2];
sx q[2];
rz(-3.0053896) q[2];
rz(-0.57701941) q[3];
sx q[3];
rz(-2.1885625) q[3];
sx q[3];
rz(2.267005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
