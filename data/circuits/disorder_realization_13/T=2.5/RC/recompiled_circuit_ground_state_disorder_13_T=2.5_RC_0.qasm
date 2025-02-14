OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.90091997) q[0];
sx q[0];
rz(-0.14578851) q[0];
sx q[0];
rz(-1.8426275) q[0];
rz(-0.19620148) q[1];
sx q[1];
rz(-1.3407522) q[1];
sx q[1];
rz(-0.10032108) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21151982) q[0];
sx q[0];
rz(-2.2317108) q[0];
sx q[0];
rz(-1.5050864) q[0];
rz(2.5415129) q[2];
sx q[2];
rz(-2.588495) q[2];
sx q[2];
rz(-1.2002522) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.92857526) q[1];
sx q[1];
rz(-1.2963821) q[1];
sx q[1];
rz(1.3540512) q[1];
x q[2];
rz(-2.5514929) q[3];
sx q[3];
rz(-1.3076772) q[3];
sx q[3];
rz(-0.89917574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4002865) q[2];
sx q[2];
rz(-2.9897959) q[2];
sx q[2];
rz(-2.3347704) q[2];
rz(-2.412879) q[3];
sx q[3];
rz(-0.7586793) q[3];
sx q[3];
rz(-1.2046643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62419409) q[0];
sx q[0];
rz(-2.875858) q[0];
sx q[0];
rz(-0.30701315) q[0];
rz(1.864805) q[1];
sx q[1];
rz(-1.1390319) q[1];
sx q[1];
rz(2.0397287) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80230882) q[0];
sx q[0];
rz(-1.7442987) q[0];
sx q[0];
rz(1.0109148) q[0];
rz(-pi) q[1];
rz(0.33874933) q[2];
sx q[2];
rz(-0.18998665) q[2];
sx q[2];
rz(-2.7471586) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.77473536) q[1];
sx q[1];
rz(-2.9178019) q[1];
sx q[1];
rz(-2.2246477) q[1];
x q[2];
rz(-0.83294932) q[3];
sx q[3];
rz(-2.2112339) q[3];
sx q[3];
rz(0.73329496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.028194204) q[2];
sx q[2];
rz(-2.5598309) q[2];
sx q[2];
rz(3.0451575) q[2];
rz(-2.9891369) q[3];
sx q[3];
rz(-1.6343445) q[3];
sx q[3];
rz(-0.69795394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.089712791) q[0];
sx q[0];
rz(-2.0413601) q[0];
sx q[0];
rz(-0.40019792) q[0];
rz(-1.5125037) q[1];
sx q[1];
rz(-0.17833231) q[1];
sx q[1];
rz(0.2581183) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5982957) q[0];
sx q[0];
rz(-0.27320293) q[0];
sx q[0];
rz(-2.1092723) q[0];
rz(-pi) q[1];
rz(1.4908954) q[2];
sx q[2];
rz(-1.4712442) q[2];
sx q[2];
rz(0.71074206) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5096881) q[1];
sx q[1];
rz(-3.1277165) q[1];
sx q[1];
rz(-2.7751776) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7430923) q[3];
sx q[3];
rz(-1.8086284) q[3];
sx q[3];
rz(-2.3034277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.021598024) q[2];
sx q[2];
rz(-1.9436516) q[2];
sx q[2];
rz(-3.1094587) q[2];
rz(-0.26767996) q[3];
sx q[3];
rz(-1.6240424) q[3];
sx q[3];
rz(1.5599498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5699919) q[0];
sx q[0];
rz(-1.5084234) q[0];
sx q[0];
rz(-0.084240325) q[0];
rz(3.1056504) q[1];
sx q[1];
rz(-0.033999559) q[1];
sx q[1];
rz(2.8004004) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55732176) q[0];
sx q[0];
rz(-2.1958618) q[0];
sx q[0];
rz(-1.8950426) q[0];
rz(-pi) q[1];
rz(2.7643599) q[2];
sx q[2];
rz(-2.5209941) q[2];
sx q[2];
rz(2.008174) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.045121) q[1];
sx q[1];
rz(-2.0092138) q[1];
sx q[1];
rz(-2.3271266) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.59935948) q[3];
sx q[3];
rz(-1.9512358) q[3];
sx q[3];
rz(-2.3886556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42698947) q[2];
sx q[2];
rz(-1.0726856) q[2];
sx q[2];
rz(-1.6061456) q[2];
rz(-0.26220194) q[3];
sx q[3];
rz(-1.6180792) q[3];
sx q[3];
rz(-0.95482993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75761211) q[0];
sx q[0];
rz(-2.7181427) q[0];
sx q[0];
rz(-2.0641932) q[0];
rz(2.7491838) q[1];
sx q[1];
rz(-0.078374021) q[1];
sx q[1];
rz(-2.1108625) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2868256) q[0];
sx q[0];
rz(-2.1041901) q[0];
sx q[0];
rz(2.0598165) q[0];
x q[1];
rz(2.9610277) q[2];
sx q[2];
rz(-1.079165) q[2];
sx q[2];
rz(-0.59300834) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3001316) q[1];
sx q[1];
rz(-2.8634954) q[1];
sx q[1];
rz(-2.2113731) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0038336) q[3];
sx q[3];
rz(-0.58708411) q[3];
sx q[3];
rz(2.7239012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3173759) q[2];
sx q[2];
rz(-2.4874096) q[2];
sx q[2];
rz(-0.83924323) q[2];
rz(-0.9134891) q[3];
sx q[3];
rz(-1.8279816) q[3];
sx q[3];
rz(2.998108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4588673) q[0];
sx q[0];
rz(-0.23696466) q[0];
sx q[0];
rz(-1.6656026) q[0];
rz(0.39235517) q[1];
sx q[1];
rz(-2.0456435) q[1];
sx q[1];
rz(2.5501693) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77350835) q[0];
sx q[0];
rz(-2.4249833) q[0];
sx q[0];
rz(-0.98978292) q[0];
rz(-1.3900969) q[2];
sx q[2];
rz(-1.2998253) q[2];
sx q[2];
rz(-2.9621552) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6900031) q[1];
sx q[1];
rz(-1.065548) q[1];
sx q[1];
rz(3.0662159) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0143842) q[3];
sx q[3];
rz(-0.74723703) q[3];
sx q[3];
rz(-2.0153449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.300294) q[2];
sx q[2];
rz(-0.66606194) q[2];
sx q[2];
rz(-2.5692614) q[2];
rz(-2.9193997) q[3];
sx q[3];
rz(-2.7100345) q[3];
sx q[3];
rz(0.64479327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38729024) q[0];
sx q[0];
rz(-3.001725) q[0];
sx q[0];
rz(-2.733316) q[0];
rz(-0.73221842) q[1];
sx q[1];
rz(-0.1258985) q[1];
sx q[1];
rz(0.29762038) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0669516) q[0];
sx q[0];
rz(-0.56109259) q[0];
sx q[0];
rz(-2.448161) q[0];
rz(3.1213644) q[2];
sx q[2];
rz(-1.0761522) q[2];
sx q[2];
rz(-1.9488364) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0113161) q[1];
sx q[1];
rz(-1.374048) q[1];
sx q[1];
rz(-1.2230243) q[1];
x q[2];
rz(-0.17642658) q[3];
sx q[3];
rz(-0.71143141) q[3];
sx q[3];
rz(-2.3030858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1709661) q[2];
sx q[2];
rz(-1.8793224) q[2];
sx q[2];
rz(2.4453898) q[2];
rz(2.2512186) q[3];
sx q[3];
rz(-1.9647157) q[3];
sx q[3];
rz(1.6740929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17906976) q[0];
sx q[0];
rz(-3.1130377) q[0];
sx q[0];
rz(2.9288375) q[0];
rz(0.46956024) q[1];
sx q[1];
rz(-0.9649562) q[1];
sx q[1];
rz(0.75417095) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5306204) q[0];
sx q[0];
rz(-1.7585131) q[0];
sx q[0];
rz(-1.8646445) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55108549) q[2];
sx q[2];
rz(-2.7646825) q[2];
sx q[2];
rz(2.2152679) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8629145) q[1];
sx q[1];
rz(-2.2154741) q[1];
sx q[1];
rz(1.7134604) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.76393083) q[3];
sx q[3];
rz(-2.3124933) q[3];
sx q[3];
rz(2.130485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.135005) q[2];
sx q[2];
rz(-0.90234119) q[2];
sx q[2];
rz(2.3593486) q[2];
rz(1.6953281) q[3];
sx q[3];
rz(-0.54636121) q[3];
sx q[3];
rz(2.7915891) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0924454) q[0];
sx q[0];
rz(-0.47098422) q[0];
sx q[0];
rz(-2.1771722) q[0];
rz(-1.2760705) q[1];
sx q[1];
rz(-1.4141021) q[1];
sx q[1];
rz(-1.5020348) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0014711) q[0];
sx q[0];
rz(-1.7709416) q[0];
sx q[0];
rz(-1.3339304) q[0];
rz(-pi) q[1];
rz(-1.1310857) q[2];
sx q[2];
rz(-0.88958101) q[2];
sx q[2];
rz(-2.0063673) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6148541) q[1];
sx q[1];
rz(-2.6236218) q[1];
sx q[1];
rz(1.5244085) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.718091) q[3];
sx q[3];
rz(-1.637666) q[3];
sx q[3];
rz(-1.221958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8883349) q[2];
sx q[2];
rz(-1.8586681) q[2];
sx q[2];
rz(-0.9453195) q[2];
rz(2.3156598) q[3];
sx q[3];
rz(-1.632894) q[3];
sx q[3];
rz(0.53502423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.076040529) q[0];
sx q[0];
rz(-1.9726418) q[0];
sx q[0];
rz(0.8031351) q[0];
rz(1.5777292) q[1];
sx q[1];
rz(-1.6604661) q[1];
sx q[1];
rz(-2.8520083) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6824376) q[0];
sx q[0];
rz(-1.4652325) q[0];
sx q[0];
rz(-3.1292874) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8044462) q[2];
sx q[2];
rz(-0.44155332) q[2];
sx q[2];
rz(1.785977) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3037422) q[1];
sx q[1];
rz(-1.9281862) q[1];
sx q[1];
rz(-0.2405432) q[1];
x q[2];
rz(1.2373459) q[3];
sx q[3];
rz(-1.8419918) q[3];
sx q[3];
rz(-1.1734133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3760066) q[2];
sx q[2];
rz(-3.0209318) q[2];
sx q[2];
rz(2.181459) q[2];
rz(2.5586186) q[3];
sx q[3];
rz(-0.65810242) q[3];
sx q[3];
rz(-0.8647024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4509907) q[0];
sx q[0];
rz(-1.7091746) q[0];
sx q[0];
rz(-1.5102392) q[0];
rz(0.040738978) q[1];
sx q[1];
rz(-0.67650411) q[1];
sx q[1];
rz(0.13112851) q[1];
rz(-2.1292674) q[2];
sx q[2];
rz(-0.53999117) q[2];
sx q[2];
rz(-2.7593031) q[2];
rz(-1.4996281) q[3];
sx q[3];
rz(-1.5200079) q[3];
sx q[3];
rz(-0.11848371) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
