OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.83710837) q[0];
sx q[0];
rz(4.8298782) q[0];
sx q[0];
rz(9.7363135) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(4.9649927) q[1];
sx q[1];
rz(8.8658219) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11249609) q[0];
sx q[0];
rz(-1.7100428) q[0];
sx q[0];
rz(1.1514366) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5842701) q[2];
sx q[2];
rz(-1.5601336) q[2];
sx q[2];
rz(0.23920857) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.47145876) q[1];
sx q[1];
rz(-1.0308627) q[1];
sx q[1];
rz(2.259841) q[1];
rz(-pi) q[2];
rz(0.43879126) q[3];
sx q[3];
rz(-1.8306797) q[3];
sx q[3];
rz(-2.4364542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3922334) q[2];
sx q[2];
rz(-1.8584676) q[2];
sx q[2];
rz(-0.63670811) q[2];
rz(0.84896815) q[3];
sx q[3];
rz(-2.520112) q[3];
sx q[3];
rz(2.9076715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7648776) q[0];
sx q[0];
rz(-2.8945518) q[0];
sx q[0];
rz(-0.15287457) q[0];
rz(0.75694594) q[1];
sx q[1];
rz(-1.5544954) q[1];
sx q[1];
rz(0.98639948) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61030412) q[0];
sx q[0];
rz(-0.072040759) q[0];
sx q[0];
rz(-2.5487367) q[0];
rz(-pi) q[1];
rz(-1.0436922) q[2];
sx q[2];
rz(-0.88000789) q[2];
sx q[2];
rz(-1.8156798) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1101148) q[1];
sx q[1];
rz(-1.8077393) q[1];
sx q[1];
rz(-2.7538339) q[1];
rz(-pi) q[2];
rz(-1.1902477) q[3];
sx q[3];
rz(-2.7893587) q[3];
sx q[3];
rz(0.16570839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5793005) q[2];
sx q[2];
rz(-1.9220756) q[2];
sx q[2];
rz(-0.78312773) q[2];
rz(0.018571818) q[3];
sx q[3];
rz(-1.637807) q[3];
sx q[3];
rz(-0.40772453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0884393) q[0];
sx q[0];
rz(-0.29215559) q[0];
sx q[0];
rz(0.96167481) q[0];
rz(0.36034521) q[1];
sx q[1];
rz(-1.1018437) q[1];
sx q[1];
rz(-0.12869421) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1291217) q[0];
sx q[0];
rz(-3.0462628) q[0];
sx q[0];
rz(2.096813) q[0];
rz(-pi) q[1];
rz(-2.6016597) q[2];
sx q[2];
rz(-2.2579102) q[2];
sx q[2];
rz(2.5643258) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3330831) q[1];
sx q[1];
rz(-1.2915478) q[1];
sx q[1];
rz(-0.46344325) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4314753) q[3];
sx q[3];
rz(-2.3686667) q[3];
sx q[3];
rz(-1.7804002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.06015691) q[2];
sx q[2];
rz(-1.9443941) q[2];
sx q[2];
rz(1.241768) q[2];
rz(2.5545819) q[3];
sx q[3];
rz(-0.95654064) q[3];
sx q[3];
rz(2.164042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.509165) q[0];
sx q[0];
rz(-0.88212633) q[0];
sx q[0];
rz(1.084491) q[0];
rz(1.4831316) q[1];
sx q[1];
rz(-0.56743923) q[1];
sx q[1];
rz(-0.09253563) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4400892) q[0];
sx q[0];
rz(-1.8532231) q[0];
sx q[0];
rz(-0.28042067) q[0];
rz(-pi) q[1];
rz(2.0501577) q[2];
sx q[2];
rz(-1.2617246) q[2];
sx q[2];
rz(1.1190337) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9738237) q[1];
sx q[1];
rz(-1.2281706) q[1];
sx q[1];
rz(1.509036) q[1];
x q[2];
rz(0.91589215) q[3];
sx q[3];
rz(-1.4426784) q[3];
sx q[3];
rz(1.7803943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3080421) q[2];
sx q[2];
rz(-1.4414859) q[2];
sx q[2];
rz(-2.8095424) q[2];
rz(2.0856693) q[3];
sx q[3];
rz(-2.8639586) q[3];
sx q[3];
rz(2.5312996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(-2.8191391) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(-2.9300368) q[0];
rz(1.8353204) q[1];
sx q[1];
rz(-1.897656) q[1];
sx q[1];
rz(-0.64770118) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54007285) q[0];
sx q[0];
rz(-1.411502) q[0];
sx q[0];
rz(-1.457731) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7366016) q[2];
sx q[2];
rz(-2.7654123) q[2];
sx q[2];
rz(-2.1567547) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.60407818) q[1];
sx q[1];
rz(-2.0734097) q[1];
sx q[1];
rz(-2.9167049) q[1];
rz(-pi) q[2];
rz(1.5557489) q[3];
sx q[3];
rz(-1.589937) q[3];
sx q[3];
rz(-1.0552989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.56090474) q[2];
sx q[2];
rz(-2.7320392) q[2];
sx q[2];
rz(0.69331759) q[2];
rz(2.4723315) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(-0.0049237331) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3174021) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(-1.90907) q[0];
rz(-1.0725853) q[1];
sx q[1];
rz(-2.0697846) q[1];
sx q[1];
rz(-2.9673064) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5151383) q[0];
sx q[0];
rz(-0.8695375) q[0];
sx q[0];
rz(0.40909543) q[0];
rz(-0.66138791) q[2];
sx q[2];
rz(-0.97860133) q[2];
sx q[2];
rz(-2.4370898) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4982521) q[1];
sx q[1];
rz(-2.2689515) q[1];
sx q[1];
rz(2.6367285) q[1];
x q[2];
rz(1.2586081) q[3];
sx q[3];
rz(-1.6815261) q[3];
sx q[3];
rz(-0.98715106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3198513) q[2];
sx q[2];
rz(-2.3345626) q[2];
sx q[2];
rz(2.9439587) q[2];
rz(-2.8526784) q[3];
sx q[3];
rz(-0.87696004) q[3];
sx q[3];
rz(1.7355841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0777247) q[0];
sx q[0];
rz(-0.48982319) q[0];
sx q[0];
rz(0.20859627) q[0];
rz(-0.96616191) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(1.5055515) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9150328) q[0];
sx q[0];
rz(-1.0462927) q[0];
sx q[0];
rz(0.71136186) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43039544) q[2];
sx q[2];
rz(-1.5376523) q[2];
sx q[2];
rz(2.7660649) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4475496) q[1];
sx q[1];
rz(-2.4195237) q[1];
sx q[1];
rz(2.2896647) q[1];
x q[2];
rz(-0.286245) q[3];
sx q[3];
rz(-1.6702594) q[3];
sx q[3];
rz(3.0486097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2839526) q[2];
sx q[2];
rz(-0.74308926) q[2];
sx q[2];
rz(3.0440142) q[2];
rz(-1.7476667) q[3];
sx q[3];
rz(-1.3503617) q[3];
sx q[3];
rz(0.71684366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7512648) q[0];
sx q[0];
rz(-1.8126235) q[0];
sx q[0];
rz(2.6275997) q[0];
rz(3.0184074) q[1];
sx q[1];
rz(-0.24736483) q[1];
sx q[1];
rz(2.2095912) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1876353) q[0];
sx q[0];
rz(-1.7062511) q[0];
sx q[0];
rz(-2.016469) q[0];
rz(-pi) q[1];
rz(0.26884218) q[2];
sx q[2];
rz(-2.1890867) q[2];
sx q[2];
rz(0.48228574) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8646647) q[1];
sx q[1];
rz(-2.118646) q[1];
sx q[1];
rz(-2.1879556) q[1];
rz(-1.7979513) q[3];
sx q[3];
rz(-2.4121768) q[3];
sx q[3];
rz(1.2136572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0294068) q[2];
sx q[2];
rz(-1.1108578) q[2];
sx q[2];
rz(-0.4294447) q[2];
rz(-1.9321692) q[3];
sx q[3];
rz(-0.37241396) q[3];
sx q[3];
rz(2.4485574) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76686239) q[0];
sx q[0];
rz(-1.6748036) q[0];
sx q[0];
rz(1.3388348) q[0];
rz(2.4354637) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(0.95058092) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6317752) q[0];
sx q[0];
rz(-2.0777367) q[0];
sx q[0];
rz(2.1419873) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0663475) q[2];
sx q[2];
rz(-1.9855472) q[2];
sx q[2];
rz(0.097397734) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.49950019) q[1];
sx q[1];
rz(-1.7576808) q[1];
sx q[1];
rz(-1.0287813) q[1];
rz(-pi) q[2];
rz(-2.1080984) q[3];
sx q[3];
rz(-2.6712382) q[3];
sx q[3];
rz(-2.5126517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6616228) q[2];
sx q[2];
rz(-1.491549) q[2];
sx q[2];
rz(2.6573112) q[2];
rz(-0.92710036) q[3];
sx q[3];
rz(-1.3092594) q[3];
sx q[3];
rz(2.5203729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1442239) q[0];
sx q[0];
rz(-3.0529418) q[0];
sx q[0];
rz(0.22928672) q[0];
rz(-0.43481049) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(-2.4226709) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15764788) q[0];
sx q[0];
rz(-0.63417182) q[0];
sx q[0];
rz(1.4112524) q[0];
x q[1];
rz(-0.11833338) q[2];
sx q[2];
rz(-2.1552857) q[2];
sx q[2];
rz(3.0362533) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.016991888) q[1];
sx q[1];
rz(-1.4233839) q[1];
sx q[1];
rz(2.7684545) q[1];
rz(-pi) q[2];
rz(2.2273916) q[3];
sx q[3];
rz(-0.73540348) q[3];
sx q[3];
rz(-0.38287336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7583313) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(-0.74679217) q[2];
rz(-2.2693999) q[3];
sx q[3];
rz(-2.3132497) q[3];
sx q[3];
rz(-0.94223589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9185716) q[0];
sx q[0];
rz(-1.3970319) q[0];
sx q[0];
rz(-1.3402517) q[0];
rz(0.37721286) q[1];
sx q[1];
rz(-1.6709534) q[1];
sx q[1];
rz(0.51660641) q[1];
rz(-0.51858356) q[2];
sx q[2];
rz(-1.2972144) q[2];
sx q[2];
rz(-1.5757061) q[2];
rz(1.4217581) q[3];
sx q[3];
rz(-0.85902135) q[3];
sx q[3];
rz(-1.600941) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
