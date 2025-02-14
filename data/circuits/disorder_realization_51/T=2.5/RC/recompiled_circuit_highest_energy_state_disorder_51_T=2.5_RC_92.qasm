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
rz(2.3866374) q[0];
sx q[0];
rz(-0.75140262) q[0];
sx q[0];
rz(2.0928535) q[0];
rz(-2.7222848) q[1];
sx q[1];
rz(-1.7555305) q[1];
sx q[1];
rz(3.0233033) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2631419) q[0];
sx q[0];
rz(-1.0174482) q[0];
sx q[0];
rz(-0.31991495) q[0];
rz(3.0906244) q[2];
sx q[2];
rz(-1.5870924) q[2];
sx q[2];
rz(0.38557926) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6403049) q[1];
sx q[1];
rz(-3.1073862) q[1];
sx q[1];
rz(-0.51987363) q[1];
rz(1.6476326) q[3];
sx q[3];
rz(-1.8733403) q[3];
sx q[3];
rz(1.3148112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2179541) q[2];
sx q[2];
rz(-0.0035692735) q[2];
sx q[2];
rz(2.9812319) q[2];
rz(1.977836) q[3];
sx q[3];
rz(-1.0842423) q[3];
sx q[3];
rz(2.3273996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9775951) q[0];
sx q[0];
rz(-1.6544592) q[0];
sx q[0];
rz(-0.98597041) q[0];
rz(-1.5892971) q[1];
sx q[1];
rz(-0.28737107) q[1];
sx q[1];
rz(-1.5602962) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47871209) q[0];
sx q[0];
rz(-1.9528188) q[0];
sx q[0];
rz(2.0860614) q[0];
x q[1];
rz(-1.549717) q[2];
sx q[2];
rz(-0.97641845) q[2];
sx q[2];
rz(-3.0836058) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6735005) q[1];
sx q[1];
rz(-1.4080268) q[1];
sx q[1];
rz(-0.36483187) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93771622) q[3];
sx q[3];
rz(-1.670853) q[3];
sx q[3];
rz(-1.2591187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.61957773) q[2];
sx q[2];
rz(-1.8176983) q[2];
sx q[2];
rz(2.1066693) q[2];
rz(2.6400635) q[3];
sx q[3];
rz(-3.0540255) q[3];
sx q[3];
rz(-0.97626221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4185249) q[0];
sx q[0];
rz(-1.9880966) q[0];
sx q[0];
rz(-1.204741) q[0];
rz(1.0459666) q[1];
sx q[1];
rz(-0.090066411) q[1];
sx q[1];
rz(0.14030309) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2818272) q[0];
sx q[0];
rz(-2.4778296) q[0];
sx q[0];
rz(-1.1361994) q[0];
rz(-pi) q[1];
rz(-1.1029805) q[2];
sx q[2];
rz(-0.78323871) q[2];
sx q[2];
rz(-0.65708438) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.86911234) q[1];
sx q[1];
rz(-1.7971874) q[1];
sx q[1];
rz(-2.5416944) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14032538) q[3];
sx q[3];
rz(-1.9092536) q[3];
sx q[3];
rz(-2.942977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3673765) q[2];
sx q[2];
rz(-2.1420631) q[2];
sx q[2];
rz(0.2224758) q[2];
rz(3.0761278) q[3];
sx q[3];
rz(-1.2832578) q[3];
sx q[3];
rz(-2.2953667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9533933) q[0];
sx q[0];
rz(-3.1173752) q[0];
sx q[0];
rz(-0.54704332) q[0];
rz(0.25829265) q[1];
sx q[1];
rz(-3.1196085) q[1];
sx q[1];
rz(2.8000854) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19778457) q[0];
sx q[0];
rz(-0.0066692624) q[0];
sx q[0];
rz(-1.3640553) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33880587) q[2];
sx q[2];
rz(-1.9680541) q[2];
sx q[2];
rz(-2.4851929) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5568513) q[1];
sx q[1];
rz(-1.9042174) q[1];
sx q[1];
rz(0.78704301) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82587256) q[3];
sx q[3];
rz(-2.4950728) q[3];
sx q[3];
rz(-2.7310128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.42939886) q[2];
sx q[2];
rz(-1.2960351) q[2];
sx q[2];
rz(2.1923501) q[2];
rz(-2.3927169) q[3];
sx q[3];
rz(-1.8604934) q[3];
sx q[3];
rz(-0.062189814) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6147989) q[0];
sx q[0];
rz(-0.034788046) q[0];
sx q[0];
rz(-1.590796) q[0];
rz(-1.3560449) q[1];
sx q[1];
rz(-0.0043914774) q[1];
sx q[1];
rz(-0.063025085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3857128) q[0];
sx q[0];
rz(-1.5405419) q[0];
sx q[0];
rz(0.066019375) q[0];
rz(-2.3895415) q[2];
sx q[2];
rz(-1.1630485) q[2];
sx q[2];
rz(-2.4501981) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5668727) q[1];
sx q[1];
rz(-1.5441455) q[1];
sx q[1];
rz(1.3625366) q[1];
rz(2.359777) q[3];
sx q[3];
rz(-1.2226579) q[3];
sx q[3];
rz(-2.0962417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.65681347) q[2];
sx q[2];
rz(-1.3098837) q[2];
sx q[2];
rz(0.64378929) q[2];
rz(2.2667609) q[3];
sx q[3];
rz(-0.30431408) q[3];
sx q[3];
rz(-2.3410102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0952045) q[0];
sx q[0];
rz(-3.0859741) q[0];
sx q[0];
rz(-0.355542) q[0];
rz(-2.9429759) q[1];
sx q[1];
rz(-0.0067409975) q[1];
sx q[1];
rz(-0.14828646) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9614811) q[0];
sx q[0];
rz(-1.5732871) q[0];
sx q[0];
rz(0.18302304) q[0];
rz(-1.6464433) q[2];
sx q[2];
rz(-1.160607) q[2];
sx q[2];
rz(0.73034053) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7050138) q[1];
sx q[1];
rz(-1.8039898) q[1];
sx q[1];
rz(1.4326653) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1790696) q[3];
sx q[3];
rz(-1.4279162) q[3];
sx q[3];
rz(-0.30331984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4699012) q[2];
sx q[2];
rz(-0.24158676) q[2];
sx q[2];
rz(0.12413231) q[2];
rz(0.5747253) q[3];
sx q[3];
rz(-2.997213) q[3];
sx q[3];
rz(0.1709443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8827051) q[0];
sx q[0];
rz(-0.12508617) q[0];
sx q[0];
rz(2.4001154) q[0];
rz(2.8575836) q[1];
sx q[1];
rz(-3.1378855) q[1];
sx q[1];
rz(-2.8264118) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18703546) q[0];
sx q[0];
rz(-1.5439568) q[0];
sx q[0];
rz(-1.5052273) q[0];
rz(-pi) q[1];
rz(-1.9839331) q[2];
sx q[2];
rz(-2.6584932) q[2];
sx q[2];
rz(-0.60483067) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8410346) q[1];
sx q[1];
rz(-1.729064) q[1];
sx q[1];
rz(0.59945089) q[1];
rz(1.3471425) q[3];
sx q[3];
rz(-1.5937433) q[3];
sx q[3];
rz(-0.060088559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.86889851) q[2];
sx q[2];
rz(-2.0468476) q[2];
sx q[2];
rz(-0.72186738) q[2];
rz(2.7591211) q[3];
sx q[3];
rz(-1.9964652) q[3];
sx q[3];
rz(-1.0684048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.566399) q[0];
sx q[0];
rz(-0.02481758) q[0];
sx q[0];
rz(1.5750634) q[0];
rz(-0.20340915) q[1];
sx q[1];
rz(-1.8433488) q[1];
sx q[1];
rz(0.64483109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0087850182) q[0];
sx q[0];
rz(-1.5624579) q[0];
sx q[0];
rz(1.7783283) q[0];
rz(1.5956085) q[2];
sx q[2];
rz(-0.96867079) q[2];
sx q[2];
rz(0.40222049) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.37721021) q[1];
sx q[1];
rz(-1.4971716) q[1];
sx q[1];
rz(-2.1097357) q[1];
rz(-pi) q[2];
rz(2.0243778) q[3];
sx q[3];
rz(-2.0189813) q[3];
sx q[3];
rz(0.77806015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5440392) q[2];
sx q[2];
rz(-2.7880703) q[2];
sx q[2];
rz(0.37975797) q[2];
rz(2.0960268) q[3];
sx q[3];
rz(-1.2328204) q[3];
sx q[3];
rz(-1.9426965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7757292) q[0];
sx q[0];
rz(-3.1079223) q[0];
sx q[0];
rz(-1.3566383) q[0];
rz(2.7011073) q[1];
sx q[1];
rz(-2.0511274) q[1];
sx q[1];
rz(-2.4408565) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7788575) q[0];
sx q[0];
rz(-2.252251) q[0];
sx q[0];
rz(-2.2454041) q[0];
rz(0.88057554) q[2];
sx q[2];
rz(-1.4172366) q[2];
sx q[2];
rz(2.409006) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8776013) q[1];
sx q[1];
rz(-0.7670247) q[1];
sx q[1];
rz(-0.33543889) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0088720284) q[3];
sx q[3];
rz(-0.43264593) q[3];
sx q[3];
rz(-0.052730058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35357722) q[2];
sx q[2];
rz(-0.37166301) q[2];
sx q[2];
rz(1.2865944) q[2];
rz(-2.6364117) q[3];
sx q[3];
rz(-0.44613871) q[3];
sx q[3];
rz(1.2222458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6219567) q[0];
sx q[0];
rz(-0.049787909) q[0];
sx q[0];
rz(-1.5420445) q[0];
rz(-2.385335) q[1];
sx q[1];
rz(-3.1344963) q[1];
sx q[1];
rz(-2.8047628) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8076421) q[0];
sx q[0];
rz(-1.5686638) q[0];
sx q[0];
rz(2.8239408) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8518108) q[2];
sx q[2];
rz(-1.769067) q[2];
sx q[2];
rz(-2.2137583) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5975128) q[1];
sx q[1];
rz(-1.4855396) q[1];
sx q[1];
rz(1.4592749) q[1];
rz(-pi) q[2];
rz(0.41542094) q[3];
sx q[3];
rz(-1.0491228) q[3];
sx q[3];
rz(0.085062438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.51356641) q[2];
sx q[2];
rz(-2.1921373) q[2];
sx q[2];
rz(0.27964082) q[2];
rz(2.5102992) q[3];
sx q[3];
rz(-2.203233) q[3];
sx q[3];
rz(0.62425557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8414128) q[0];
sx q[0];
rz(-1.5914088) q[0];
sx q[0];
rz(1.7803022) q[0];
rz(0.77371669) q[1];
sx q[1];
rz(-0.63540375) q[1];
sx q[1];
rz(0.2159963) q[1];
rz(2.2778289) q[2];
sx q[2];
rz(-0.84767271) q[2];
sx q[2];
rz(1.5172971) q[2];
rz(1.5935043) q[3];
sx q[3];
rz(-1.1946214) q[3];
sx q[3];
rz(-0.0047636845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
