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
rz(0.1641195) q[0];
sx q[0];
rz(-2.1174705) q[0];
sx q[0];
rz(0.48409387) q[0];
rz(2.3789499) q[1];
sx q[1];
rz(4.7197309) q[1];
sx q[1];
rz(9.3698256) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8380175) q[0];
sx q[0];
rz(-2.2406881) q[0];
sx q[0];
rz(-0.039867) q[0];
rz(1.7968494) q[2];
sx q[2];
rz(-0.25616562) q[2];
sx q[2];
rz(-2.0992172) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.37147) q[1];
sx q[1];
rz(-2.9767163) q[1];
sx q[1];
rz(-0.6368963) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0255768) q[3];
sx q[3];
rz(-2.4686128) q[3];
sx q[3];
rz(-2.2448886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.83076465) q[2];
sx q[2];
rz(-2.1711633) q[2];
sx q[2];
rz(-0.27466276) q[2];
rz(-2.2667387) q[3];
sx q[3];
rz(-1.091205) q[3];
sx q[3];
rz(0.27354512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65983588) q[0];
sx q[0];
rz(-2.4517224) q[0];
sx q[0];
rz(2.7428395) q[0];
rz(-0.2440456) q[1];
sx q[1];
rz(-1.9052541) q[1];
sx q[1];
rz(-2.0358548) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.669848) q[0];
sx q[0];
rz(-2.6378184) q[0];
sx q[0];
rz(0.48724799) q[0];
rz(2.5697487) q[2];
sx q[2];
rz(-2.1430121) q[2];
sx q[2];
rz(-2.4586611) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2011021) q[1];
sx q[1];
rz(-0.0041714287) q[1];
sx q[1];
rz(0.68912403) q[1];
rz(-1.1532743) q[3];
sx q[3];
rz(-0.29051775) q[3];
sx q[3];
rz(2.6444496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4411321) q[2];
sx q[2];
rz(-1.1819785) q[2];
sx q[2];
rz(1.0750809) q[2];
rz(-0.69617802) q[3];
sx q[3];
rz(-1.2735561) q[3];
sx q[3];
rz(-2.9785494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27580801) q[0];
sx q[0];
rz(-1.2514665) q[0];
sx q[0];
rz(2.9248917) q[0];
rz(1.3801344) q[1];
sx q[1];
rz(-2.7237027) q[1];
sx q[1];
rz(-0.61947852) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19720896) q[0];
sx q[0];
rz(-1.608662) q[0];
sx q[0];
rz(-3.1119431) q[0];
rz(-0.046929788) q[2];
sx q[2];
rz(-1.4478193) q[2];
sx q[2];
rz(2.8356981) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.8113831) q[1];
sx q[1];
rz(-0.37322497) q[1];
sx q[1];
rz(0.95287816) q[1];
rz(-pi) q[2];
rz(-0.9914753) q[3];
sx q[3];
rz(-2.1001171) q[3];
sx q[3];
rz(-2.4059699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8850024) q[2];
sx q[2];
rz(-1.2991178) q[2];
sx q[2];
rz(3.0510862) q[2];
rz(0.45804405) q[3];
sx q[3];
rz(-0.48726714) q[3];
sx q[3];
rz(1.1080144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3312382) q[0];
sx q[0];
rz(-2.2585223) q[0];
sx q[0];
rz(-2.2099387) q[0];
rz(2.9725507) q[1];
sx q[1];
rz(-2.5773498) q[1];
sx q[1];
rz(-2.1021252) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4783206) q[0];
sx q[0];
rz(-2.2815858) q[0];
sx q[0];
rz(1.8904314) q[0];
rz(-0.2257963) q[2];
sx q[2];
rz(-1.8183299) q[2];
sx q[2];
rz(-0.087866656) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66450602) q[1];
sx q[1];
rz(-1.5110104) q[1];
sx q[1];
rz(-1.3633481) q[1];
rz(-pi) q[2];
x q[2];
rz(0.28431953) q[3];
sx q[3];
rz(-1.7756724) q[3];
sx q[3];
rz(-2.7903872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.52183759) q[2];
sx q[2];
rz(-1.1014742) q[2];
sx q[2];
rz(2.6452765) q[2];
rz(-2.3450092) q[3];
sx q[3];
rz(-2.8556672) q[3];
sx q[3];
rz(2.0485785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26688823) q[0];
sx q[0];
rz(-0.58458352) q[0];
sx q[0];
rz(-0.97187483) q[0];
rz(-3.019849) q[1];
sx q[1];
rz(-1.5309445) q[1];
sx q[1];
rz(1.8121388) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9864101) q[0];
sx q[0];
rz(-1.8768132) q[0];
sx q[0];
rz(3.0772692) q[0];
rz(2.2586063) q[2];
sx q[2];
rz(-2.5341883) q[2];
sx q[2];
rz(0.6907874) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54534528) q[1];
sx q[1];
rz(-1.9831428) q[1];
sx q[1];
rz(-2.418787) q[1];
x q[2];
rz(1.4917144) q[3];
sx q[3];
rz(-1.5800161) q[3];
sx q[3];
rz(3.0371333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.10151265) q[2];
sx q[2];
rz(-1.2206581) q[2];
sx q[2];
rz(-0.15138781) q[2];
rz(0.76465145) q[3];
sx q[3];
rz(-0.83573666) q[3];
sx q[3];
rz(-1.5966655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0720035) q[0];
sx q[0];
rz(-1.4729426) q[0];
sx q[0];
rz(-1.0714916) q[0];
rz(-2.6103861) q[1];
sx q[1];
rz(-1.6516282) q[1];
sx q[1];
rz(0.62613097) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1086639) q[0];
sx q[0];
rz(-1.5538006) q[0];
sx q[0];
rz(0.0047608382) q[0];
x q[1];
rz(-0.59711908) q[2];
sx q[2];
rz(-1.7749987) q[2];
sx q[2];
rz(2.0300421) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5007361) q[1];
sx q[1];
rz(-1.2656414) q[1];
sx q[1];
rz(-0.63668294) q[1];
rz(-pi) q[2];
rz(2.4247384) q[3];
sx q[3];
rz(-0.62590137) q[3];
sx q[3];
rz(-2.7519873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7834187) q[2];
sx q[2];
rz(-0.61966115) q[2];
sx q[2];
rz(2.1273071) q[2];
rz(-1.2554393) q[3];
sx q[3];
rz(-1.7858601) q[3];
sx q[3];
rz(2.0881418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96100539) q[0];
sx q[0];
rz(-2.2659681) q[0];
sx q[0];
rz(-2.2210806) q[0];
rz(-3.0275184) q[1];
sx q[1];
rz(-1.4055077) q[1];
sx q[1];
rz(-1.1309518) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75926542) q[0];
sx q[0];
rz(-2.5438519) q[0];
sx q[0];
rz(-0.010058479) q[0];
x q[1];
rz(-0.80318309) q[2];
sx q[2];
rz(-0.72942643) q[2];
sx q[2];
rz(-0.27952295) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9662094) q[1];
sx q[1];
rz(-1.4384603) q[1];
sx q[1];
rz(0.90465178) q[1];
x q[2];
rz(2.6252296) q[3];
sx q[3];
rz(-1.5011906) q[3];
sx q[3];
rz(-2.7854837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.34746927) q[2];
sx q[2];
rz(-2.4994714) q[2];
sx q[2];
rz(0.40880173) q[2];
rz(-0.18320228) q[3];
sx q[3];
rz(-1.5513523) q[3];
sx q[3];
rz(1.1387775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.113753) q[0];
sx q[0];
rz(-3.0945859) q[0];
sx q[0];
rz(-0.29079944) q[0];
rz(-2.193702) q[1];
sx q[1];
rz(-2.6746076) q[1];
sx q[1];
rz(-1.9416169) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79340807) q[0];
sx q[0];
rz(-0.68773341) q[0];
sx q[0];
rz(2.4853112) q[0];
x q[1];
rz(-2.3340204) q[2];
sx q[2];
rz(-2.117273) q[2];
sx q[2];
rz(-1.9131017) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0702182) q[1];
sx q[1];
rz(-1.5244251) q[1];
sx q[1];
rz(1.8646311) q[1];
rz(-pi) q[2];
rz(-3.0407716) q[3];
sx q[3];
rz(-1.3885048) q[3];
sx q[3];
rz(1.3829447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.36134186) q[2];
sx q[2];
rz(-1.4486518) q[2];
sx q[2];
rz(1.7928436) q[2];
rz(-1.5032984) q[3];
sx q[3];
rz(-2.2595451) q[3];
sx q[3];
rz(-2.9564814) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1143484) q[0];
sx q[0];
rz(-2.769727) q[0];
sx q[0];
rz(-1.0741023) q[0];
rz(-0.4298003) q[1];
sx q[1];
rz(-2.0975515) q[1];
sx q[1];
rz(2.6780186) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060424711) q[0];
sx q[0];
rz(-1.7588108) q[0];
sx q[0];
rz(-0.096676143) q[0];
x q[1];
rz(1.7724636) q[2];
sx q[2];
rz(-1.8275765) q[2];
sx q[2];
rz(-0.078291206) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9755755) q[1];
sx q[1];
rz(-0.75812997) q[1];
sx q[1];
rz(1.0075955) q[1];
rz(-pi) q[2];
rz(-0.14318569) q[3];
sx q[3];
rz(-2.5783263) q[3];
sx q[3];
rz(1.3384502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8672436) q[2];
sx q[2];
rz(-0.63259071) q[2];
sx q[2];
rz(-2.3331433) q[2];
rz(0.48458734) q[3];
sx q[3];
rz(-2.7633568) q[3];
sx q[3];
rz(2.7073879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7782068) q[0];
sx q[0];
rz(-0.45553842) q[0];
sx q[0];
rz(-1.1325915) q[0];
rz(0.038381902) q[1];
sx q[1];
rz(-1.8033586) q[1];
sx q[1];
rz(-2.8448232) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39958056) q[0];
sx q[0];
rz(-2.098408) q[0];
sx q[0];
rz(0.3216089) q[0];
rz(-2.692524) q[2];
sx q[2];
rz(-1.1658323) q[2];
sx q[2];
rz(1.8835889) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3546875) q[1];
sx q[1];
rz(-2.9070435) q[1];
sx q[1];
rz(2.812318) q[1];
rz(-pi) q[2];
rz(1.1007916) q[3];
sx q[3];
rz(-0.74277011) q[3];
sx q[3];
rz(-0.41760412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.177629) q[2];
sx q[2];
rz(-1.12135) q[2];
sx q[2];
rz(2.5753944) q[2];
rz(-1.1244134) q[3];
sx q[3];
rz(-3.0163613) q[3];
sx q[3];
rz(-1.9521149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.88846702) q[0];
sx q[0];
rz(-0.69857004) q[0];
sx q[0];
rz(0.10733124) q[0];
rz(2.879907) q[1];
sx q[1];
rz(-1.9410004) q[1];
sx q[1];
rz(0.95071361) q[1];
rz(-1.5952806) q[2];
sx q[2];
rz(-1.4888121) q[2];
sx q[2];
rz(-0.76406995) q[2];
rz(-2.2903743) q[3];
sx q[3];
rz(-2.2891582) q[3];
sx q[3];
rz(-1.4386406) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
