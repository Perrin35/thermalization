OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6899941) q[0];
sx q[0];
rz(-2.8328083) q[0];
sx q[0];
rz(0.2398332) q[0];
rz(-3.7026703) q[1];
sx q[1];
rz(2.3993888) q[1];
sx q[1];
rz(11.053434) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37555239) q[0];
sx q[0];
rz(-0.39012018) q[0];
sx q[0];
rz(1.9128156) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0885029) q[2];
sx q[2];
rz(-2.6183207) q[2];
sx q[2];
rz(2.5918617) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.35084942) q[1];
sx q[1];
rz(-2.6294957) q[1];
sx q[1];
rz(0.30330203) q[1];
x q[2];
rz(-1.2414819) q[3];
sx q[3];
rz(-0.17767492) q[3];
sx q[3];
rz(2.2840281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.98770398) q[2];
sx q[2];
rz(-2.14553) q[2];
sx q[2];
rz(1.055701) q[2];
rz(-0.40134564) q[3];
sx q[3];
rz(-1.6258806) q[3];
sx q[3];
rz(-1.6373985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.5098679) q[0];
sx q[0];
rz(-2.8724176) q[0];
sx q[0];
rz(-1.4058231) q[0];
rz(0.74613219) q[1];
sx q[1];
rz(-1.965062) q[1];
sx q[1];
rz(3.0768118) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21532741) q[0];
sx q[0];
rz(-2.522922) q[0];
sx q[0];
rz(-0.72665213) q[0];
rz(-pi) q[1];
rz(1.4006813) q[2];
sx q[2];
rz(-0.89808849) q[2];
sx q[2];
rz(1.9667039) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.244811) q[1];
sx q[1];
rz(-2.9771388) q[1];
sx q[1];
rz(1.0599778) q[1];
rz(1.4218036) q[3];
sx q[3];
rz(-2.9284366) q[3];
sx q[3];
rz(-3.0557291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.259321) q[2];
sx q[2];
rz(-0.52053014) q[2];
sx q[2];
rz(-0.041291324) q[2];
rz(2.0222372) q[3];
sx q[3];
rz(-1.3675523) q[3];
sx q[3];
rz(1.1411427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3248046) q[0];
sx q[0];
rz(-2.7140706) q[0];
sx q[0];
rz(-0.69931716) q[0];
rz(-2.518867) q[1];
sx q[1];
rz(-0.8131578) q[1];
sx q[1];
rz(2.8831388) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39614933) q[0];
sx q[0];
rz(-0.93702261) q[0];
sx q[0];
rz(1.5859912) q[0];
x q[1];
rz(2.223143) q[2];
sx q[2];
rz(-1.47965) q[2];
sx q[2];
rz(-0.80207588) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.099951) q[1];
sx q[1];
rz(-1.1219684) q[1];
sx q[1];
rz(-1.5484518) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.650189) q[3];
sx q[3];
rz(-0.43118048) q[3];
sx q[3];
rz(3.0357547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.93566018) q[2];
sx q[2];
rz(-1.6814597) q[2];
sx q[2];
rz(0.70651954) q[2];
rz(2.5126854) q[3];
sx q[3];
rz(-2.3157412) q[3];
sx q[3];
rz(2.0188873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.127447) q[0];
sx q[0];
rz(-1.4268459) q[0];
sx q[0];
rz(-2.1674147) q[0];
rz(1.4840508) q[1];
sx q[1];
rz(-1.0933135) q[1];
sx q[1];
rz(-2.1407703) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7324686) q[0];
sx q[0];
rz(-1.7268469) q[0];
sx q[0];
rz(2.1740995) q[0];
rz(-pi) q[1];
rz(-2.328726) q[2];
sx q[2];
rz(-2.1846601) q[2];
sx q[2];
rz(0.46068207) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.82881935) q[1];
sx q[1];
rz(-0.76550325) q[1];
sx q[1];
rz(2.440052) q[1];
rz(0.3206887) q[3];
sx q[3];
rz(-1.5237892) q[3];
sx q[3];
rz(0.21002029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4546844) q[2];
sx q[2];
rz(-1.3409706) q[2];
sx q[2];
rz(-2.3243813) q[2];
rz(-2.5456083) q[3];
sx q[3];
rz(-1.417336) q[3];
sx q[3];
rz(1.7433085) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3572094) q[0];
sx q[0];
rz(-1.7359808) q[0];
sx q[0];
rz(-1.0719517) q[0];
rz(-1.1373854) q[1];
sx q[1];
rz(-1.6268566) q[1];
sx q[1];
rz(0.17328182) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77861518) q[0];
sx q[0];
rz(-2.5234875) q[0];
sx q[0];
rz(0.72160665) q[0];
rz(-pi) q[1];
rz(-0.38727792) q[2];
sx q[2];
rz(-1.1114632) q[2];
sx q[2];
rz(2.0044553) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2318423) q[1];
sx q[1];
rz(-0.76059231) q[1];
sx q[1];
rz(1.4666345) q[1];
rz(-0.16290476) q[3];
sx q[3];
rz(-1.6626193) q[3];
sx q[3];
rz(3.0260928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5002354) q[2];
sx q[2];
rz(-0.45318979) q[2];
sx q[2];
rz(0.87453169) q[2];
rz(-0.66679653) q[3];
sx q[3];
rz(-1.4614481) q[3];
sx q[3];
rz(-0.82505208) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85103971) q[0];
sx q[0];
rz(-0.14851004) q[0];
sx q[0];
rz(0.17459757) q[0];
rz(-0.54706508) q[1];
sx q[1];
rz(-0.51536307) q[1];
sx q[1];
rz(1.2778767) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39659835) q[0];
sx q[0];
rz(-1.5088668) q[0];
sx q[0];
rz(1.2838191) q[0];
x q[1];
rz(-0.49680423) q[2];
sx q[2];
rz(-0.81648177) q[2];
sx q[2];
rz(1.7330488) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3681952) q[1];
sx q[1];
rz(-2.2822579) q[1];
sx q[1];
rz(-1.0991286) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2383591) q[3];
sx q[3];
rz(-1.0683016) q[3];
sx q[3];
rz(0.63695217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7586907) q[2];
sx q[2];
rz(-2.8768657) q[2];
sx q[2];
rz(-0.36941377) q[2];
rz(-0.82768011) q[3];
sx q[3];
rz(-1.8281432) q[3];
sx q[3];
rz(-2.9874492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5740042) q[0];
sx q[0];
rz(-2.2249157) q[0];
sx q[0];
rz(2.9045203) q[0];
rz(1.7489307) q[1];
sx q[1];
rz(-1.9475513) q[1];
sx q[1];
rz(-1.0038092) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0802409) q[0];
sx q[0];
rz(-2.6727242) q[0];
sx q[0];
rz(-2.7659225) q[0];
rz(-pi) q[1];
rz(1.1091484) q[2];
sx q[2];
rz(-2.2301815) q[2];
sx q[2];
rz(2.387799) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8449515) q[1];
sx q[1];
rz(-1.9012723) q[1];
sx q[1];
rz(0.56092324) q[1];
x q[2];
rz(0.69606651) q[3];
sx q[3];
rz(-2.4579774) q[3];
sx q[3];
rz(0.044134951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6134593) q[2];
sx q[2];
rz(-1.5789092) q[2];
sx q[2];
rz(-2.6252739) q[2];
rz(-2.3033219) q[3];
sx q[3];
rz(-1.6603989) q[3];
sx q[3];
rz(2.9576438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2876005) q[0];
sx q[0];
rz(-2.8599399) q[0];
sx q[0];
rz(-0.91424346) q[0];
rz(-2.5573348) q[1];
sx q[1];
rz(-1.7353568) q[1];
sx q[1];
rz(-0.90726888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7210684) q[0];
sx q[0];
rz(-2.4926909) q[0];
sx q[0];
rz(1.6008928) q[0];
rz(-1.3970988) q[2];
sx q[2];
rz(-0.96250421) q[2];
sx q[2];
rz(-2.4243958) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8311892) q[1];
sx q[1];
rz(-1.1205648) q[1];
sx q[1];
rz(0.72527171) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97213718) q[3];
sx q[3];
rz(-1.7936417) q[3];
sx q[3];
rz(3.0226662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3153136) q[2];
sx q[2];
rz(-1.7118914) q[2];
sx q[2];
rz(-2.2798174) q[2];
rz(1.9050542) q[3];
sx q[3];
rz(-3.0573513) q[3];
sx q[3];
rz(1.9305362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5009907) q[0];
sx q[0];
rz(-0.7826829) q[0];
sx q[0];
rz(2.1114517) q[0];
rz(-0.31632272) q[1];
sx q[1];
rz(-1.7681237) q[1];
sx q[1];
rz(0.28824678) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7310627) q[0];
sx q[0];
rz(-1.1708492) q[0];
sx q[0];
rz(1.7284786) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.64916237) q[2];
sx q[2];
rz(-1.9948975) q[2];
sx q[2];
rz(2.8555388) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.94142524) q[1];
sx q[1];
rz(-1.1516478) q[1];
sx q[1];
rz(2.2302385) q[1];
x q[2];
rz(-2.8602198) q[3];
sx q[3];
rz(-0.69123879) q[3];
sx q[3];
rz(-0.465525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.85252243) q[2];
sx q[2];
rz(-1.4129637) q[2];
sx q[2];
rz(2.9150325) q[2];
rz(-1.6864927) q[3];
sx q[3];
rz(-2.2454567) q[3];
sx q[3];
rz(0.69211012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15014547) q[0];
sx q[0];
rz(-1.2757855) q[0];
sx q[0];
rz(2.5164497) q[0];
rz(-1.9236247) q[1];
sx q[1];
rz(-0.70911276) q[1];
sx q[1];
rz(-1.2120754) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0481794) q[0];
sx q[0];
rz(-2.163889) q[0];
sx q[0];
rz(-2.1593447) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3122968) q[2];
sx q[2];
rz(-0.91583672) q[2];
sx q[2];
rz(-2.0531246) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4752848) q[1];
sx q[1];
rz(-1.5130677) q[1];
sx q[1];
rz(-1.1294731) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2792222) q[3];
sx q[3];
rz(-1.8599038) q[3];
sx q[3];
rz(1.5559199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43442279) q[2];
sx q[2];
rz(-1.7662798) q[2];
sx q[2];
rz(-2.0909069) q[2];
rz(2.5896942) q[3];
sx q[3];
rz(-2.4514908) q[3];
sx q[3];
rz(-1.8570073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62591775) q[0];
sx q[0];
rz(-0.82763012) q[0];
sx q[0];
rz(-0.81830842) q[0];
rz(0.7863518) q[1];
sx q[1];
rz(-2.4813589) q[1];
sx q[1];
rz(-2.9207041) q[1];
rz(-1.2073866) q[2];
sx q[2];
rz(-1.5354234) q[2];
sx q[2];
rz(-1.0101049) q[2];
rz(2.4554853) q[3];
sx q[3];
rz(-1.6391564) q[3];
sx q[3];
rz(2.2166722) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
