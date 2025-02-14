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
rz(-1.6260835) q[0];
sx q[0];
rz(-0.27684394) q[0];
sx q[0];
rz(0.11685637) q[0];
rz(-0.46861831) q[1];
sx q[1];
rz(-2.7873971) q[1];
sx q[1];
rz(2.6503704) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011176414) q[0];
sx q[0];
rz(-0.49228406) q[0];
sx q[0];
rz(1.4468132) q[0];
x q[1];
rz(-0.38552706) q[2];
sx q[2];
rz(-2.6221199) q[2];
sx q[2];
rz(-0.2310209) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.95211103) q[1];
sx q[1];
rz(-2.5006388) q[1];
sx q[1];
rz(2.9753995) q[1];
x q[2];
rz(-3.1163636) q[3];
sx q[3];
rz(-1.9140035) q[3];
sx q[3];
rz(1.1017088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.7626071) q[2];
sx q[2];
rz(-2.8540322) q[2];
sx q[2];
rz(-0.22988764) q[2];
rz(-3.0843132) q[3];
sx q[3];
rz(-2.4725437) q[3];
sx q[3];
rz(0.91276401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7951935) q[0];
sx q[0];
rz(-0.38723463) q[0];
sx q[0];
rz(-0.21695319) q[0];
rz(0.03288658) q[1];
sx q[1];
rz(-2.2078728) q[1];
sx q[1];
rz(0.65509534) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0467619) q[0];
sx q[0];
rz(-2.3944986) q[0];
sx q[0];
rz(1.2741035) q[0];
rz(2.7883165) q[2];
sx q[2];
rz(-0.52918079) q[2];
sx q[2];
rz(-0.25111094) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.078121878) q[1];
sx q[1];
rz(-0.95032108) q[1];
sx q[1];
rz(-0.93777754) q[1];
rz(-0.83637303) q[3];
sx q[3];
rz(-1.4640088) q[3];
sx q[3];
rz(-1.9143145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.81755) q[2];
sx q[2];
rz(-0.62411672) q[2];
sx q[2];
rz(-1.8435271) q[2];
rz(1.921418) q[3];
sx q[3];
rz(-1.7036567) q[3];
sx q[3];
rz(-0.65742457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(2.2636181) q[0];
sx q[0];
rz(-2.2442696) q[0];
sx q[0];
rz(-2.5275912) q[0];
rz(-1.3677431) q[1];
sx q[1];
rz(-0.65679336) q[1];
sx q[1];
rz(2.6460389) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6959636) q[0];
sx q[0];
rz(-0.82484748) q[0];
sx q[0];
rz(-0.64795171) q[0];
x q[1];
rz(2.6265385) q[2];
sx q[2];
rz(-0.94991377) q[2];
sx q[2];
rz(0.99109288) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0125791) q[1];
sx q[1];
rz(-1.8099752) q[1];
sx q[1];
rz(-0.19197477) q[1];
rz(-2.6747781) q[3];
sx q[3];
rz(-1.2520077) q[3];
sx q[3];
rz(-2.7960645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3765748) q[2];
sx q[2];
rz(-2.0707264) q[2];
sx q[2];
rz(0.59256727) q[2];
rz(-2.9060034) q[3];
sx q[3];
rz(-2.3841136) q[3];
sx q[3];
rz(0.45430115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(0.6465103) q[0];
sx q[0];
rz(-2.3401234) q[0];
sx q[0];
rz(-3.1357646) q[0];
rz(0.5799154) q[1];
sx q[1];
rz(-0.35570759) q[1];
sx q[1];
rz(-0.52111202) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65893764) q[0];
sx q[0];
rz(-1.5636496) q[0];
sx q[0];
rz(0.004645017) q[0];
rz(0.96984152) q[2];
sx q[2];
rz(-1.5950987) q[2];
sx q[2];
rz(-1.5497335) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.191491) q[1];
sx q[1];
rz(-1.1707715) q[1];
sx q[1];
rz(0.49329503) q[1];
rz(2.829414) q[3];
sx q[3];
rz(-1.494656) q[3];
sx q[3];
rz(2.5654441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2516605) q[2];
sx q[2];
rz(-2.6106788) q[2];
sx q[2];
rz(-2.8233675) q[2];
rz(-0.71277726) q[3];
sx q[3];
rz(-0.30064279) q[3];
sx q[3];
rz(1.6596863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37079674) q[0];
sx q[0];
rz(-0.21368055) q[0];
sx q[0];
rz(-0.036291432) q[0];
rz(0.52641422) q[1];
sx q[1];
rz(-1.4585835) q[1];
sx q[1];
rz(2.859419) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68809915) q[0];
sx q[0];
rz(-2.4465313) q[0];
sx q[0];
rz(2.009925) q[0];
rz(-pi) q[1];
rz(-2.2065978) q[2];
sx q[2];
rz(-0.80342573) q[2];
sx q[2];
rz(-0.73551169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.020318839) q[1];
sx q[1];
rz(-1.6409546) q[1];
sx q[1];
rz(-2.4362045) q[1];
x q[2];
rz(0.75936486) q[3];
sx q[3];
rz(-2.1753484) q[3];
sx q[3];
rz(0.60455632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.37461773) q[2];
sx q[2];
rz(-0.66156113) q[2];
sx q[2];
rz(-3.0401163) q[2];
rz(-2.1756419) q[3];
sx q[3];
rz(-0.92065293) q[3];
sx q[3];
rz(0.11400338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5517956) q[0];
sx q[0];
rz(-2.9647201) q[0];
sx q[0];
rz(1.1740603) q[0];
rz(2.7575098) q[1];
sx q[1];
rz(-1.1840772) q[1];
sx q[1];
rz(-3.1016268) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55520457) q[0];
sx q[0];
rz(-0.38643906) q[0];
sx q[0];
rz(0.9903592) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33127516) q[2];
sx q[2];
rz(-2.185198) q[2];
sx q[2];
rz(-1.5671052) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0607073) q[1];
sx q[1];
rz(-0.97700143) q[1];
sx q[1];
rz(1.1546561) q[1];
rz(-pi) q[2];
rz(-2.3153855) q[3];
sx q[3];
rz(-1.4109269) q[3];
sx q[3];
rz(1.9217971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2018955) q[2];
sx q[2];
rz(-2.5494718) q[2];
sx q[2];
rz(0.8197909) q[2];
rz(2.36006) q[3];
sx q[3];
rz(-0.90618366) q[3];
sx q[3];
rz(-1.7614822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49067295) q[0];
sx q[0];
rz(-1.0395721) q[0];
sx q[0];
rz(-1.8248722) q[0];
rz(-1.7735749) q[1];
sx q[1];
rz(-1.2317692) q[1];
sx q[1];
rz(-2.7046611) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36501554) q[0];
sx q[0];
rz(-1.3804589) q[0];
sx q[0];
rz(0.037419407) q[0];
x q[1];
rz(2.390449) q[2];
sx q[2];
rz(-2.9808384) q[2];
sx q[2];
rz(-0.017226177) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2059784) q[1];
sx q[1];
rz(-2.2956478) q[1];
sx q[1];
rz(2.0589252) q[1];
rz(-pi) q[2];
rz(-2.715492) q[3];
sx q[3];
rz(-1.0747194) q[3];
sx q[3];
rz(2.8471198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4009319) q[2];
sx q[2];
rz(-1.0793945) q[2];
sx q[2];
rz(-0.13460049) q[2];
rz(-1.918119) q[3];
sx q[3];
rz(-1.384602) q[3];
sx q[3];
rz(-1.1659762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80158919) q[0];
sx q[0];
rz(-2.8497301) q[0];
sx q[0];
rz(-2.959751) q[0];
rz(-2.8002411) q[1];
sx q[1];
rz(-1.1266484) q[1];
sx q[1];
rz(2.9827319) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0358064) q[0];
sx q[0];
rz(-1.6393568) q[0];
sx q[0];
rz(-2.6935656) q[0];
x q[1];
rz(2.7188825) q[2];
sx q[2];
rz(-1.8313932) q[2];
sx q[2];
rz(-2.8963331) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78217172) q[1];
sx q[1];
rz(-1.3119932) q[1];
sx q[1];
rz(-1.7951629) q[1];
rz(-pi) q[2];
rz(-2.8980632) q[3];
sx q[3];
rz(-2.5741044) q[3];
sx q[3];
rz(-1.5195888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86026704) q[2];
sx q[2];
rz(-2.211326) q[2];
sx q[2];
rz(-0.88480985) q[2];
rz(1.4191215) q[3];
sx q[3];
rz(-1.0249187) q[3];
sx q[3];
rz(-2.6917246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.849843) q[0];
sx q[0];
rz(-1.869864) q[0];
sx q[0];
rz(-2.1480609) q[0];
rz(1.3604856) q[1];
sx q[1];
rz(-2.2527756) q[1];
sx q[1];
rz(0.55307585) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6926413) q[0];
sx q[0];
rz(-1.3425437) q[0];
sx q[0];
rz(0.3897764) q[0];
rz(-pi) q[1];
rz(-2.4031248) q[2];
sx q[2];
rz(-1.5212125) q[2];
sx q[2];
rz(-0.91623437) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7874541) q[1];
sx q[1];
rz(-2.1864357) q[1];
sx q[1];
rz(1.8880185) q[1];
x q[2];
rz(-1.9286128) q[3];
sx q[3];
rz(-0.67104895) q[3];
sx q[3];
rz(0.34085654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.20455827) q[2];
sx q[2];
rz(-2.5356346) q[2];
sx q[2];
rz(-0.23703144) q[2];
rz(-1.5406746) q[3];
sx q[3];
rz(-0.76434869) q[3];
sx q[3];
rz(1.9717533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071526214) q[0];
sx q[0];
rz(-1.5109753) q[0];
sx q[0];
rz(0.1189098) q[0];
rz(2.7023244) q[1];
sx q[1];
rz(-1.8426789) q[1];
sx q[1];
rz(-3.0781504) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7164511) q[0];
sx q[0];
rz(-0.91380318) q[0];
sx q[0];
rz(2.444066) q[0];
rz(-pi) q[1];
rz(0.25524885) q[2];
sx q[2];
rz(-0.63997696) q[2];
sx q[2];
rz(-0.6682804) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6193004) q[1];
sx q[1];
rz(-1.5792773) q[1];
sx q[1];
rz(1.897476) q[1];
rz(0.48096913) q[3];
sx q[3];
rz(-1.0578553) q[3];
sx q[3];
rz(0.45403593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0110317) q[2];
sx q[2];
rz(-1.0728027) q[2];
sx q[2];
rz(-1.6434742) q[2];
rz(0.9324075) q[3];
sx q[3];
rz(-2.0216209) q[3];
sx q[3];
rz(-1.9278661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7508004) q[0];
sx q[0];
rz(-1.1815429) q[0];
sx q[0];
rz(-1.0829157) q[0];
rz(0.87286585) q[1];
sx q[1];
rz(-1.083359) q[1];
sx q[1];
rz(-0.92360003) q[1];
rz(-2.4058061) q[2];
sx q[2];
rz(-2.308261) q[2];
sx q[2];
rz(-2.8125662) q[2];
rz(0.11029966) q[3];
sx q[3];
rz(-0.58021373) q[3];
sx q[3];
rz(1.8662682) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
