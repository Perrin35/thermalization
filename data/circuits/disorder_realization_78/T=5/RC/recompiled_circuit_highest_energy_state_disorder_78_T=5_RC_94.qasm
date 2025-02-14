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
rz(1.5155091) q[0];
sx q[0];
rz(-2.8647487) q[0];
sx q[0];
rz(-0.11685637) q[0];
rz(-0.46861831) q[1];
sx q[1];
rz(-2.7873971) q[1];
sx q[1];
rz(2.6503704) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.011176414) q[0];
sx q[0];
rz(-2.6493086) q[0];
sx q[0];
rz(-1.6947794) q[0];
rz(-pi) q[1];
rz(0.38552706) q[2];
sx q[2];
rz(-2.6221199) q[2];
sx q[2];
rz(0.2310209) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.74581214) q[1];
sx q[1];
rz(-2.2015101) q[1];
sx q[1];
rz(-1.4480026) q[1];
rz(3.1163636) q[3];
sx q[3];
rz(-1.9140035) q[3];
sx q[3];
rz(2.0398839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3789856) q[2];
sx q[2];
rz(-2.8540322) q[2];
sx q[2];
rz(-2.911705) q[2];
rz(-0.057279438) q[3];
sx q[3];
rz(-2.4725437) q[3];
sx q[3];
rz(-0.91276401) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7951935) q[0];
sx q[0];
rz(-0.38723463) q[0];
sx q[0];
rz(0.21695319) q[0];
rz(0.03288658) q[1];
sx q[1];
rz(-2.2078728) q[1];
sx q[1];
rz(0.65509534) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.396916) q[0];
sx q[0];
rz(-1.370805) q[0];
sx q[0];
rz(-2.295664) q[0];
rz(-2.7883165) q[2];
sx q[2];
rz(-2.6124119) q[2];
sx q[2];
rz(2.8904817) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8959143) q[1];
sx q[1];
rz(-2.0729471) q[1];
sx q[1];
rz(-2.4163626) q[1];
rz(-0.14343393) q[3];
sx q[3];
rz(-2.3000882) q[3];
sx q[3];
rz(-2.7021331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.32404262) q[2];
sx q[2];
rz(-2.5174759) q[2];
sx q[2];
rz(-1.8435271) q[2];
rz(-1.921418) q[3];
sx q[3];
rz(-1.7036567) q[3];
sx q[3];
rz(-2.4841681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2636181) q[0];
sx q[0];
rz(-2.2442696) q[0];
sx q[0];
rz(0.61400145) q[0];
rz(-1.3677431) q[1];
sx q[1];
rz(-2.4847993) q[1];
sx q[1];
rz(0.49555379) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6959636) q[0];
sx q[0];
rz(-0.82484748) q[0];
sx q[0];
rz(-2.4936409) q[0];
x q[1];
rz(2.2587288) q[2];
sx q[2];
rz(-1.1585711) q[2];
sx q[2];
rz(2.8799873) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4877971) q[1];
sx q[1];
rz(-1.7572409) q[1];
sx q[1];
rz(-1.814278) q[1];
x q[2];
rz(-1.2167778) q[3];
sx q[3];
rz(-1.1292158) q[3];
sx q[3];
rz(-1.7596678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.76501781) q[2];
sx q[2];
rz(-1.0708662) q[2];
sx q[2];
rz(-0.59256727) q[2];
rz(2.9060034) q[3];
sx q[3];
rz(-2.3841136) q[3];
sx q[3];
rz(-0.45430115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(2.4950824) q[0];
sx q[0];
rz(-2.3401234) q[0];
sx q[0];
rz(3.1357646) q[0];
rz(-2.5616772) q[1];
sx q[1];
rz(-2.7858851) q[1];
sx q[1];
rz(0.52111202) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2297008) q[0];
sx q[0];
rz(-1.5661514) q[0];
sx q[0];
rz(1.5636496) q[0];
rz(-pi) q[1];
x q[1];
rz(0.02946202) q[2];
sx q[2];
rz(-0.97004393) q[2];
sx q[2];
rz(-0.0044057152) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.191491) q[1];
sx q[1];
rz(-1.9708212) q[1];
sx q[1];
rz(-2.6482976) q[1];
rz(-pi) q[2];
rz(1.6507876) q[3];
sx q[3];
rz(-1.2595525) q[3];
sx q[3];
rz(-2.1224006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8899322) q[2];
sx q[2];
rz(-2.6106788) q[2];
sx q[2];
rz(2.8233675) q[2];
rz(2.4288154) q[3];
sx q[3];
rz(-2.8409499) q[3];
sx q[3];
rz(-1.6596863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37079674) q[0];
sx q[0];
rz(-0.21368055) q[0];
sx q[0];
rz(-3.1053012) q[0];
rz(-0.52641422) q[1];
sx q[1];
rz(-1.6830091) q[1];
sx q[1];
rz(2.859419) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68809915) q[0];
sx q[0];
rz(-0.69506139) q[0];
sx q[0];
rz(2.009925) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5897609) q[2];
sx q[2];
rz(-0.95317344) q[2];
sx q[2];
rz(-0.080452563) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6732521) q[1];
sx q[1];
rz(-0.70827228) q[1];
sx q[1];
rz(3.0336223) q[1];
x q[2];
rz(2.3822278) q[3];
sx q[3];
rz(-2.1753484) q[3];
sx q[3];
rz(-0.60455632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.37461773) q[2];
sx q[2];
rz(-2.4800315) q[2];
sx q[2];
rz(3.0401163) q[2];
rz(2.1756419) q[3];
sx q[3];
rz(-2.2209397) q[3];
sx q[3];
rz(0.11400338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58979708) q[0];
sx q[0];
rz(-0.17687251) q[0];
sx q[0];
rz(1.9675323) q[0];
rz(-0.38408285) q[1];
sx q[1];
rz(-1.9575155) q[1];
sx q[1];
rz(-0.039965872) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5863881) q[0];
sx q[0];
rz(-0.38643906) q[0];
sx q[0];
rz(-2.1512335) q[0];
x q[1];
rz(0.33127516) q[2];
sx q[2];
rz(-2.185198) q[2];
sx q[2];
rz(1.5744874) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.080885305) q[1];
sx q[1];
rz(-2.1645912) q[1];
sx q[1];
rz(1.1546561) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8043936) q[3];
sx q[3];
rz(-0.75832483) q[3];
sx q[3];
rz(-0.17994623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.93969718) q[2];
sx q[2];
rz(-0.59212089) q[2];
sx q[2];
rz(2.3218018) q[2];
rz(0.78153265) q[3];
sx q[3];
rz(-2.235409) q[3];
sx q[3];
rz(1.3801105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.49067295) q[0];
sx q[0];
rz(-2.1020205) q[0];
sx q[0];
rz(1.8248722) q[0];
rz(1.7735749) q[1];
sx q[1];
rz(-1.2317692) q[1];
sx q[1];
rz(2.7046611) q[1];
rz(-pi) q[2];
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
rz(3.1243665) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8812411) q[1];
sx q[1];
rz(-0.8484183) q[1];
sx q[1];
rz(-0.48697014) q[1];
rz(-pi) q[2];
rz(-2.223001) q[3];
sx q[3];
rz(-0.64213412) q[3];
sx q[3];
rz(-1.0561424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4009319) q[2];
sx q[2];
rz(-1.0793945) q[2];
sx q[2];
rz(-0.13460049) q[2];
rz(1.2234737) q[3];
sx q[3];
rz(-1.7569907) q[3];
sx q[3];
rz(-1.9756165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3400035) q[0];
sx q[0];
rz(-0.29186258) q[0];
sx q[0];
rz(2.959751) q[0];
rz(2.8002411) q[1];
sx q[1];
rz(-1.1266484) q[1];
sx q[1];
rz(-2.9827319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10578622) q[0];
sx q[0];
rz(-1.5022359) q[0];
sx q[0];
rz(-0.44802702) q[0];
x q[1];
rz(-2.7188825) q[2];
sx q[2];
rz(-1.8313932) q[2];
sx q[2];
rz(-0.24525951) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.73028683) q[1];
sx q[1];
rz(-1.3540233) q[1];
sx q[1];
rz(-2.8764399) q[1];
x q[2];
rz(2.8980632) q[3];
sx q[3];
rz(-0.56748828) q[3];
sx q[3];
rz(1.6220038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2813256) q[2];
sx q[2];
rz(-0.93026668) q[2];
sx q[2];
rz(0.88480985) q[2];
rz(1.4191215) q[3];
sx q[3];
rz(-1.0249187) q[3];
sx q[3];
rz(-2.6917246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29174969) q[0];
sx q[0];
rz(-1.869864) q[0];
sx q[0];
rz(2.1480609) q[0];
rz(1.3604856) q[1];
sx q[1];
rz(-2.2527756) q[1];
sx q[1];
rz(0.55307585) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6926413) q[0];
sx q[0];
rz(-1.3425437) q[0];
sx q[0];
rz(2.7518163) q[0];
x q[1];
rz(-3.0680066) q[2];
sx q[2];
rz(-0.739817) q[2];
sx q[2];
rz(0.70895665) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7375664) q[1];
sx q[1];
rz(-1.3133004) q[1];
sx q[1];
rz(2.5016258) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9286128) q[3];
sx q[3];
rz(-2.4705437) q[3];
sx q[3];
rz(-0.34085654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.20455827) q[2];
sx q[2];
rz(-2.5356346) q[2];
sx q[2];
rz(-2.9045612) q[2];
rz(-1.6009181) q[3];
sx q[3];
rz(-2.377244) q[3];
sx q[3];
rz(1.9717533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071526214) q[0];
sx q[0];
rz(-1.6306174) q[0];
sx q[0];
rz(-0.1189098) q[0];
rz(-2.7023244) q[1];
sx q[1];
rz(-1.8426789) q[1];
sx q[1];
rz(-0.063442245) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6574888) q[0];
sx q[0];
rz(-2.2229337) q[0];
sx q[0];
rz(0.87638292) q[0];
x q[1];
rz(0.62427743) q[2];
sx q[2];
rz(-1.7221525) q[2];
sx q[2];
rz(-0.69619149) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6193004) q[1];
sx q[1];
rz(-1.5792773) q[1];
sx q[1];
rz(-1.897476) q[1];
rz(0.48096913) q[3];
sx q[3];
rz(-2.0837373) q[3];
sx q[3];
rz(-0.45403593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0110317) q[2];
sx q[2];
rz(-1.0728027) q[2];
sx q[2];
rz(1.4981184) q[2];
rz(2.2091852) q[3];
sx q[3];
rz(-2.0216209) q[3];
sx q[3];
rz(1.9278661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7508004) q[0];
sx q[0];
rz(-1.9600497) q[0];
sx q[0];
rz(2.0586769) q[0];
rz(-2.2687268) q[1];
sx q[1];
rz(-1.083359) q[1];
sx q[1];
rz(-0.92360003) q[1];
rz(-0.68442576) q[2];
sx q[2];
rz(-1.050907) q[2];
sx q[2];
rz(2.4466865) q[2];
rz(-0.57742248) q[3];
sx q[3];
rz(-1.631177) q[3];
sx q[3];
rz(-2.9384818) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
