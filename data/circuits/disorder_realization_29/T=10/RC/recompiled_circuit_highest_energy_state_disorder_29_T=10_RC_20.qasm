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
rz(-2.2588377) q[0];
sx q[0];
rz(-0.14482276) q[0];
sx q[0];
rz(-0.75135279) q[0];
rz(-1.6169647) q[1];
sx q[1];
rz(5.3137988) q[1];
sx q[1];
rz(9.245524) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67231945) q[0];
sx q[0];
rz(-1.0814965) q[0];
sx q[0];
rz(-2.8239522) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1044311) q[2];
sx q[2];
rz(-2.4183309) q[2];
sx q[2];
rz(0.0497555) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9589267) q[1];
sx q[1];
rz(-1.2616871) q[1];
sx q[1];
rz(1.6462417) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4901913) q[3];
sx q[3];
rz(-0.63524073) q[3];
sx q[3];
rz(1.0280392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2396607) q[2];
sx q[2];
rz(-1.8472981) q[2];
sx q[2];
rz(-2.530976) q[2];
rz(-2.2527952) q[3];
sx q[3];
rz(-2.461268) q[3];
sx q[3];
rz(2.4077267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69773847) q[0];
sx q[0];
rz(-2.839851) q[0];
sx q[0];
rz(-1.8119716) q[0];
rz(-1.656172) q[1];
sx q[1];
rz(-1.4621567) q[1];
sx q[1];
rz(-2.2775473) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3812019) q[0];
sx q[0];
rz(-1.4684033) q[0];
sx q[0];
rz(-1.824462) q[0];
rz(-2.0338221) q[2];
sx q[2];
rz(-1.290375) q[2];
sx q[2];
rz(-0.011649557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0985579) q[1];
sx q[1];
rz(-1.7435257) q[1];
sx q[1];
rz(0.14202001) q[1];
rz(0.95324272) q[3];
sx q[3];
rz(-0.88897486) q[3];
sx q[3];
rz(0.98495959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0580505) q[2];
sx q[2];
rz(-2.4594049) q[2];
sx q[2];
rz(-0.19134276) q[2];
rz(-0.30609104) q[3];
sx q[3];
rz(-1.6813262) q[3];
sx q[3];
rz(0.93650854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8323583) q[0];
sx q[0];
rz(-1.8620055) q[0];
sx q[0];
rz(2.7171296) q[0];
rz(-1.8008495) q[1];
sx q[1];
rz(-0.91573358) q[1];
sx q[1];
rz(-0.65139604) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0754335) q[0];
sx q[0];
rz(-1.4071583) q[0];
sx q[0];
rz(3.1343824) q[0];
x q[1];
rz(0.38369757) q[2];
sx q[2];
rz(-1.7406775) q[2];
sx q[2];
rz(-3.1184514) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0733903) q[1];
sx q[1];
rz(-1.7531839) q[1];
sx q[1];
rz(-2.6368111) q[1];
rz(0.15768361) q[3];
sx q[3];
rz(-1.7947949) q[3];
sx q[3];
rz(1.7832613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1059619) q[2];
sx q[2];
rz(-2.0023846) q[2];
sx q[2];
rz(-0.6967217) q[2];
rz(-0.68909711) q[3];
sx q[3];
rz(-1.1545811) q[3];
sx q[3];
rz(-2.885163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7032787) q[0];
sx q[0];
rz(-0.2916446) q[0];
sx q[0];
rz(2.1412204) q[0];
rz(1.6814303) q[1];
sx q[1];
rz(-1.6078452) q[1];
sx q[1];
rz(-1.9283074) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1532261) q[0];
sx q[0];
rz(-1.2906162) q[0];
sx q[0];
rz(3.0154254) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51474173) q[2];
sx q[2];
rz(-0.79823433) q[2];
sx q[2];
rz(-1.9412184) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.91649969) q[1];
sx q[1];
rz(-1.2429534) q[1];
sx q[1];
rz(-1.6986548) q[1];
rz(-2.1093879) q[3];
sx q[3];
rz(-1.4929139) q[3];
sx q[3];
rz(-0.36161446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1022243) q[2];
sx q[2];
rz(-0.82304707) q[2];
sx q[2];
rz(0.064112045) q[2];
rz(-0.72470775) q[3];
sx q[3];
rz(-1.9331845) q[3];
sx q[3];
rz(2.7418315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940014) q[0];
sx q[0];
rz(-1.8444909) q[0];
sx q[0];
rz(-2.3498348) q[0];
rz(1.1119615) q[1];
sx q[1];
rz(-2.0826191) q[1];
sx q[1];
rz(1.3349104) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4348818) q[0];
sx q[0];
rz(-1.1010873) q[0];
sx q[0];
rz(-2.241894) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58081268) q[2];
sx q[2];
rz(-1.6585095) q[2];
sx q[2];
rz(2.650819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.47274703) q[1];
sx q[1];
rz(-0.18660523) q[1];
sx q[1];
rz(-1.9884459) q[1];
rz(-pi) q[2];
rz(0.99467268) q[3];
sx q[3];
rz(-1.398456) q[3];
sx q[3];
rz(1.4927499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2935334) q[2];
sx q[2];
rz(-0.89911014) q[2];
sx q[2];
rz(1.1406356) q[2];
rz(-0.10661495) q[3];
sx q[3];
rz(-1.5299608) q[3];
sx q[3];
rz(-0.30985668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3265729) q[0];
sx q[0];
rz(-2.6928379) q[0];
sx q[0];
rz(-3.1383681) q[0];
rz(-0.047317304) q[1];
sx q[1];
rz(-1.4553757) q[1];
sx q[1];
rz(1.7914194) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2605476) q[0];
sx q[0];
rz(-2.8821917) q[0];
sx q[0];
rz(1.3400643) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5642942) q[2];
sx q[2];
rz(-0.98492981) q[2];
sx q[2];
rz(-0.68638869) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9474831) q[1];
sx q[1];
rz(-1.0110185) q[1];
sx q[1];
rz(-2.4653788) q[1];
x q[2];
rz(-1.4357225) q[3];
sx q[3];
rz(-0.22264847) q[3];
sx q[3];
rz(0.50375736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.99001592) q[2];
sx q[2];
rz(-0.185597) q[2];
sx q[2];
rz(-0.081175096) q[2];
rz(0.89933991) q[3];
sx q[3];
rz(-2.3266561) q[3];
sx q[3];
rz(-1.2871294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7154295) q[0];
sx q[0];
rz(-1.3112712) q[0];
sx q[0];
rz(0.50450605) q[0];
rz(-2.1465178) q[1];
sx q[1];
rz(-1.0202531) q[1];
sx q[1];
rz(3.1212433) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58344719) q[0];
sx q[0];
rz(-2.7648395) q[0];
sx q[0];
rz(1.8280297) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0145524) q[2];
sx q[2];
rz(-1.3641225) q[2];
sx q[2];
rz(-2.5673696) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.099906057) q[1];
sx q[1];
rz(-1.761224) q[1];
sx q[1];
rz(1.8407525) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.048681569) q[3];
sx q[3];
rz(-1.657007) q[3];
sx q[3];
rz(1.8415368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4703579) q[2];
sx q[2];
rz(-1.8378259) q[2];
sx q[2];
rz(1.7669558) q[2];
rz(-1.6124407) q[3];
sx q[3];
rz(-2.0917454) q[3];
sx q[3];
rz(3.0488739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90877157) q[0];
sx q[0];
rz(-0.11255539) q[0];
sx q[0];
rz(1.0821279) q[0];
rz(0.96083653) q[1];
sx q[1];
rz(-1.4832393) q[1];
sx q[1];
rz(2.198641) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0884224) q[0];
sx q[0];
rz(-2.5665847) q[0];
sx q[0];
rz(1.0002329) q[0];
x q[1];
rz(-2.2029938) q[2];
sx q[2];
rz(-1.7603121) q[2];
sx q[2];
rz(2.0507484) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3228938) q[1];
sx q[1];
rz(-2.2320691) q[1];
sx q[1];
rz(-0.7612919) q[1];
rz(-3.037463) q[3];
sx q[3];
rz(-1.2648598) q[3];
sx q[3];
rz(2.1776074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1711787) q[2];
sx q[2];
rz(-1.8381939) q[2];
sx q[2];
rz(-1.9937493) q[2];
rz(-2.7796699) q[3];
sx q[3];
rz(-1.7636834) q[3];
sx q[3];
rz(1.743478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4107133) q[0];
sx q[0];
rz(-2.78237) q[0];
sx q[0];
rz(3.0391589) q[0];
rz(2.5275285) q[1];
sx q[1];
rz(-1.026261) q[1];
sx q[1];
rz(-2.3238497) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9303904) q[0];
sx q[0];
rz(-1.0387392) q[0];
sx q[0];
rz(-2.3389057) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2417415) q[2];
sx q[2];
rz(-1.6616115) q[2];
sx q[2];
rz(-2.1911774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.39108155) q[1];
sx q[1];
rz(-0.72946786) q[1];
sx q[1];
rz(-1.5991648) q[1];
x q[2];
rz(-2.6800167) q[3];
sx q[3];
rz(-0.99050039) q[3];
sx q[3];
rz(0.020315276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3488591) q[2];
sx q[2];
rz(-1.0772971) q[2];
sx q[2];
rz(0.31361541) q[2];
rz(2.6775728) q[3];
sx q[3];
rz(-0.84657621) q[3];
sx q[3];
rz(-2.2217506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87026507) q[0];
sx q[0];
rz(-2.3509404) q[0];
sx q[0];
rz(-2.9050997) q[0];
rz(0.21367167) q[1];
sx q[1];
rz(-2.2387319) q[1];
sx q[1];
rz(0.22458354) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7296627) q[0];
sx q[0];
rz(-0.87422919) q[0];
sx q[0];
rz(2.4921472) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6170391) q[2];
sx q[2];
rz(-0.65215014) q[2];
sx q[2];
rz(-0.99936501) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1285291) q[1];
sx q[1];
rz(-2.6010102) q[1];
sx q[1];
rz(-2.9180727) q[1];
rz(-pi) q[2];
rz(2.7216689) q[3];
sx q[3];
rz(-2.2444199) q[3];
sx q[3];
rz(-2.6633584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1647722) q[2];
sx q[2];
rz(-0.87146622) q[2];
sx q[2];
rz(0.15178794) q[2];
rz(-3.1101036) q[3];
sx q[3];
rz(-1.7446691) q[3];
sx q[3];
rz(0.12846863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0781773) q[0];
sx q[0];
rz(-1.5237533) q[0];
sx q[0];
rz(-0.40405003) q[0];
rz(-1.0833441) q[1];
sx q[1];
rz(-1.5647519) q[1];
sx q[1];
rz(-1.5819989) q[1];
rz(1.622621) q[2];
sx q[2];
rz(-1.3695649) q[2];
sx q[2];
rz(2.8893378) q[2];
rz(2.8491889) q[3];
sx q[3];
rz(-0.3429827) q[3];
sx q[3];
rz(-0.61957785) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
