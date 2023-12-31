OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2330115) q[0];
sx q[0];
rz(-1.1865948) q[0];
sx q[0];
rz(-1.9957805) q[0];
rz(0.97283483) q[1];
sx q[1];
rz(-1.4714779) q[1];
sx q[1];
rz(0.29247984) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777893) q[0];
sx q[0];
rz(-1.5153432) q[0];
sx q[0];
rz(1.1763014) q[0];
rz(-pi) q[1];
rz(-0.77941676) q[2];
sx q[2];
rz(-1.6072818) q[2];
sx q[2];
rz(-1.905029) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.18947345) q[1];
sx q[1];
rz(-1.6166302) q[1];
sx q[1];
rz(2.1650251) q[1];
rz(-pi) q[2];
rz(-3.1396477) q[3];
sx q[3];
rz(-1.0969775) q[3];
sx q[3];
rz(2.4634944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.73074377) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(0.4494108) q[2];
rz(2.4959026) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(-1.2908363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9280424) q[0];
sx q[0];
rz(-2.2556861) q[0];
sx q[0];
rz(0.74203062) q[0];
rz(1.4713326) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(1.3630294) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33107685) q[0];
sx q[0];
rz(-1.1805981) q[0];
sx q[0];
rz(0.66731989) q[0];
rz(0.46911932) q[2];
sx q[2];
rz(-1.763478) q[2];
sx q[2];
rz(-2.0872781) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0568309) q[1];
sx q[1];
rz(-2.8014604) q[1];
sx q[1];
rz(0.057573307) q[1];
rz(-pi) q[2];
rz(2.1915073) q[3];
sx q[3];
rz(-0.32334298) q[3];
sx q[3];
rz(-0.33085535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8403975) q[2];
sx q[2];
rz(-2.6448554) q[2];
sx q[2];
rz(-1.2724686) q[2];
rz(-2.8043591) q[3];
sx q[3];
rz(-1.971258) q[3];
sx q[3];
rz(3.1315394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(3.0740046) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(2.8787676) q[0];
rz(2.7858531) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(1.9062818) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0316186) q[0];
sx q[0];
rz(-1.5556766) q[0];
sx q[0];
rz(-1.5605687) q[0];
rz(-pi) q[1];
rz(-3.1149408) q[2];
sx q[2];
rz(-1.3304552) q[2];
sx q[2];
rz(1.196256) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1097818) q[1];
sx q[1];
rz(-0.33966741) q[1];
sx q[1];
rz(2.8207645) q[1];
rz(0.34706195) q[3];
sx q[3];
rz(-2.2357781) q[3];
sx q[3];
rz(0.60225981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0170903) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(2.9612605) q[2];
rz(0.63594937) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(0.37160555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76698774) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(-2.4225127) q[0];
rz(-0.51302296) q[1];
sx q[1];
rz(-1.2027556) q[1];
sx q[1];
rz(2.1069353) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35095222) q[0];
sx q[0];
rz(-1.6296903) q[0];
sx q[0];
rz(-1.5509997) q[0];
rz(-pi) q[1];
rz(3.0940975) q[2];
sx q[2];
rz(-1.7058027) q[2];
sx q[2];
rz(2.6035655) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3116341) q[1];
sx q[1];
rz(-1.2882075) q[1];
sx q[1];
rz(0.71210536) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47311584) q[3];
sx q[3];
rz(-1.72662) q[3];
sx q[3];
rz(1.5161878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93306142) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(2.9053524) q[2];
rz(-1.158372) q[3];
sx q[3];
rz(-2.1393496) q[3];
sx q[3];
rz(-2.2896144) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1119969) q[0];
sx q[0];
rz(-1.010226) q[0];
sx q[0];
rz(2.239256) q[0];
rz(-2.2205655) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(0.62228084) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.92684) q[0];
sx q[0];
rz(-1.8536957) q[0];
sx q[0];
rz(0.19006417) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12139608) q[2];
sx q[2];
rz(-2.4897794) q[2];
sx q[2];
rz(2.7898942) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2942549) q[1];
sx q[1];
rz(-1.6604742) q[1];
sx q[1];
rz(1.4011821) q[1];
rz(-pi) q[2];
rz(1.9641818) q[3];
sx q[3];
rz(-1.1227566) q[3];
sx q[3];
rz(-3.0569227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2174012) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(3.1203111) q[2];
rz(-1.4422669) q[3];
sx q[3];
rz(-0.71443021) q[3];
sx q[3];
rz(-2.2309979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845881) q[0];
sx q[0];
rz(-2.5452884) q[0];
sx q[0];
rz(3.0969627) q[0];
rz(-1.0021098) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(-0.034428509) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1761985) q[0];
sx q[0];
rz(-1.5571496) q[0];
sx q[0];
rz(-1.3854881) q[0];
rz(-pi) q[1];
rz(0.12182932) q[2];
sx q[2];
rz(-1.0708772) q[2];
sx q[2];
rz(0.027539754) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.039636314) q[1];
sx q[1];
rz(-0.64412457) q[1];
sx q[1];
rz(2.1307751) q[1];
x q[2];
rz(-1.4061635) q[3];
sx q[3];
rz(-2.2693686) q[3];
sx q[3];
rz(0.64100953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.78952152) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(-2.1748523) q[2];
rz(0.012185193) q[3];
sx q[3];
rz(-0.92958486) q[3];
sx q[3];
rz(3.0325082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1352585) q[0];
sx q[0];
rz(-0.2671347) q[0];
sx q[0];
rz(2.1475041) q[0];
rz(-3.0864691) q[1];
sx q[1];
rz(-1.8693285) q[1];
sx q[1];
rz(0.73928839) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8871043) q[0];
sx q[0];
rz(-1.4758037) q[0];
sx q[0];
rz(1.4332921) q[0];
rz(-pi) q[1];
rz(-1.0871068) q[2];
sx q[2];
rz(-1.0488044) q[2];
sx q[2];
rz(-0.99624485) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1623605) q[1];
sx q[1];
rz(-1.7602966) q[1];
sx q[1];
rz(0.53342553) q[1];
x q[2];
rz(1.2823295) q[3];
sx q[3];
rz(-0.97086421) q[3];
sx q[3];
rz(2.695042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.58549515) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(-0.56376702) q[2];
rz(2.901315) q[3];
sx q[3];
rz(-0.53301817) q[3];
sx q[3];
rz(1.7668004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4450842) q[0];
sx q[0];
rz(-0.17906469) q[0];
sx q[0];
rz(-0.58404303) q[0];
rz(0.42075992) q[1];
sx q[1];
rz(-1.0412443) q[1];
sx q[1];
rz(-0.27841321) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95094341) q[0];
sx q[0];
rz(-1.6961432) q[0];
sx q[0];
rz(-0.9992674) q[0];
x q[1];
rz(-0.43012302) q[2];
sx q[2];
rz(-0.42818907) q[2];
sx q[2];
rz(0.45454121) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.003493017) q[1];
sx q[1];
rz(-1.3521191) q[1];
sx q[1];
rz(1.2117282) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1182546) q[3];
sx q[3];
rz(-1.7806782) q[3];
sx q[3];
rz(-1.7367712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.86928308) q[2];
sx q[2];
rz(-2.3217106) q[2];
sx q[2];
rz(0.38796866) q[2];
rz(-1.7913473) q[3];
sx q[3];
rz(-0.46348059) q[3];
sx q[3];
rz(-0.2750245) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85912722) q[0];
sx q[0];
rz(-2.9282741) q[0];
sx q[0];
rz(-0.7318837) q[0];
rz(2.9751119) q[1];
sx q[1];
rz(-2.6961168) q[1];
sx q[1];
rz(-3.0678715) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026422231) q[0];
sx q[0];
rz(-2.0048884) q[0];
sx q[0];
rz(3.0300481) q[0];
rz(0.080967112) q[2];
sx q[2];
rz(-1.4514187) q[2];
sx q[2];
rz(2.2031684) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8970866) q[1];
sx q[1];
rz(-0.58468854) q[1];
sx q[1];
rz(0.067824407) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1707465) q[3];
sx q[3];
rz(-0.57983825) q[3];
sx q[3];
rz(-0.29614007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.59447294) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(0.42868844) q[2];
rz(1.3171014) q[3];
sx q[3];
rz(-1.7842112) q[3];
sx q[3];
rz(1.930442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73228943) q[0];
sx q[0];
rz(-1.6560873) q[0];
sx q[0];
rz(2.0987341) q[0];
rz(1.5111142) q[1];
sx q[1];
rz(-1.7621367) q[1];
sx q[1];
rz(-1.3483378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3955584) q[0];
sx q[0];
rz(-1.6840044) q[0];
sx q[0];
rz(3.1157225) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.83421631) q[2];
sx q[2];
rz(-1.5110821) q[2];
sx q[2];
rz(0.63876736) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.18689449) q[1];
sx q[1];
rz(-0.3317301) q[1];
sx q[1];
rz(-1.7897254) q[1];
x q[2];
rz(2.4753307) q[3];
sx q[3];
rz(-2.8386335) q[3];
sx q[3];
rz(-3.016824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3412791) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(0.094816118) q[2];
rz(-1.9684277) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(0.15869424) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.512758) q[0];
sx q[0];
rz(-0.47369581) q[0];
sx q[0];
rz(1.0439903) q[0];
rz(-1.6013153) q[1];
sx q[1];
rz(-1.5581144) q[1];
sx q[1];
rz(1.5099572) q[1];
rz(-3.0586254) q[2];
sx q[2];
rz(-1.7832179) q[2];
sx q[2];
rz(-0.28945343) q[2];
rz(-2.1223162) q[3];
sx q[3];
rz(-0.84001361) q[3];
sx q[3];
rz(2.7915814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
