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
rz(-2.1036086) q[0];
sx q[0];
rz(-1.684364) q[0];
sx q[0];
rz(1.3527704) q[0];
rz(-5.0985131) q[1];
sx q[1];
rz(2.4689622) q[1];
sx q[1];
rz(11.050635) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0294581) q[0];
sx q[0];
rz(-1.7687651) q[0];
sx q[0];
rz(-0.1898808) q[0];
x q[1];
rz(0.47907655) q[2];
sx q[2];
rz(-2.276439) q[2];
sx q[2];
rz(2.3912663) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3799073) q[1];
sx q[1];
rz(-2.5406079) q[1];
sx q[1];
rz(0.65324776) q[1];
x q[2];
rz(1.2836841) q[3];
sx q[3];
rz(-1.1653882) q[3];
sx q[3];
rz(1.4472345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0036014) q[2];
sx q[2];
rz(-0.10819745) q[2];
sx q[2];
rz(-0.18599621) q[2];
rz(-0.018639175) q[3];
sx q[3];
rz(-1.2942856) q[3];
sx q[3];
rz(-1.0703465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5828534) q[0];
sx q[0];
rz(-2.5753729) q[0];
sx q[0];
rz(2.8232316) q[0];
rz(-2.5045577) q[1];
sx q[1];
rz(-2.36167) q[1];
sx q[1];
rz(2.6723518) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56457389) q[0];
sx q[0];
rz(-1.7697608) q[0];
sx q[0];
rz(-1.7451203) q[0];
x q[1];
rz(-1.7809243) q[2];
sx q[2];
rz(-0.38345018) q[2];
sx q[2];
rz(2.0166778) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7982895) q[1];
sx q[1];
rz(-1.251601) q[1];
sx q[1];
rz(-1.809824) q[1];
x q[2];
rz(0.0057159609) q[3];
sx q[3];
rz(-2.2607231) q[3];
sx q[3];
rz(-1.5100513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3176754) q[2];
sx q[2];
rz(-2.9313512) q[2];
sx q[2];
rz(-2.4483185) q[2];
rz(-1.3200101) q[3];
sx q[3];
rz(-1.3258679) q[3];
sx q[3];
rz(-1.0953974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1410809) q[0];
sx q[0];
rz(-2.5040099) q[0];
sx q[0];
rz(-0.47955036) q[0];
rz(-2.6374822) q[1];
sx q[1];
rz(-1.9934374) q[1];
sx q[1];
rz(-3.0832916) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6710799) q[0];
sx q[0];
rz(-1.8288178) q[0];
sx q[0];
rz(-0.51768984) q[0];
rz(-pi) q[1];
rz(2.53251) q[2];
sx q[2];
rz(-0.78819617) q[2];
sx q[2];
rz(-2.4969375) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.41861288) q[1];
sx q[1];
rz(-1.449391) q[1];
sx q[1];
rz(-2.1239807) q[1];
rz(-pi) q[2];
rz(-0.4543484) q[3];
sx q[3];
rz(-0.55255167) q[3];
sx q[3];
rz(2.1228028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7551859) q[2];
sx q[2];
rz(-1.9599954) q[2];
sx q[2];
rz(3.1043261) q[2];
rz(1.5197598) q[3];
sx q[3];
rz(-2.683679) q[3];
sx q[3];
rz(0.38145414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.195049) q[0];
sx q[0];
rz(-2.6730972) q[0];
sx q[0];
rz(-0.066548912) q[0];
rz(1.6171803) q[1];
sx q[1];
rz(-1.8981372) q[1];
sx q[1];
rz(2.4066511) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83072216) q[0];
sx q[0];
rz(-1.3390216) q[0];
sx q[0];
rz(0.59344296) q[0];
x q[1];
rz(-0.73012107) q[2];
sx q[2];
rz(-1.7459622) q[2];
sx q[2];
rz(-2.1736682) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3614703) q[1];
sx q[1];
rz(-1.2166942) q[1];
sx q[1];
rz(1.2722871) q[1];
x q[2];
rz(-0.14566497) q[3];
sx q[3];
rz(-2.2861135) q[3];
sx q[3];
rz(-0.13038929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8526326) q[2];
sx q[2];
rz(-1.5776878) q[2];
sx q[2];
rz(0.10870474) q[2];
rz(-2.2147801) q[3];
sx q[3];
rz(-0.97065297) q[3];
sx q[3];
rz(1.577781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0842593) q[0];
sx q[0];
rz(-0.63705343) q[0];
sx q[0];
rz(-2.0507226) q[0];
rz(-2.5690761) q[1];
sx q[1];
rz(-1.1895836) q[1];
sx q[1];
rz(-2.3172839) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.123468) q[0];
sx q[0];
rz(-1.9459007) q[0];
sx q[0];
rz(0.25183046) q[0];
x q[1];
rz(0.9277497) q[2];
sx q[2];
rz(-2.0111736) q[2];
sx q[2];
rz(-1.7945031) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.5368176) q[1];
sx q[1];
rz(-2.1120434) q[1];
sx q[1];
rz(2.6763335) q[1];
x q[2];
rz(2.0044672) q[3];
sx q[3];
rz(-1.1992362) q[3];
sx q[3];
rz(0.56321689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.480964) q[2];
sx q[2];
rz(-2.169256) q[2];
sx q[2];
rz(0.29120293) q[2];
rz(2.5045942) q[3];
sx q[3];
rz(-1.6773418) q[3];
sx q[3];
rz(-1.8891107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7630735) q[0];
sx q[0];
rz(-0.28994361) q[0];
sx q[0];
rz(0.13105233) q[0];
rz(-1.1669) q[1];
sx q[1];
rz(-0.98375541) q[1];
sx q[1];
rz(-0.022620591) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0221828) q[0];
sx q[0];
rz(-1.9227322) q[0];
sx q[0];
rz(-2.5608313) q[0];
rz(-2.9282369) q[2];
sx q[2];
rz(-1.6913026) q[2];
sx q[2];
rz(1.0671187) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2056668) q[1];
sx q[1];
rz(-0.64595157) q[1];
sx q[1];
rz(-1.0420858) q[1];
rz(-pi) q[2];
rz(0.87225391) q[3];
sx q[3];
rz(-2.8917851) q[3];
sx q[3];
rz(1.2605309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.268198) q[2];
sx q[2];
rz(-0.70316535) q[2];
sx q[2];
rz(-2.0568636) q[2];
rz(1.9966513) q[3];
sx q[3];
rz(-1.8549253) q[3];
sx q[3];
rz(-1.0219319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5893843) q[0];
sx q[0];
rz(-1.5188058) q[0];
sx q[0];
rz(2.6680706) q[0];
rz(-1.2208968) q[1];
sx q[1];
rz(-2.7291606) q[1];
sx q[1];
rz(-1.8730877) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3800156) q[0];
sx q[0];
rz(-2.3748739) q[0];
sx q[0];
rz(2.3469427) q[0];
rz(-0.94077295) q[2];
sx q[2];
rz(-2.2029999) q[2];
sx q[2];
rz(0.8559627) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2655609) q[1];
sx q[1];
rz(-0.68194333) q[1];
sx q[1];
rz(2.4443786) q[1];
rz(-1.1348261) q[3];
sx q[3];
rz(-1.4957168) q[3];
sx q[3];
rz(1.0385385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1641757) q[2];
sx q[2];
rz(-1.5722534) q[2];
sx q[2];
rz(-2.8951728) q[2];
rz(0.94270802) q[3];
sx q[3];
rz(-1.0859414) q[3];
sx q[3];
rz(-2.3781618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3281658) q[0];
sx q[0];
rz(-1.270371) q[0];
sx q[0];
rz(3.0810007) q[0];
rz(-1.6391485) q[1];
sx q[1];
rz(-1.7785347) q[1];
sx q[1];
rz(1.5477808) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12822026) q[0];
sx q[0];
rz(-3.1376079) q[0];
sx q[0];
rz(0.26040034) q[0];
x q[1];
rz(-0.28569371) q[2];
sx q[2];
rz(-0.59425747) q[2];
sx q[2];
rz(-0.93462925) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2912078) q[1];
sx q[1];
rz(-2.3602848) q[1];
sx q[1];
rz(-2.2154097) q[1];
rz(-2.4259858) q[3];
sx q[3];
rz(-1.9247492) q[3];
sx q[3];
rz(-2.7277814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8760406) q[2];
sx q[2];
rz(-2.2206842) q[2];
sx q[2];
rz(1.6843686) q[2];
rz(0.11988104) q[3];
sx q[3];
rz(-2.9370152) q[3];
sx q[3];
rz(1.0286819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2940755) q[0];
sx q[0];
rz(-2.1471922) q[0];
sx q[0];
rz(-2.0062398) q[0];
rz(0.10336939) q[1];
sx q[1];
rz(-1.7366948) q[1];
sx q[1];
rz(-2.1939383) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97449873) q[0];
sx q[0];
rz(-1.5038953) q[0];
sx q[0];
rz(-0.6386021) q[0];
rz(-pi) q[1];
rz(0.88236188) q[2];
sx q[2];
rz(-1.8970006) q[2];
sx q[2];
rz(0.38089124) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7458742) q[1];
sx q[1];
rz(-1.304879) q[1];
sx q[1];
rz(2.2661792) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2470948) q[3];
sx q[3];
rz(-0.94493659) q[3];
sx q[3];
rz(-1.8602399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.4244298) q[2];
sx q[2];
rz(-2.5134176) q[2];
sx q[2];
rz(-0.90432811) q[2];
rz(1.2633911) q[3];
sx q[3];
rz(-1.2369913) q[3];
sx q[3];
rz(0.5164856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.271027) q[0];
sx q[0];
rz(-0.75815433) q[0];
sx q[0];
rz(-2.1562321) q[0];
rz(-0.35161463) q[1];
sx q[1];
rz(-0.96013394) q[1];
sx q[1];
rz(0.34686372) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0766633) q[0];
sx q[0];
rz(-0.67333013) q[0];
sx q[0];
rz(-1.2162186) q[0];
rz(-1.1248671) q[2];
sx q[2];
rz(-1.6724148) q[2];
sx q[2];
rz(-1.9856688) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5849196) q[1];
sx q[1];
rz(-0.69603633) q[1];
sx q[1];
rz(3.0026376) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7010491) q[3];
sx q[3];
rz(-1.9048573) q[3];
sx q[3];
rz(-1.7831217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2615307) q[2];
sx q[2];
rz(-0.71892771) q[2];
sx q[2];
rz(-0.13492179) q[2];
rz(-0.12923446) q[3];
sx q[3];
rz(-2.7415469) q[3];
sx q[3];
rz(-1.038704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1226817) q[0];
sx q[0];
rz(-1.2553348) q[0];
sx q[0];
rz(-3.0154764) q[0];
rz(-2.6799754) q[1];
sx q[1];
rz(-1.4001662) q[1];
sx q[1];
rz(-2.9495159) q[1];
rz(0.89712894) q[2];
sx q[2];
rz(-1.9403602) q[2];
sx q[2];
rz(-0.42862949) q[2];
rz(2.5712874) q[3];
sx q[3];
rz(-0.30443301) q[3];
sx q[3];
rz(-2.6058578) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
