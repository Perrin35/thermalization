OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3044843) q[0];
sx q[0];
rz(-1.6882856) q[0];
sx q[0];
rz(-0.31153554) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(-1.3181926) q[1];
sx q[1];
rz(-0.55895609) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1564724) q[0];
sx q[0];
rz(-0.44056842) q[0];
sx q[0];
rz(1.23929) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5582325) q[2];
sx q[2];
rz(-2.1280834) q[2];
sx q[2];
rz(-1.81665) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6701339) q[1];
sx q[1];
rz(-1.0308627) q[1];
sx q[1];
rz(0.88175168) q[1];
rz(-pi) q[2];
rz(1.2851089) q[3];
sx q[3];
rz(-1.1477071) q[3];
sx q[3];
rz(0.74564122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7493593) q[2];
sx q[2];
rz(-1.2831251) q[2];
sx q[2];
rz(0.63670811) q[2];
rz(-0.84896815) q[3];
sx q[3];
rz(-0.62148062) q[3];
sx q[3];
rz(-0.23392114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.7648776) q[0];
sx q[0];
rz(-2.8945518) q[0];
sx q[0];
rz(2.9887181) q[0];
rz(0.75694594) q[1];
sx q[1];
rz(-1.5544954) q[1];
sx q[1];
rz(-2.1551932) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7727535) q[0];
sx q[0];
rz(-1.6110238) q[0];
sx q[0];
rz(-0.059779151) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.37843) q[2];
sx q[2];
rz(-1.9689416) q[2];
sx q[2];
rz(-0.60000186) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6978554) q[1];
sx q[1];
rz(-1.9471696) q[1];
sx q[1];
rz(1.8259551) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9513449) q[3];
sx q[3];
rz(-2.7893587) q[3];
sx q[3];
rz(-0.16570839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5793005) q[2];
sx q[2];
rz(-1.9220756) q[2];
sx q[2];
rz(2.3584649) q[2];
rz(-0.018571818) q[3];
sx q[3];
rz(-1.637807) q[3];
sx q[3];
rz(0.40772453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0531533) q[0];
sx q[0];
rz(-2.8494371) q[0];
sx q[0];
rz(-2.1799178) q[0];
rz(-2.7812474) q[1];
sx q[1];
rz(-1.1018437) q[1];
sx q[1];
rz(-0.12869421) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51552396) q[0];
sx q[0];
rz(-1.6532073) q[0];
sx q[0];
rz(-0.047973085) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53993291) q[2];
sx q[2];
rz(-2.2579102) q[2];
sx q[2];
rz(0.5772669) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.62540185) q[1];
sx q[1];
rz(-1.1266202) q[1];
sx q[1];
rz(1.8810012) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4314753) q[3];
sx q[3];
rz(-0.77292597) q[3];
sx q[3];
rz(-1.7804002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0814357) q[2];
sx q[2];
rz(-1.9443941) q[2];
sx q[2];
rz(-1.241768) q[2];
rz(-0.5870108) q[3];
sx q[3];
rz(-0.95654064) q[3];
sx q[3];
rz(-0.97755066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.509165) q[0];
sx q[0];
rz(-0.88212633) q[0];
sx q[0];
rz(-2.0571016) q[0];
rz(-1.658461) q[1];
sx q[1];
rz(-2.5741534) q[1];
sx q[1];
rz(0.09253563) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10034427) q[0];
sx q[0];
rz(-2.7462602) q[0];
sx q[0];
rz(0.80907099) q[0];
rz(-pi) q[1];
rz(-0.34543583) q[2];
sx q[2];
rz(-1.1159117) q[2];
sx q[2];
rz(2.5330184) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.98574084) q[1];
sx q[1];
rz(-0.34793138) q[1];
sx q[1];
rz(2.9702529) q[1];
rz(-pi) q[2];
rz(0.16102287) q[3];
sx q[3];
rz(-2.2194214) q[3];
sx q[3];
rz(-0.30740689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3080421) q[2];
sx q[2];
rz(-1.7001067) q[2];
sx q[2];
rz(2.8095424) q[2];
rz(1.0559233) q[3];
sx q[3];
rz(-2.8639586) q[3];
sx q[3];
rz(0.61029303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32245359) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(0.21155587) q[0];
rz(1.8353204) q[1];
sx q[1];
rz(-1.2439367) q[1];
sx q[1];
rz(-2.4938915) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0928597) q[0];
sx q[0];
rz(-1.6824241) q[0];
sx q[0];
rz(2.9812921) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7251882) q[2];
sx q[2];
rz(-1.9152181) q[2];
sx q[2];
rz(-1.7248578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9808637) q[1];
sx q[1];
rz(-0.54667066) q[1];
sx q[1];
rz(1.1854118) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1224498) q[3];
sx q[3];
rz(-1.585841) q[3];
sx q[3];
rz(2.6263833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.56090474) q[2];
sx q[2];
rz(-0.40955341) q[2];
sx q[2];
rz(2.4482751) q[2];
rz(0.66926113) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(-3.1366689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82419056) q[0];
sx q[0];
rz(-2.9142002) q[0];
sx q[0];
rz(-1.90907) q[0];
rz(-1.0725853) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(-0.17428621) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5151383) q[0];
sx q[0];
rz(-0.8695375) q[0];
sx q[0];
rz(-0.40909543) q[0];
rz(-pi) q[1];
rz(-2.3107489) q[2];
sx q[2];
rz(-2.2846662) q[2];
sx q[2];
rz(1.6531528) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7277158) q[1];
sx q[1];
rz(-1.1912279) q[1];
sx q[1];
rz(-0.80645251) q[1];
rz(1.8829846) q[3];
sx q[3];
rz(-1.6815261) q[3];
sx q[3];
rz(-2.1544416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3198513) q[2];
sx q[2];
rz(-2.3345626) q[2];
sx q[2];
rz(2.9439587) q[2];
rz(2.8526784) q[3];
sx q[3];
rz(-2.2646326) q[3];
sx q[3];
rz(-1.4060085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.063868) q[0];
sx q[0];
rz(-0.48982319) q[0];
sx q[0];
rz(2.9329964) q[0];
rz(-2.1754307) q[1];
sx q[1];
rz(-1.1599133) q[1];
sx q[1];
rz(-1.6360412) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9150328) q[0];
sx q[0];
rz(-2.0953) q[0];
sx q[0];
rz(-2.4302308) q[0];
rz(-pi) q[1];
rz(-2.7111972) q[2];
sx q[2];
rz(-1.6039404) q[2];
sx q[2];
rz(2.7660649) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70430763) q[1];
sx q[1];
rz(-1.1204801) q[1];
sx q[1];
rz(-0.98547658) q[1];
rz(2.8553477) q[3];
sx q[3];
rz(-1.6702594) q[3];
sx q[3];
rz(3.0486097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.85764) q[2];
sx q[2];
rz(-2.3985034) q[2];
sx q[2];
rz(3.0440142) q[2];
rz(1.3939259) q[3];
sx q[3];
rz(-1.3503617) q[3];
sx q[3];
rz(-2.424749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7512648) q[0];
sx q[0];
rz(-1.3289691) q[0];
sx q[0];
rz(0.51399291) q[0];
rz(-0.12318525) q[1];
sx q[1];
rz(-2.8942278) q[1];
sx q[1];
rz(0.93200144) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4829464) q[0];
sx q[0];
rz(-2.6771149) q[0];
sx q[0];
rz(-1.8770201) q[0];
x q[1];
rz(-1.9281689) q[2];
sx q[2];
rz(-2.4744611) q[2];
sx q[2];
rz(0.92598976) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0814221) q[1];
sx q[1];
rz(-2.0875071) q[1];
sx q[1];
rz(0.64232773) q[1];
x q[2];
rz(0.85429116) q[3];
sx q[3];
rz(-1.4201418) q[3];
sx q[3];
rz(0.52779576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1121858) q[2];
sx q[2];
rz(-2.0307348) q[2];
sx q[2];
rz(-2.712148) q[2];
rz(-1.2094234) q[3];
sx q[3];
rz(-0.37241396) q[3];
sx q[3];
rz(-2.4485574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3747303) q[0];
sx q[0];
rz(-1.466789) q[0];
sx q[0];
rz(1.8027579) q[0];
rz(-2.4354637) q[1];
sx q[1];
rz(-1.8920205) q[1];
sx q[1];
rz(0.95058092) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7079733) q[0];
sx q[0];
rz(-0.74445671) q[0];
sx q[0];
rz(-0.77197335) q[0];
rz(1.9865932) q[2];
sx q[2];
rz(-1.6396513) q[2];
sx q[2];
rz(1.6378251) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3705759) q[1];
sx q[1];
rz(-2.5713213) q[1];
sx q[1];
rz(-1.9221406) q[1];
rz(-pi) q[2];
rz(-1.1589963) q[3];
sx q[3];
rz(-1.8048865) q[3];
sx q[3];
rz(1.4300508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6616228) q[2];
sx q[2];
rz(-1.6500436) q[2];
sx q[2];
rz(2.6573112) q[2];
rz(0.92710036) q[3];
sx q[3];
rz(-1.8323332) q[3];
sx q[3];
rz(-0.62121975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1442239) q[0];
sx q[0];
rz(-0.08865083) q[0];
sx q[0];
rz(2.9123059) q[0];
rz(2.7067822) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(0.71892175) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9839448) q[0];
sx q[0];
rz(-2.5074208) q[0];
sx q[0];
rz(1.4112524) q[0];
rz(-3.0232593) q[2];
sx q[2];
rz(-2.1552857) q[2];
sx q[2];
rz(-3.0362533) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6112411) q[1];
sx q[1];
rz(-1.2019005) q[1];
sx q[1];
rz(-1.4126652) q[1];
rz(-pi) q[2];
rz(-0.94902456) q[3];
sx q[3];
rz(-1.9927295) q[3];
sx q[3];
rz(2.4728647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3832613) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(2.3948005) q[2];
rz(-0.87219277) q[3];
sx q[3];
rz(-0.82834297) q[3];
sx q[3];
rz(-0.94223589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(2.7643798) q[1];
sx q[1];
rz(-1.4706392) q[1];
sx q[1];
rz(-2.6249862) q[1];
rz(0.51858356) q[2];
sx q[2];
rz(-1.8443783) q[2];
sx q[2];
rz(1.5658866) q[2];
rz(-1.4217581) q[3];
sx q[3];
rz(-2.2825713) q[3];
sx q[3];
rz(1.5406516) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
