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
rz(1.4662161) q[0];
sx q[0];
rz(-0.94453064) q[0];
sx q[0];
rz(0.20139995) q[0];
rz(-2.0693076) q[1];
sx q[1];
rz(-2.3787002) q[1];
sx q[1];
rz(1.720517) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.759255) q[0];
sx q[0];
rz(-1.6256285) q[0];
sx q[0];
rz(-2.1817529) q[0];
rz(-pi) q[1];
rz(-2.6236963) q[2];
sx q[2];
rz(-0.96518789) q[2];
sx q[2];
rz(-2.6014858) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.046467543) q[1];
sx q[1];
rz(-2.3541942) q[1];
sx q[1];
rz(2.5745029) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2675377) q[3];
sx q[3];
rz(-0.49421453) q[3];
sx q[3];
rz(2.6074802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1066771) q[2];
sx q[2];
rz(-0.80433977) q[2];
sx q[2];
rz(0.2790645) q[2];
rz(-1.8883102) q[3];
sx q[3];
rz(-0.24975714) q[3];
sx q[3];
rz(-2.9899924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6107553) q[0];
sx q[0];
rz(-0.84722561) q[0];
sx q[0];
rz(2.8952059) q[0];
rz(1.9465744) q[1];
sx q[1];
rz(-0.89391005) q[1];
sx q[1];
rz(-1.6349207) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89013571) q[0];
sx q[0];
rz(-1.1201753) q[0];
sx q[0];
rz(-2.8681953) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8737239) q[2];
sx q[2];
rz(-1.611735) q[2];
sx q[2];
rz(-2.8589279) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.40827258) q[1];
sx q[1];
rz(-0.36767861) q[1];
sx q[1];
rz(-1.4094844) q[1];
x q[2];
rz(1.0492658) q[3];
sx q[3];
rz(-1.9535011) q[3];
sx q[3];
rz(2.692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.40111497) q[2];
sx q[2];
rz(-0.97492188) q[2];
sx q[2];
rz(1.231989) q[2];
rz(-1.1022107) q[3];
sx q[3];
rz(-0.004318459) q[3];
sx q[3];
rz(-1.8050885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6983637) q[0];
sx q[0];
rz(-2.4446428) q[0];
sx q[0];
rz(0.90782905) q[0];
rz(-0.56697956) q[1];
sx q[1];
rz(-2.1588529) q[1];
sx q[1];
rz(-3.0388015) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11353569) q[0];
sx q[0];
rz(-0.60657078) q[0];
sx q[0];
rz(2.7675178) q[0];
rz(-pi) q[1];
rz(0.081016531) q[2];
sx q[2];
rz(-1.4668494) q[2];
sx q[2];
rz(-0.31652094) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56270617) q[1];
sx q[1];
rz(-1.9745312) q[1];
sx q[1];
rz(-2.5078234) q[1];
rz(-pi) q[2];
rz(-1.0166753) q[3];
sx q[3];
rz(-1.301487) q[3];
sx q[3];
rz(2.179972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8583782) q[2];
sx q[2];
rz(-2.3329222) q[2];
sx q[2];
rz(1.7041448) q[2];
rz(0.39250675) q[3];
sx q[3];
rz(-1.1350574) q[3];
sx q[3];
rz(-0.2984305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0933541) q[0];
sx q[0];
rz(-1.7788576) q[0];
sx q[0];
rz(1.1887953) q[0];
rz(-2.8554754) q[1];
sx q[1];
rz(-2.5240099) q[1];
sx q[1];
rz(1.5123222) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6962357) q[0];
sx q[0];
rz(-1.0508448) q[0];
sx q[0];
rz(-0.32132863) q[0];
rz(-2.3740923) q[2];
sx q[2];
rz(-1.8165117) q[2];
sx q[2];
rz(1.2188639) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0864073) q[1];
sx q[1];
rz(-0.98208222) q[1];
sx q[1];
rz(2.9179395) q[1];
rz(-pi) q[2];
x q[2];
rz(1.063594) q[3];
sx q[3];
rz(-1.7559057) q[3];
sx q[3];
rz(1.1066268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.35147038) q[2];
sx q[2];
rz(-2.0581547) q[2];
sx q[2];
rz(-0.67698395) q[2];
rz(-1.2466768) q[3];
sx q[3];
rz(-1.2723943) q[3];
sx q[3];
rz(1.5164794) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2077654) q[0];
sx q[0];
rz(-2.8849869) q[0];
sx q[0];
rz(3.1166792) q[0];
rz(-1.0554396) q[1];
sx q[1];
rz(-0.98506227) q[1];
sx q[1];
rz(0.28344646) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058539778) q[0];
sx q[0];
rz(-1.3136567) q[0];
sx q[0];
rz(1.8399946) q[0];
rz(-pi) q[1];
rz(-0.13155664) q[2];
sx q[2];
rz(-1.3765928) q[2];
sx q[2];
rz(2.921791) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1920212) q[1];
sx q[1];
rz(-1.5402435) q[1];
sx q[1];
rz(2.0716689) q[1];
rz(-pi) q[2];
rz(-2.2437566) q[3];
sx q[3];
rz(-1.3011609) q[3];
sx q[3];
rz(1.2854888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8289566) q[2];
sx q[2];
rz(-0.41746155) q[2];
sx q[2];
rz(0.32583315) q[2];
rz(-1.2827986) q[3];
sx q[3];
rz(-1.157434) q[3];
sx q[3];
rz(1.5939943) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7376937) q[0];
sx q[0];
rz(-1.4827381) q[0];
sx q[0];
rz(-0.37495908) q[0];
rz(2.6629579) q[1];
sx q[1];
rz(-2.4691983) q[1];
sx q[1];
rz(2.7714444) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70765342) q[0];
sx q[0];
rz(-1.5824727) q[0];
sx q[0];
rz(3.1184831) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0022718) q[2];
sx q[2];
rz(-1.3874386) q[2];
sx q[2];
rz(-2.574711) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1620336) q[1];
sx q[1];
rz(-1.1104854) q[1];
sx q[1];
rz(-2.3737337) q[1];
rz(-pi) q[2];
rz(-0.010995098) q[3];
sx q[3];
rz(-2.7137626) q[3];
sx q[3];
rz(-1.4233936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7438573) q[2];
sx q[2];
rz(-2.9501259) q[2];
sx q[2];
rz(-1.7605304) q[2];
rz(2.0928275) q[3];
sx q[3];
rz(-1.7966813) q[3];
sx q[3];
rz(2.5019808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25310707) q[0];
sx q[0];
rz(-2.3024004) q[0];
sx q[0];
rz(0.92415586) q[0];
rz(-0.91668516) q[1];
sx q[1];
rz(-1.0662096) q[1];
sx q[1];
rz(1.4013269) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9144672) q[0];
sx q[0];
rz(-1.0773398) q[0];
sx q[0];
rz(-1.6363793) q[0];
x q[1];
rz(-0.22221128) q[2];
sx q[2];
rz(-0.66721081) q[2];
sx q[2];
rz(1.3363584) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8055011) q[1];
sx q[1];
rz(-2.3004043) q[1];
sx q[1];
rz(1.355624) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0823233) q[3];
sx q[3];
rz(-1.3783749) q[3];
sx q[3];
rz(-0.50979641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9132729) q[2];
sx q[2];
rz(-1.4618123) q[2];
sx q[2];
rz(1.224996) q[2];
rz(-3.1332968) q[3];
sx q[3];
rz(-3.1305997) q[3];
sx q[3];
rz(-0.86709658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5834354) q[0];
sx q[0];
rz(-2.3112516) q[0];
sx q[0];
rz(2.661327) q[0];
rz(2.9971314) q[1];
sx q[1];
rz(-2.2339349) q[1];
sx q[1];
rz(0.86437782) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1594161) q[0];
sx q[0];
rz(-1.5242531) q[0];
sx q[0];
rz(-1.4070562) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5369253) q[2];
sx q[2];
rz(-1.4761756) q[2];
sx q[2];
rz(-0.44924863) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0740114) q[1];
sx q[1];
rz(-1.0339104) q[1];
sx q[1];
rz(-2.8606961) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8671667) q[3];
sx q[3];
rz(-1.84019) q[3];
sx q[3];
rz(-1.5116949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9353443) q[2];
sx q[2];
rz(-2.1312921) q[2];
sx q[2];
rz(2.5065191) q[2];
rz(0.57279974) q[3];
sx q[3];
rz(-1.7856995) q[3];
sx q[3];
rz(2.2662381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0315345) q[0];
sx q[0];
rz(-2.3681971) q[0];
sx q[0];
rz(-3.0173259) q[0];
rz(-2.7763413) q[1];
sx q[1];
rz(-1.4271913) q[1];
sx q[1];
rz(-1.431538) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54473684) q[0];
sx q[0];
rz(-1.3787621) q[0];
sx q[0];
rz(1.3903244) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1373491) q[2];
sx q[2];
rz(-0.99663094) q[2];
sx q[2];
rz(-3.1176709) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4540951) q[1];
sx q[1];
rz(-0.75088596) q[1];
sx q[1];
rz(2.6862957) q[1];
rz(-pi) q[2];
rz(1.133119) q[3];
sx q[3];
rz(-1.0854377) q[3];
sx q[3];
rz(2.7818668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8753836) q[2];
sx q[2];
rz(-2.2201846) q[2];
sx q[2];
rz(1.9688152) q[2];
rz(-1.734599) q[3];
sx q[3];
rz(-1.809779) q[3];
sx q[3];
rz(-1.2146568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59239546) q[0];
sx q[0];
rz(-2.0068491) q[0];
sx q[0];
rz(0.59463516) q[0];
rz(1.0058962) q[1];
sx q[1];
rz(-1.9391831) q[1];
sx q[1];
rz(2.2524021) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.206736) q[0];
sx q[0];
rz(-1.6214658) q[0];
sx q[0];
rz(-1.7913179) q[0];
rz(1.4290733) q[2];
sx q[2];
rz(-2.6772039) q[2];
sx q[2];
rz(1.022867) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2500221) q[1];
sx q[1];
rz(-1.6305939) q[1];
sx q[1];
rz(1.476261) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51639207) q[3];
sx q[3];
rz(-1.603873) q[3];
sx q[3];
rz(-0.75504485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1903926) q[2];
sx q[2];
rz(-1.1756281) q[2];
sx q[2];
rz(2.5962043) q[2];
rz(-0.96327463) q[3];
sx q[3];
rz(-1.9656209) q[3];
sx q[3];
rz(1.6776599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.6496898) q[0];
sx q[0];
rz(-2.2408673) q[0];
sx q[0];
rz(1.3229205) q[0];
rz(0.84359618) q[1];
sx q[1];
rz(-2.0633162) q[1];
sx q[1];
rz(1.8819527) q[1];
rz(1.1314992) q[2];
sx q[2];
rz(-1.3901303) q[2];
sx q[2];
rz(0.11700665) q[2];
rz(0.58335494) q[3];
sx q[3];
rz(-2.0045378) q[3];
sx q[3];
rz(-1.0160673) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
