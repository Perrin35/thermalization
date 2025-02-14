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
rz(0.67233664) q[0];
sx q[0];
rz(5.0007102) q[0];
sx q[0];
rz(11.034427) q[0];
rz(1.2338282) q[1];
sx q[1];
rz(-1.9246074) q[1];
sx q[1];
rz(1.6977065) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1579698) q[0];
sx q[0];
rz(-2.2394272) q[0];
sx q[0];
rz(3.0374531) q[0];
rz(-2.2488689) q[2];
sx q[2];
rz(-0.44995445) q[2];
sx q[2];
rz(1.8261248) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.083084785) q[1];
sx q[1];
rz(-2.1821686) q[1];
sx q[1];
rz(-1.0669003) q[1];
x q[2];
rz(-1.9099376) q[3];
sx q[3];
rz(-1.5445441) q[3];
sx q[3];
rz(1.3052125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3621138) q[2];
sx q[2];
rz(-1.9356091) q[2];
sx q[2];
rz(-0.77902478) q[2];
rz(-1.3075167) q[3];
sx q[3];
rz(-2.8179759) q[3];
sx q[3];
rz(1.754508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5883412) q[0];
sx q[0];
rz(-1.3655751) q[0];
sx q[0];
rz(2.8643082) q[0];
rz(-0.39018997) q[1];
sx q[1];
rz(-2.0398102) q[1];
sx q[1];
rz(2.0139093) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71806112) q[0];
sx q[0];
rz(-3.09457) q[0];
sx q[0];
rz(2.7014947) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28553005) q[2];
sx q[2];
rz(-1.4737616) q[2];
sx q[2];
rz(2.3842528) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8133328) q[1];
sx q[1];
rz(-1.6333123) q[1];
sx q[1];
rz(-1.8429788) q[1];
rz(-pi) q[2];
rz(0.66172285) q[3];
sx q[3];
rz(-1.8764917) q[3];
sx q[3];
rz(2.2591285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4358501) q[2];
sx q[2];
rz(-0.59429344) q[2];
sx q[2];
rz(1.1951813) q[2];
rz(-2.4863906) q[3];
sx q[3];
rz(-1.6569258) q[3];
sx q[3];
rz(-0.5932194) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6708267) q[0];
sx q[0];
rz(-2.4556181) q[0];
sx q[0];
rz(2.2210448) q[0];
rz(1.2819194) q[1];
sx q[1];
rz(-2.0125407) q[1];
sx q[1];
rz(1.4528073) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1143415) q[0];
sx q[0];
rz(-0.14273164) q[0];
sx q[0];
rz(0.65988512) q[0];
rz(-pi) q[1];
rz(1.0017582) q[2];
sx q[2];
rz(-1.3033221) q[2];
sx q[2];
rz(-0.75050577) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6026028) q[1];
sx q[1];
rz(-0.33947152) q[1];
sx q[1];
rz(0.15442876) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7887855) q[3];
sx q[3];
rz(-3.0485287) q[3];
sx q[3];
rz(-2.399202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.61536962) q[2];
sx q[2];
rz(-1.0542032) q[2];
sx q[2];
rz(0.39895454) q[2];
rz(0.51359549) q[3];
sx q[3];
rz(-0.31702888) q[3];
sx q[3];
rz(1.1494466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32617635) q[0];
sx q[0];
rz(-0.57891095) q[0];
sx q[0];
rz(-1.2166566) q[0];
rz(0.51219621) q[1];
sx q[1];
rz(-1.1108111) q[1];
sx q[1];
rz(1.1042575) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51795635) q[0];
sx q[0];
rz(-0.41455634) q[0];
sx q[0];
rz(2.3974933) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3338564) q[2];
sx q[2];
rz(-2.8931354) q[2];
sx q[2];
rz(-0.31794869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7514551) q[1];
sx q[1];
rz(-1.5105167) q[1];
sx q[1];
rz(1.3740609) q[1];
x q[2];
rz(-1.8301303) q[3];
sx q[3];
rz(-1.7899198) q[3];
sx q[3];
rz(-2.863209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5964552) q[2];
sx q[2];
rz(-1.9497654) q[2];
sx q[2];
rz(-0.54984251) q[2];
rz(-0.95692974) q[3];
sx q[3];
rz(-2.3423829) q[3];
sx q[3];
rz(-1.6526615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6560219) q[0];
sx q[0];
rz(-0.85000426) q[0];
sx q[0];
rz(-2.3824084) q[0];
rz(-2.1929269) q[1];
sx q[1];
rz(-1.9976839) q[1];
sx q[1];
rz(-0.33160981) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7973855) q[0];
sx q[0];
rz(-1.8669882) q[0];
sx q[0];
rz(-2.3621817) q[0];
x q[1];
rz(-2.6668915) q[2];
sx q[2];
rz(-1.5384288) q[2];
sx q[2];
rz(-2.2970478) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3892203) q[1];
sx q[1];
rz(-1.3942317) q[1];
sx q[1];
rz(-2.9851488) q[1];
rz(-2.7301627) q[3];
sx q[3];
rz(-1.6959018) q[3];
sx q[3];
rz(2.8730212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9552976) q[2];
sx q[2];
rz(-2.1897327) q[2];
sx q[2];
rz(-2.1985506) q[2];
rz(2.9694563) q[3];
sx q[3];
rz(-0.94477263) q[3];
sx q[3];
rz(-2.9183563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73582369) q[0];
sx q[0];
rz(-2.3961234) q[0];
sx q[0];
rz(-1.5752342) q[0];
rz(2.9834421) q[1];
sx q[1];
rz(-0.76845303) q[1];
sx q[1];
rz(-2.2134773) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0029308) q[0];
sx q[0];
rz(-1.2671906) q[0];
sx q[0];
rz(-1.4360484) q[0];
x q[1];
rz(-2.1989397) q[2];
sx q[2];
rz(-2.4898875) q[2];
sx q[2];
rz(-2.556837) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5852768) q[1];
sx q[1];
rz(-2.2158874) q[1];
sx q[1];
rz(0.32151244) q[1];
rz(-0.32676231) q[3];
sx q[3];
rz(-2.770677) q[3];
sx q[3];
rz(-2.4075631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1930262) q[2];
sx q[2];
rz(-2.8552449) q[2];
sx q[2];
rz(2.7890653) q[2];
rz(-2.7726717) q[3];
sx q[3];
rz(-2.5589294) q[3];
sx q[3];
rz(-2.2841456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7848659) q[0];
sx q[0];
rz(-3.1251188) q[0];
sx q[0];
rz(2.4705868) q[0];
rz(-0.10637936) q[1];
sx q[1];
rz(-0.89768592) q[1];
sx q[1];
rz(-2.5118714) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8170057) q[0];
sx q[0];
rz(-1.3793077) q[0];
sx q[0];
rz(-1.5517439) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3307569) q[2];
sx q[2];
rz(-2.5784284) q[2];
sx q[2];
rz(2.882618) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2196673) q[1];
sx q[1];
rz(-1.6451854) q[1];
sx q[1];
rz(-0.32385918) q[1];
x q[2];
rz(-1.2604146) q[3];
sx q[3];
rz(-0.93573007) q[3];
sx q[3];
rz(2.6657651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.070039198) q[2];
sx q[2];
rz(-0.89954251) q[2];
sx q[2];
rz(1.0413991) q[2];
rz(1.533016) q[3];
sx q[3];
rz(-0.93419111) q[3];
sx q[3];
rz(-1.1478724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94569412) q[0];
sx q[0];
rz(-2.4905289) q[0];
sx q[0];
rz(2.661929) q[0];
rz(3.116963) q[1];
sx q[1];
rz(-1.004091) q[1];
sx q[1];
rz(-3.0020795) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5845156) q[0];
sx q[0];
rz(-0.13894835) q[0];
sx q[0];
rz(-0.90451972) q[0];
x q[1];
rz(-1.0647667) q[2];
sx q[2];
rz(-1.2373018) q[2];
sx q[2];
rz(0.34766742) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.3132726) q[1];
sx q[1];
rz(-1.7370151) q[1];
sx q[1];
rz(2.4795824) q[1];
x q[2];
rz(-2.471227) q[3];
sx q[3];
rz(-1.6135525) q[3];
sx q[3];
rz(-1.6898516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2554243) q[2];
sx q[2];
rz(-1.5972127) q[2];
sx q[2];
rz(1.1523979) q[2];
rz(-1.7216916) q[3];
sx q[3];
rz(-2.4072188) q[3];
sx q[3];
rz(-1.3091807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.060147978) q[0];
sx q[0];
rz(-1.6763433) q[0];
sx q[0];
rz(1.4981221) q[0];
rz(0.76088798) q[1];
sx q[1];
rz(-2.1983169) q[1];
sx q[1];
rz(1.6852185) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24600226) q[0];
sx q[0];
rz(-2.1293679) q[0];
sx q[0];
rz(1.3590525) q[0];
rz(-pi) q[1];
rz(2.6713702) q[2];
sx q[2];
rz(-2.6986327) q[2];
sx q[2];
rz(2.265919) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.93865004) q[1];
sx q[1];
rz(-0.65697815) q[1];
sx q[1];
rz(0.88969661) q[1];
rz(-0.96181576) q[3];
sx q[3];
rz(-0.66003335) q[3];
sx q[3];
rz(0.48840085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3219519) q[2];
sx q[2];
rz(-2.7072622) q[2];
sx q[2];
rz(-2.1257909) q[2];
rz(0.89456931) q[3];
sx q[3];
rz(-1.2457448) q[3];
sx q[3];
rz(2.3999124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9858518) q[0];
sx q[0];
rz(-2.4116801) q[0];
sx q[0];
rz(0.99005449) q[0];
rz(2.2186225) q[1];
sx q[1];
rz(-1.3786517) q[1];
sx q[1];
rz(-2.1754481) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9618398) q[0];
sx q[0];
rz(-1.5003164) q[0];
sx q[0];
rz(-1.671179) q[0];
rz(-pi) q[1];
rz(-0.82426333) q[2];
sx q[2];
rz(-2.5474862) q[2];
sx q[2];
rz(2.4115393) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6539972) q[1];
sx q[1];
rz(-1.8593809) q[1];
sx q[1];
rz(1.1638132) q[1];
rz(3.0279492) q[3];
sx q[3];
rz(-1.3727649) q[3];
sx q[3];
rz(0.5403434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3155589) q[2];
sx q[2];
rz(-0.78353271) q[2];
sx q[2];
rz(-2.3890736) q[2];
rz(0.13713947) q[3];
sx q[3];
rz(-2.5004041) q[3];
sx q[3];
rz(2.7916059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.242545) q[0];
sx q[0];
rz(-0.62472961) q[0];
sx q[0];
rz(-0.20150264) q[0];
rz(1.3626199) q[1];
sx q[1];
rz(-2.2687804) q[1];
sx q[1];
rz(-0.16558095) q[1];
rz(1.4819862) q[2];
sx q[2];
rz(-1.6203151) q[2];
sx q[2];
rz(2.0697894) q[2];
rz(2.7794502) q[3];
sx q[3];
rz(-1.1884304) q[3];
sx q[3];
rz(-1.8588369) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
