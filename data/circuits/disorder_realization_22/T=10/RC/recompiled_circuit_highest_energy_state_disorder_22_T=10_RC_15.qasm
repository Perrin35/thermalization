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
rz(-2.469256) q[0];
sx q[0];
rz(-1.8591175) q[0];
sx q[0];
rz(-1.6096492) q[0];
rz(1.2338282) q[1];
sx q[1];
rz(-1.9246074) q[1];
sx q[1];
rz(1.6977065) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47752781) q[0];
sx q[0];
rz(-1.6524547) q[0];
sx q[0];
rz(0.89952138) q[0];
rz(-0.89272372) q[2];
sx q[2];
rz(-2.6916382) q[2];
sx q[2];
rz(1.8261248) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4593764) q[1];
sx q[1];
rz(-0.77118783) q[1];
sx q[1];
rz(0.60318635) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6495632) q[3];
sx q[3];
rz(-0.34011671) q[3];
sx q[3];
rz(-0.3398557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3621138) q[2];
sx q[2];
rz(-1.2059836) q[2];
sx q[2];
rz(0.77902478) q[2];
rz(1.834076) q[3];
sx q[3];
rz(-2.8179759) q[3];
sx q[3];
rz(-1.3870846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55325145) q[0];
sx q[0];
rz(-1.7760176) q[0];
sx q[0];
rz(0.27728444) q[0];
rz(-2.7514027) q[1];
sx q[1];
rz(-2.0398102) q[1];
sx q[1];
rz(-2.0139093) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1585856) q[0];
sx q[0];
rz(-1.6133353) q[0];
sx q[0];
rz(-1.5507513) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28553005) q[2];
sx q[2];
rz(-1.6678311) q[2];
sx q[2];
rz(-2.3842528) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8816205) q[1];
sx q[1];
rz(-1.8424336) q[1];
sx q[1];
rz(0.064898811) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9513177) q[3];
sx q[3];
rz(-2.1968958) q[3];
sx q[3];
rz(-2.2230119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7057425) q[2];
sx q[2];
rz(-2.5472992) q[2];
sx q[2];
rz(-1.1951813) q[2];
rz(-2.4863906) q[3];
sx q[3];
rz(-1.6569258) q[3];
sx q[3];
rz(2.5483733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6708267) q[0];
sx q[0];
rz(-2.4556181) q[0];
sx q[0];
rz(0.92054787) q[0];
rz(-1.2819194) q[1];
sx q[1];
rz(-2.0125407) q[1];
sx q[1];
rz(1.6887853) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11139599) q[0];
sx q[0];
rz(-1.483484) q[0];
sx q[0];
rz(-3.0285378) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0412943) q[2];
sx q[2];
rz(-0.62242939) q[2];
sx q[2];
rz(-0.4284455) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7662168) q[1];
sx q[1];
rz(-1.2355243) q[1];
sx q[1];
rz(1.6250618) q[1];
rz(-1.6369197) q[3];
sx q[3];
rz(-1.6363314) q[3];
sx q[3];
rz(1.5333444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.526223) q[2];
sx q[2];
rz(-1.0542032) q[2];
sx q[2];
rz(0.39895454) q[2];
rz(-2.6279972) q[3];
sx q[3];
rz(-0.31702888) q[3];
sx q[3];
rz(-1.992146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32617635) q[0];
sx q[0];
rz(-0.57891095) q[0];
sx q[0];
rz(-1.9249361) q[0];
rz(-2.6293964) q[1];
sx q[1];
rz(-1.1108111) q[1];
sx q[1];
rz(1.1042575) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3062631) q[0];
sx q[0];
rz(-1.2699513) q[0];
sx q[0];
rz(1.2811238) q[0];
rz(-pi) q[1];
rz(-0.80773621) q[2];
sx q[2];
rz(-0.24845727) q[2];
sx q[2];
rz(-0.31794869) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.19266549) q[1];
sx q[1];
rz(-1.3744229) q[1];
sx q[1];
rz(-0.061462359) q[1];
rz(-pi) q[2];
rz(-0.22645183) q[3];
sx q[3];
rz(-1.8237916) q[3];
sx q[3];
rz(1.3500201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5451374) q[2];
sx q[2];
rz(-1.1918273) q[2];
sx q[2];
rz(-2.5917501) q[2];
rz(0.95692974) q[3];
sx q[3];
rz(-0.79920971) q[3];
sx q[3];
rz(1.4889312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4855708) q[0];
sx q[0];
rz(-0.85000426) q[0];
sx q[0];
rz(-0.75918424) q[0];
rz(-0.94866577) q[1];
sx q[1];
rz(-1.1439088) q[1];
sx q[1];
rz(2.8099828) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0873957) q[0];
sx q[0];
rz(-2.3080462) q[0];
sx q[0];
rz(-1.1655318) q[0];
rz(0.47470113) q[2];
sx q[2];
rz(-1.6031638) q[2];
sx q[2];
rz(-0.84454483) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3508671) q[1];
sx q[1];
rz(-1.7247883) q[1];
sx q[1];
rz(-1.3920946) q[1];
rz(0.41142996) q[3];
sx q[3];
rz(-1.4456909) q[3];
sx q[3];
rz(0.26857146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.18629508) q[2];
sx q[2];
rz(-2.1897327) q[2];
sx q[2];
rz(-2.1985506) q[2];
rz(-2.9694563) q[3];
sx q[3];
rz(-2.19682) q[3];
sx q[3];
rz(0.22323639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.405769) q[0];
sx q[0];
rz(-2.3961234) q[0];
sx q[0];
rz(-1.5663585) q[0];
rz(0.15815059) q[1];
sx q[1];
rz(-2.3731396) q[1];
sx q[1];
rz(0.92811531) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52735746) q[0];
sx q[0];
rz(-1.4422461) q[0];
sx q[0];
rz(-2.8353793) q[0];
x q[1];
rz(-2.1989397) q[2];
sx q[2];
rz(-0.65170519) q[2];
sx q[2];
rz(2.556837) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.050452) q[1];
sx q[1];
rz(-0.71040043) q[1];
sx q[1];
rz(-1.968376) q[1];
rz(1.6949878) q[3];
sx q[3];
rz(-1.9212133) q[3];
sx q[3];
rz(-2.7563376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1930262) q[2];
sx q[2];
rz(-0.28634772) q[2];
sx q[2];
rz(2.7890653) q[2];
rz(-2.7726717) q[3];
sx q[3];
rz(-0.58266321) q[3];
sx q[3];
rz(-0.85744706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3567268) q[0];
sx q[0];
rz(-3.1251188) q[0];
sx q[0];
rz(-0.67100588) q[0];
rz(-3.0352133) q[1];
sx q[1];
rz(-2.2439067) q[1];
sx q[1];
rz(-2.5118714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9167921) q[0];
sx q[0];
rz(-0.19242254) q[0];
sx q[0];
rz(3.0436361) q[0];
x q[1];
rz(-1.9999973) q[2];
sx q[2];
rz(-1.1941807) q[2];
sx q[2];
rz(-1.1528328) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2196673) q[1];
sx q[1];
rz(-1.4964073) q[1];
sx q[1];
rz(0.32385918) q[1];
x q[2];
rz(-0.39291556) q[3];
sx q[3];
rz(-0.69732058) q[3];
sx q[3];
rz(-0.97148767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0715535) q[2];
sx q[2];
rz(-2.2420501) q[2];
sx q[2];
rz(2.1001935) q[2];
rz(-1.6085767) q[3];
sx q[3];
rz(-2.2074015) q[3];
sx q[3];
rz(-1.9937203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1958985) q[0];
sx q[0];
rz(-2.4905289) q[0];
sx q[0];
rz(-2.661929) q[0];
rz(3.116963) q[1];
sx q[1];
rz(-1.004091) q[1];
sx q[1];
rz(0.13951313) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91352816) q[0];
sx q[0];
rz(-1.6798927) q[0];
sx q[0];
rz(3.0553711) q[0];
rz(-pi) q[1];
x q[1];
rz(2.191338) q[2];
sx q[2];
rz(-0.59796158) q[2];
sx q[2];
rz(-0.68956748) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0483267) q[1];
sx q[1];
rz(-2.4620926) q[1];
sx q[1];
rz(2.8751586) q[1];
x q[2];
rz(-1.6253396) q[3];
sx q[3];
rz(-2.2404376) q[3];
sx q[3];
rz(-0.085179335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2554243) q[2];
sx q[2];
rz(-1.5972127) q[2];
sx q[2];
rz(-1.9891948) q[2];
rz(-1.7216916) q[3];
sx q[3];
rz(-2.4072188) q[3];
sx q[3];
rz(1.8324119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8955904) q[0];
sx q[0];
rz(-1.0122247) q[0];
sx q[0];
rz(1.7825401) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3590737) q[2];
sx q[2];
rz(-1.9628606) q[2];
sx q[2];
rz(1.388092) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.93865004) q[1];
sx q[1];
rz(-0.65697815) q[1];
sx q[1];
rz(-2.251896) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41784235) q[3];
sx q[3];
rz(-1.0438232) q[3];
sx q[3];
rz(2.9067519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8196408) q[2];
sx q[2];
rz(-0.4343304) q[2];
sx q[2];
rz(2.1257909) q[2];
rz(0.89456931) q[3];
sx q[3];
rz(-1.2457448) q[3];
sx q[3];
rz(-0.7416803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9858518) q[0];
sx q[0];
rz(-0.72991252) q[0];
sx q[0];
rz(0.99005449) q[0];
rz(-2.2186225) q[1];
sx q[1];
rz(-1.3786517) q[1];
sx q[1];
rz(-0.96614456) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9618398) q[0];
sx q[0];
rz(-1.5003164) q[0];
sx q[0];
rz(-1.4704136) q[0];
rz(1.1104634) q[2];
sx q[2];
rz(-1.9607414) q[2];
sx q[2];
rz(-1.4947948) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6412539) q[1];
sx q[1];
rz(-0.49420824) q[1];
sx q[1];
rz(2.2142835) q[1];
rz(-1.3715128) q[3];
sx q[3];
rz(-1.6822094) q[3];
sx q[3];
rz(-2.1335909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.82603377) q[2];
sx q[2];
rz(-2.3580599) q[2];
sx q[2];
rz(-0.75251904) q[2];
rz(-0.13713947) q[3];
sx q[3];
rz(-0.64118853) q[3];
sx q[3];
rz(-0.34998676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8990477) q[0];
sx q[0];
rz(-2.516863) q[0];
sx q[0];
rz(2.94009) q[0];
rz(1.7789727) q[1];
sx q[1];
rz(-0.87281223) q[1];
sx q[1];
rz(2.9760117) q[1];
rz(0.049714391) q[2];
sx q[2];
rz(-1.6594973) q[2];
sx q[2];
rz(0.49458557) q[2];
rz(-1.1646545) q[3];
sx q[3];
rz(-1.2358758) q[3];
sx q[3];
rz(-0.14761543) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
