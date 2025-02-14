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
rz(-1.9077644) q[1];
sx q[1];
rz(5.0662) q[1];
sx q[1];
rz(10.868664) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47752781) q[0];
sx q[0];
rz(-1.489138) q[0];
sx q[0];
rz(0.89952138) q[0];
x q[1];
rz(2.2488689) q[2];
sx q[2];
rz(-0.44995445) q[2];
sx q[2];
rz(-1.8261248) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7942209) q[1];
sx q[1];
rz(-1.9772775) q[1];
sx q[1];
rz(-2.4665753) q[1];
rz(-1.9099376) q[3];
sx q[3];
rz(-1.5970486) q[3];
sx q[3];
rz(1.8363801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7794789) q[2];
sx q[2];
rz(-1.9356091) q[2];
sx q[2];
rz(0.77902478) q[2];
rz(1.3075167) q[3];
sx q[3];
rz(-2.8179759) q[3];
sx q[3];
rz(-1.754508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5883412) q[0];
sx q[0];
rz(-1.7760176) q[0];
sx q[0];
rz(-0.27728444) q[0];
rz(0.39018997) q[1];
sx q[1];
rz(-2.0398102) q[1];
sx q[1];
rz(1.1276833) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9830071) q[0];
sx q[0];
rz(-1.5282573) q[0];
sx q[0];
rz(1.5908414) q[0];
rz(-1.4696944) q[2];
sx q[2];
rz(-1.2866469) q[2];
sx q[2];
rz(0.84188879) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1192644) q[1];
sx q[1];
rz(-2.8624967) q[1];
sx q[1];
rz(-1.7995681) q[1];
x q[2];
rz(0.66172285) q[3];
sx q[3];
rz(-1.265101) q[3];
sx q[3];
rz(0.88246417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7057425) q[2];
sx q[2];
rz(-0.59429344) q[2];
sx q[2];
rz(-1.9464114) q[2];
rz(-0.65520203) q[3];
sx q[3];
rz(-1.6569258) q[3];
sx q[3];
rz(0.5932194) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4707659) q[0];
sx q[0];
rz(-0.68597454) q[0];
sx q[0];
rz(0.92054787) q[0];
rz(-1.8596733) q[1];
sx q[1];
rz(-1.1290519) q[1];
sx q[1];
rz(-1.4528073) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4494999) q[0];
sx q[0];
rz(-1.458174) q[0];
sx q[0];
rz(-1.6586668) q[0];
rz(-2.0412943) q[2];
sx q[2];
rz(-0.62242939) q[2];
sx q[2];
rz(0.4284455) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7662168) q[1];
sx q[1];
rz(-1.2355243) q[1];
sx q[1];
rz(-1.5165308) q[1];
x q[2];
rz(3.0759144) q[3];
sx q[3];
rz(-1.6367776) q[3];
sx q[3];
rz(-3.0998041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.61536962) q[2];
sx q[2];
rz(-2.0873895) q[2];
sx q[2];
rz(0.39895454) q[2];
rz(-0.51359549) q[3];
sx q[3];
rz(-2.8245638) q[3];
sx q[3];
rz(-1.992146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8154163) q[0];
sx q[0];
rz(-0.57891095) q[0];
sx q[0];
rz(-1.2166566) q[0];
rz(-0.51219621) q[1];
sx q[1];
rz(-2.0307816) q[1];
sx q[1];
rz(-2.0373352) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51795635) q[0];
sx q[0];
rz(-0.41455634) q[0];
sx q[0];
rz(-0.74409938) q[0];
rz(-pi) q[1];
rz(-0.17357628) q[2];
sx q[2];
rz(-1.7494698) q[2];
sx q[2];
rz(1.0965958) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3901375) q[1];
sx q[1];
rz(-1.631076) q[1];
sx q[1];
rz(-1.3740609) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85569546) q[3];
sx q[3];
rz(-0.33791204) q[3];
sx q[3];
rz(-2.5355946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5964552) q[2];
sx q[2];
rz(-1.1918273) q[2];
sx q[2];
rz(0.54984251) q[2];
rz(0.95692974) q[3];
sx q[3];
rz(-0.79920971) q[3];
sx q[3];
rz(-1.6526615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4855708) q[0];
sx q[0];
rz(-0.85000426) q[0];
sx q[0];
rz(-0.75918424) q[0];
rz(0.94866577) q[1];
sx q[1];
rz(-1.1439088) q[1];
sx q[1];
rz(-2.8099828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.054197) q[0];
sx q[0];
rz(-2.3080462) q[0];
sx q[0];
rz(1.1655318) q[0];
rz(-0.47470113) q[2];
sx q[2];
rz(-1.5384288) q[2];
sx q[2];
rz(2.2970478) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3508671) q[1];
sx q[1];
rz(-1.7247883) q[1];
sx q[1];
rz(-1.749498) q[1];
x q[2];
rz(2.7301627) q[3];
sx q[3];
rz(-1.6959018) q[3];
sx q[3];
rz(0.26857146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9552976) q[2];
sx q[2];
rz(-0.95185995) q[2];
sx q[2];
rz(-2.1985506) q[2];
rz(0.17213639) q[3];
sx q[3];
rz(-0.94477263) q[3];
sx q[3];
rz(2.9183563) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.405769) q[0];
sx q[0];
rz(-0.74546927) q[0];
sx q[0];
rz(1.5752342) q[0];
rz(0.15815059) q[1];
sx q[1];
rz(-2.3731396) q[1];
sx q[1];
rz(0.92811531) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6142352) q[0];
sx q[0];
rz(-1.6993465) q[0];
sx q[0];
rz(-0.30621333) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.720143) q[2];
sx q[2];
rz(-2.0837651) q[2];
sx q[2];
rz(2.9861272) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5852768) q[1];
sx q[1];
rz(-0.92570526) q[1];
sx q[1];
rz(-0.32151244) q[1];
rz(-1.4466049) q[3];
sx q[3];
rz(-1.9212133) q[3];
sx q[3];
rz(-2.7563376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.9485665) q[2];
sx q[2];
rz(-0.28634772) q[2];
sx q[2];
rz(-2.7890653) q[2];
rz(-2.7726717) q[3];
sx q[3];
rz(-0.58266321) q[3];
sx q[3];
rz(2.2841456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3567268) q[0];
sx q[0];
rz(-0.016473869) q[0];
sx q[0];
rz(2.4705868) q[0];
rz(0.10637936) q[1];
sx q[1];
rz(-2.2439067) q[1];
sx q[1];
rz(-2.5118714) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9167921) q[0];
sx q[0];
rz(-0.19242254) q[0];
sx q[0];
rz(-0.097956603) q[0];
rz(-pi) q[1];
rz(0.4102629) q[2];
sx q[2];
rz(-1.9681491) q[2];
sx q[2];
rz(-0.58471459) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5725745) q[1];
sx q[1];
rz(-0.33199939) q[1];
sx q[1];
rz(0.23004679) q[1];
rz(-pi) q[2];
rz(-1.2604146) q[3];
sx q[3];
rz(-0.93573007) q[3];
sx q[3];
rz(2.6657651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0715535) q[2];
sx q[2];
rz(-2.2420501) q[2];
sx q[2];
rz(2.1001935) q[2];
rz(1.533016) q[3];
sx q[3];
rz(-0.93419111) q[3];
sx q[3];
rz(1.9937203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.94569412) q[0];
sx q[0];
rz(-0.65106374) q[0];
sx q[0];
rz(2.661929) q[0];
rz(-0.02462968) q[1];
sx q[1];
rz(-1.004091) q[1];
sx q[1];
rz(-3.0020795) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64785731) q[0];
sx q[0];
rz(-1.656504) q[0];
sx q[0];
rz(-1.6802963) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37711511) q[2];
sx q[2];
rz(-2.0465436) q[2];
sx q[2];
rz(-1.7390133) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7558328) q[1];
sx q[1];
rz(-2.2221098) q[1];
sx q[1];
rz(-1.7803704) q[1];
x q[2];
rz(-3.0728389) q[3];
sx q[3];
rz(-0.67151755) q[3];
sx q[3];
rz(-2.9686787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8861683) q[2];
sx q[2];
rz(-1.5972127) q[2];
sx q[2];
rz(1.1523979) q[2];
rz(1.7216916) q[3];
sx q[3];
rz(-2.4072188) q[3];
sx q[3];
rz(-1.8324119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0814447) q[0];
sx q[0];
rz(-1.4652493) q[0];
sx q[0];
rz(-1.4981221) q[0];
rz(-0.76088798) q[1];
sx q[1];
rz(-2.1983169) q[1];
sx q[1];
rz(1.4563742) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4382317) q[0];
sx q[0];
rz(-1.7499763) q[0];
sx q[0];
rz(0.56877259) q[0];
rz(-pi) q[1];
rz(1.782519) q[2];
sx q[2];
rz(-1.178732) q[2];
sx q[2];
rz(1.388092) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2029426) q[1];
sx q[1];
rz(-2.4846145) q[1];
sx q[1];
rz(-0.88969661) q[1];
rz(-pi) q[2];
rz(2.1797769) q[3];
sx q[3];
rz(-0.66003335) q[3];
sx q[3];
rz(0.48840085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3219519) q[2];
sx q[2];
rz(-2.7072622) q[2];
sx q[2];
rz(2.1257909) q[2];
rz(-2.2470233) q[3];
sx q[3];
rz(-1.2457448) q[3];
sx q[3];
rz(-0.7416803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9858518) q[0];
sx q[0];
rz(-2.4116801) q[0];
sx q[0];
rz(-2.1515382) q[0];
rz(-2.2186225) q[1];
sx q[1];
rz(-1.3786517) q[1];
sx q[1];
rz(-0.96614456) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7434563) q[0];
sx q[0];
rz(-1.4706637) q[0];
sx q[0];
rz(3.0707573) q[0];
rz(2.0311293) q[2];
sx q[2];
rz(-1.9607414) q[2];
sx q[2];
rz(1.4947948) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.50033874) q[1];
sx q[1];
rz(-0.49420824) q[1];
sx q[1];
rz(-2.2142835) q[1];
x q[2];
rz(1.0564141) q[3];
sx q[3];
rz(-2.9136411) q[3];
sx q[3];
rz(-2.0755656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3155589) q[2];
sx q[2];
rz(-0.78353271) q[2];
sx q[2];
rz(-0.75251904) q[2];
rz(0.13713947) q[3];
sx q[3];
rz(-2.5004041) q[3];
sx q[3];
rz(2.7916059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8990477) q[0];
sx q[0];
rz(-0.62472961) q[0];
sx q[0];
rz(-0.20150264) q[0];
rz(1.3626199) q[1];
sx q[1];
rz(-2.2687804) q[1];
sx q[1];
rz(-0.16558095) q[1];
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
