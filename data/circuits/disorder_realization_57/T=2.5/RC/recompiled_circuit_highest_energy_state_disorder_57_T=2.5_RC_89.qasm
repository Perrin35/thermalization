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
rz(-1.7419185) q[0];
sx q[0];
rz(-1.4631441) q[0];
sx q[0];
rz(1.4611257) q[0];
rz(2.6598016) q[1];
sx q[1];
rz(-2.3269589) q[1];
sx q[1];
rz(2.2021267) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0622913) q[0];
sx q[0];
rz(-1.5271375) q[0];
sx q[0];
rz(-0.0096304499) q[0];
x q[1];
rz(-2.6109735) q[2];
sx q[2];
rz(-1.8920027) q[2];
sx q[2];
rz(1.2752334) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.11096) q[1];
sx q[1];
rz(-0.041555066) q[1];
sx q[1];
rz(-1.3135733) q[1];
rz(-2.7339805) q[3];
sx q[3];
rz(-1.3702754) q[3];
sx q[3];
rz(1.149857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7751094) q[2];
sx q[2];
rz(-1.6747687) q[2];
sx q[2];
rz(0.62466204) q[2];
rz(0.96056739) q[3];
sx q[3];
rz(-1.4550236) q[3];
sx q[3];
rz(-0.87489405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7333154) q[0];
sx q[0];
rz(-0.69077078) q[0];
sx q[0];
rz(-2.9498003) q[0];
rz(-2.5674112) q[1];
sx q[1];
rz(-1.5230813) q[1];
sx q[1];
rz(0.40723732) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8005368) q[0];
sx q[0];
rz(-1.558702) q[0];
sx q[0];
rz(1.3950134) q[0];
x q[1];
rz(-3.0767976) q[2];
sx q[2];
rz(-1.7470834) q[2];
sx q[2];
rz(0.48688146) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5246689) q[1];
sx q[1];
rz(-1.7546931) q[1];
sx q[1];
rz(2.6016629) q[1];
rz(-0.089669946) q[3];
sx q[3];
rz(-1.0700135) q[3];
sx q[3];
rz(-1.9865685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.189397) q[2];
sx q[2];
rz(-1.1789097) q[2];
sx q[2];
rz(0.65323812) q[2];
rz(2.2630528) q[3];
sx q[3];
rz(-2.0432751) q[3];
sx q[3];
rz(-3.0285335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99881309) q[0];
sx q[0];
rz(-1.4343364) q[0];
sx q[0];
rz(-1.0502195) q[0];
rz(-2.5255919) q[1];
sx q[1];
rz(-1.417825) q[1];
sx q[1];
rz(0.91316694) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28024689) q[0];
sx q[0];
rz(-1.4644091) q[0];
sx q[0];
rz(1.4004571) q[0];
rz(-pi) q[1];
rz(-1.2840604) q[2];
sx q[2];
rz(-1.456919) q[2];
sx q[2];
rz(1.8658474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2985067) q[1];
sx q[1];
rz(-1.8131327) q[1];
sx q[1];
rz(1.0964811) q[1];
rz(-0.88190474) q[3];
sx q[3];
rz(-1.2744474) q[3];
sx q[3];
rz(-2.0810082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.22127613) q[2];
sx q[2];
rz(-2.2344799) q[2];
sx q[2];
rz(-3.082412) q[2];
rz(2.3024998) q[3];
sx q[3];
rz(-1.9234761) q[3];
sx q[3];
rz(0.86735094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7729823) q[0];
sx q[0];
rz(-2.1324069) q[0];
sx q[0];
rz(2.7918145) q[0];
rz(-1.7971669) q[1];
sx q[1];
rz(-0.58660048) q[1];
sx q[1];
rz(3.0288568) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37575133) q[0];
sx q[0];
rz(-2.021702) q[0];
sx q[0];
rz(-1.9519873) q[0];
rz(1.2069682) q[2];
sx q[2];
rz(-2.4487961) q[2];
sx q[2];
rz(2.7464904) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4055914) q[1];
sx q[1];
rz(-1.2697769) q[1];
sx q[1];
rz(2.244414) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2318683) q[3];
sx q[3];
rz(-2.1436757) q[3];
sx q[3];
rz(2.6343144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3939487) q[2];
sx q[2];
rz(-1.5564206) q[2];
sx q[2];
rz(3.0843688) q[2];
rz(1.8852437) q[3];
sx q[3];
rz(-0.57259721) q[3];
sx q[3];
rz(0.73067874) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2644065) q[0];
sx q[0];
rz(-1.974396) q[0];
sx q[0];
rz(0.71006829) q[0];
rz(-0.26142985) q[1];
sx q[1];
rz(-1.556004) q[1];
sx q[1];
rz(0.95103055) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20999058) q[0];
sx q[0];
rz(-1.8977802) q[0];
sx q[0];
rz(0.47173421) q[0];
x q[1];
rz(-2.5158993) q[2];
sx q[2];
rz(-0.73391418) q[2];
sx q[2];
rz(-0.87482051) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3465865) q[1];
sx q[1];
rz(-1.3685431) q[1];
sx q[1];
rz(-2.4319885) q[1];
x q[2];
rz(-0.099154648) q[3];
sx q[3];
rz(-2.2014849) q[3];
sx q[3];
rz(2.658228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49932536) q[2];
sx q[2];
rz(-0.32411164) q[2];
sx q[2];
rz(-0.30936852) q[2];
rz(-2.046509) q[3];
sx q[3];
rz(-1.1284004) q[3];
sx q[3];
rz(2.7373718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22769044) q[0];
sx q[0];
rz(-0.18276437) q[0];
sx q[0];
rz(2.2392654) q[0];
rz(3.0790216) q[1];
sx q[1];
rz(-1.2689271) q[1];
sx q[1];
rz(-1.2455469) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7880873) q[0];
sx q[0];
rz(-0.89998945) q[0];
sx q[0];
rz(-1.7680579) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4772513) q[2];
sx q[2];
rz(-0.54160833) q[2];
sx q[2];
rz(-2.6111185) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3712103) q[1];
sx q[1];
rz(-2.4744445) q[1];
sx q[1];
rz(2.4333111) q[1];
rz(-pi) q[2];
rz(2.5596722) q[3];
sx q[3];
rz(-0.76184638) q[3];
sx q[3];
rz(-0.91530734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.41783276) q[2];
sx q[2];
rz(-1.87169) q[2];
sx q[2];
rz(-1.7977686) q[2];
rz(1.8580681) q[3];
sx q[3];
rz(-2.696974) q[3];
sx q[3];
rz(0.72492391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1062687) q[0];
sx q[0];
rz(-3.1043053) q[0];
sx q[0];
rz(-0.37621793) q[0];
rz(-0.47209921) q[1];
sx q[1];
rz(-2.4596228) q[1];
sx q[1];
rz(0.45904407) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5482785) q[0];
sx q[0];
rz(-0.86515831) q[0];
sx q[0];
rz(-1.0039019) q[0];
x q[1];
rz(-2.6920175) q[2];
sx q[2];
rz(-1.7284942) q[2];
sx q[2];
rz(2.0248991) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9320786) q[1];
sx q[1];
rz(-2.9649295) q[1];
sx q[1];
rz(2.7448369) q[1];
rz(-pi) q[2];
rz(0.59765069) q[3];
sx q[3];
rz(-1.5603746) q[3];
sx q[3];
rz(2.8437993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0692733) q[2];
sx q[2];
rz(-2.4299712) q[2];
sx q[2];
rz(-2.7834328) q[2];
rz(-2.2530796) q[3];
sx q[3];
rz(-1.0212967) q[3];
sx q[3];
rz(-0.93069881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096980378) q[0];
sx q[0];
rz(-0.29878578) q[0];
sx q[0];
rz(0.97691798) q[0];
rz(0.75374976) q[1];
sx q[1];
rz(-1.4434554) q[1];
sx q[1];
rz(-1.5026198) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3328505) q[0];
sx q[0];
rz(-0.70510222) q[0];
sx q[0];
rz(-0.1211959) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6303275) q[2];
sx q[2];
rz(-2.3663372) q[2];
sx q[2];
rz(-2.023319) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1267008) q[1];
sx q[1];
rz(-2.3897244) q[1];
sx q[1];
rz(-2.6905462) q[1];
rz(-pi) q[2];
rz(0.23076337) q[3];
sx q[3];
rz(-2.181859) q[3];
sx q[3];
rz(-2.3045674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.94540191) q[2];
sx q[2];
rz(-1.6795029) q[2];
sx q[2];
rz(-0.84575829) q[2];
rz(-0.78019199) q[3];
sx q[3];
rz(-1.3581685) q[3];
sx q[3];
rz(-2.5030524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82275259) q[0];
sx q[0];
rz(-2.0697116) q[0];
sx q[0];
rz(0.34291294) q[0];
rz(-2.137939) q[1];
sx q[1];
rz(-1.2316848) q[1];
sx q[1];
rz(2.4249605) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0239923) q[0];
sx q[0];
rz(-0.85451925) q[0];
sx q[0];
rz(2.0407659) q[0];
rz(-pi) q[1];
rz(1.253867) q[2];
sx q[2];
rz(-1.0980596) q[2];
sx q[2];
rz(2.6726618) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9616507) q[1];
sx q[1];
rz(-1.3030717) q[1];
sx q[1];
rz(-0.60795345) q[1];
rz(-pi) q[2];
rz(1.1716722) q[3];
sx q[3];
rz(-2.2033327) q[3];
sx q[3];
rz(-0.61033953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0922682) q[2];
sx q[2];
rz(-1.5169787) q[2];
sx q[2];
rz(-2.1410904) q[2];
rz(0.81691027) q[3];
sx q[3];
rz(-1.2714081) q[3];
sx q[3];
rz(2.7774155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3889121) q[0];
sx q[0];
rz(-2.2081544) q[0];
sx q[0];
rz(-2.3626784) q[0];
rz(2.3498416) q[1];
sx q[1];
rz(-1.2261032) q[1];
sx q[1];
rz(2.7955999) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2591922) q[0];
sx q[0];
rz(-1.1244052) q[0];
sx q[0];
rz(-1.4784228) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1295165) q[2];
sx q[2];
rz(-1.8740956) q[2];
sx q[2];
rz(2.8331918) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.71042594) q[1];
sx q[1];
rz(-1.0252608) q[1];
sx q[1];
rz(-0.04108508) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84684144) q[3];
sx q[3];
rz(-1.2961642) q[3];
sx q[3];
rz(-1.7881029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9927696) q[2];
sx q[2];
rz(-1.537348) q[2];
sx q[2];
rz(1.0619304) q[2];
rz(0.81260931) q[3];
sx q[3];
rz(-1.9989697) q[3];
sx q[3];
rz(0.44212166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84520311) q[0];
sx q[0];
rz(-2.3136105) q[0];
sx q[0];
rz(-1.4985341) q[0];
rz(-2.8614112) q[1];
sx q[1];
rz(-1.9277086) q[1];
sx q[1];
rz(2.3452506) q[1];
rz(-2.7325999) q[2];
sx q[2];
rz(-1.5341285) q[2];
sx q[2];
rz(-1.6899012) q[2];
rz(2.3128458) q[3];
sx q[3];
rz(-1.256905) q[3];
sx q[3];
rz(-1.4245413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
