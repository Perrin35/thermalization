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
rz(-0.91322652) q[0];
sx q[0];
rz(-0.93728596) q[0];
sx q[0];
rz(-0.32567853) q[0];
rz(0.20618662) q[1];
sx q[1];
rz(-0.34178692) q[1];
sx q[1];
rz(0.49965247) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2316811) q[0];
sx q[0];
rz(-1.9107816) q[0];
sx q[0];
rz(2.489734) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0023053447) q[2];
sx q[2];
rz(-2.1426149) q[2];
sx q[2];
rz(-0.44849685) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.25741577) q[1];
sx q[1];
rz(-1.9828772) q[1];
sx q[1];
rz(0.33452371) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1127171) q[3];
sx q[3];
rz(-1.5124953) q[3];
sx q[3];
rz(-1.3227303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.25834045) q[2];
sx q[2];
rz(-0.84104937) q[2];
sx q[2];
rz(0.19634761) q[2];
rz(-1.9718862) q[3];
sx q[3];
rz(-1.615808) q[3];
sx q[3];
rz(-1.660478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.779988) q[0];
sx q[0];
rz(-1.1971645) q[0];
sx q[0];
rz(-2.8874604) q[0];
rz(0.77468553) q[1];
sx q[1];
rz(-1.4137555) q[1];
sx q[1];
rz(-0.69001251) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7609576) q[0];
sx q[0];
rz(-1.1948144) q[0];
sx q[0];
rz(1.4248821) q[0];
rz(1.5845845) q[2];
sx q[2];
rz(-2.1796527) q[2];
sx q[2];
rz(-1.4964455) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.96448) q[1];
sx q[1];
rz(-1.2487673) q[1];
sx q[1];
rz(0.42735512) q[1];
rz(-pi) q[2];
rz(-1.0496185) q[3];
sx q[3];
rz(-1.3173977) q[3];
sx q[3];
rz(0.99072841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0320354) q[2];
sx q[2];
rz(-0.037470128) q[2];
sx q[2];
rz(2.1407342) q[2];
rz(-1.9199269) q[3];
sx q[3];
rz(-1.6358717) q[3];
sx q[3];
rz(-3.1083623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.074653305) q[0];
sx q[0];
rz(-0.017843857) q[0];
sx q[0];
rz(-0.53027207) q[0];
rz(2.6182134) q[1];
sx q[1];
rz(-1.7019848) q[1];
sx q[1];
rz(-1.9022) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8187934) q[0];
sx q[0];
rz(-1.8086495) q[0];
sx q[0];
rz(-2.4616694) q[0];
x q[1];
rz(1.3472413) q[2];
sx q[2];
rz(-1.8164779) q[2];
sx q[2];
rz(-1.9839191) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7433333) q[1];
sx q[1];
rz(-1.6930237) q[1];
sx q[1];
rz(1.2714703) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1887458) q[3];
sx q[3];
rz(-1.1955559) q[3];
sx q[3];
rz(0.98179152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8221028) q[2];
sx q[2];
rz(-0.73188657) q[2];
sx q[2];
rz(-0.11719318) q[2];
rz(2.9486837) q[3];
sx q[3];
rz(-1.1567876) q[3];
sx q[3];
rz(-2.9624654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8837226) q[0];
sx q[0];
rz(-1.1494881) q[0];
sx q[0];
rz(-3.1174739) q[0];
rz(0.78512496) q[1];
sx q[1];
rz(-1.5631915) q[1];
sx q[1];
rz(1.1611353) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28306357) q[0];
sx q[0];
rz(-1.8743321) q[0];
sx q[0];
rz(0.86743939) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5490521) q[2];
sx q[2];
rz(-1.7090343) q[2];
sx q[2];
rz(-2.5858226) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9332446) q[1];
sx q[1];
rz(-1.0050432) q[1];
sx q[1];
rz(1.1858334) q[1];
x q[2];
rz(0.9735017) q[3];
sx q[3];
rz(-1.4184059) q[3];
sx q[3];
rz(0.77578557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44196096) q[2];
sx q[2];
rz(-0.72579759) q[2];
sx q[2];
rz(-0.3886784) q[2];
rz(-1.9937438) q[3];
sx q[3];
rz(-2.0129222) q[3];
sx q[3];
rz(-1.8584937) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.024014274) q[0];
sx q[0];
rz(-0.83947459) q[0];
sx q[0];
rz(2.1003387) q[0];
rz(-2.4196692) q[1];
sx q[1];
rz(-1.8903774) q[1];
sx q[1];
rz(1.3232683) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6206317) q[0];
sx q[0];
rz(-2.8537321) q[0];
sx q[0];
rz(2.5858712) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2385246) q[2];
sx q[2];
rz(-0.67004824) q[2];
sx q[2];
rz(2.5243896) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8787707) q[1];
sx q[1];
rz(-0.77057099) q[1];
sx q[1];
rz(-2.7380323) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3106662) q[3];
sx q[3];
rz(-1.8046265) q[3];
sx q[3];
rz(1.1375858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2920275) q[2];
sx q[2];
rz(-0.592507) q[2];
sx q[2];
rz(2.3805857) q[2];
rz(-1.546953) q[3];
sx q[3];
rz(-1.2706815) q[3];
sx q[3];
rz(2.7775619) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7008692) q[0];
sx q[0];
rz(-2.4969164) q[0];
sx q[0];
rz(-2.0194637) q[0];
rz(-0.912965) q[1];
sx q[1];
rz(-2.2367621) q[1];
sx q[1];
rz(-0.95776552) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4816613) q[0];
sx q[0];
rz(-2.5536116) q[0];
sx q[0];
rz(1.6350782) q[0];
x q[1];
rz(1.945044) q[2];
sx q[2];
rz(-0.23434445) q[2];
sx q[2];
rz(-1.0398911) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1407048) q[1];
sx q[1];
rz(-1.1240214) q[1];
sx q[1];
rz(1.3088777) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9676858) q[3];
sx q[3];
rz(-1.4507796) q[3];
sx q[3];
rz(0.83156026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9271348) q[2];
sx q[2];
rz(-0.99481499) q[2];
sx q[2];
rz(-2.4289995) q[2];
rz(1.384895) q[3];
sx q[3];
rz(-1.2777998) q[3];
sx q[3];
rz(-0.28759292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49823847) q[0];
sx q[0];
rz(-2.9297628) q[0];
sx q[0];
rz(2.8969452) q[0];
rz(-0.08617607) q[1];
sx q[1];
rz(-1.1174322) q[1];
sx q[1];
rz(2.2138331) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13138765) q[0];
sx q[0];
rz(-0.82021111) q[0];
sx q[0];
rz(1.6774235) q[0];
x q[1];
rz(2.5439569) q[2];
sx q[2];
rz(-0.25271395) q[2];
sx q[2];
rz(-1.4022577) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.40585247) q[1];
sx q[1];
rz(-2.3946792) q[1];
sx q[1];
rz(-2.4290415) q[1];
x q[2];
rz(-0.71998511) q[3];
sx q[3];
rz(-0.84331028) q[3];
sx q[3];
rz(2.2184531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1722395) q[2];
sx q[2];
rz(-0.66102207) q[2];
sx q[2];
rz(-1.6924525) q[2];
rz(-1.0714072) q[3];
sx q[3];
rz(-1.7457733) q[3];
sx q[3];
rz(-1.2034108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093056738) q[0];
sx q[0];
rz(-3.0003248) q[0];
sx q[0];
rz(2.5988044) q[0];
rz(-3.0154748) q[1];
sx q[1];
rz(-1.7593971) q[1];
sx q[1];
rz(0.19475591) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4493748) q[0];
sx q[0];
rz(-2.0701773) q[0];
sx q[0];
rz(-1.7591018) q[0];
rz(1.0967238) q[2];
sx q[2];
rz(-1.1685089) q[2];
sx q[2];
rz(1.2685217) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.69079) q[1];
sx q[1];
rz(-2.6581785) q[1];
sx q[1];
rz(-0.40287896) q[1];
rz(-pi) q[2];
rz(-0.1925999) q[3];
sx q[3];
rz(-2.3705774) q[3];
sx q[3];
rz(-0.23480496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18437082) q[2];
sx q[2];
rz(-1.3536644) q[2];
sx q[2];
rz(-1.703519) q[2];
rz(-2.461869) q[3];
sx q[3];
rz(-0.34346911) q[3];
sx q[3];
rz(0.22469416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1723802) q[0];
sx q[0];
rz(-0.42177105) q[0];
sx q[0];
rz(-2.1286185) q[0];
rz(-0.29533932) q[1];
sx q[1];
rz(-1.4238009) q[1];
sx q[1];
rz(-1.0257592) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8689629) q[0];
sx q[0];
rz(-0.40237793) q[0];
sx q[0];
rz(2.8182882) q[0];
rz(-pi) q[1];
rz(2.1842833) q[2];
sx q[2];
rz(-2.0762091) q[2];
sx q[2];
rz(-1.2135722) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6597627) q[1];
sx q[1];
rz(-0.24630486) q[1];
sx q[1];
rz(-0.15683163) q[1];
rz(0.72297289) q[3];
sx q[3];
rz(-1.3275258) q[3];
sx q[3];
rz(1.9169501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.50461304) q[2];
sx q[2];
rz(-0.29256088) q[2];
sx q[2];
rz(2.2620849) q[2];
rz(-2.8579936) q[3];
sx q[3];
rz(-1.2181543) q[3];
sx q[3];
rz(-1.2951736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2844236) q[0];
sx q[0];
rz(-1.5791945) q[0];
sx q[0];
rz(-2.5403585) q[0];
rz(-0.55016905) q[1];
sx q[1];
rz(-1.0160805) q[1];
sx q[1];
rz(0.30204958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45363891) q[0];
sx q[0];
rz(-2.1778202) q[0];
sx q[0];
rz(-2.2974854) q[0];
rz(-1.0112125) q[2];
sx q[2];
rz(-2.4200984) q[2];
sx q[2];
rz(-1.341429) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2419521) q[1];
sx q[1];
rz(-1.5664827) q[1];
sx q[1];
rz(-0.14546995) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9317231) q[3];
sx q[3];
rz(-0.59912938) q[3];
sx q[3];
rz(-0.15219469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.067513116) q[2];
sx q[2];
rz(-0.16177598) q[2];
sx q[2];
rz(-1.3241241) q[2];
rz(-2.6586804) q[3];
sx q[3];
rz(-0.78871471) q[3];
sx q[3];
rz(-2.8632792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0637689) q[0];
sx q[0];
rz(-1.9763197) q[0];
sx q[0];
rz(1.9015953) q[0];
rz(2.3700312) q[1];
sx q[1];
rz(-2.0004708) q[1];
sx q[1];
rz(-3.0761459) q[1];
rz(-2.9289519) q[2];
sx q[2];
rz(-2.2937951) q[2];
sx q[2];
rz(-0.19611025) q[2];
rz(2.2167233) q[3];
sx q[3];
rz(-0.51565167) q[3];
sx q[3];
rz(-2.6272163) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
