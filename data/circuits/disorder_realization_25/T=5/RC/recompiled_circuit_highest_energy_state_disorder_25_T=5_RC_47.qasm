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
rz(-0.93936062) q[0];
sx q[0];
rz(-2.0225749) q[0];
sx q[0];
rz(0.033666704) q[0];
rz(1.6839924) q[1];
sx q[1];
rz(4.9479244) q[1];
sx q[1];
rz(7.91628) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38070883) q[0];
sx q[0];
rz(-0.11822001) q[0];
sx q[0];
rz(1.3143172) q[0];
rz(-1.326697) q[2];
sx q[2];
rz(-2.8389922) q[2];
sx q[2];
rz(-1.5426334) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0284217) q[1];
sx q[1];
rz(-1.426371) q[1];
sx q[1];
rz(2.444205) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63920984) q[3];
sx q[3];
rz(-0.91851252) q[3];
sx q[3];
rz(-0.88310421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17026751) q[2];
sx q[2];
rz(-2.2496932) q[2];
sx q[2];
rz(-2.3049221) q[2];
rz(0.80638805) q[3];
sx q[3];
rz(-2.216335) q[3];
sx q[3];
rz(2.9562318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(3.0760913) q[0];
sx q[0];
rz(-2.0747023) q[0];
sx q[0];
rz(1.9245032) q[0];
rz(-2.941046) q[1];
sx q[1];
rz(-0.79981083) q[1];
sx q[1];
rz(-0.54380551) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6652619) q[0];
sx q[0];
rz(-1.7996739) q[0];
sx q[0];
rz(0.32461597) q[0];
rz(-pi) q[1];
rz(-0.5949623) q[2];
sx q[2];
rz(-1.5820855) q[2];
sx q[2];
rz(1.6807229) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4098379) q[1];
sx q[1];
rz(-1.8588762) q[1];
sx q[1];
rz(-1.9695915) q[1];
x q[2];
rz(-2.2460098) q[3];
sx q[3];
rz(-1.5434268) q[3];
sx q[3];
rz(2.7252008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.66889888) q[2];
sx q[2];
rz(-2.3987179) q[2];
sx q[2];
rz(-1.5354068) q[2];
rz(1.499048) q[3];
sx q[3];
rz(-2.5583772) q[3];
sx q[3];
rz(-1.8892052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1010308) q[0];
sx q[0];
rz(-2.5767548) q[0];
sx q[0];
rz(-3.132013) q[0];
rz(1.0293695) q[1];
sx q[1];
rz(-1.5681489) q[1];
sx q[1];
rz(-1.8052489) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1642417) q[0];
sx q[0];
rz(-1.0625496) q[0];
sx q[0];
rz(-0.68911106) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32149489) q[2];
sx q[2];
rz(-1.090637) q[2];
sx q[2];
rz(-0.26115044) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.0071466999) q[1];
sx q[1];
rz(-1.2255417) q[1];
sx q[1];
rz(0.18642872) q[1];
rz(-1.1747178) q[3];
sx q[3];
rz(-1.4142591) q[3];
sx q[3];
rz(-0.96352947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0336527) q[2];
sx q[2];
rz(-0.5113217) q[2];
sx q[2];
rz(1.8281724) q[2];
rz(-3.1016453) q[3];
sx q[3];
rz(-0.91244709) q[3];
sx q[3];
rz(2.0136755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0082598) q[0];
sx q[0];
rz(-1.2902211) q[0];
sx q[0];
rz(-2.7449352) q[0];
rz(2.2465514) q[1];
sx q[1];
rz(-1.4245234) q[1];
sx q[1];
rz(-2.3779714) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1930924) q[0];
sx q[0];
rz(-1.4903127) q[0];
sx q[0];
rz(1.3250687) q[0];
rz(1.5537797) q[2];
sx q[2];
rz(-1.6913171) q[2];
sx q[2];
rz(-1.8486384) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9634526) q[1];
sx q[1];
rz(-1.3333433) q[1];
sx q[1];
rz(-3.1141539) q[1];
x q[2];
rz(-2.850789) q[3];
sx q[3];
rz(-1.606308) q[3];
sx q[3];
rz(2.5574656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7274196) q[2];
sx q[2];
rz(-1.5199993) q[2];
sx q[2];
rz(-1.1302036) q[2];
rz(2.0027509) q[3];
sx q[3];
rz(-1.9500705) q[3];
sx q[3];
rz(-1.9738919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8324757) q[0];
sx q[0];
rz(-2.583355) q[0];
sx q[0];
rz(0.48156893) q[0];
rz(-2.8140977) q[1];
sx q[1];
rz(-1.4188674) q[1];
sx q[1];
rz(2.17735) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5835579) q[0];
sx q[0];
rz(-2.6917771) q[0];
sx q[0];
rz(-2.0925453) q[0];
x q[1];
rz(-2.1523508) q[2];
sx q[2];
rz(-1.4099401) q[2];
sx q[2];
rz(0.040078321) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3754978) q[1];
sx q[1];
rz(-1.8154486) q[1];
sx q[1];
rz(-3.0159113) q[1];
rz(-1.7751415) q[3];
sx q[3];
rz(-1.1875523) q[3];
sx q[3];
rz(-2.1658773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.97011906) q[2];
sx q[2];
rz(-1.118266) q[2];
sx q[2];
rz(-1.4946651) q[2];
rz(-2.1286879) q[3];
sx q[3];
rz(-0.92618147) q[3];
sx q[3];
rz(-2.6265898) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2038912) q[0];
sx q[0];
rz(-1.0045811) q[0];
sx q[0];
rz(1.0585693) q[0];
rz(0.92264289) q[1];
sx q[1];
rz(-1.0170931) q[1];
sx q[1];
rz(-0.10291084) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.424868) q[0];
sx q[0];
rz(-2.1916336) q[0];
sx q[0];
rz(-2.5607361) q[0];
x q[1];
rz(-1.2936959) q[2];
sx q[2];
rz(-0.99969802) q[2];
sx q[2];
rz(1.7383611) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.52572374) q[1];
sx q[1];
rz(-1.3868887) q[1];
sx q[1];
rz(-1.0332107) q[1];
rz(-pi) q[2];
rz(-0.0062557022) q[3];
sx q[3];
rz(-2.2872074) q[3];
sx q[3];
rz(-0.34230086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.70204488) q[2];
sx q[2];
rz(-0.44946686) q[2];
sx q[2];
rz(-0.31583819) q[2];
rz(2.1465007) q[3];
sx q[3];
rz(-1.2732048) q[3];
sx q[3];
rz(2.4163213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(2.134326) q[0];
sx q[0];
rz(-0.26679978) q[0];
sx q[0];
rz(0.63359487) q[0];
rz(0.38465706) q[1];
sx q[1];
rz(-1.1793914) q[1];
sx q[1];
rz(0.42919174) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68394008) q[0];
sx q[0];
rz(-2.5999617) q[0];
sx q[0];
rz(2.0787879) q[0];
rz(3.09893) q[2];
sx q[2];
rz(-0.68603078) q[2];
sx q[2];
rz(1.7614642) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8365578) q[1];
sx q[1];
rz(-1.9227131) q[1];
sx q[1];
rz(2.5266527) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6771992) q[3];
sx q[3];
rz(-2.0804318) q[3];
sx q[3];
rz(1.2023752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0523494) q[2];
sx q[2];
rz(-2.3640335) q[2];
sx q[2];
rz(0.26697978) q[2];
rz(-0.071768196) q[3];
sx q[3];
rz(-1.0177344) q[3];
sx q[3];
rz(1.4192386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4778022) q[0];
sx q[0];
rz(-0.1013805) q[0];
sx q[0];
rz(-1.7410393) q[0];
rz(1.1446704) q[1];
sx q[1];
rz(-2.0796227) q[1];
sx q[1];
rz(0.58715075) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42204866) q[0];
sx q[0];
rz(-1.9912155) q[0];
sx q[0];
rz(0.69164135) q[0];
rz(-pi) q[1];
rz(2.8835587) q[2];
sx q[2];
rz(-0.43083015) q[2];
sx q[2];
rz(1.8547386) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.78899164) q[1];
sx q[1];
rz(-1.7956211) q[1];
sx q[1];
rz(-2.4487752) q[1];
rz(0.31717698) q[3];
sx q[3];
rz(-0.97517386) q[3];
sx q[3];
rz(0.45444876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8553541) q[2];
sx q[2];
rz(-1.1796395) q[2];
sx q[2];
rz(-2.722495) q[2];
rz(1.6176443) q[3];
sx q[3];
rz(-0.91846275) q[3];
sx q[3];
rz(-1.5036478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.5986901) q[0];
sx q[0];
rz(-1.221523) q[0];
sx q[0];
rz(0.30938095) q[0];
rz(1.7912553) q[1];
sx q[1];
rz(-0.58809909) q[1];
sx q[1];
rz(-2.5877171) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0623035) q[0];
sx q[0];
rz(-1.8895747) q[0];
sx q[0];
rz(2.7026524) q[0];
rz(-pi) q[1];
rz(-1.9543086) q[2];
sx q[2];
rz(-1.4330136) q[2];
sx q[2];
rz(0.99064186) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.22928897) q[1];
sx q[1];
rz(-0.69432753) q[1];
sx q[1];
rz(-2.0635384) q[1];
rz(-pi) q[2];
rz(1.8819412) q[3];
sx q[3];
rz(-2.2486945) q[3];
sx q[3];
rz(2.4571612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89173335) q[2];
sx q[2];
rz(-1.6102108) q[2];
sx q[2];
rz(-1.3472793) q[2];
rz(2.3740718) q[3];
sx q[3];
rz(-1.1859505) q[3];
sx q[3];
rz(-0.84848136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7193741) q[0];
sx q[0];
rz(-0.30000559) q[0];
sx q[0];
rz(-1.2147709) q[0];
rz(-2.1396554) q[1];
sx q[1];
rz(-1.6286214) q[1];
sx q[1];
rz(2.217206) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0262605) q[0];
sx q[0];
rz(-0.9965653) q[0];
sx q[0];
rz(2.6747245) q[0];
x q[1];
rz(-0.027410313) q[2];
sx q[2];
rz(-2.1152585) q[2];
sx q[2];
rz(-1.4278864) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1781285) q[1];
sx q[1];
rz(-0.22451065) q[1];
sx q[1];
rz(0.1046532) q[1];
rz(2.2460292) q[3];
sx q[3];
rz(-2.319284) q[3];
sx q[3];
rz(-3.0117717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7944916) q[2];
sx q[2];
rz(-2.2091986) q[2];
sx q[2];
rz(-2.7744228) q[2];
rz(-1.9764887) q[3];
sx q[3];
rz(-0.96787435) q[3];
sx q[3];
rz(-0.86618209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1871056) q[0];
sx q[0];
rz(-2.4024873) q[0];
sx q[0];
rz(-3.0229229) q[0];
rz(2.8754996) q[1];
sx q[1];
rz(-2.2292744) q[1];
sx q[1];
rz(1.6899065) q[1];
rz(2.8546147) q[2];
sx q[2];
rz(-1.3084044) q[2];
sx q[2];
rz(-0.85592782) q[2];
rz(0.20826505) q[3];
sx q[3];
rz(-0.97934813) q[3];
sx q[3];
rz(-1.4772268) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
