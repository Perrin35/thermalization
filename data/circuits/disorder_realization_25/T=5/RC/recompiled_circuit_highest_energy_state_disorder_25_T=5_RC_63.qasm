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
rz(4.2606104) q[0];
sx q[0];
rz(9.4584447) q[0];
rz(-1.4576003) q[1];
sx q[1];
rz(-1.8063318) q[1];
sx q[1];
rz(1.508498) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7608838) q[0];
sx q[0];
rz(-3.0233726) q[0];
sx q[0];
rz(-1.8272754) q[0];
rz(-1.8148957) q[2];
sx q[2];
rz(-0.30260049) q[2];
sx q[2];
rz(-1.5426334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7693387) q[1];
sx q[1];
rz(-2.4318691) q[1];
sx q[1];
rz(-0.22270568) q[1];
x q[2];
rz(-0.90773031) q[3];
sx q[3];
rz(-2.2624348) q[3];
sx q[3];
rz(-1.7691106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.17026751) q[2];
sx q[2];
rz(-2.2496932) q[2];
sx q[2];
rz(0.83667052) q[2];
rz(0.80638805) q[3];
sx q[3];
rz(-2.216335) q[3];
sx q[3];
rz(2.9562318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.065501325) q[0];
sx q[0];
rz(-2.0747023) q[0];
sx q[0];
rz(-1.2170894) q[0];
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
rz(1.4763308) q[0];
sx q[0];
rz(-1.7996739) q[0];
sx q[0];
rz(0.32461597) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1214525) q[2];
sx q[2];
rz(-2.5465362) q[2];
sx q[2];
rz(3.014987) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.041808279) q[1];
sx q[1];
rz(-1.9522894) q[1];
sx q[1];
rz(0.31111335) q[1];
rz(-pi) q[2];
rz(1.614566) q[3];
sx q[3];
rz(-2.4659116) q[3];
sx q[3];
rz(-2.0213493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4726938) q[2];
sx q[2];
rz(-2.3987179) q[2];
sx q[2];
rz(-1.5354068) q[2];
rz(1.499048) q[3];
sx q[3];
rz(-0.58321548) q[3];
sx q[3];
rz(1.8892052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0405618) q[0];
sx q[0];
rz(-0.56483785) q[0];
sx q[0];
rz(3.132013) q[0];
rz(-1.0293695) q[1];
sx q[1];
rz(-1.5681489) q[1];
sx q[1];
rz(1.8052489) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1270008) q[0];
sx q[0];
rz(-2.3107502) q[0];
sx q[0];
rz(0.71944351) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.025454) q[2];
sx q[2];
rz(-0.57078123) q[2];
sx q[2];
rz(-2.2557543) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6402898) q[1];
sx q[1];
rz(-2.7510018) q[1];
sx q[1];
rz(2.0466481) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9668749) q[3];
sx q[3];
rz(-1.4142591) q[3];
sx q[3];
rz(2.1780632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.10793992) q[2];
sx q[2];
rz(-0.5113217) q[2];
sx q[2];
rz(-1.3134202) q[2];
rz(0.039947346) q[3];
sx q[3];
rz(-2.2291456) q[3];
sx q[3];
rz(-2.0136755) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1333328) q[0];
sx q[0];
rz(-1.8513716) q[0];
sx q[0];
rz(-2.7449352) q[0];
rz(0.89504129) q[1];
sx q[1];
rz(-1.4245234) q[1];
sx q[1];
rz(-0.76362124) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9485003) q[0];
sx q[0];
rz(-1.6512799) q[0];
sx q[0];
rz(1.3250687) q[0];
rz(-pi) q[1];
rz(-3.0020048) q[2];
sx q[2];
rz(-0.12171037) q[2];
sx q[2];
rz(-1.1523397) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29429193) q[1];
sx q[1];
rz(-0.2390034) q[1];
sx q[1];
rz(-1.6836749) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12327776) q[3];
sx q[3];
rz(-2.8486898) q[3];
sx q[3];
rz(-0.86859217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7274196) q[2];
sx q[2];
rz(-1.6215934) q[2];
sx q[2];
rz(-2.011389) q[2];
rz(-2.0027509) q[3];
sx q[3];
rz(-1.1915221) q[3];
sx q[3];
rz(1.1677008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8324757) q[0];
sx q[0];
rz(-2.583355) q[0];
sx q[0];
rz(-0.48156893) q[0];
rz(-2.8140977) q[1];
sx q[1];
rz(-1.4188674) q[1];
sx q[1];
rz(-0.96424261) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5580348) q[0];
sx q[0];
rz(-0.44981558) q[0];
sx q[0];
rz(1.0490473) q[0];
rz(0.19179253) q[2];
sx q[2];
rz(-0.99770498) q[2];
sx q[2];
rz(1.7157784) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.89470071) q[1];
sx q[1];
rz(-2.8671226) q[1];
sx q[1];
rz(2.0361221) q[1];
rz(-pi) q[2];
rz(-0.39059095) q[3];
sx q[3];
rz(-1.3814622) q[3];
sx q[3];
rz(0.51774294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1714736) q[2];
sx q[2];
rz(-2.0233266) q[2];
sx q[2];
rz(-1.4946651) q[2];
rz(-1.0129048) q[3];
sx q[3];
rz(-0.92618147) q[3];
sx q[3];
rz(-0.51500285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2038912) q[0];
sx q[0];
rz(-2.1370115) q[0];
sx q[0];
rz(1.0585693) q[0];
rz(-0.92264289) q[1];
sx q[1];
rz(-2.1244996) q[1];
sx q[1];
rz(-0.10291084) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12891586) q[0];
sx q[0];
rz(-0.8230477) q[0];
sx q[0];
rz(-0.91632592) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58892624) q[2];
sx q[2];
rz(-1.803033) q[2];
sx q[2];
rz(-3.1265772) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3427149) q[1];
sx q[1];
rz(-2.5763576) q[1];
sx q[1];
rz(1.2223627) q[1];
rz(0.0062557022) q[3];
sx q[3];
rz(-0.85438529) q[3];
sx q[3];
rz(-0.34230086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4395478) q[2];
sx q[2];
rz(-2.6921258) q[2];
sx q[2];
rz(0.31583819) q[2];
rz(2.1465007) q[3];
sx q[3];
rz(-1.8683878) q[3];
sx q[3];
rz(0.72527138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.0072667) q[0];
sx q[0];
rz(-2.8747929) q[0];
sx q[0];
rz(-2.5079978) q[0];
rz(-2.7569356) q[1];
sx q[1];
rz(-1.1793914) q[1];
sx q[1];
rz(0.42919174) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68394008) q[0];
sx q[0];
rz(-2.5999617) q[0];
sx q[0];
rz(-1.0628048) q[0];
x q[1];
rz(-1.5358938) q[2];
sx q[2];
rz(-2.2560824) q[2];
sx q[2];
rz(-1.8165782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.30503482) q[1];
sx q[1];
rz(-1.2188796) q[1];
sx q[1];
rz(-2.5266527) q[1];
rz(-pi) q[2];
rz(0.89520888) q[3];
sx q[3];
rz(-2.4662202) q[3];
sx q[3];
rz(2.7378367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0523494) q[2];
sx q[2];
rz(-2.3640335) q[2];
sx q[2];
rz(-0.26697978) q[2];
rz(-3.0698245) q[3];
sx q[3];
rz(-1.0177344) q[3];
sx q[3];
rz(1.7223541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.4778022) q[0];
sx q[0];
rz(-0.1013805) q[0];
sx q[0];
rz(1.4005533) q[0];
rz(-1.1446704) q[1];
sx q[1];
rz(-1.0619699) q[1];
sx q[1];
rz(-2.5544419) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42204866) q[0];
sx q[0];
rz(-1.1503771) q[0];
sx q[0];
rz(2.4499513) q[0];
rz(-2.7233974) q[2];
sx q[2];
rz(-1.4640239) q[2];
sx q[2];
rz(-0.51929861) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0441832) q[1];
sx q[1];
rz(-0.72260586) q[1];
sx q[1];
rz(-2.797762) q[1];
x q[2];
rz(-2.1903992) q[3];
sx q[3];
rz(-1.3096598) q[3];
sx q[3];
rz(-1.2984683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8553541) q[2];
sx q[2];
rz(-1.1796395) q[2];
sx q[2];
rz(-2.722495) q[2];
rz(-1.5239483) q[3];
sx q[3];
rz(-2.2231299) q[3];
sx q[3];
rz(-1.6379448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5986901) q[0];
sx q[0];
rz(-1.221523) q[0];
sx q[0];
rz(-0.30938095) q[0];
rz(1.3503374) q[1];
sx q[1];
rz(-0.58809909) q[1];
sx q[1];
rz(-0.5538756) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2387008) q[0];
sx q[0];
rz(-2.6052778) q[0];
sx q[0];
rz(2.4812921) q[0];
x q[1];
rz(-1.925681) q[2];
sx q[2];
rz(-2.7352374) q[2];
sx q[2];
rz(0.90824897) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1914158) q[1];
sx q[1];
rz(-1.8783057) q[1];
sx q[1];
rz(0.93788304) q[1];
rz(-pi) q[2];
rz(-2.43954) q[3];
sx q[3];
rz(-1.8115731) q[3];
sx q[3];
rz(-1.0853827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.89173335) q[2];
sx q[2];
rz(-1.6102108) q[2];
sx q[2];
rz(1.7943133) q[2];
rz(0.76752082) q[3];
sx q[3];
rz(-1.1859505) q[3];
sx q[3];
rz(0.84848136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42221853) q[0];
sx q[0];
rz(-2.8415871) q[0];
sx q[0];
rz(-1.2147709) q[0];
rz(2.1396554) q[1];
sx q[1];
rz(-1.6286214) q[1];
sx q[1];
rz(0.92438662) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11533212) q[0];
sx q[0];
rz(-2.1450274) q[0];
sx q[0];
rz(-0.46686812) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1154249) q[2];
sx q[2];
rz(-1.5473502) q[2];
sx q[2];
rz(-0.12870994) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.9634642) q[1];
sx q[1];
rz(-2.917082) q[1];
sx q[1];
rz(3.0369395) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2460292) q[3];
sx q[3];
rz(-0.82230869) q[3];
sx q[3];
rz(-3.0117717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7944916) q[2];
sx q[2];
rz(-2.2091986) q[2];
sx q[2];
rz(0.36716983) q[2];
rz(1.1651039) q[3];
sx q[3];
rz(-0.96787435) q[3];
sx q[3];
rz(2.2754106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9544871) q[0];
sx q[0];
rz(-0.73910537) q[0];
sx q[0];
rz(0.11866971) q[0];
rz(-0.26609303) q[1];
sx q[1];
rz(-2.2292744) q[1];
sx q[1];
rz(1.6899065) q[1];
rz(-2.8546147) q[2];
sx q[2];
rz(-1.8331883) q[2];
sx q[2];
rz(2.2856648) q[2];
rz(-0.20826505) q[3];
sx q[3];
rz(-2.1622445) q[3];
sx q[3];
rz(1.6643659) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
