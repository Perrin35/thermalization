OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3888336) q[0];
sx q[0];
rz(-1.692481) q[0];
sx q[0];
rz(-1.4504855) q[0];
rz(-5.3095498) q[1];
sx q[1];
rz(7.7205478) q[1];
sx q[1];
rz(7.2024495) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1796003) q[0];
sx q[0];
rz(-1.5536947) q[0];
sx q[0];
rz(-3.1110065) q[0];
x q[1];
rz(-2.4953142) q[2];
sx q[2];
rz(-2.9913104) q[2];
sx q[2];
rz(2.8719547) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31499681) q[1];
sx q[1];
rz(-2.2758141) q[1];
sx q[1];
rz(-0.12048529) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2926225) q[3];
sx q[3];
rz(-2.5458286) q[3];
sx q[3];
rz(1.9879544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3322525) q[2];
sx q[2];
rz(-1.4490178) q[2];
sx q[2];
rz(1.4130886) q[2];
rz(0.20279065) q[3];
sx q[3];
rz(-1.7594124) q[3];
sx q[3];
rz(-3.0387759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24421144) q[0];
sx q[0];
rz(-2.057071) q[0];
sx q[0];
rz(0.71075034) q[0];
rz(2.3731025) q[1];
sx q[1];
rz(-1.0667421) q[1];
sx q[1];
rz(-1.01952) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7026414) q[0];
sx q[0];
rz(-2.0176689) q[0];
sx q[0];
rz(-1.6146445) q[0];
rz(-pi) q[1];
rz(-0.076093535) q[2];
sx q[2];
rz(-2.00092) q[2];
sx q[2];
rz(-0.96904749) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2166357) q[1];
sx q[1];
rz(-1.5397738) q[1];
sx q[1];
rz(1.9418632) q[1];
rz(-pi) q[2];
rz(-3.0583303) q[3];
sx q[3];
rz(-1.0181277) q[3];
sx q[3];
rz(3.0577615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3780313) q[2];
sx q[2];
rz(-1.2445933) q[2];
sx q[2];
rz(-1.9192609) q[2];
rz(1.2126806) q[3];
sx q[3];
rz(-2.1247532) q[3];
sx q[3];
rz(0.71162629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0843622) q[0];
sx q[0];
rz(-2.1570692) q[0];
sx q[0];
rz(1.2179751) q[0];
rz(0.51586622) q[1];
sx q[1];
rz(-0.54793826) q[1];
sx q[1];
rz(2.2191494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24127125) q[0];
sx q[0];
rz(-1.5301955) q[0];
sx q[0];
rz(-0.057828219) q[0];
rz(-pi) q[1];
rz(0.5567906) q[2];
sx q[2];
rz(-2.3344731) q[2];
sx q[2];
rz(3.0592164) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6092012) q[1];
sx q[1];
rz(-0.89840404) q[1];
sx q[1];
rz(2.3228541) q[1];
rz(-pi) q[2];
rz(-2.6785128) q[3];
sx q[3];
rz(-2.7768917) q[3];
sx q[3];
rz(1.1187584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1232274) q[2];
sx q[2];
rz(-1.0845228) q[2];
sx q[2];
rz(1.8017192) q[2];
rz(-2.5943622) q[3];
sx q[3];
rz(-2.7180143) q[3];
sx q[3];
rz(0.71119285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.055534) q[0];
sx q[0];
rz(-0.14500293) q[0];
sx q[0];
rz(-2.9823629) q[0];
rz(3.1314462) q[1];
sx q[1];
rz(-2.1187014) q[1];
sx q[1];
rz(-0.25746447) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5437357) q[0];
sx q[0];
rz(-0.86559767) q[0];
sx q[0];
rz(-1.0969093) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5045549) q[2];
sx q[2];
rz(-0.77755723) q[2];
sx q[2];
rz(0.8929127) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60002335) q[1];
sx q[1];
rz(-2.5065055) q[1];
sx q[1];
rz(2.9734008) q[1];
rz(-pi) q[2];
rz(2.5268528) q[3];
sx q[3];
rz(-2.2343544) q[3];
sx q[3];
rz(-0.033294423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3782392) q[2];
sx q[2];
rz(-1.2946318) q[2];
sx q[2];
rz(-0.47284687) q[2];
rz(2.4371448) q[3];
sx q[3];
rz(-1.364578) q[3];
sx q[3];
rz(0.64594597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6351629) q[0];
sx q[0];
rz(-1.9671054) q[0];
sx q[0];
rz(-0.065486431) q[0];
rz(2.7323515) q[1];
sx q[1];
rz(-1.1301273) q[1];
sx q[1];
rz(-1.7154891) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55520362) q[0];
sx q[0];
rz(-0.23815933) q[0];
sx q[0];
rz(0.40443964) q[0];
rz(2.9272396) q[2];
sx q[2];
rz(-0.39696908) q[2];
sx q[2];
rz(-3.1052542) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3196484) q[1];
sx q[1];
rz(-1.2927755) q[1];
sx q[1];
rz(1.8026428) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17284837) q[3];
sx q[3];
rz(-1.2842872) q[3];
sx q[3];
rz(0.89927538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19491974) q[2];
sx q[2];
rz(-2.0544923) q[2];
sx q[2];
rz(1.3067513) q[2];
rz(1.9715747) q[3];
sx q[3];
rz(-2.7368059) q[3];
sx q[3];
rz(-3.034333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4829247) q[0];
sx q[0];
rz(-1.985745) q[0];
sx q[0];
rz(-2.7401127) q[0];
rz(0.67277706) q[1];
sx q[1];
rz(-2.3428226) q[1];
sx q[1];
rz(-0.85404095) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5082762) q[0];
sx q[0];
rz(-0.22249732) q[0];
sx q[0];
rz(-1.359324) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9561192) q[2];
sx q[2];
rz(-1.7181686) q[2];
sx q[2];
rz(0.81306785) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2931765) q[1];
sx q[1];
rz(-0.59976116) q[1];
sx q[1];
rz(0.17677115) q[1];
x q[2];
rz(2.5307114) q[3];
sx q[3];
rz(-1.1958836) q[3];
sx q[3];
rz(1.5883816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4868769) q[2];
sx q[2];
rz(-2.2701023) q[2];
sx q[2];
rz(1.1775449) q[2];
rz(-0.55142895) q[3];
sx q[3];
rz(-1.5827554) q[3];
sx q[3];
rz(2.3059755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7198782) q[0];
sx q[0];
rz(-0.28110176) q[0];
sx q[0];
rz(-1.1623435) q[0];
rz(-1.989919) q[1];
sx q[1];
rz(-1.3538227) q[1];
sx q[1];
rz(-2.2231359) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8480523) q[0];
sx q[0];
rz(-1.2402724) q[0];
sx q[0];
rz(2.8577515) q[0];
rz(-pi) q[1];
rz(0.8753885) q[2];
sx q[2];
rz(-1.3686485) q[2];
sx q[2];
rz(3.1234158) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7734997) q[1];
sx q[1];
rz(-1.5640904) q[1];
sx q[1];
rz(-0.0073324629) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9145394) q[3];
sx q[3];
rz(-1.731428) q[3];
sx q[3];
rz(0.70987046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.034417001) q[2];
sx q[2];
rz(-1.3849266) q[2];
sx q[2];
rz(2.6386063) q[2];
rz(0.77477396) q[3];
sx q[3];
rz(-2.6796902) q[3];
sx q[3];
rz(-1.3431965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4485432) q[0];
sx q[0];
rz(-0.22668426) q[0];
sx q[0];
rz(-1.6402798) q[0];
rz(-0.34128183) q[1];
sx q[1];
rz(-1.4309859) q[1];
sx q[1];
rz(2.0223845) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2687896) q[0];
sx q[0];
rz(-0.69376341) q[0];
sx q[0];
rz(2.540178) q[0];
rz(-2.881024) q[2];
sx q[2];
rz(-1.1036647) q[2];
sx q[2];
rz(1.9595944) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.70485605) q[1];
sx q[1];
rz(-2.4697127) q[1];
sx q[1];
rz(-0.19499548) q[1];
rz(-pi) q[2];
rz(1.7057034) q[3];
sx q[3];
rz(-1.9622318) q[3];
sx q[3];
rz(1.6329444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.98562733) q[2];
sx q[2];
rz(-2.1817744) q[2];
sx q[2];
rz(0.36925527) q[2];
rz(-2.8333832) q[3];
sx q[3];
rz(-2.1741368) q[3];
sx q[3];
rz(-1.5306028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8592598) q[0];
sx q[0];
rz(-1.3066602) q[0];
sx q[0];
rz(-0.34307137) q[0];
rz(1.9725017) q[1];
sx q[1];
rz(-1.7770551) q[1];
sx q[1];
rz(1.7128568) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.13641) q[0];
sx q[0];
rz(-0.47629582) q[0];
sx q[0];
rz(0.57203697) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83447225) q[2];
sx q[2];
rz(-1.021046) q[2];
sx q[2];
rz(0.85981926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.27738849) q[1];
sx q[1];
rz(-2.7608747) q[1];
sx q[1];
rz(-1.4053132) q[1];
x q[2];
rz(-0.87782209) q[3];
sx q[3];
rz(-1.9435427) q[3];
sx q[3];
rz(2.8543548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2009361) q[2];
sx q[2];
rz(-0.20874615) q[2];
sx q[2];
rz(1.3915871) q[2];
rz(2.7348147) q[3];
sx q[3];
rz(-1.4429561) q[3];
sx q[3];
rz(-2.2627635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0402891) q[0];
sx q[0];
rz(-0.60281301) q[0];
sx q[0];
rz(1.3579177) q[0];
rz(-1.2376002) q[1];
sx q[1];
rz(-1.0126746) q[1];
sx q[1];
rz(-2.3416669) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023066085) q[0];
sx q[0];
rz(-0.79341054) q[0];
sx q[0];
rz(1.8961468) q[0];
rz(2.5913057) q[2];
sx q[2];
rz(-1.3105416) q[2];
sx q[2];
rz(-1.8837007) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.057803072) q[1];
sx q[1];
rz(-2.4531595) q[1];
sx q[1];
rz(-2.8958984) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8915845) q[3];
sx q[3];
rz(-2.0431314) q[3];
sx q[3];
rz(-1.8369305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.37821975) q[2];
sx q[2];
rz(-1.7174481) q[2];
sx q[2];
rz(-0.073089449) q[2];
rz(2.8857005) q[3];
sx q[3];
rz(-0.81232324) q[3];
sx q[3];
rz(-1.488744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8681317) q[0];
sx q[0];
rz(-1.68597) q[0];
sx q[0];
rz(-1.2706533) q[0];
rz(-1.801626) q[1];
sx q[1];
rz(-1.3506964) q[1];
sx q[1];
rz(0.66257308) q[1];
rz(0.70941464) q[2];
sx q[2];
rz(-1.5723036) q[2];
sx q[2];
rz(0.2607762) q[2];
rz(-0.64107278) q[3];
sx q[3];
rz(-1.0985634) q[3];
sx q[3];
rz(-2.6177277) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
