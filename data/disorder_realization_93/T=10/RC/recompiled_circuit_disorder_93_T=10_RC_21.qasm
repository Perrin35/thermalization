OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4424326) q[0];
sx q[0];
rz(-1.3843098) q[0];
sx q[0];
rz(-1.260489) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(-1.7927875) q[1];
sx q[1];
rz(-0.92372149) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7605654) q[0];
sx q[0];
rz(-1.1228704) q[0];
sx q[0];
rz(2.7678124) q[0];
rz(-1.2878296) q[2];
sx q[2];
rz(-1.7093061) q[2];
sx q[2];
rz(1.1793009) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0337861) q[1];
sx q[1];
rz(-0.68192712) q[1];
sx q[1];
rz(2.4297907) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.190891) q[3];
sx q[3];
rz(-2.0348747) q[3];
sx q[3];
rz(-2.3604148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.91360056) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(-2.9795734) q[2];
rz(2.2062733) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(-2.4285765) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0682003) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(-1.1967999) q[0];
rz(-0.67990047) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(1.4555567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.992422) q[0];
sx q[0];
rz(-0.99827168) q[0];
sx q[0];
rz(2.9426129) q[0];
rz(-pi) q[1];
rz(0.37462072) q[2];
sx q[2];
rz(-1.4943559) q[2];
sx q[2];
rz(-0.74707109) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.93745366) q[1];
sx q[1];
rz(-2.0935645) q[1];
sx q[1];
rz(0.015603113) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5973813) q[3];
sx q[3];
rz(-1.9823325) q[3];
sx q[3];
rz(-2.6708024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.42852795) q[2];
sx q[2];
rz(-1.6899127) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(0.18243608) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3882554) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(2.341111) q[0];
rz(3.1128186) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(-1.9690537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5867509) q[0];
sx q[0];
rz(-1.8790073) q[0];
sx q[0];
rz(0.59535938) q[0];
x q[1];
rz(1.8981947) q[2];
sx q[2];
rz(-2.2990169) q[2];
sx q[2];
rz(-2.3936405) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.4634358) q[1];
sx q[1];
rz(-1.0532182) q[1];
sx q[1];
rz(0.28097681) q[1];
rz(-pi) q[2];
rz(-0.058733744) q[3];
sx q[3];
rz(-1.8211094) q[3];
sx q[3];
rz(-2.5854923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0744434) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(-2.2303936) q[2];
rz(-2.1905812) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(-0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7610385) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(0.12810853) q[0];
rz(0.076106636) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(0.52350837) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0353161) q[0];
sx q[0];
rz(-1.9395394) q[0];
sx q[0];
rz(-0.62555255) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8115494) q[2];
sx q[2];
rz(-1.4903307) q[2];
sx q[2];
rz(-0.45804322) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7586786) q[1];
sx q[1];
rz(-0.74712336) q[1];
sx q[1];
rz(1.1927356) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13682271) q[3];
sx q[3];
rz(-0.98383437) q[3];
sx q[3];
rz(0.50859355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52544242) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(0.564044) q[2];
rz(0.28856746) q[3];
sx q[3];
rz(-0.42268649) q[3];
sx q[3];
rz(-0.55571663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600835) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(-1.6500641) q[0];
rz(-0.87961698) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(-0.99194828) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.49682) q[0];
sx q[0];
rz(-2.9368375) q[0];
sx q[0];
rz(0.55069189) q[0];
rz(-3.132658) q[2];
sx q[2];
rz(-1.8736471) q[2];
sx q[2];
rz(2.6645899) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.93047749) q[1];
sx q[1];
rz(-2.8301297) q[1];
sx q[1];
rz(0.13179563) q[1];
rz(-pi) q[2];
rz(0.12985142) q[3];
sx q[3];
rz(-0.81749812) q[3];
sx q[3];
rz(-1.070147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1296967) q[2];
sx q[2];
rz(-0.36964881) q[2];
sx q[2];
rz(2.8707855) q[2];
rz(-2.9233542) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(2.9158084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72702423) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(1.7927992) q[0];
rz(-2.7596966) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(-1.4250925) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2857367) q[0];
sx q[0];
rz(-0.55186134) q[0];
sx q[0];
rz(0.43666552) q[0];
rz(-pi) q[1];
rz(3.0503057) q[2];
sx q[2];
rz(-1.8249776) q[2];
sx q[2];
rz(-0.19677481) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1849991) q[1];
sx q[1];
rz(-1.5041659) q[1];
sx q[1];
rz(1.1207629) q[1];
rz(-2.1720042) q[3];
sx q[3];
rz(-2.0242656) q[3];
sx q[3];
rz(2.3525402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1340593) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(-0.56387222) q[2];
rz(2.9610736) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6329704) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(-2.7222743) q[0];
rz(-1.5527027) q[1];
sx q[1];
rz(-1.8808552) q[1];
sx q[1];
rz(0.82180506) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2215295) q[0];
sx q[0];
rz(-0.95525817) q[0];
sx q[0];
rz(-2.2277742) q[0];
x q[1];
rz(2.6452438) q[2];
sx q[2];
rz(-2.2349572) q[2];
sx q[2];
rz(-1.4585444) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2591178) q[1];
sx q[1];
rz(-1.7977409) q[1];
sx q[1];
rz(-2.4005753) q[1];
rz(-pi) q[2];
rz(0.52845593) q[3];
sx q[3];
rz(-0.65675694) q[3];
sx q[3];
rz(0.17880759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3372779) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(-0.24469963) q[2];
rz(0.129536) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(-1.6285508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.725175) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(2.3186671) q[0];
rz(2.8322463) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(1.8364505) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96852036) q[0];
sx q[0];
rz(-2.1497796) q[0];
sx q[0];
rz(1.9645683) q[0];
rz(2.4370286) q[2];
sx q[2];
rz(-1.8578055) q[2];
sx q[2];
rz(-1.2003843) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0481938) q[1];
sx q[1];
rz(-1.1953925) q[1];
sx q[1];
rz(-2.4673389) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5042138) q[3];
sx q[3];
rz(-0.3399907) q[3];
sx q[3];
rz(0.41223994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(1.4902327) q[2];
rz(-1.0772609) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(3.0100477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3867144) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(-0.3219147) q[0];
rz(1.5362668) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(-0.70294356) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2962869) q[0];
sx q[0];
rz(-1.1830813) q[0];
sx q[0];
rz(1.182343) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1019194) q[2];
sx q[2];
rz(-1.0618321) q[2];
sx q[2];
rz(-0.85658011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1785537) q[1];
sx q[1];
rz(-0.56561618) q[1];
sx q[1];
rz(-1.2600684) q[1];
rz(-pi) q[2];
rz(1.1144936) q[3];
sx q[3];
rz(-2.4745686) q[3];
sx q[3];
rz(-1.7175355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2408509) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(-1.8019603) q[2];
rz(-2.83589) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(-1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9777578) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(1.9158069) q[0];
rz(-0.90351358) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(0.46863619) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43384957) q[0];
sx q[0];
rz(-1.187547) q[0];
sx q[0];
rz(1.7953403) q[0];
rz(-pi) q[1];
rz(0.87601985) q[2];
sx q[2];
rz(-2.6555736) q[2];
sx q[2];
rz(1.8234058) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.74433078) q[1];
sx q[1];
rz(-1.8363734) q[1];
sx q[1];
rz(0.070752146) q[1];
rz(-pi) q[2];
rz(2.7077984) q[3];
sx q[3];
rz(-1.1489831) q[3];
sx q[3];
rz(2.1123561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.48352155) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(-1.1432077) q[2];
rz(-3.0269567) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(-1.5293998) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951915) q[0];
sx q[0];
rz(-1.9468745) q[0];
sx q[0];
rz(-0.68328802) q[0];
rz(-0.62190965) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(-2.4516104) q[2];
sx q[2];
rz(-0.98914115) q[2];
sx q[2];
rz(3.0422899) q[2];
rz(-1.6735531) q[3];
sx q[3];
rz(-2.4874874) q[3];
sx q[3];
rz(-2.7174674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
