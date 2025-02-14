OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8671829) q[0];
sx q[0];
rz(-2.6743968) q[0];
sx q[0];
rz(0.89611563) q[0];
rz(2.3236302) q[1];
sx q[1];
rz(-1.5939413) q[1];
sx q[1];
rz(-2.154248) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5817669) q[0];
sx q[0];
rz(-2.2004413) q[0];
sx q[0];
rz(-0.47578509) q[0];
x q[1];
rz(0.9573663) q[2];
sx q[2];
rz(-2.2920458) q[2];
sx q[2];
rz(2.3867875) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3707177) q[1];
sx q[1];
rz(-0.32568078) q[1];
sx q[1];
rz(-0.54906719) q[1];
x q[2];
rz(-0.39405502) q[3];
sx q[3];
rz(-0.39356227) q[3];
sx q[3];
rz(-1.2979899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2200615) q[2];
sx q[2];
rz(-2.107373) q[2];
sx q[2];
rz(0.49911487) q[2];
rz(-0.81579298) q[3];
sx q[3];
rz(-2.54839) q[3];
sx q[3];
rz(2.7844586) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4895184) q[0];
sx q[0];
rz(-2.7777785) q[0];
sx q[0];
rz(2.9079085) q[0];
rz(0.09985996) q[1];
sx q[1];
rz(-1.654511) q[1];
sx q[1];
rz(1.608009) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31029345) q[0];
sx q[0];
rz(-1.6673894) q[0];
sx q[0];
rz(-1.4486758) q[0];
x q[1];
rz(-0.43849053) q[2];
sx q[2];
rz(-1.4027565) q[2];
sx q[2];
rz(-2.4166256) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.93615426) q[1];
sx q[1];
rz(-0.9673276) q[1];
sx q[1];
rz(-3.0304069) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8321974) q[3];
sx q[3];
rz(-2.0338879) q[3];
sx q[3];
rz(0.17406305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0443772) q[2];
sx q[2];
rz(-2.2558687) q[2];
sx q[2];
rz(-3.0866747) q[2];
rz(0.6461668) q[3];
sx q[3];
rz(-1.61444) q[3];
sx q[3];
rz(3.0687148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9570479) q[0];
sx q[0];
rz(-2.8962729) q[0];
sx q[0];
rz(2.3967337) q[0];
rz(1.8648196) q[1];
sx q[1];
rz(-1.0294015) q[1];
sx q[1];
rz(0.863711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.058216393) q[0];
sx q[0];
rz(-1.6291219) q[0];
sx q[0];
rz(0.44022473) q[0];
rz(-pi) q[1];
rz(-2.3529007) q[2];
sx q[2];
rz(-1.1967107) q[2];
sx q[2];
rz(1.3918592) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59434592) q[1];
sx q[1];
rz(-0.9556831) q[1];
sx q[1];
rz(-2.1662461) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66074547) q[3];
sx q[3];
rz(-1.163394) q[3];
sx q[3];
rz(-1.416172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6700217) q[2];
sx q[2];
rz(-1.5315285) q[2];
sx q[2];
rz(-2.7460597) q[2];
rz(1.0223355) q[3];
sx q[3];
rz(-2.6774355) q[3];
sx q[3];
rz(2.6497604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96524298) q[0];
sx q[0];
rz(-1.1818161) q[0];
sx q[0];
rz(-0.055971948) q[0];
rz(-0.61198992) q[1];
sx q[1];
rz(-1.5925708) q[1];
sx q[1];
rz(2.2956119) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50091417) q[0];
sx q[0];
rz(-2.0862824) q[0];
sx q[0];
rz(0.2978931) q[0];
x q[1];
rz(1.0998639) q[2];
sx q[2];
rz(-2.467017) q[2];
sx q[2];
rz(-1.2686307) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8236134) q[1];
sx q[1];
rz(-2.7732964) q[1];
sx q[1];
rz(-0.55693926) q[1];
x q[2];
rz(-1.0958798) q[3];
sx q[3];
rz(-1.2327415) q[3];
sx q[3];
rz(-0.89382225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0116288) q[2];
sx q[2];
rz(-1.6946946) q[2];
sx q[2];
rz(-0.72032991) q[2];
rz(0.35308853) q[3];
sx q[3];
rz(-1.1938286) q[3];
sx q[3];
rz(-0.24875719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7406152) q[0];
sx q[0];
rz(-1.4530797) q[0];
sx q[0];
rz(-2.154696) q[0];
rz(-1.1445716) q[1];
sx q[1];
rz(-0.83895504) q[1];
sx q[1];
rz(-0.95485895) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78518822) q[0];
sx q[0];
rz(-2.4739728) q[0];
sx q[0];
rz(1.6130527) q[0];
rz(-pi) q[1];
rz(-2.0667162) q[2];
sx q[2];
rz(-2.5212147) q[2];
sx q[2];
rz(1.4324783) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.25453) q[1];
sx q[1];
rz(-1.1963486) q[1];
sx q[1];
rz(-2.6545054) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14820672) q[3];
sx q[3];
rz(-1.3199521) q[3];
sx q[3];
rz(2.776752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5130875) q[2];
sx q[2];
rz(-2.2272765) q[2];
sx q[2];
rz(2.0950192) q[2];
rz(-0.70932499) q[3];
sx q[3];
rz(-1.1449292) q[3];
sx q[3];
rz(2.5078702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3777622) q[0];
sx q[0];
rz(-1.7338294) q[0];
sx q[0];
rz(-2.8705226) q[0];
rz(0.079004869) q[1];
sx q[1];
rz(-0.43402356) q[1];
sx q[1];
rz(1.7105506) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1880497) q[0];
sx q[0];
rz(-0.255628) q[0];
sx q[0];
rz(0.60582692) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1151401) q[2];
sx q[2];
rz(-2.0664795) q[2];
sx q[2];
rz(-3.0718671) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.14306127) q[1];
sx q[1];
rz(-2.240671) q[1];
sx q[1];
rz(1.2550485) q[1];
x q[2];
rz(-1.1523139) q[3];
sx q[3];
rz(-1.7177204) q[3];
sx q[3];
rz(2.4148586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0926823) q[2];
sx q[2];
rz(-1.9581257) q[2];
sx q[2];
rz(0.023822039) q[2];
rz(1.228099) q[3];
sx q[3];
rz(-2.7619599) q[3];
sx q[3];
rz(-0.14321271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.989885) q[0];
sx q[0];
rz(-2.8360974) q[0];
sx q[0];
rz(3.049343) q[0];
rz(-0.30091885) q[1];
sx q[1];
rz(-1.2291127) q[1];
sx q[1];
rz(-2.0613964) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0940822) q[0];
sx q[0];
rz(-1.2258397) q[0];
sx q[0];
rz(3.0408919) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2713234) q[2];
sx q[2];
rz(-0.53127161) q[2];
sx q[2];
rz(0.16933717) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9321334) q[1];
sx q[1];
rz(-1.569063) q[1];
sx q[1];
rz(-2.9936547) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8870513) q[3];
sx q[3];
rz(-0.90083921) q[3];
sx q[3];
rz(-3.0681821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0098308) q[2];
sx q[2];
rz(-0.35085446) q[2];
sx q[2];
rz(0.62823137) q[2];
rz(0.96588165) q[3];
sx q[3];
rz(-1.5768257) q[3];
sx q[3];
rz(0.32111827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.0308762) q[0];
sx q[0];
rz(-2.4346011) q[0];
sx q[0];
rz(1.4068756) q[0];
rz(-3.0137317) q[1];
sx q[1];
rz(-1.7117932) q[1];
sx q[1];
rz(-2.0157287) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.073407452) q[0];
sx q[0];
rz(-2.3940384) q[0];
sx q[0];
rz(-2.8697467) q[0];
x q[1];
rz(3.1229805) q[2];
sx q[2];
rz(-2.0224704) q[2];
sx q[2];
rz(-1.9960595) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5347157) q[1];
sx q[1];
rz(-1.770853) q[1];
sx q[1];
rz(1.6354531) q[1];
rz(-pi) q[2];
rz(-2.5223658) q[3];
sx q[3];
rz(-1.9398749) q[3];
sx q[3];
rz(0.31821966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13059482) q[2];
sx q[2];
rz(-1.1329634) q[2];
sx q[2];
rz(3.0478743) q[2];
rz(1.7081918) q[3];
sx q[3];
rz(-2.8580557) q[3];
sx q[3];
rz(1.3933498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.691064) q[0];
sx q[0];
rz(-2.0557623) q[0];
sx q[0];
rz(-1.0040671) q[0];
rz(2.0380691) q[1];
sx q[1];
rz(-2.6595778) q[1];
sx q[1];
rz(1.7335266) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.69312) q[0];
sx q[0];
rz(-1.5768484) q[0];
sx q[0];
rz(-0.5151185) q[0];
rz(-pi) q[1];
rz(-1.4595217) q[2];
sx q[2];
rz(-1.8812875) q[2];
sx q[2];
rz(-1.7686012) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.50251657) q[1];
sx q[1];
rz(-2.6682819) q[1];
sx q[1];
rz(1.972354) q[1];
x q[2];
rz(0.086661913) q[3];
sx q[3];
rz(-0.52059697) q[3];
sx q[3];
rz(0.39256061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.39190009) q[2];
sx q[2];
rz(-1.6204648) q[2];
sx q[2];
rz(2.5938972) q[2];
rz(2.8294166) q[3];
sx q[3];
rz(-1.4145989) q[3];
sx q[3];
rz(-1.7267905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.3373229) q[0];
sx q[0];
rz(-1.6707358) q[0];
sx q[0];
rz(-0.7789337) q[0];
rz(-0.3745105) q[1];
sx q[1];
rz(-1.51314) q[1];
sx q[1];
rz(2.3096854) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13369416) q[0];
sx q[0];
rz(-1.3092894) q[0];
sx q[0];
rz(-0.71870872) q[0];
rz(0.86082117) q[2];
sx q[2];
rz(-0.92433017) q[2];
sx q[2];
rz(2.8614028) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.191997) q[1];
sx q[1];
rz(-2.1131385) q[1];
sx q[1];
rz(-0.80710141) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79672265) q[3];
sx q[3];
rz(-1.5490388) q[3];
sx q[3];
rz(-0.12783229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1067918) q[2];
sx q[2];
rz(-2.782244) q[2];
sx q[2];
rz(-3.1331151) q[2];
rz(0.40555412) q[3];
sx q[3];
rz(-1.6885992) q[3];
sx q[3];
rz(-1.4439553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1494898) q[0];
sx q[0];
rz(-1.533951) q[0];
sx q[0];
rz(0.78716192) q[0];
rz(-1.0943195) q[1];
sx q[1];
rz(-0.37847395) q[1];
sx q[1];
rz(2.7056221) q[1];
rz(1.0395861) q[2];
sx q[2];
rz(-0.24460228) q[2];
sx q[2];
rz(-1.2096418) q[2];
rz(1.5658436) q[3];
sx q[3];
rz(-1.4765783) q[3];
sx q[3];
rz(0.12259132) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
