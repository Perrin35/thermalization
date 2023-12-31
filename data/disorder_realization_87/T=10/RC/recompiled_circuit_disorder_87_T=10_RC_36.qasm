OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9053469) q[0];
sx q[0];
rz(-0.72609225) q[0];
sx q[0];
rz(-0.2015764) q[0];
rz(-2.6456614) q[1];
sx q[1];
rz(-2.6013241) q[1];
sx q[1];
rz(0.93710605) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0694323) q[0];
sx q[0];
rz(-1.9919112) q[0];
sx q[0];
rz(-2.518597) q[0];
rz(-0.65242282) q[2];
sx q[2];
rz(-2.4180275) q[2];
sx q[2];
rz(3.0384118) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.01552445) q[1];
sx q[1];
rz(-1.4926732) q[1];
sx q[1];
rz(-1.3874345) q[1];
x q[2];
rz(1.2308146) q[3];
sx q[3];
rz(-2.8570606) q[3];
sx q[3];
rz(-1.9576548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15930882) q[2];
sx q[2];
rz(-2.8757877) q[2];
sx q[2];
rz(-0.086159555) q[2];
rz(-0.75749767) q[3];
sx q[3];
rz(-2.3870654) q[3];
sx q[3];
rz(1.0936201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57698292) q[0];
sx q[0];
rz(-1.6596721) q[0];
sx q[0];
rz(1.8151059) q[0];
rz(-1.8857229) q[1];
sx q[1];
rz(-1.5763667) q[1];
sx q[1];
rz(-0.27145162) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87646987) q[0];
sx q[0];
rz(-2.070283) q[0];
sx q[0];
rz(0.45544099) q[0];
rz(-1.8797305) q[2];
sx q[2];
rz(-2.0663107) q[2];
sx q[2];
rz(-0.77906424) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6352306) q[1];
sx q[1];
rz(-2.2681384) q[1];
sx q[1];
rz(0.49763775) q[1];
rz(-0.80941697) q[3];
sx q[3];
rz(-1.3076233) q[3];
sx q[3];
rz(2.0078878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0945956) q[2];
sx q[2];
rz(-1.3826028) q[2];
sx q[2];
rz(2.5644152) q[2];
rz(-2.2180637) q[3];
sx q[3];
rz(-0.90463343) q[3];
sx q[3];
rz(-1.7318168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.8975824) q[0];
sx q[0];
rz(-1.0616466) q[0];
sx q[0];
rz(2.999021) q[0];
rz(1.3525195) q[1];
sx q[1];
rz(-2.1058154) q[1];
sx q[1];
rz(-2.9325063) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3554879) q[0];
sx q[0];
rz(-2.3007563) q[0];
sx q[0];
rz(-0.79746042) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4019037) q[2];
sx q[2];
rz(-1.2581173) q[2];
sx q[2];
rz(-2.3198421) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1413893) q[1];
sx q[1];
rz(-0.63265002) q[1];
sx q[1];
rz(0.81694095) q[1];
rz(-1.3554058) q[3];
sx q[3];
rz(-1.8672018) q[3];
sx q[3];
rz(2.4733558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3645939) q[2];
sx q[2];
rz(-0.89407388) q[2];
sx q[2];
rz(0.40412942) q[2];
rz(-1.2858307) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(0.16734853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4362815) q[0];
sx q[0];
rz(-3.0349338) q[0];
sx q[0];
rz(0.96281111) q[0];
rz(0.46936938) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(3.1406291) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92361802) q[0];
sx q[0];
rz(-1.0110564) q[0];
sx q[0];
rz(-2.2773507) q[0];
x q[1];
rz(0.12400603) q[2];
sx q[2];
rz(-1.5713072) q[2];
sx q[2];
rz(2.309547) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8648659) q[1];
sx q[1];
rz(-0.58502561) q[1];
sx q[1];
rz(-1.8086955) q[1];
x q[2];
rz(2.5721781) q[3];
sx q[3];
rz(-1.4953519) q[3];
sx q[3];
rz(-2.8358592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0659539) q[2];
sx q[2];
rz(-1.6299738) q[2];
sx q[2];
rz(-0.11165079) q[2];
rz(0.81104898) q[3];
sx q[3];
rz(-2.6871197) q[3];
sx q[3];
rz(3.1276935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1902996) q[0];
sx q[0];
rz(-0.90068978) q[0];
sx q[0];
rz(0.35650373) q[0];
rz(0.50645343) q[1];
sx q[1];
rz(-1.7849779) q[1];
sx q[1];
rz(2.8809663) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.168419) q[0];
sx q[0];
rz(-1.5289613) q[0];
sx q[0];
rz(2.9996458) q[0];
rz(2.8904301) q[2];
sx q[2];
rz(-1.8166944) q[2];
sx q[2];
rz(-2.9159301) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4958447) q[1];
sx q[1];
rz(-0.9965047) q[1];
sx q[1];
rz(1.5498284) q[1];
rz(-pi) q[2];
rz(-2.3652606) q[3];
sx q[3];
rz(-2.9388802) q[3];
sx q[3];
rz(2.4622038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2300718) q[2];
sx q[2];
rz(-1.6610961) q[2];
sx q[2];
rz(-0.22932209) q[2];
rz(-0.54245943) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(-2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3689573) q[0];
sx q[0];
rz(-2.2326523) q[0];
sx q[0];
rz(1.4916346) q[0];
rz(1.0391327) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(1.7274436) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4947858) q[0];
sx q[0];
rz(-1.8456869) q[0];
sx q[0];
rz(-1.6199714) q[0];
rz(-pi) q[1];
rz(-1.3041777) q[2];
sx q[2];
rz(-2.0732217) q[2];
sx q[2];
rz(-2.7555639) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1097088) q[1];
sx q[1];
rz(-1.110853) q[1];
sx q[1];
rz(-1.0201395) q[1];
x q[2];
rz(1.2392483) q[3];
sx q[3];
rz(-2.9132531) q[3];
sx q[3];
rz(2.6339298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1386537) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(-0.027739851) q[2];
rz(-2.6489143) q[3];
sx q[3];
rz(-1.6511107) q[3];
sx q[3];
rz(-1.3940575) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7572927) q[0];
sx q[0];
rz(-1.5526271) q[0];
sx q[0];
rz(0.12167715) q[0];
rz(-1.9901468) q[1];
sx q[1];
rz(-2.6897488) q[1];
sx q[1];
rz(-2.7391403) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2936195) q[0];
sx q[0];
rz(-1.3012039) q[0];
sx q[0];
rz(0.86975354) q[0];
rz(0.046026344) q[2];
sx q[2];
rz(-1.5791025) q[2];
sx q[2];
rz(-1.9866895) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.23558815) q[1];
sx q[1];
rz(-1.2885805) q[1];
sx q[1];
rz(1.3868335) q[1];
rz(1.4943069) q[3];
sx q[3];
rz(-2.7923931) q[3];
sx q[3];
rz(0.66222144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1272614) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(0.86223117) q[2];
rz(-2.6640653) q[3];
sx q[3];
rz(-1.9600441) q[3];
sx q[3];
rz(0.91807085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7858793) q[0];
sx q[0];
rz(-0.15545758) q[0];
sx q[0];
rz(0.73295897) q[0];
rz(-3.0015302) q[1];
sx q[1];
rz(-2.1439794) q[1];
sx q[1];
rz(-1.0345116) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93787545) q[0];
sx q[0];
rz(-2.8563742) q[0];
sx q[0];
rz(-0.9271778) q[0];
rz(-pi) q[1];
rz(-0.60136749) q[2];
sx q[2];
rz(-0.12430087) q[2];
sx q[2];
rz(2.3483495) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.58971436) q[1];
sx q[1];
rz(-1.0052181) q[1];
sx q[1];
rz(-0.043885529) q[1];
rz(-pi) q[2];
rz(-1.4626059) q[3];
sx q[3];
rz(-1.9868317) q[3];
sx q[3];
rz(-3.0707404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.90116477) q[2];
sx q[2];
rz(-1.9463836) q[2];
sx q[2];
rz(2.365716) q[2];
rz(2.4173229) q[3];
sx q[3];
rz(-2.3359719) q[3];
sx q[3];
rz(-1.4340713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4416606) q[0];
sx q[0];
rz(-2.3775546) q[0];
sx q[0];
rz(-3.124776) q[0];
rz(-0.018521221) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(-2.3628078) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9176365) q[0];
sx q[0];
rz(-2.7178239) q[0];
sx q[0];
rz(-2.5013322) q[0];
x q[1];
rz(-2.8674556) q[2];
sx q[2];
rz(-2.0275653) q[2];
sx q[2];
rz(2.0704839) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15021819) q[1];
sx q[1];
rz(-2.7017936) q[1];
sx q[1];
rz(-0.99779731) q[1];
rz(-1.7917463) q[3];
sx q[3];
rz(-1.5077935) q[3];
sx q[3];
rz(0.27730478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6918216) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(-0.71643913) q[2];
rz(1.6843494) q[3];
sx q[3];
rz(-1.4102035) q[3];
sx q[3];
rz(-1.8523857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.107782) q[0];
sx q[0];
rz(-0.53492117) q[0];
sx q[0];
rz(-0.15144908) q[0];
rz(0.17627136) q[1];
sx q[1];
rz(-1.1947894) q[1];
sx q[1];
rz(0.72296468) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4484946) q[0];
sx q[0];
rz(-1.7211595) q[0];
sx q[0];
rz(-3.0340414) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1935913) q[2];
sx q[2];
rz(-1.4133487) q[2];
sx q[2];
rz(2.8709656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1180229) q[1];
sx q[1];
rz(-1.6260864) q[1];
sx q[1];
rz(-2.0748595) q[1];
rz(-pi) q[2];
rz(2.74182) q[3];
sx q[3];
rz(-0.82953605) q[3];
sx q[3];
rz(-1.6403891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4621949) q[2];
sx q[2];
rz(-1.858819) q[2];
sx q[2];
rz(-2.6160713) q[2];
rz(-0.28371352) q[3];
sx q[3];
rz(-2.0339537) q[3];
sx q[3];
rz(-2.7450558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5491966) q[0];
sx q[0];
rz(-2.017673) q[0];
sx q[0];
rz(2.8155433) q[0];
rz(-1.0992959) q[1];
sx q[1];
rz(-2.9922843) q[1];
sx q[1];
rz(-0.86984787) q[1];
rz(2.1500258) q[2];
sx q[2];
rz(-0.75800037) q[2];
sx q[2];
rz(-2.6108685) q[2];
rz(-0.59393926) q[3];
sx q[3];
rz(-0.38729061) q[3];
sx q[3];
rz(-0.95872986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
