OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.23624578) q[0];
sx q[0];
rz(-2.4155004) q[0];
sx q[0];
rz(-2.9400163) q[0];
rz(-2.6456614) q[1];
sx q[1];
rz(-2.6013241) q[1];
sx q[1];
rz(0.93710605) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0694323) q[0];
sx q[0];
rz(-1.9919112) q[0];
sx q[0];
rz(-0.62299563) q[0];
rz(-1.0785157) q[2];
sx q[2];
rz(-1.0168076) q[2];
sx q[2];
rz(0.89821494) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1569251) q[1];
sx q[1];
rz(-2.942454) q[1];
sx q[1];
rz(-1.1652633) q[1];
rz(-pi) q[2];
rz(-1.2308146) q[3];
sx q[3];
rz(-2.8570606) q[3];
sx q[3];
rz(1.9576548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9822838) q[2];
sx q[2];
rz(-2.8757877) q[2];
sx q[2];
rz(-3.0554331) q[2];
rz(0.75749767) q[3];
sx q[3];
rz(-2.3870654) q[3];
sx q[3];
rz(-1.0936201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.5646097) q[0];
sx q[0];
rz(-1.4819205) q[0];
sx q[0];
rz(1.3264867) q[0];
rz(-1.8857229) q[1];
sx q[1];
rz(-1.5763667) q[1];
sx q[1];
rz(2.870141) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6730246) q[0];
sx q[0];
rz(-0.66267555) q[0];
sx q[0];
rz(-2.2492692) q[0];
rz(1.2618622) q[2];
sx q[2];
rz(-2.0663107) q[2];
sx q[2];
rz(-0.77906424) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6352306) q[1];
sx q[1];
rz(-2.2681384) q[1];
sx q[1];
rz(0.49763775) q[1];
rz(1.1984962) q[3];
sx q[3];
rz(-0.79685235) q[3];
sx q[3];
rz(-0.70355319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0945956) q[2];
sx q[2];
rz(-1.7589898) q[2];
sx q[2];
rz(2.5644152) q[2];
rz(2.2180637) q[3];
sx q[3];
rz(-2.2369592) q[3];
sx q[3];
rz(-1.7318168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8975824) q[0];
sx q[0];
rz(-1.0616466) q[0];
sx q[0];
rz(0.14257167) q[0];
rz(1.7890731) q[1];
sx q[1];
rz(-2.1058154) q[1];
sx q[1];
rz(-0.20908633) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36172983) q[0];
sx q[0];
rz(-1.0233101) q[0];
sx q[0];
rz(-2.2453383) q[0];
rz(-pi) q[1];
rz(0.44720165) q[2];
sx q[2];
rz(-2.3502091) q[2];
sx q[2];
rz(-2.7176822) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0641388) q[1];
sx q[1];
rz(-1.1541379) q[1];
sx q[1];
rz(-1.0799079) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.3029996) q[3];
sx q[3];
rz(-1.7766561) q[3];
sx q[3];
rz(0.96637615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3645939) q[2];
sx q[2];
rz(-0.89407388) q[2];
sx q[2];
rz(-2.7374632) q[2];
rz(1.8557619) q[3];
sx q[3];
rz(-2.0127681) q[3];
sx q[3];
rz(-2.9742441) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7053112) q[0];
sx q[0];
rz(-0.10665882) q[0];
sx q[0];
rz(-2.1787815) q[0];
rz(-0.46936938) q[1];
sx q[1];
rz(-0.58987394) q[1];
sx q[1];
rz(0.00096360047) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2036635) q[0];
sx q[0];
rz(-2.2712049) q[0];
sx q[0];
rz(0.80313375) q[0];
rz(-pi) q[1];
rz(0.12400603) q[2];
sx q[2];
rz(-1.5702855) q[2];
sx q[2];
rz(0.83204568) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.094593781) q[1];
sx q[1];
rz(-1.440289) q[1];
sx q[1];
rz(-0.99884896) q[1];
rz(-0.56941454) q[3];
sx q[3];
rz(-1.6462407) q[3];
sx q[3];
rz(2.8358592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0756388) q[2];
sx q[2];
rz(-1.6299738) q[2];
sx q[2];
rz(0.11165079) q[2];
rz(-0.81104898) q[3];
sx q[3];
rz(-0.45447293) q[3];
sx q[3];
rz(-0.013899175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.951293) q[0];
sx q[0];
rz(-0.90068978) q[0];
sx q[0];
rz(0.35650373) q[0];
rz(0.50645343) q[1];
sx q[1];
rz(-1.7849779) q[1];
sx q[1];
rz(2.8809663) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97317364) q[0];
sx q[0];
rz(-1.6126313) q[0];
sx q[0];
rz(0.14194685) q[0];
x q[1];
rz(1.8243276) q[2];
sx q[2];
rz(-1.8142482) q[2];
sx q[2];
rz(-1.4075116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0780328) q[1];
sx q[1];
rz(-1.5884001) q[1];
sx q[1];
rz(-2.5672008) q[1];
x q[2];
rz(-0.77633206) q[3];
sx q[3];
rz(-2.9388802) q[3];
sx q[3];
rz(0.67938882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2300718) q[2];
sx q[2];
rz(-1.6610961) q[2];
sx q[2];
rz(-0.22932209) q[2];
rz(-2.5991332) q[3];
sx q[3];
rz(-0.31033236) q[3];
sx q[3];
rz(2.1054161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3689573) q[0];
sx q[0];
rz(-0.90894037) q[0];
sx q[0];
rz(1.4916346) q[0];
rz(-2.1024599) q[1];
sx q[1];
rz(-1.2960478) q[1];
sx q[1];
rz(-1.414149) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4947858) q[0];
sx q[0];
rz(-1.2959058) q[0];
sx q[0];
rz(1.6199714) q[0];
rz(-pi) q[1];
rz(2.6238407) q[2];
sx q[2];
rz(-1.8038097) q[2];
sx q[2];
rz(2.0875967) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8687739) q[1];
sx q[1];
rz(-2.0588074) q[1];
sx q[1];
rz(-2.6150319) q[1];
rz(-pi) q[2];
rz(-3.066091) q[3];
sx q[3];
rz(-1.3550948) q[3];
sx q[3];
rz(0.84738934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1386537) q[2];
sx q[2];
rz(-2.5230375) q[2];
sx q[2];
rz(3.1138528) q[2];
rz(-0.49267832) q[3];
sx q[3];
rz(-1.6511107) q[3];
sx q[3];
rz(-1.7475351) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3843) q[0];
sx q[0];
rz(-1.5526271) q[0];
sx q[0];
rz(-3.0199155) q[0];
rz(-1.9901468) q[1];
sx q[1];
rz(-0.45184389) q[1];
sx q[1];
rz(-0.40245232) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6399022) q[0];
sx q[0];
rz(-2.2417289) q[0];
sx q[0];
rz(0.34696607) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5624814) q[2];
sx q[2];
rz(-1.6168211) q[2];
sx q[2];
rz(2.7260821) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.35343364) q[1];
sx q[1];
rz(-0.33553365) q[1];
sx q[1];
rz(0.56281705) q[1];
x q[2];
rz(3.1137755) q[3];
sx q[3];
rz(-1.9189315) q[3];
sx q[3];
rz(-0.58084014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1272614) q[2];
sx q[2];
rz(-1.1703337) q[2];
sx q[2];
rz(-2.2793615) q[2];
rz(2.6640653) q[3];
sx q[3];
rz(-1.1815485) q[3];
sx q[3];
rz(0.91807085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35571337) q[0];
sx q[0];
rz(-2.9861351) q[0];
sx q[0];
rz(2.4086337) q[0];
rz(-3.0015302) q[1];
sx q[1];
rz(-2.1439794) q[1];
sx q[1];
rz(2.1070811) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0089793423) q[0];
sx q[0];
rz(-1.4011369) q[0];
sx q[0];
rz(1.8011814) q[0];
rz(-1.5002235) q[2];
sx q[2];
rz(-1.6732054) q[2];
sx q[2];
rz(0.18825738) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1840399) q[1];
sx q[1];
rz(-1.5337481) q[1];
sx q[1];
rz(-2.1368105) q[1];
rz(-pi) q[2];
rz(1.6789867) q[3];
sx q[3];
rz(-1.1547609) q[3];
sx q[3];
rz(3.0707404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.90116477) q[2];
sx q[2];
rz(-1.1952091) q[2];
sx q[2];
rz(0.77587664) q[2];
rz(-2.4173229) q[3];
sx q[3];
rz(-2.3359719) q[3];
sx q[3];
rz(-1.7075214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6999321) q[0];
sx q[0];
rz(-0.76403809) q[0];
sx q[0];
rz(-3.124776) q[0];
rz(3.1230714) q[1];
sx q[1];
rz(-1.8073945) q[1];
sx q[1];
rz(0.7787849) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6803857) q[0];
sx q[0];
rz(-1.9068423) q[0];
sx q[0];
rz(1.3075605) q[0];
x q[1];
rz(2.8674556) q[2];
sx q[2];
rz(-1.1140274) q[2];
sx q[2];
rz(2.0704839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15021819) q[1];
sx q[1];
rz(-2.7017936) q[1];
sx q[1];
rz(2.1437953) q[1];
rz(-pi) q[2];
rz(-0.064568297) q[3];
sx q[3];
rz(-1.350292) q[3];
sx q[3];
rz(1.2793503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4497711) q[2];
sx q[2];
rz(-1.8372953) q[2];
sx q[2];
rz(-2.4251535) q[2];
rz(-1.6843494) q[3];
sx q[3];
rz(-1.7313892) q[3];
sx q[3];
rz(1.289207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.107782) q[0];
sx q[0];
rz(-2.6066715) q[0];
sx q[0];
rz(0.15144908) q[0];
rz(2.9653213) q[1];
sx q[1];
rz(-1.1947894) q[1];
sx q[1];
rz(-0.72296468) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6930981) q[0];
sx q[0];
rz(-1.4204331) q[0];
sx q[0];
rz(-0.10755121) q[0];
x q[1];
rz(1.1638219) q[2];
sx q[2];
rz(-2.7343035) q[2];
sx q[2];
rz(-0.92330698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7192917) q[1];
sx q[1];
rz(-2.0740168) q[1];
sx q[1];
rz(-3.0784688) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9728174) q[3];
sx q[3];
rz(-2.3178181) q[3];
sx q[3];
rz(1.0812425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.67939776) q[2];
sx q[2];
rz(-1.858819) q[2];
sx q[2];
rz(-0.52552137) q[2];
rz(0.28371352) q[3];
sx q[3];
rz(-2.0339537) q[3];
sx q[3];
rz(-0.39653683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59239607) q[0];
sx q[0];
rz(-2.017673) q[0];
sx q[0];
rz(2.8155433) q[0];
rz(-2.0422968) q[1];
sx q[1];
rz(-0.14930832) q[1];
sx q[1];
rz(2.2717448) q[1];
rz(2.2407871) q[2];
sx q[2];
rz(-1.1849891) q[2];
sx q[2];
rz(1.6580788) q[2];
rz(2.8156149) q[3];
sx q[3];
rz(-1.3578284) q[3];
sx q[3];
rz(1.1708543) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
