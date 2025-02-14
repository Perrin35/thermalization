OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.355964) q[0];
sx q[0];
rz(-2.845919) q[0];
sx q[0];
rz(2.7276738) q[0];
rz(3.4105372) q[1];
sx q[1];
rz(6.6528448) q[1];
sx q[1];
rz(11.289968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39525014) q[0];
sx q[0];
rz(-1.6621886) q[0];
sx q[0];
rz(-1.4797828) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90246713) q[2];
sx q[2];
rz(-1.5215708) q[2];
sx q[2];
rz(-2.7572981) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1943097) q[1];
sx q[1];
rz(-0.81856809) q[1];
sx q[1];
rz(-1.1771997) q[1];
rz(-pi) q[2];
rz(-2.5475575) q[3];
sx q[3];
rz(-1.4184991) q[3];
sx q[3];
rz(1.9826979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.74497491) q[2];
sx q[2];
rz(-1.2367542) q[2];
sx q[2];
rz(-1.940894) q[2];
rz(2.4602304) q[3];
sx q[3];
rz(-1.6187637) q[3];
sx q[3];
rz(-2.7549226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.50614) q[0];
sx q[0];
rz(-0.30158392) q[0];
sx q[0];
rz(-2.8801081) q[0];
rz(1.6188949) q[1];
sx q[1];
rz(-0.50917429) q[1];
sx q[1];
rz(-1.8656628) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2842002) q[0];
sx q[0];
rz(-2.5821857) q[0];
sx q[0];
rz(0.18226923) q[0];
rz(-pi) q[1];
rz(-2.475937) q[2];
sx q[2];
rz(-0.7705858) q[2];
sx q[2];
rz(1.6536755) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2535227) q[1];
sx q[1];
rz(-1.8948529) q[1];
sx q[1];
rz(0.56323207) q[1];
x q[2];
rz(1.4218296) q[3];
sx q[3];
rz(-1.5368122) q[3];
sx q[3];
rz(-2.8016718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9246989) q[2];
sx q[2];
rz(-1.0328707) q[2];
sx q[2];
rz(-0.72803289) q[2];
rz(-1.3701471) q[3];
sx q[3];
rz(-2.5244505) q[3];
sx q[3];
rz(0.035577687) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.034123357) q[0];
sx q[0];
rz(-1.1119482) q[0];
sx q[0];
rz(-2.4566101) q[0];
rz(-2.0383539) q[1];
sx q[1];
rz(-0.61956844) q[1];
sx q[1];
rz(2.129295) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023803614) q[0];
sx q[0];
rz(-0.53142457) q[0];
sx q[0];
rz(0.84367911) q[0];
rz(-pi) q[1];
rz(1.3755765) q[2];
sx q[2];
rz(-1.7417656) q[2];
sx q[2];
rz(0.43087105) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6382513) q[1];
sx q[1];
rz(-1.4443399) q[1];
sx q[1];
rz(2.9328037) q[1];
rz(-0.7638974) q[3];
sx q[3];
rz(-1.1925081) q[3];
sx q[3];
rz(-2.0495142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.181902) q[2];
sx q[2];
rz(-1.5896553) q[2];
sx q[2];
rz(2.336179) q[2];
rz(0.81625932) q[3];
sx q[3];
rz(-0.62556848) q[3];
sx q[3];
rz(-0.094154112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1212946) q[0];
sx q[0];
rz(-2.2327406) q[0];
sx q[0];
rz(-2.5138309) q[0];
rz(0.10920564) q[1];
sx q[1];
rz(-0.86995482) q[1];
sx q[1];
rz(2.3039718) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9484884) q[0];
sx q[0];
rz(-0.49608663) q[0];
sx q[0];
rz(-0.14004616) q[0];
rz(-1.8960647) q[2];
sx q[2];
rz(-2.5211689) q[2];
sx q[2];
rz(-1.9577946) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.301103) q[1];
sx q[1];
rz(-0.87574848) q[1];
sx q[1];
rz(1.7453421) q[1];
rz(-pi) q[2];
rz(-2.7660699) q[3];
sx q[3];
rz(-1.2831472) q[3];
sx q[3];
rz(-0.2660397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8236905) q[2];
sx q[2];
rz(-2.6573942) q[2];
sx q[2];
rz(-2.79706) q[2];
rz(-1.7349998) q[3];
sx q[3];
rz(-2.78648) q[3];
sx q[3];
rz(-3.0936892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.542881) q[0];
sx q[0];
rz(-0.6364091) q[0];
sx q[0];
rz(0.5058381) q[0];
rz(2.089031) q[1];
sx q[1];
rz(-0.86741766) q[1];
sx q[1];
rz(-0.94598407) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82362433) q[0];
sx q[0];
rz(-1.8232947) q[0];
sx q[0];
rz(-2.859145) q[0];
x q[1];
rz(2.7840678) q[2];
sx q[2];
rz(-1.7058183) q[2];
sx q[2];
rz(-2.6808878) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2767503) q[1];
sx q[1];
rz(-2.4227243) q[1];
sx q[1];
rz(-0.5867153) q[1];
rz(-1.9993773) q[3];
sx q[3];
rz(-0.16332291) q[3];
sx q[3];
rz(1.6960916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0537009) q[2];
sx q[2];
rz(-1.0975857) q[2];
sx q[2];
rz(-1.5637448) q[2];
rz(1.0846694) q[3];
sx q[3];
rz(-2.2659149) q[3];
sx q[3];
rz(-2.5904371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56187335) q[0];
sx q[0];
rz(-0.60573524) q[0];
sx q[0];
rz(-0.24000034) q[0];
rz(-2.2198246) q[1];
sx q[1];
rz(-2.0926937) q[1];
sx q[1];
rz(1.1704495) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3791948) q[0];
sx q[0];
rz(-2.1403011) q[0];
sx q[0];
rz(-1.5194917) q[0];
rz(1.4438774) q[2];
sx q[2];
rz(-1.9245714) q[2];
sx q[2];
rz(0.60777174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7237548) q[1];
sx q[1];
rz(-1.0772632) q[1];
sx q[1];
rz(0.51512169) q[1];
rz(-1.7975054) q[3];
sx q[3];
rz(-2.7313958) q[3];
sx q[3];
rz(2.6997379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74448284) q[2];
sx q[2];
rz(-0.8195256) q[2];
sx q[2];
rz(2.5605555) q[2];
rz(2.2033612) q[3];
sx q[3];
rz(-0.56648985) q[3];
sx q[3];
rz(-0.70403045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0051603) q[0];
sx q[0];
rz(-2.4649824) q[0];
sx q[0];
rz(0.50049472) q[0];
rz(2.9035134) q[1];
sx q[1];
rz(-1.6903189) q[1];
sx q[1];
rz(-0.64661017) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80386418) q[0];
sx q[0];
rz(-3.1067305) q[0];
sx q[0];
rz(3.0906223) q[0];
x q[1];
rz(-3.1076381) q[2];
sx q[2];
rz(-1.6852549) q[2];
sx q[2];
rz(1.7781228) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.096488086) q[1];
sx q[1];
rz(-0.50042168) q[1];
sx q[1];
rz(1.3184271) q[1];
rz(-pi) q[2];
rz(-0.061731652) q[3];
sx q[3];
rz(-1.928707) q[3];
sx q[3];
rz(-2.6883812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2639192) q[2];
sx q[2];
rz(-2.4223902) q[2];
sx q[2];
rz(-0.70362299) q[2];
rz(1.2879114) q[3];
sx q[3];
rz(-2.0813007) q[3];
sx q[3];
rz(2.6357486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8877761) q[0];
sx q[0];
rz(-2.0075338) q[0];
sx q[0];
rz(-1.8164841) q[0];
rz(-0.77249402) q[1];
sx q[1];
rz(-1.1224116) q[1];
sx q[1];
rz(2.966029) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52161509) q[0];
sx q[0];
rz(-0.78536638) q[0];
sx q[0];
rz(-2.8761441) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5042261) q[2];
sx q[2];
rz(-2.2295579) q[2];
sx q[2];
rz(-0.8085685) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.404322) q[1];
sx q[1];
rz(-2.0312738) q[1];
sx q[1];
rz(0.15562017) q[1];
rz(-pi) q[2];
rz(-1.6962849) q[3];
sx q[3];
rz(-0.85832046) q[3];
sx q[3];
rz(3.0355796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0385711) q[2];
sx q[2];
rz(-1.0744289) q[2];
sx q[2];
rz(1.0938905) q[2];
rz(1.1755747) q[3];
sx q[3];
rz(-2.3837377) q[3];
sx q[3];
rz(-2.8567543) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4232101) q[0];
sx q[0];
rz(-3.0960313) q[0];
sx q[0];
rz(1.3257931) q[0];
rz(-0.52628738) q[1];
sx q[1];
rz(-2.278639) q[1];
sx q[1];
rz(0.78904271) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51368749) q[0];
sx q[0];
rz(-0.18842489) q[0];
sx q[0];
rz(2.1112676) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5297935) q[2];
sx q[2];
rz(-0.42610301) q[2];
sx q[2];
rz(-0.95624051) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10875094) q[1];
sx q[1];
rz(-0.95966731) q[1];
sx q[1];
rz(-2.717589) q[1];
rz(-0.6700683) q[3];
sx q[3];
rz(-0.89691554) q[3];
sx q[3];
rz(-1.4953227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9290756) q[2];
sx q[2];
rz(-0.93364659) q[2];
sx q[2];
rz(2.2361501) q[2];
rz(0.91600156) q[3];
sx q[3];
rz(-1.6986877) q[3];
sx q[3];
rz(2.0974832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7020029) q[0];
sx q[0];
rz(-2.3083394) q[0];
sx q[0];
rz(-2.2138017) q[0];
rz(-2.0480428) q[1];
sx q[1];
rz(-0.55892006) q[1];
sx q[1];
rz(0.13339001) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5800574) q[0];
sx q[0];
rz(-1.8609338) q[0];
sx q[0];
rz(-1.7878637) q[0];
rz(-pi) q[1];
rz(2.9355644) q[2];
sx q[2];
rz(-2.1325754) q[2];
sx q[2];
rz(0.59625193) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9815326) q[1];
sx q[1];
rz(-2.0813749) q[1];
sx q[1];
rz(2.4909944) q[1];
x q[2];
rz(-2.0285151) q[3];
sx q[3];
rz(-1.7594171) q[3];
sx q[3];
rz(2.6831804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5043958) q[2];
sx q[2];
rz(-1.7375676) q[2];
sx q[2];
rz(-1.2633598) q[2];
rz(-2.0718306) q[3];
sx q[3];
rz(-2.1061149) q[3];
sx q[3];
rz(3.0647762) q[3];
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
rz(-pi) q[0];
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
rz(0.99061154) q[0];
sx q[0];
rz(-2.4438416) q[0];
sx q[0];
rz(-1.7846815) q[0];
rz(1.4487343) q[1];
sx q[1];
rz(-2.3285463) q[1];
sx q[1];
rz(2.9978233) q[1];
rz(2.0781381) q[2];
sx q[2];
rz(-2.5567071) q[2];
sx q[2];
rz(-1.1699789) q[2];
rz(1.9142022) q[3];
sx q[3];
rz(-1.3385942) q[3];
sx q[3];
rz(-0.75029324) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
