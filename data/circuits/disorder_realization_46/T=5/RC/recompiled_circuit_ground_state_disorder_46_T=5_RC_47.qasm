OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.159654) q[0];
sx q[0];
rz(-1.5647793) q[0];
sx q[0];
rz(1.2793581) q[0];
rz(0.98969412) q[1];
sx q[1];
rz(-1.1969748) q[1];
sx q[1];
rz(2.9762414) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5419614) q[0];
sx q[0];
rz(-1.8867869) q[0];
sx q[0];
rz(-0.91008474) q[0];
rz(-pi) q[1];
rz(2.3897091) q[2];
sx q[2];
rz(-0.57752973) q[2];
sx q[2];
rz(2.5381388) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0667116) q[1];
sx q[1];
rz(-2.0600658) q[1];
sx q[1];
rz(0.46164767) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8836423) q[3];
sx q[3];
rz(-1.2688302) q[3];
sx q[3];
rz(-1.0416958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2692261) q[2];
sx q[2];
rz(-1.1507611) q[2];
sx q[2];
rz(1.9744138) q[2];
rz(-0.40927467) q[3];
sx q[3];
rz(-0.27547488) q[3];
sx q[3];
rz(-0.95612139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(3.069365) q[0];
sx q[0];
rz(-0.15412155) q[0];
sx q[0];
rz(-1.3519721) q[0];
rz(1.6701472) q[1];
sx q[1];
rz(-0.92266005) q[1];
sx q[1];
rz(-0.89554375) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35947511) q[0];
sx q[0];
rz(-1.5687587) q[0];
sx q[0];
rz(-0.26709132) q[0];
x q[1];
rz(-1.6487213) q[2];
sx q[2];
rz(-1.0977626) q[2];
sx q[2];
rz(-0.89240197) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1481759) q[1];
sx q[1];
rz(-1.9636781) q[1];
sx q[1];
rz(2.3200289) q[1];
x q[2];
rz(2.9053365) q[3];
sx q[3];
rz(-1.1614262) q[3];
sx q[3];
rz(-0.3796814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.46513778) q[2];
sx q[2];
rz(-2.319591) q[2];
sx q[2];
rz(1.6949534) q[2];
rz(-1.0195352) q[3];
sx q[3];
rz(-1.3186224) q[3];
sx q[3];
rz(-0.29729602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.4296221) q[0];
sx q[0];
rz(-1.1638887) q[0];
sx q[0];
rz(-0.0042560552) q[0];
rz(-0.70107067) q[1];
sx q[1];
rz(-2.6316167) q[1];
sx q[1];
rz(3.1140936) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4002747) q[0];
sx q[0];
rz(-1.7257462) q[0];
sx q[0];
rz(1.1688031) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84553366) q[2];
sx q[2];
rz(-1.4730244) q[2];
sx q[2];
rz(-2.1883983) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.356952) q[1];
sx q[1];
rz(-2.7773547) q[1];
sx q[1];
rz(1.0491153) q[1];
x q[2];
rz(3.1276567) q[3];
sx q[3];
rz(-1.2346141) q[3];
sx q[3];
rz(2.5232836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.49715257) q[2];
sx q[2];
rz(-1.3500328) q[2];
sx q[2];
rz(-0.39786097) q[2];
rz(-1.0128939) q[3];
sx q[3];
rz(-1.0254859) q[3];
sx q[3];
rz(2.7329172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6560646) q[0];
sx q[0];
rz(-1.8714454) q[0];
sx q[0];
rz(-2.2960466) q[0];
rz(0.0088508765) q[1];
sx q[1];
rz(-2.2933941) q[1];
sx q[1];
rz(-1.9890076) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76149135) q[0];
sx q[0];
rz(-1.9753755) q[0];
sx q[0];
rz(-0.92104162) q[0];
rz(-pi) q[1];
rz(-2.8834226) q[2];
sx q[2];
rz(-2.008778) q[2];
sx q[2];
rz(0.5817619) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2680451) q[1];
sx q[1];
rz(-1.5427886) q[1];
sx q[1];
rz(-2.5491879) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5043855) q[3];
sx q[3];
rz(-1.2577269) q[3];
sx q[3];
rz(0.84740365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74959603) q[2];
sx q[2];
rz(-2.1139202) q[2];
sx q[2];
rz(1.2006302) q[2];
rz(2.6642753) q[3];
sx q[3];
rz(-2.6715607) q[3];
sx q[3];
rz(0.70394713) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93477997) q[0];
sx q[0];
rz(-0.88345695) q[0];
sx q[0];
rz(-0.90428895) q[0];
rz(-0.9710871) q[1];
sx q[1];
rz(-1.1957518) q[1];
sx q[1];
rz(2.8579874) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6672517) q[0];
sx q[0];
rz(-1.5548238) q[0];
sx q[0];
rz(-3.0706997) q[0];
x q[1];
rz(0.798337) q[2];
sx q[2];
rz(-2.5588648) q[2];
sx q[2];
rz(-2.5495364) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7452104) q[1];
sx q[1];
rz(-1.3644427) q[1];
sx q[1];
rz(2.2554805) q[1];
x q[2];
rz(1.0002076) q[3];
sx q[3];
rz(-0.772744) q[3];
sx q[3];
rz(0.34602133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7859555) q[2];
sx q[2];
rz(-1.2355618) q[2];
sx q[2];
rz(0.23814417) q[2];
rz(2.1613878) q[3];
sx q[3];
rz(-1.9832059) q[3];
sx q[3];
rz(-2.3908206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7439483) q[0];
sx q[0];
rz(-1.605796) q[0];
sx q[0];
rz(2.7342947) q[0];
rz(2.7382355) q[1];
sx q[1];
rz(-2.401001) q[1];
sx q[1];
rz(-0.86722803) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7776398) q[0];
sx q[0];
rz(-2.1356076) q[0];
sx q[0];
rz(-0.070567957) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8132956) q[2];
sx q[2];
rz(-0.5813501) q[2];
sx q[2];
rz(0.7996847) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2867409) q[1];
sx q[1];
rz(-0.84567243) q[1];
sx q[1];
rz(-1.185536) q[1];
rz(-pi) q[2];
rz(-1.1759472) q[3];
sx q[3];
rz(-2.8215652) q[3];
sx q[3];
rz(-1.4744028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0200218) q[2];
sx q[2];
rz(-0.82841221) q[2];
sx q[2];
rz(1.7106445) q[2];
rz(-2.6089148) q[3];
sx q[3];
rz(-1.0360274) q[3];
sx q[3];
rz(-1.7238341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1581887) q[0];
sx q[0];
rz(-0.15022763) q[0];
sx q[0];
rz(-2.0576553) q[0];
rz(0.051941959) q[1];
sx q[1];
rz(-0.40819326) q[1];
sx q[1];
rz(3.1165677) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1442102) q[0];
sx q[0];
rz(-2.6591427) q[0];
sx q[0];
rz(0.061003322) q[0];
rz(1.0919326) q[2];
sx q[2];
rz(-0.81170481) q[2];
sx q[2];
rz(-2.1213437) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6083303) q[1];
sx q[1];
rz(-0.46583781) q[1];
sx q[1];
rz(-1.9598258) q[1];
rz(-pi) q[2];
rz(-2.8353416) q[3];
sx q[3];
rz(-2.6970377) q[3];
sx q[3];
rz(-1.8646835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0673125) q[2];
sx q[2];
rz(-1.4558027) q[2];
sx q[2];
rz(-2.32453) q[2];
rz(2.1909511) q[3];
sx q[3];
rz(-1.0167511) q[3];
sx q[3];
rz(-1.954621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0364712) q[0];
sx q[0];
rz(-0.28577411) q[0];
sx q[0];
rz(-0.60086077) q[0];
rz(-3.1071682) q[1];
sx q[1];
rz(-1.333678) q[1];
sx q[1];
rz(-1.6011802) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87266028) q[0];
sx q[0];
rz(-1.0300228) q[0];
sx q[0];
rz(-2.2898104) q[0];
rz(2.41018) q[2];
sx q[2];
rz(-2.3225975) q[2];
sx q[2];
rz(-0.42289823) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.29112962) q[1];
sx q[1];
rz(-0.4402059) q[1];
sx q[1];
rz(3.090211) q[1];
rz(-pi) q[2];
rz(-0.65798379) q[3];
sx q[3];
rz(-0.65449981) q[3];
sx q[3];
rz(2.2006048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8704661) q[2];
sx q[2];
rz(-1.2721964) q[2];
sx q[2];
rz(-1.6477443) q[2];
rz(-0.77914023) q[3];
sx q[3];
rz(-1.4804163) q[3];
sx q[3];
rz(-2.486865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7036024) q[0];
sx q[0];
rz(-1.0810738) q[0];
sx q[0];
rz(-2.3533452) q[0];
rz(1.2429271) q[1];
sx q[1];
rz(-1.5595167) q[1];
sx q[1];
rz(2.8020614) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96974715) q[0];
sx q[0];
rz(-1.3783499) q[0];
sx q[0];
rz(0.11809531) q[0];
rz(-0.56216424) q[2];
sx q[2];
rz(-1.8149281) q[2];
sx q[2];
rz(0.31329271) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1052983) q[1];
sx q[1];
rz(-0.52043167) q[1];
sx q[1];
rz(-2.4776513) q[1];
x q[2];
rz(-0.92206436) q[3];
sx q[3];
rz(-0.90708238) q[3];
sx q[3];
rz(1.0459378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.01123151) q[2];
sx q[2];
rz(-1.1684544) q[2];
sx q[2];
rz(-1.8846903) q[2];
rz(-2.0508749) q[3];
sx q[3];
rz(-1.2033477) q[3];
sx q[3];
rz(0.91712657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(-2.286769) q[0];
sx q[0];
rz(-1.6463065) q[0];
sx q[0];
rz(-2.0538034) q[0];
rz(2.4913359) q[1];
sx q[1];
rz(-1.6842664) q[1];
sx q[1];
rz(3.1204209) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8842554) q[0];
sx q[0];
rz(-1.7790964) q[0];
sx q[0];
rz(-0.45680372) q[0];
x q[1];
rz(0.4707321) q[2];
sx q[2];
rz(-1.3952878) q[2];
sx q[2];
rz(2.144475) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11243992) q[1];
sx q[1];
rz(-2.314056) q[1];
sx q[1];
rz(-0.2874916) q[1];
x q[2];
rz(-0.043280799) q[3];
sx q[3];
rz(-1.0387522) q[3];
sx q[3];
rz(-2.8325641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0040969) q[2];
sx q[2];
rz(-1.7986412) q[2];
sx q[2];
rz(-2.1583978) q[2];
rz(-2.148597) q[3];
sx q[3];
rz(-2.0296622) q[3];
sx q[3];
rz(0.23428169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8773593) q[0];
sx q[0];
rz(-0.79240427) q[0];
sx q[0];
rz(-2.0808676) q[0];
rz(-2.0296774) q[1];
sx q[1];
rz(-0.30525515) q[1];
sx q[1];
rz(0.12259604) q[1];
rz(1.0884825) q[2];
sx q[2];
rz(-1.7822722) q[2];
sx q[2];
rz(1.7581802) q[2];
rz(1.5817011) q[3];
sx q[3];
rz(-2.0937216) q[3];
sx q[3];
rz(1.0918319) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
