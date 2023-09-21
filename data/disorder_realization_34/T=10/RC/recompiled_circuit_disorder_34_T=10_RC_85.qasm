OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.37378398) q[0];
sx q[0];
rz(-2.7019579) q[0];
sx q[0];
rz(0.08131942) q[0];
rz(0.65027872) q[1];
sx q[1];
rz(-1.283409) q[1];
sx q[1];
rz(0.78279701) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5188132) q[0];
sx q[0];
rz(-2.8840027) q[0];
sx q[0];
rz(-1.8151087) q[0];
x q[1];
rz(-1.9185669) q[2];
sx q[2];
rz(-0.43453056) q[2];
sx q[2];
rz(2.6968616) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4408735) q[1];
sx q[1];
rz(-2.6881725) q[1];
sx q[1];
rz(-1.2749626) q[1];
x q[2];
rz(-1.6601895) q[3];
sx q[3];
rz(-1.2561744) q[3];
sx q[3];
rz(-2.6320626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.89447442) q[2];
sx q[2];
rz(-1.0049745) q[2];
sx q[2];
rz(-0.11581126) q[2];
rz(1.5420906) q[3];
sx q[3];
rz(-3.0452947) q[3];
sx q[3];
rz(-2.0882864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88749921) q[0];
sx q[0];
rz(-0.54953456) q[0];
sx q[0];
rz(-0.19533531) q[0];
rz(2.7665566) q[1];
sx q[1];
rz(-1.6655567) q[1];
sx q[1];
rz(0.23981747) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6759684) q[0];
sx q[0];
rz(-1.8166607) q[0];
sx q[0];
rz(1.4093536) q[0];
x q[1];
rz(0.96599483) q[2];
sx q[2];
rz(-0.36075337) q[2];
sx q[2];
rz(0.15168562) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5259243) q[1];
sx q[1];
rz(-2.3500372) q[1];
sx q[1];
rz(-0.073992373) q[1];
x q[2];
rz(0.014157045) q[3];
sx q[3];
rz(-1.5870968) q[3];
sx q[3];
rz(-2.1099159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5043162) q[2];
sx q[2];
rz(-2.3392623) q[2];
sx q[2];
rz(1.3298539) q[2];
rz(1.7999533) q[3];
sx q[3];
rz(-1.642671) q[3];
sx q[3];
rz(-1.8168861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.864569) q[0];
sx q[0];
rz(-1.8947911) q[0];
sx q[0];
rz(-0.66816107) q[0];
rz(1.4913303) q[1];
sx q[1];
rz(-0.69258339) q[1];
sx q[1];
rz(2.0756762) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0695755) q[0];
sx q[0];
rz(-1.6117275) q[0];
sx q[0];
rz(1.3224365) q[0];
x q[1];
rz(-0.22914825) q[2];
sx q[2];
rz(-0.56729588) q[2];
sx q[2];
rz(1.0457525) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4855811) q[1];
sx q[1];
rz(-2.1335568) q[1];
sx q[1];
rz(2.807711) q[1];
rz(-pi) q[2];
rz(1.7245674) q[3];
sx q[3];
rz(-0.66699281) q[3];
sx q[3];
rz(2.1949777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.187591) q[2];
sx q[2];
rz(-1.9090586) q[2];
sx q[2];
rz(-3.1090453) q[2];
rz(0.35999808) q[3];
sx q[3];
rz(-2.0149752) q[3];
sx q[3];
rz(-0.79157296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2739094) q[0];
sx q[0];
rz(-1.5554579) q[0];
sx q[0];
rz(-2.4348863) q[0];
rz(1.2061521) q[1];
sx q[1];
rz(-2.8001092) q[1];
sx q[1];
rz(1.4867841) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33071163) q[0];
sx q[0];
rz(-1.8510305) q[0];
sx q[0];
rz(2.7964554) q[0];
x q[1];
rz(-1.9183667) q[2];
sx q[2];
rz(-0.98993694) q[2];
sx q[2];
rz(0.52085224) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.023991) q[1];
sx q[1];
rz(-1.4211434) q[1];
sx q[1];
rz(2.8531122) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5634469) q[3];
sx q[3];
rz(-2.3971933) q[3];
sx q[3];
rz(0.053754427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5965745) q[2];
sx q[2];
rz(-0.89349616) q[2];
sx q[2];
rz(-1.0085227) q[2];
rz(1.0962983) q[3];
sx q[3];
rz(-1.228628) q[3];
sx q[3];
rz(-3.0116459) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1059882) q[0];
sx q[0];
rz(-2.2898219) q[0];
sx q[0];
rz(1.6812356) q[0];
rz(1.5530855) q[1];
sx q[1];
rz(-1.9057143) q[1];
sx q[1];
rz(3.1255186) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8668629) q[0];
sx q[0];
rz(-2.1511973) q[0];
sx q[0];
rz(-0.6443364) q[0];
x q[1];
rz(2.0015012) q[2];
sx q[2];
rz(-1.5783731) q[2];
sx q[2];
rz(1.8097144) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6962471) q[1];
sx q[1];
rz(-1.0187341) q[1];
sx q[1];
rz(0.50260431) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5752951) q[3];
sx q[3];
rz(-1.2084949) q[3];
sx q[3];
rz(-3.0294861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1092704) q[2];
sx q[2];
rz(-2.1123999) q[2];
sx q[2];
rz(-0.62189046) q[2];
rz(2.0444929) q[3];
sx q[3];
rz(-2.3648839) q[3];
sx q[3];
rz(2.9212852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9428403) q[0];
sx q[0];
rz(-3.1382914) q[0];
sx q[0];
rz(0.90676701) q[0];
rz(2.3268907) q[1];
sx q[1];
rz(-0.68836132) q[1];
sx q[1];
rz(-1.2247359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82299267) q[0];
sx q[0];
rz(-2.4791414) q[0];
sx q[0];
rz(1.9255161) q[0];
rz(-pi) q[1];
rz(1.1282863) q[2];
sx q[2];
rz(-1.3820717) q[2];
sx q[2];
rz(-2.92885) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8606557) q[1];
sx q[1];
rz(-2.5334364) q[1];
sx q[1];
rz(1.6673191) q[1];
x q[2];
rz(-0.25572689) q[3];
sx q[3];
rz(-2.7284107) q[3];
sx q[3];
rz(0.43858389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1427052) q[2];
sx q[2];
rz(-2.3110516) q[2];
sx q[2];
rz(2.2018946) q[2];
rz(2.9283004) q[3];
sx q[3];
rz(-0.34049884) q[3];
sx q[3];
rz(-1.7512158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1619103) q[0];
sx q[0];
rz(-0.96452159) q[0];
sx q[0];
rz(-2.5612223) q[0];
rz(-2.0866108) q[1];
sx q[1];
rz(-1.4529198) q[1];
sx q[1];
rz(2.4408128) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.128292) q[0];
sx q[0];
rz(-1.6830781) q[0];
sx q[0];
rz(-2.0237472) q[0];
rz(-0.3406616) q[2];
sx q[2];
rz(-0.33764687) q[2];
sx q[2];
rz(-1.9549191) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.339387) q[1];
sx q[1];
rz(-1.4805668) q[1];
sx q[1];
rz(1.0046563) q[1];
rz(0.084214597) q[3];
sx q[3];
rz(-0.84184064) q[3];
sx q[3];
rz(-2.6654411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7769527) q[2];
sx q[2];
rz(-2.8312455) q[2];
sx q[2];
rz(-0.022162612) q[2];
rz(-2.3968905) q[3];
sx q[3];
rz(-1.1170758) q[3];
sx q[3];
rz(-2.7409592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.3547524) q[0];
sx q[0];
rz(-2.1180034) q[0];
sx q[0];
rz(-1.7250852) q[0];
rz(-1.3757061) q[1];
sx q[1];
rz(-1.4027275) q[1];
sx q[1];
rz(0.8917121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.007908) q[0];
sx q[0];
rz(-0.79151756) q[0];
sx q[0];
rz(-1.3484811) q[0];
rz(-pi) q[1];
rz(1.1197929) q[2];
sx q[2];
rz(-1.0378671) q[2];
sx q[2];
rz(1.9110796) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8430431) q[1];
sx q[1];
rz(-1.1617359) q[1];
sx q[1];
rz(0.72960735) q[1];
rz(-1.1710839) q[3];
sx q[3];
rz(-2.4418695) q[3];
sx q[3];
rz(1.8613929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4884168) q[2];
sx q[2];
rz(-1.9501053) q[2];
sx q[2];
rz(2.3573504) q[2];
rz(2.6358321) q[3];
sx q[3];
rz(-2.2890746) q[3];
sx q[3];
rz(0.23323664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6655675) q[0];
sx q[0];
rz(-1.5371756) q[0];
sx q[0];
rz(2.419557) q[0];
rz(-2.8083535) q[1];
sx q[1];
rz(-1.1958586) q[1];
sx q[1];
rz(-1.3649712) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6745973) q[0];
sx q[0];
rz(-0.82775263) q[0];
sx q[0];
rz(1.6270431) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1599837) q[2];
sx q[2];
rz(-1.4940726) q[2];
sx q[2];
rz(-3.0551747) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8553798) q[1];
sx q[1];
rz(-2.3533822) q[1];
sx q[1];
rz(2.6469995) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2737208) q[3];
sx q[3];
rz(-1.3798514) q[3];
sx q[3];
rz(1.2862658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1086796) q[2];
sx q[2];
rz(-1.7620757) q[2];
sx q[2];
rz(1.2314679) q[2];
rz(-3.1075297) q[3];
sx q[3];
rz(-1.8647727) q[3];
sx q[3];
rz(-0.6347707) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0868527) q[0];
sx q[0];
rz(-2.5755136) q[0];
sx q[0];
rz(-1.6636794) q[0];
rz(-2.058303) q[1];
sx q[1];
rz(-1.3996841) q[1];
sx q[1];
rz(2.1733984) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9111239) q[0];
sx q[0];
rz(-1.231912) q[0];
sx q[0];
rz(1.1662657) q[0];
rz(-pi) q[1];
rz(2.7231611) q[2];
sx q[2];
rz(-1.1088587) q[2];
sx q[2];
rz(-0.07428169) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.248368) q[1];
sx q[1];
rz(-1.4239422) q[1];
sx q[1];
rz(2.8180772) q[1];
rz(-pi) q[2];
rz(-1.2486542) q[3];
sx q[3];
rz(-1.570574) q[3];
sx q[3];
rz(-2.145888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0782464) q[2];
sx q[2];
rz(-0.72786704) q[2];
sx q[2];
rz(0.001312288) q[2];
rz(2.0007658) q[3];
sx q[3];
rz(-1.8246548) q[3];
sx q[3];
rz(1.7989981) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.39682) q[0];
sx q[0];
rz(-2.0347432) q[0];
sx q[0];
rz(1.9532935) q[0];
rz(-2.7753579) q[1];
sx q[1];
rz(-1.9402505) q[1];
sx q[1];
rz(-1.8016626) q[1];
rz(2.3920849) q[2];
sx q[2];
rz(-1.3386249) q[2];
sx q[2];
rz(-2.7809536) q[2];
rz(-2.1063741) q[3];
sx q[3];
rz(-2.513701) q[3];
sx q[3];
rz(2.7660478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];