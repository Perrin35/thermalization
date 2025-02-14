OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7689826) q[0];
sx q[0];
rz(3.1867653) q[0];
sx q[0];
rz(10.09633) q[0];
rz(2.1454732) q[1];
sx q[1];
rz(-0.78875199) q[1];
sx q[1];
rz(-2.8532343) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0147284) q[0];
sx q[0];
rz(-1.5343976) q[0];
sx q[0];
rz(1.3728736) q[0];
rz(2.8003807) q[2];
sx q[2];
rz(-2.4483969) q[2];
sx q[2];
rz(-2.5643519) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8446286) q[1];
sx q[1];
rz(-1.1229672) q[1];
sx q[1];
rz(-1.5289246) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7156473) q[3];
sx q[3];
rz(-1.1205691) q[3];
sx q[3];
rz(-1.0369155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70197376) q[2];
sx q[2];
rz(-0.34586033) q[2];
sx q[2];
rz(-2.4227179) q[2];
rz(-1.6950722) q[3];
sx q[3];
rz(-1.4868163) q[3];
sx q[3];
rz(2.1046624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0816536) q[0];
sx q[0];
rz(-0.43123284) q[0];
sx q[0];
rz(-2.8531139) q[0];
rz(-2.5892995) q[1];
sx q[1];
rz(-1.0918795) q[1];
sx q[1];
rz(-1.1757895) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2644076) q[0];
sx q[0];
rz(-2.7938936) q[0];
sx q[0];
rz(1.1176137) q[0];
rz(3.0251558) q[2];
sx q[2];
rz(-1.7431362) q[2];
sx q[2];
rz(2.4209173) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.9811444) q[1];
sx q[1];
rz(-1.4466043) q[1];
sx q[1];
rz(-0.9915413) q[1];
rz(3.1404488) q[3];
sx q[3];
rz(-1.5387156) q[3];
sx q[3];
rz(2.6828121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4484078) q[2];
sx q[2];
rz(-1.2625445) q[2];
sx q[2];
rz(3.1191678) q[2];
rz(-1.7153995) q[3];
sx q[3];
rz(-1.2778927) q[3];
sx q[3];
rz(0.9001596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65396032) q[0];
sx q[0];
rz(-1.910169) q[0];
sx q[0];
rz(-2.3768429) q[0];
rz(-1.7817616) q[1];
sx q[1];
rz(-1.2218852) q[1];
sx q[1];
rz(-1.8870032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5391961) q[0];
sx q[0];
rz(-2.975292) q[0];
sx q[0];
rz(-2.2434851) q[0];
rz(-pi) q[1];
rz(-0.41608918) q[2];
sx q[2];
rz(-0.35230428) q[2];
sx q[2];
rz(-1.8606297) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4855328) q[1];
sx q[1];
rz(-2.0295709) q[1];
sx q[1];
rz(-1.2478254) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7929501) q[3];
sx q[3];
rz(-0.36262929) q[3];
sx q[3];
rz(-1.8809821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7356871) q[2];
sx q[2];
rz(-0.76852208) q[2];
sx q[2];
rz(-0.020966919) q[2];
rz(1.5603125) q[3];
sx q[3];
rz(-2.0798648) q[3];
sx q[3];
rz(2.235152) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86879325) q[0];
sx q[0];
rz(-0.11196207) q[0];
sx q[0];
rz(1.772076) q[0];
rz(0.70961332) q[1];
sx q[1];
rz(-1.2779002) q[1];
sx q[1];
rz(2.4443464) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0727834) q[0];
sx q[0];
rz(-1.261473) q[0];
sx q[0];
rz(-0.42865045) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3241858) q[2];
sx q[2];
rz(-2.1073034) q[2];
sx q[2];
rz(2.4101933) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.972447) q[1];
sx q[1];
rz(-0.94408542) q[1];
sx q[1];
rz(2.7657126) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3342821) q[3];
sx q[3];
rz(-1.7380889) q[3];
sx q[3];
rz(3.0574556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9559481) q[2];
sx q[2];
rz(-2.1944025) q[2];
sx q[2];
rz(0.67898018) q[2];
rz(1.125157) q[3];
sx q[3];
rz(-1.9966639) q[3];
sx q[3];
rz(-2.9185435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0515902) q[0];
sx q[0];
rz(-1.9380049) q[0];
sx q[0];
rz(-0.077127174) q[0];
rz(1.9902825) q[1];
sx q[1];
rz(-1.5682033) q[1];
sx q[1];
rz(2.3627949) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1098233) q[0];
sx q[0];
rz(-1.5228049) q[0];
sx q[0];
rz(-0.035295156) q[0];
x q[1];
rz(-2.9735317) q[2];
sx q[2];
rz(-1.8521143) q[2];
sx q[2];
rz(2.9945053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59757876) q[1];
sx q[1];
rz(-0.42002267) q[1];
sx q[1];
rz(2.2060664) q[1];
rz(-pi) q[2];
rz(-0.058249931) q[3];
sx q[3];
rz(-0.16213972) q[3];
sx q[3];
rz(1.3221962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6246346) q[2];
sx q[2];
rz(-1.7004852) q[2];
sx q[2];
rz(-1.9514294) q[2];
rz(-2.5879351) q[3];
sx q[3];
rz(-2.5167969) q[3];
sx q[3];
rz(0.73738086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6051642) q[0];
sx q[0];
rz(-1.9909415) q[0];
sx q[0];
rz(2.8705257) q[0];
rz(-2.1412663) q[1];
sx q[1];
rz(-1.2157636) q[1];
sx q[1];
rz(-2.2727374) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9199182) q[0];
sx q[0];
rz(-0.98548792) q[0];
sx q[0];
rz(-2.1409537) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.087681282) q[2];
sx q[2];
rz(-0.73747915) q[2];
sx q[2];
rz(-2.7840441) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4255526) q[1];
sx q[1];
rz(-1.2231266) q[1];
sx q[1];
rz(0.1802318) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0398851) q[3];
sx q[3];
rz(-0.84381754) q[3];
sx q[3];
rz(-0.25694822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.102999) q[2];
sx q[2];
rz(-2.8391892) q[2];
sx q[2];
rz(-2.0043066) q[2];
rz(0.31050995) q[3];
sx q[3];
rz(-2.2313084) q[3];
sx q[3];
rz(2.212132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0248658) q[0];
sx q[0];
rz(-1.4288582) q[0];
sx q[0];
rz(2.0615935) q[0];
rz(2.179821) q[1];
sx q[1];
rz(-0.56517833) q[1];
sx q[1];
rz(2.9186509) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5436542) q[0];
sx q[0];
rz(-3.1059382) q[0];
sx q[0];
rz(-1.7810506) q[0];
rz(-0.30855631) q[2];
sx q[2];
rz(-0.94538222) q[2];
sx q[2];
rz(0.61185123) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2126951) q[1];
sx q[1];
rz(-1.846781) q[1];
sx q[1];
rz(2.3174965) q[1];
rz(1.5925528) q[3];
sx q[3];
rz(-0.59594107) q[3];
sx q[3];
rz(0.21778743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96233931) q[2];
sx q[2];
rz(-2.8620359) q[2];
sx q[2];
rz(-1.1780098) q[2];
rz(-2.8152605) q[3];
sx q[3];
rz(-2.3309565) q[3];
sx q[3];
rz(1.1403722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2278263) q[0];
sx q[0];
rz(-2.1837809) q[0];
sx q[0];
rz(-2.3334099) q[0];
rz(2.783964) q[1];
sx q[1];
rz(-1.459815) q[1];
sx q[1];
rz(-2.0590032) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6480761) q[0];
sx q[0];
rz(-2.2449623) q[0];
sx q[0];
rz(-2.9229259) q[0];
x q[1];
rz(1.9679734) q[2];
sx q[2];
rz(-0.589314) q[2];
sx q[2];
rz(2.9795424) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.030176) q[1];
sx q[1];
rz(-1.7716265) q[1];
sx q[1];
rz(-0.80441956) q[1];
x q[2];
rz(-2.6851366) q[3];
sx q[3];
rz(-2.0051867) q[3];
sx q[3];
rz(2.1214157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0766582) q[2];
sx q[2];
rz(-2.6532463) q[2];
sx q[2];
rz(1.2154382) q[2];
rz(0.91228929) q[3];
sx q[3];
rz(-2.2513697) q[3];
sx q[3];
rz(0.54273763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8488345) q[0];
sx q[0];
rz(-0.64105761) q[0];
sx q[0];
rz(-0.42375281) q[0];
rz(1.1774225) q[1];
sx q[1];
rz(-2.2260428) q[1];
sx q[1];
rz(0.37568572) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3553414) q[0];
sx q[0];
rz(-1.6415941) q[0];
sx q[0];
rz(1.7485445) q[0];
x q[1];
rz(1.5974974) q[2];
sx q[2];
rz(-1.809568) q[2];
sx q[2];
rz(1.8516061) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2089952) q[1];
sx q[1];
rz(-1.0699341) q[1];
sx q[1];
rz(0.074525699) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2311835) q[3];
sx q[3];
rz(-2.2582128) q[3];
sx q[3];
rz(-0.31654762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.11289135) q[2];
sx q[2];
rz(-1.4914923) q[2];
sx q[2];
rz(2.4493307) q[2];
rz(-0.68474692) q[3];
sx q[3];
rz(-0.94963494) q[3];
sx q[3];
rz(0.14028604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.5570062) q[0];
sx q[0];
rz(-2.1742915) q[0];
sx q[0];
rz(-1.6242356) q[0];
rz(-1.5380305) q[1];
sx q[1];
rz(-1.3471194) q[1];
sx q[1];
rz(-1.2204407) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.206911) q[0];
sx q[0];
rz(-2.5162272) q[0];
sx q[0];
rz(-1.4689133) q[0];
rz(-2.7292029) q[2];
sx q[2];
rz(-2.6021716) q[2];
sx q[2];
rz(-2.256765) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0076589) q[1];
sx q[1];
rz(-0.87532212) q[1];
sx q[1];
rz(-1.9401866) q[1];
rz(-pi) q[2];
rz(-1.8655928) q[3];
sx q[3];
rz(-0.85569564) q[3];
sx q[3];
rz(-3.0218647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5753691) q[2];
sx q[2];
rz(-2.1341925) q[2];
sx q[2];
rz(-2.2929906) q[2];
rz(0.97950116) q[3];
sx q[3];
rz(-1.5749911) q[3];
sx q[3];
rz(-0.037467329) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.033087) q[0];
sx q[0];
rz(-1.6642234) q[0];
sx q[0];
rz(-2.4432175) q[0];
rz(-2.5820844) q[1];
sx q[1];
rz(-2.5026176) q[1];
sx q[1];
rz(2.7284596) q[1];
rz(-2.8200061) q[2];
sx q[2];
rz(-1.3272616) q[2];
sx q[2];
rz(2.7250901) q[2];
rz(-1.63089) q[3];
sx q[3];
rz(-0.99219764) q[3];
sx q[3];
rz(0.037430684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
