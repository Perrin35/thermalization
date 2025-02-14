OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0903836) q[0];
sx q[0];
rz(-2.3530355) q[0];
sx q[0];
rz(0.30858827) q[0];
rz(0.42449549) q[1];
sx q[1];
rz(-2.0574044) q[1];
sx q[1];
rz(-0.37556136) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6294) q[0];
sx q[0];
rz(-1.4951271) q[0];
sx q[0];
rz(1.8528432) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0388366) q[2];
sx q[2];
rz(-0.30347541) q[2];
sx q[2];
rz(3.0927401) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8164067) q[1];
sx q[1];
rz(-2.6734787) q[1];
sx q[1];
rz(0.31030841) q[1];
rz(-pi) q[2];
rz(-1.7143024) q[3];
sx q[3];
rz(-1.5670027) q[3];
sx q[3];
rz(-0.11303813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3382807) q[2];
sx q[2];
rz(-0.25140005) q[2];
sx q[2];
rz(-0.54331642) q[2];
rz(-1.4079037) q[3];
sx q[3];
rz(-1.4145773) q[3];
sx q[3];
rz(0.89731115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5498085) q[0];
sx q[0];
rz(-1.9928638) q[0];
sx q[0];
rz(0.59709221) q[0];
rz(1.0626556) q[1];
sx q[1];
rz(-2.5025044) q[1];
sx q[1];
rz(-2.7426681) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0380533) q[0];
sx q[0];
rz(-1.000043) q[0];
sx q[0];
rz(0.45989238) q[0];
x q[1];
rz(1.609732) q[2];
sx q[2];
rz(-0.557347) q[2];
sx q[2];
rz(-0.61902905) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0831757) q[1];
sx q[1];
rz(-2.1444369) q[1];
sx q[1];
rz(-0.87223335) q[1];
rz(-pi) q[2];
rz(-0.53360231) q[3];
sx q[3];
rz(-2.1416365) q[3];
sx q[3];
rz(0.10122964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.67542934) q[2];
sx q[2];
rz(-1.2389946) q[2];
sx q[2];
rz(-1.1413525) q[2];
rz(-0.91836786) q[3];
sx q[3];
rz(-2.4310302) q[3];
sx q[3];
rz(-2.5520798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68701768) q[0];
sx q[0];
rz(-0.37746754) q[0];
sx q[0];
rz(-2.9662509) q[0];
rz(1.31458) q[1];
sx q[1];
rz(-1.250123) q[1];
sx q[1];
rz(-0.66741991) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6771362) q[0];
sx q[0];
rz(-0.51119294) q[0];
sx q[0];
rz(1.808046) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0834789) q[2];
sx q[2];
rz(-2.1082775) q[2];
sx q[2];
rz(1.7318685) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2167425) q[1];
sx q[1];
rz(-1.5702254) q[1];
sx q[1];
rz(-2.8257829) q[1];
rz(-1.3046566) q[3];
sx q[3];
rz(-2.168306) q[3];
sx q[3];
rz(1.430703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0257001) q[2];
sx q[2];
rz(-2.1860055) q[2];
sx q[2];
rz(-1.0828177) q[2];
rz(1.3057905) q[3];
sx q[3];
rz(-0.81597733) q[3];
sx q[3];
rz(-0.79735565) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0037271518) q[0];
sx q[0];
rz(-0.19345134) q[0];
sx q[0];
rz(2.5306012) q[0];
rz(-0.56426471) q[1];
sx q[1];
rz(-2.1128426) q[1];
sx q[1];
rz(0.70704031) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90908289) q[0];
sx q[0];
rz(-2.3444026) q[0];
sx q[0];
rz(0.042122201) q[0];
rz(-pi) q[1];
rz(0.66989278) q[2];
sx q[2];
rz(-2.6351053) q[2];
sx q[2];
rz(-2.2746925) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8312861) q[1];
sx q[1];
rz(-1.093068) q[1];
sx q[1];
rz(-2.0654486) q[1];
rz(-0.50102542) q[3];
sx q[3];
rz(-1.820537) q[3];
sx q[3];
rz(1.6660894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.45346144) q[2];
sx q[2];
rz(-0.42669272) q[2];
sx q[2];
rz(2.1235662) q[2];
rz(-2.7536143) q[3];
sx q[3];
rz(-1.4464902) q[3];
sx q[3];
rz(0.33897266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6729048) q[0];
sx q[0];
rz(-2.0141116) q[0];
sx q[0];
rz(-0.52184033) q[0];
rz(-2.8537967) q[1];
sx q[1];
rz(-2.5591873) q[1];
sx q[1];
rz(0.60307455) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1408242) q[0];
sx q[0];
rz(-0.10061564) q[0];
sx q[0];
rz(-0.510143) q[0];
rz(-pi) q[1];
rz(-0.65400161) q[2];
sx q[2];
rz(-0.79833657) q[2];
sx q[2];
rz(2.4545074) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.52023639) q[1];
sx q[1];
rz(-1.6545611) q[1];
sx q[1];
rz(-2.7158974) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8766635) q[3];
sx q[3];
rz(-1.2390803) q[3];
sx q[3];
rz(-1.1379366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7120984) q[2];
sx q[2];
rz(-0.79271972) q[2];
sx q[2];
rz(0.73954868) q[2];
rz(-2.126501) q[3];
sx q[3];
rz(-2.6138217) q[3];
sx q[3];
rz(-0.019006193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5571112) q[0];
sx q[0];
rz(-0.33753532) q[0];
sx q[0];
rz(0.8648411) q[0];
rz(2.5014014) q[1];
sx q[1];
rz(-0.65018153) q[1];
sx q[1];
rz(0.4943628) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1304754) q[0];
sx q[0];
rz(-1.1251161) q[0];
sx q[0];
rz(-2.9167239) q[0];
rz(-pi) q[1];
rz(0.041448822) q[2];
sx q[2];
rz(-2.1728656) q[2];
sx q[2];
rz(-2.1727501) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.84029217) q[1];
sx q[1];
rz(-1.6052142) q[1];
sx q[1];
rz(0.057386668) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5036656) q[3];
sx q[3];
rz(-2.0450838) q[3];
sx q[3];
rz(2.4573457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2690027) q[2];
sx q[2];
rz(-2.5073017) q[2];
sx q[2];
rz(-0.96160257) q[2];
rz(-1.6327935) q[3];
sx q[3];
rz(-1.5140994) q[3];
sx q[3];
rz(0.076920286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85840571) q[0];
sx q[0];
rz(-0.34342331) q[0];
sx q[0];
rz(2.0984233) q[0];
rz(-0.6768325) q[1];
sx q[1];
rz(-0.64462858) q[1];
sx q[1];
rz(2.4136995) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4235509) q[0];
sx q[0];
rz(-1.4598993) q[0];
sx q[0];
rz(-1.9865722) q[0];
x q[1];
rz(-0.0038996242) q[2];
sx q[2];
rz(-0.74206381) q[2];
sx q[2];
rz(0.33483792) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.7866061) q[1];
sx q[1];
rz(-2.1713272) q[1];
sx q[1];
rz(0.40734603) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9123825) q[3];
sx q[3];
rz(-2.1331969) q[3];
sx q[3];
rz(-1.7096568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4595043) q[2];
sx q[2];
rz(-2.3900718) q[2];
sx q[2];
rz(-0.53306836) q[2];
rz(1.3047949) q[3];
sx q[3];
rz(-2.6663836) q[3];
sx q[3];
rz(-2.3615725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.099139) q[0];
sx q[0];
rz(-3.070153) q[0];
sx q[0];
rz(-0.024600994) q[0];
rz(3.091231) q[1];
sx q[1];
rz(-0.55614007) q[1];
sx q[1];
rz(0.76739001) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3022223) q[0];
sx q[0];
rz(-1.5451533) q[0];
sx q[0];
rz(-0.54021949) q[0];
rz(-1.3774302) q[2];
sx q[2];
rz(-1.2467209) q[2];
sx q[2];
rz(1.1149158) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.394262) q[1];
sx q[1];
rz(-0.74968265) q[1];
sx q[1];
rz(-0.59632991) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1786091) q[3];
sx q[3];
rz(-1.0309891) q[3];
sx q[3];
rz(-1.3572049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.47160992) q[2];
sx q[2];
rz(-2.6675384) q[2];
sx q[2];
rz(1.6101884) q[2];
rz(-2.1286185) q[3];
sx q[3];
rz(-1.4644724) q[3];
sx q[3];
rz(0.60232919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7811964) q[0];
sx q[0];
rz(-0.66806) q[0];
sx q[0];
rz(2.6511731) q[0];
rz(2.7410653) q[1];
sx q[1];
rz(-2.782395) q[1];
sx q[1];
rz(1.5708956) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84874451) q[0];
sx q[0];
rz(-1.2146753) q[0];
sx q[0];
rz(0.57909052) q[0];
x q[1];
rz(1.2068858) q[2];
sx q[2];
rz(-1.8779953) q[2];
sx q[2];
rz(-2.4138434) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8388673) q[1];
sx q[1];
rz(-1.6714393) q[1];
sx q[1];
rz(-1.9221646) q[1];
rz(-pi) q[2];
rz(0.79244198) q[3];
sx q[3];
rz(-1.0002478) q[3];
sx q[3];
rz(-1.5287409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5639497) q[2];
sx q[2];
rz(-0.22259139) q[2];
sx q[2];
rz(2.8655748) q[2];
rz(2.4769619) q[3];
sx q[3];
rz(-0.97739995) q[3];
sx q[3];
rz(0.19653921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68263245) q[0];
sx q[0];
rz(-0.17550547) q[0];
sx q[0];
rz(1.578414) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-2.1140607) q[1];
sx q[1];
rz(-0.064090699) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10326021) q[0];
sx q[0];
rz(-1.1128367) q[0];
sx q[0];
rz(3.1151014) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59696609) q[2];
sx q[2];
rz(-0.709049) q[2];
sx q[2];
rz(0.89151357) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0371263) q[1];
sx q[1];
rz(-2.0580225) q[1];
sx q[1];
rz(-3.0477316) q[1];
rz(-0.80497165) q[3];
sx q[3];
rz(-2.0197778) q[3];
sx q[3];
rz(1.4624553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2854707) q[2];
sx q[2];
rz(-2.8557114) q[2];
sx q[2];
rz(0.90606436) q[2];
rz(1.2178577) q[3];
sx q[3];
rz(-1.2713615) q[3];
sx q[3];
rz(-1.0539894) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8080407) q[0];
sx q[0];
rz(-1.6186436) q[0];
sx q[0];
rz(2.9343395) q[0];
rz(-2.8605657) q[1];
sx q[1];
rz(-1.749975) q[1];
sx q[1];
rz(-0.71798807) q[1];
rz(0.40194527) q[2];
sx q[2];
rz(-1.6464169) q[2];
sx q[2];
rz(0.34435549) q[2];
rz(0.47697502) q[3];
sx q[3];
rz(-1.5581445) q[3];
sx q[3];
rz(2.1731044) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
