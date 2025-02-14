OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11197055) q[0];
sx q[0];
rz(-2.114871) q[0];
sx q[0];
rz(-1.2929027) q[0];
rz(-1.0678043) q[1];
sx q[1];
rz(-0.81967241) q[1];
sx q[1];
rz(-1.9424865) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0584058) q[0];
sx q[0];
rz(-1.633344) q[0];
sx q[0];
rz(-3.0330171) q[0];
rz(-pi) q[1];
rz(-2.4680696) q[2];
sx q[2];
rz(-1.4269514) q[2];
sx q[2];
rz(1.0811445) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.73429062) q[1];
sx q[1];
rz(-1.3774187) q[1];
sx q[1];
rz(0.9949245) q[1];
x q[2];
rz(-2.2584003) q[3];
sx q[3];
rz(-0.84005648) q[3];
sx q[3];
rz(-0.20945621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6370411) q[2];
sx q[2];
rz(-1.7402288) q[2];
sx q[2];
rz(-1.6999647) q[2];
rz(1.0124892) q[3];
sx q[3];
rz(-1.9952521) q[3];
sx q[3];
rz(2.6700524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83990324) q[0];
sx q[0];
rz(-2.58044) q[0];
sx q[0];
rz(-1.0294234) q[0];
rz(-1.3599716) q[1];
sx q[1];
rz(-1.9554892) q[1];
sx q[1];
rz(-0.50672466) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27924867) q[0];
sx q[0];
rz(-1.5515258) q[0];
sx q[0];
rz(1.6750245) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0015745) q[2];
sx q[2];
rz(-2.3990941) q[2];
sx q[2];
rz(-3.1378821) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3054333) q[1];
sx q[1];
rz(-2.1483114) q[1];
sx q[1];
rz(-3.0633395) q[1];
x q[2];
rz(1.2348753) q[3];
sx q[3];
rz(-0.75390076) q[3];
sx q[3];
rz(0.17673211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3158675) q[2];
sx q[2];
rz(-0.59810144) q[2];
sx q[2];
rz(2.9244847) q[2];
rz(-2.121296) q[3];
sx q[3];
rz(-1.3291357) q[3];
sx q[3];
rz(0.98062688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26647907) q[0];
sx q[0];
rz(-1.7784235) q[0];
sx q[0];
rz(-0.76817051) q[0];
rz(-2.8504596) q[1];
sx q[1];
rz(-1.7765287) q[1];
sx q[1];
rz(-1.9545782) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465837) q[0];
sx q[0];
rz(-1.4872941) q[0];
sx q[0];
rz(0.21159391) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0303346) q[2];
sx q[2];
rz(-3.0560962) q[2];
sx q[2];
rz(0.70921113) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95988454) q[1];
sx q[1];
rz(-1.2165637) q[1];
sx q[1];
rz(-1.6003057) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95245016) q[3];
sx q[3];
rz(-0.72126167) q[3];
sx q[3];
rz(-2.9968468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27176303) q[2];
sx q[2];
rz(-0.84438476) q[2];
sx q[2];
rz(-0.96770206) q[2];
rz(-0.23396954) q[3];
sx q[3];
rz(-0.7015737) q[3];
sx q[3];
rz(-0.35681891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02803239) q[0];
sx q[0];
rz(-1.5505294) q[0];
sx q[0];
rz(-1.0796984) q[0];
rz(-2.8375541) q[1];
sx q[1];
rz(-1.1540776) q[1];
sx q[1];
rz(-1.3160204) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42342351) q[0];
sx q[0];
rz(-1.4470513) q[0];
sx q[0];
rz(2.607671) q[0];
rz(-pi) q[1];
rz(-1.448668) q[2];
sx q[2];
rz(-2.5425362) q[2];
sx q[2];
rz(1.0775104) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58319891) q[1];
sx q[1];
rz(-2.2474562) q[1];
sx q[1];
rz(-0.79460245) q[1];
rz(-pi) q[2];
rz(-1.7075482) q[3];
sx q[3];
rz(-2.7117808) q[3];
sx q[3];
rz(-0.97625247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6198081) q[2];
sx q[2];
rz(-1.2336171) q[2];
sx q[2];
rz(0.141315) q[2];
rz(-2.5610793) q[3];
sx q[3];
rz(-1.4760735) q[3];
sx q[3];
rz(2.9132402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88437951) q[0];
sx q[0];
rz(-2.3452106) q[0];
sx q[0];
rz(-1.4759395) q[0];
rz(-2.2085704) q[1];
sx q[1];
rz(-2.4304183) q[1];
sx q[1];
rz(1.7500386) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0071268) q[0];
sx q[0];
rz(-0.8392121) q[0];
sx q[0];
rz(-2.9842019) q[0];
rz(2.8093908) q[2];
sx q[2];
rz(-2.6966288) q[2];
sx q[2];
rz(-1.5277758) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.4310225) q[1];
sx q[1];
rz(-0.52174134) q[1];
sx q[1];
rz(1.6304593) q[1];
x q[2];
rz(1.0443893) q[3];
sx q[3];
rz(-0.59848753) q[3];
sx q[3];
rz(0.55803669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.89620227) q[2];
sx q[2];
rz(-0.37084493) q[2];
sx q[2];
rz(2.89213) q[2];
rz(-1.6438515) q[3];
sx q[3];
rz(-1.7353053) q[3];
sx q[3];
rz(2.7264285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6849218) q[0];
sx q[0];
rz(-0.49600729) q[0];
sx q[0];
rz(-1.7556835) q[0];
rz(-1.8798401) q[1];
sx q[1];
rz(-1.2077786) q[1];
sx q[1];
rz(-2.6127846) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.38976) q[0];
sx q[0];
rz(-2.6963137) q[0];
sx q[0];
rz(-1.5543429) q[0];
rz(2.4086921) q[2];
sx q[2];
rz(-0.7846047) q[2];
sx q[2];
rz(2.8200788) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4836135) q[1];
sx q[1];
rz(-2.2083862) q[1];
sx q[1];
rz(3.0581091) q[1];
x q[2];
rz(-1.7011794) q[3];
sx q[3];
rz(-0.88997872) q[3];
sx q[3];
rz(0.21475204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1888107) q[2];
sx q[2];
rz(-2.6691801) q[2];
sx q[2];
rz(-1.3713651) q[2];
rz(-2.93907) q[3];
sx q[3];
rz(-1.4684418) q[3];
sx q[3];
rz(-1.9523581) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33223575) q[0];
sx q[0];
rz(-1.3487331) q[0];
sx q[0];
rz(-2.9877544) q[0];
rz(-0.73506749) q[1];
sx q[1];
rz(-1.091833) q[1];
sx q[1];
rz(2.99248) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1510493) q[0];
sx q[0];
rz(-1.4874389) q[0];
sx q[0];
rz(0.91611422) q[0];
rz(-pi) q[1];
rz(1.1400827) q[2];
sx q[2];
rz(-0.39719279) q[2];
sx q[2];
rz(-2.6818962) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11483773) q[1];
sx q[1];
rz(-2.2828365) q[1];
sx q[1];
rz(-2.4898082) q[1];
rz(-pi) q[2];
rz(-3.0508196) q[3];
sx q[3];
rz(-2.0730264) q[3];
sx q[3];
rz(1.9050364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6814802) q[2];
sx q[2];
rz(-1.511938) q[2];
sx q[2];
rz(-2.6849449) q[2];
rz(1.7234507) q[3];
sx q[3];
rz(-2.2599615) q[3];
sx q[3];
rz(-0.68896967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5099455) q[0];
sx q[0];
rz(-2.8896285) q[0];
sx q[0];
rz(0.052635996) q[0];
rz(-1.3098199) q[1];
sx q[1];
rz(-0.24886623) q[1];
sx q[1];
rz(1.9356669) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1261189) q[0];
sx q[0];
rz(-2.0778225) q[0];
sx q[0];
rz(2.6958204) q[0];
x q[1];
rz(1.8032372) q[2];
sx q[2];
rz(-1.0224059) q[2];
sx q[2];
rz(0.91969925) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5284023) q[1];
sx q[1];
rz(-0.90189108) q[1];
sx q[1];
rz(1.7714798) q[1];
x q[2];
rz(-2.9067791) q[3];
sx q[3];
rz(-2.1800332) q[3];
sx q[3];
rz(-1.2138194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2974818) q[2];
sx q[2];
rz(-1.380912) q[2];
sx q[2];
rz(-2.5816176) q[2];
rz(2.0848134) q[3];
sx q[3];
rz(-0.70928514) q[3];
sx q[3];
rz(0.81282508) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63447222) q[0];
sx q[0];
rz(-1.7422603) q[0];
sx q[0];
rz(-1.5647474) q[0];
rz(1.5693846) q[1];
sx q[1];
rz(-1.1610616) q[1];
sx q[1];
rz(-0.80894583) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29306889) q[0];
sx q[0];
rz(-1.5331755) q[0];
sx q[0];
rz(0.49531181) q[0];
rz(-pi) q[1];
rz(0.81124146) q[2];
sx q[2];
rz(-1.8643856) q[2];
sx q[2];
rz(2.9305803) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29076838) q[1];
sx q[1];
rz(-1.5838669) q[1];
sx q[1];
rz(-0.95990659) q[1];
rz(-pi) q[2];
rz(1.7478757) q[3];
sx q[3];
rz(-1.5054323) q[3];
sx q[3];
rz(-3.1013427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6605777) q[2];
sx q[2];
rz(-0.95260859) q[2];
sx q[2];
rz(-3.0926404) q[2];
rz(-1.281721) q[3];
sx q[3];
rz(-0.47544026) q[3];
sx q[3];
rz(-1.4648645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82185164) q[0];
sx q[0];
rz(-2.8806683) q[0];
sx q[0];
rz(-1.9191746) q[0];
rz(-1.4756731) q[1];
sx q[1];
rz(-1.9404575) q[1];
sx q[1];
rz(2.6840721) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8638141) q[0];
sx q[0];
rz(-0.84054986) q[0];
sx q[0];
rz(1.6919096) q[0];
rz(-pi) q[1];
rz(2.2234688) q[2];
sx q[2];
rz(-0.61033536) q[2];
sx q[2];
rz(0.68896919) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.88120115) q[1];
sx q[1];
rz(-2.7335751) q[1];
sx q[1];
rz(2.6101006) q[1];
rz(-1.1531257) q[3];
sx q[3];
rz(-0.92781298) q[3];
sx q[3];
rz(-0.90211274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7484625) q[2];
sx q[2];
rz(-2.8905383) q[2];
sx q[2];
rz(-2.1012696) q[2];
rz(-0.48804247) q[3];
sx q[3];
rz(-1.7010242) q[3];
sx q[3];
rz(-2.603781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3671065) q[0];
sx q[0];
rz(-1.0716866) q[0];
sx q[0];
rz(-2.5197784) q[0];
rz(1.8472916) q[1];
sx q[1];
rz(-0.58302561) q[1];
sx q[1];
rz(2.2517712) q[1];
rz(-2.0001844) q[2];
sx q[2];
rz(-2.5038237) q[2];
sx q[2];
rz(0.43546168) q[2];
rz(0.46066416) q[3];
sx q[3];
rz(-0.74429269) q[3];
sx q[3];
rz(-0.85891354) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
