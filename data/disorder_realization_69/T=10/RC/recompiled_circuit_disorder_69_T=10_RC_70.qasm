OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34713137) q[0];
sx q[0];
rz(-1.0153271) q[0];
sx q[0];
rz(0.46749687) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(1.7953035) q[1];
sx q[1];
rz(10.105159) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9677833) q[0];
sx q[0];
rz(-0.73497226) q[0];
sx q[0];
rz(-1.3761139) q[0];
rz(-pi) q[1];
rz(1.8664076) q[2];
sx q[2];
rz(-1.7586853) q[2];
sx q[2];
rz(1.0653898) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.89741325) q[1];
sx q[1];
rz(-1.3379828) q[1];
sx q[1];
rz(1.0531055) q[1];
x q[2];
rz(-2.8075571) q[3];
sx q[3];
rz(-1.738027) q[3];
sx q[3];
rz(-2.2807896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.15158571) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(3.0541259) q[2];
rz(2.4123689) q[3];
sx q[3];
rz(-0.37344033) q[3];
sx q[3];
rz(2.8698486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7049578) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(-1.0013642) q[0];
rz(0.17240605) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(-0.52406812) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8952626) q[0];
sx q[0];
rz(-1.2753914) q[0];
sx q[0];
rz(-2.2569879) q[0];
x q[1];
rz(0.42627724) q[2];
sx q[2];
rz(-0.91765109) q[2];
sx q[2];
rz(-3.1003568) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8088278) q[1];
sx q[1];
rz(-1.1953224) q[1];
sx q[1];
rz(-1.3198225) q[1];
rz(-pi) q[2];
rz(1.4161795) q[3];
sx q[3];
rz(-2.1025476) q[3];
sx q[3];
rz(-1.0242517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.15741631) q[2];
sx q[2];
rz(-0.45980644) q[2];
sx q[2];
rz(1.3400419) q[2];
rz(-2.346758) q[3];
sx q[3];
rz(-1.1398311) q[3];
sx q[3];
rz(-0.036858233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3443417) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(-2.5168193) q[0];
rz(0.96427381) q[1];
sx q[1];
rz(-2.6565629) q[1];
sx q[1];
rz(-2.952081) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0441168) q[0];
sx q[0];
rz(-0.86320816) q[0];
sx q[0];
rz(-1.9301027) q[0];
rz(-pi) q[1];
rz(-1.6763716) q[2];
sx q[2];
rz(-1.1974679) q[2];
sx q[2];
rz(-0.8651274) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3526926) q[1];
sx q[1];
rz(-2.8955724) q[1];
sx q[1];
rz(-0.23060631) q[1];
rz(-pi) q[2];
rz(-2.965851) q[3];
sx q[3];
rz(-0.97676859) q[3];
sx q[3];
rz(0.16080805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1069964) q[2];
sx q[2];
rz(-0.84356374) q[2];
sx q[2];
rz(1.754388) q[2];
rz(0.43131367) q[3];
sx q[3];
rz(-1.286819) q[3];
sx q[3];
rz(2.0534024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65961924) q[0];
sx q[0];
rz(-0.94835931) q[0];
sx q[0];
rz(1.5455998) q[0];
rz(1.1384456) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(-1.0901573) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4340934) q[0];
sx q[0];
rz(-0.543569) q[0];
sx q[0];
rz(0.71732934) q[0];
rz(-1.0341187) q[2];
sx q[2];
rz(-2.3148429) q[2];
sx q[2];
rz(-0.85988753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.62413299) q[1];
sx q[1];
rz(-1.0003261) q[1];
sx q[1];
rz(2.0681417) q[1];
rz(-pi) q[2];
rz(-0.43153901) q[3];
sx q[3];
rz(-2.6772237) q[3];
sx q[3];
rz(2.4114256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99889341) q[2];
sx q[2];
rz(-2.6901851) q[2];
sx q[2];
rz(2.2857655) q[2];
rz(-1.1936197) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6511433) q[0];
sx q[0];
rz(-2.1676846) q[0];
sx q[0];
rz(2.9918616) q[0];
rz(-0.99114746) q[1];
sx q[1];
rz(-1.9350355) q[1];
sx q[1];
rz(-1.1700464) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65483246) q[0];
sx q[0];
rz(-2.1954384) q[0];
sx q[0];
rz(0.39958582) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32355002) q[2];
sx q[2];
rz(-0.91634446) q[2];
sx q[2];
rz(0.3277342) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7157117) q[1];
sx q[1];
rz(-2.5320146) q[1];
sx q[1];
rz(-0.73519911) q[1];
rz(-pi) q[2];
rz(-0.31801362) q[3];
sx q[3];
rz(-2.5807096) q[3];
sx q[3];
rz(-0.47919264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5503927) q[2];
sx q[2];
rz(3.0043547) q[2];
rz(-1.3373226) q[3];
sx q[3];
rz(-0.58115712) q[3];
sx q[3];
rz(1.8491245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9542434) q[0];
sx q[0];
rz(-1.3784778) q[0];
sx q[0];
rz(-1.3504008) q[0];
rz(2.2987135) q[1];
sx q[1];
rz(-2.402669) q[1];
sx q[1];
rz(2.4687016) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2900989) q[0];
sx q[0];
rz(-2.5585472) q[0];
sx q[0];
rz(-1.7593988) q[0];
rz(-pi) q[1];
x q[1];
rz(0.024725155) q[2];
sx q[2];
rz(-0.85660663) q[2];
sx q[2];
rz(-1.4662454) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3402965) q[1];
sx q[1];
rz(-2.166689) q[1];
sx q[1];
rz(-1.8227541) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17500413) q[3];
sx q[3];
rz(-0.41737469) q[3];
sx q[3];
rz(2.008703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.3488397) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(-0.43506452) q[2];
rz(-1.7815636) q[3];
sx q[3];
rz(-0.74917787) q[3];
sx q[3];
rz(-0.24766651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.172794) q[0];
sx q[0];
rz(-0.89247576) q[0];
sx q[0];
rz(2.0794179) q[0];
rz(2.0299714) q[1];
sx q[1];
rz(-1.9042791) q[1];
sx q[1];
rz(-1.7395082) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7715493) q[0];
sx q[0];
rz(-1.84329) q[0];
sx q[0];
rz(0.38099654) q[0];
rz(2.3767396) q[2];
sx q[2];
rz(-1.7424889) q[2];
sx q[2];
rz(-2.9219251) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3558823) q[1];
sx q[1];
rz(-1.4175698) q[1];
sx q[1];
rz(-0.36550891) q[1];
rz(-pi) q[2];
rz(-3.0038463) q[3];
sx q[3];
rz(-2.1546954) q[3];
sx q[3];
rz(0.13970845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7193675) q[2];
sx q[2];
rz(-0.56754595) q[2];
sx q[2];
rz(0.68022234) q[2];
rz(0.42823544) q[3];
sx q[3];
rz(-1.2546344) q[3];
sx q[3];
rz(1.359882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.632804) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(1.592214) q[0];
rz(-2.8920065) q[1];
sx q[1];
rz(-1.1523749) q[1];
sx q[1];
rz(-0.54135281) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66051018) q[0];
sx q[0];
rz(-0.56429243) q[0];
sx q[0];
rz(0.020945992) q[0];
rz(-1.2301684) q[2];
sx q[2];
rz(-2.0001786) q[2];
sx q[2];
rz(-2.9609749) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.036877) q[1];
sx q[1];
rz(-2.2762205) q[1];
sx q[1];
rz(-2.5887262) q[1];
rz(-pi) q[2];
x q[2];
rz(0.0046644966) q[3];
sx q[3];
rz(-2.7612673) q[3];
sx q[3];
rz(-0.62957803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0970739) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(-0.77735916) q[2];
rz(0.85123953) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(-1.3210993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.7580496) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(-0.0099649075) q[0];
rz(-1.0154356) q[1];
sx q[1];
rz(-0.76428691) q[1];
sx q[1];
rz(1.6962956) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89307071) q[0];
sx q[0];
rz(-1.0501486) q[0];
sx q[0];
rz(-2.394886) q[0];
x q[1];
rz(-2.2575349) q[2];
sx q[2];
rz(-2.3381655) q[2];
sx q[2];
rz(-0.49056177) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6362793) q[1];
sx q[1];
rz(-1.5552102) q[1];
sx q[1];
rz(0.63025766) q[1];
rz(-pi) q[2];
rz(2.290756) q[3];
sx q[3];
rz(-0.81504226) q[3];
sx q[3];
rz(-1.1834061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4426667) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(2.5342069) q[2];
rz(-1.7025042) q[3];
sx q[3];
rz(-1.7581698) q[3];
sx q[3];
rz(-0.79090345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8059175) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(0.6138531) q[0];
rz(2.0954258) q[1];
sx q[1];
rz(-2.8764953) q[1];
sx q[1];
rz(-0.39224958) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0948254) q[0];
sx q[0];
rz(-1.1724768) q[0];
sx q[0];
rz(2.4368068) q[0];
x q[1];
rz(-2.5160772) q[2];
sx q[2];
rz(-2.054347) q[2];
sx q[2];
rz(0.12550437) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2892462) q[1];
sx q[1];
rz(-0.95818555) q[1];
sx q[1];
rz(-0.5457408) q[1];
rz(-pi) q[2];
rz(-2.936593) q[3];
sx q[3];
rz(-1.5991755) q[3];
sx q[3];
rz(1.8207707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.067651) q[2];
sx q[2];
rz(-1.757688) q[2];
sx q[2];
rz(2.5411141) q[2];
rz(-2.0742119) q[3];
sx q[3];
rz(-1.821527) q[3];
sx q[3];
rz(2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7077211) q[0];
sx q[0];
rz(-0.43294551) q[0];
sx q[0];
rz(-1.659163) q[0];
rz(0.99647635) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(-2.773182) q[2];
sx q[2];
rz(-0.80130063) q[2];
sx q[2];
rz(-0.038566312) q[2];
rz(0.92580933) q[3];
sx q[3];
rz(-0.79173294) q[3];
sx q[3];
rz(-1.790495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
