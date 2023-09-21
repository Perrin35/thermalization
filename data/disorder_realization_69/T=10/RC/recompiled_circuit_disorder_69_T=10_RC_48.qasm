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
rz(-1.3462892) q[1];
sx q[1];
rz(-0.68038124) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43359783) q[0];
sx q[0];
rz(-2.2888219) q[0];
sx q[0];
rz(-0.17311592) q[0];
rz(-pi) q[1];
rz(-1.275185) q[2];
sx q[2];
rz(-1.7586853) q[2];
sx q[2];
rz(1.0653898) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3375652) q[1];
sx q[1];
rz(-1.0684039) q[1];
sx q[1];
rz(2.875209) q[1];
rz(-pi) q[2];
rz(-0.33403553) q[3];
sx q[3];
rz(-1.738027) q[3];
sx q[3];
rz(2.2807896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.15158571) q[2];
sx q[2];
rz(-2.1534584) q[2];
sx q[2];
rz(-3.0541259) q[2];
rz(-2.4123689) q[3];
sx q[3];
rz(-2.7681523) q[3];
sx q[3];
rz(2.8698486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7049578) q[0];
sx q[0];
rz(-0.74902642) q[0];
sx q[0];
rz(-1.0013642) q[0];
rz(-2.9691866) q[1];
sx q[1];
rz(-1.1162076) q[1];
sx q[1];
rz(-0.52406812) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8952626) q[0];
sx q[0];
rz(-1.2753914) q[0];
sx q[0];
rz(-2.2569879) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0662202) q[2];
sx q[2];
rz(-0.76250695) q[2];
sx q[2];
rz(-2.5410595) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9429051) q[1];
sx q[1];
rz(-0.44829475) q[1];
sx q[1];
rz(-2.5793736) q[1];
rz(-2.8855521) q[3];
sx q[3];
rz(-0.55169332) q[3];
sx q[3];
rz(1.3224758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9841763) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(-1.3400419) q[2];
rz(2.346758) q[3];
sx q[3];
rz(-1.1398311) q[3];
sx q[3];
rz(0.036858233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79725093) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(-0.62477338) q[0];
rz(0.96427381) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(2.952081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7150869) q[0];
sx q[0];
rz(-2.3622892) q[0];
sx q[0];
rz(2.7515609) q[0];
rz(1.465221) q[2];
sx q[2];
rz(-1.9441248) q[2];
sx q[2];
rz(-2.2764652) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1358007) q[1];
sx q[1];
rz(-1.5151007) q[1];
sx q[1];
rz(0.23975753) q[1];
rz(-pi) q[2];
rz(-0.96954815) q[3];
sx q[3];
rz(-1.7161955) q[3];
sx q[3];
rz(-1.5090514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0345962) q[2];
sx q[2];
rz(-0.84356374) q[2];
sx q[2];
rz(-1.754388) q[2];
rz(0.43131367) q[3];
sx q[3];
rz(-1.8547736) q[3];
sx q[3];
rz(1.0881902) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65961924) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(-1.5959928) q[0];
rz(-1.1384456) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(-2.0514354) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7074993) q[0];
sx q[0];
rz(-0.543569) q[0];
sx q[0];
rz(-0.71732934) q[0];
rz(-pi) q[1];
rz(-2.3218669) q[2];
sx q[2];
rz(-1.956454) q[2];
sx q[2];
rz(-1.0939329) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.412147) q[1];
sx q[1];
rz(-0.73819654) q[1];
sx q[1];
rz(-0.63936887) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7100536) q[3];
sx q[3];
rz(-0.46436897) q[3];
sx q[3];
rz(0.73016703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99889341) q[2];
sx q[2];
rz(-0.45140758) q[2];
sx q[2];
rz(-0.85582716) q[2];
rz(-1.9479729) q[3];
sx q[3];
rz(-1.5193628) q[3];
sx q[3];
rz(2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6511433) q[0];
sx q[0];
rz(-2.1676846) q[0];
sx q[0];
rz(-2.9918616) q[0];
rz(-2.1504452) q[1];
sx q[1];
rz(-1.9350355) q[1];
sx q[1];
rz(1.1700464) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67384185) q[0];
sx q[0];
rz(-1.8918599) q[0];
sx q[0];
rz(-2.2348316) q[0];
x q[1];
rz(1.9636376) q[2];
sx q[2];
rz(-2.422214) q[2];
sx q[2];
rz(2.9658085) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7829195) q[1];
sx q[1];
rz(-1.1766608) q[1];
sx q[1];
rz(2.6637117) q[1];
x q[2];
rz(2.823579) q[3];
sx q[3];
rz(-0.56088305) q[3];
sx q[3];
rz(0.47919264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5912) q[2];
sx q[2];
rz(0.13723792) q[2];
rz(1.3373226) q[3];
sx q[3];
rz(-0.58115712) q[3];
sx q[3];
rz(1.2924682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9542434) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(-1.3504008) q[0];
rz(0.84287914) q[1];
sx q[1];
rz(-2.402669) q[1];
sx q[1];
rz(-2.4687016) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85149375) q[0];
sx q[0];
rz(-0.58304542) q[0];
sx q[0];
rz(1.3821938) q[0];
rz(0.024725155) q[2];
sx q[2];
rz(-0.85660663) q[2];
sx q[2];
rz(1.6753472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2313378) q[1];
sx q[1];
rz(-0.64097039) q[1];
sx q[1];
rz(2.7892968) q[1];
rz(-pi) q[2];
rz(2.7298922) q[3];
sx q[3];
rz(-1.6414335) q[3];
sx q[3];
rz(-2.8639346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.3488397) q[2];
sx q[2];
rz(-1.2092084) q[2];
sx q[2];
rz(-0.43506452) q[2];
rz(1.3600291) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(-2.8939261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(0.9687987) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(-2.0794179) q[0];
rz(2.0299714) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(-1.4020845) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7715493) q[0];
sx q[0];
rz(-1.84329) q[0];
sx q[0];
rz(0.38099654) q[0];
x q[1];
rz(-0.76485302) q[2];
sx q[2];
rz(-1.7424889) q[2];
sx q[2];
rz(-2.9219251) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2981616) q[1];
sx q[1];
rz(-1.9318252) q[1];
sx q[1];
rz(-1.7346738) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0038463) q[3];
sx q[3];
rz(-0.98689729) q[3];
sx q[3];
rz(-3.0018842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7193675) q[2];
sx q[2];
rz(-0.56754595) q[2];
sx q[2];
rz(-0.68022234) q[2];
rz(0.42823544) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(-1.359882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5087886) q[0];
sx q[0];
rz(-1.4161685) q[0];
sx q[0];
rz(1.5493786) q[0];
rz(2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(-0.54135281) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66051018) q[0];
sx q[0];
rz(-0.56429243) q[0];
sx q[0];
rz(3.1206467) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9114242) q[2];
sx q[2];
rz(-2.0001786) q[2];
sx q[2];
rz(0.18061772) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1047157) q[1];
sx q[1];
rz(-0.86537213) q[1];
sx q[1];
rz(0.55286644) q[1];
rz(-1.5726611) q[3];
sx q[3];
rz(-1.9511173) q[3];
sx q[3];
rz(-0.62455458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.044518746) q[2];
sx q[2];
rz(-0.35733435) q[2];
sx q[2];
rz(-0.77735916) q[2];
rz(-2.2903531) q[3];
sx q[3];
rz(-2.080353) q[3];
sx q[3];
rz(-1.8204934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.7580496) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(-0.0099649075) q[0];
rz(2.126157) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(1.4452971) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2485219) q[0];
sx q[0];
rz(-1.0501486) q[0];
sx q[0];
rz(-0.74670665) q[0];
rz(-2.5601013) q[2];
sx q[2];
rz(-0.98052374) q[2];
sx q[2];
rz(1.3587388) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0547486) q[1];
sx q[1];
rz(-0.63042414) q[1];
sx q[1];
rz(3.1151506) q[1];
x q[2];
rz(2.2441838) q[3];
sx q[3];
rz(-2.071278) q[3];
sx q[3];
rz(-2.9874779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4426667) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(-2.5342069) q[2];
rz(1.7025042) q[3];
sx q[3];
rz(-1.7581698) q[3];
sx q[3];
rz(-2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8059175) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(-0.6138531) q[0];
rz(-1.0461668) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(0.39224958) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0961571) q[0];
sx q[0];
rz(-2.3491612) q[0];
sx q[0];
rz(-2.5655454) q[0];
x q[1];
rz(-2.1456111) q[2];
sx q[2];
rz(-2.1157584) q[2];
sx q[2];
rz(-1.121322) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1019198) q[1];
sx q[1];
rz(-2.3452248) q[1];
sx q[1];
rz(2.2069195) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0030389) q[3];
sx q[3];
rz(-2.9346653) q[3];
sx q[3];
rz(0.11433992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0739416) q[2];
sx q[2];
rz(-1.757688) q[2];
sx q[2];
rz(0.60047853) q[2];
rz(-2.0742119) q[3];
sx q[3];
rz(-1.821527) q[3];
sx q[3];
rz(2.0480806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4338715) q[0];
sx q[0];
rz(-2.7086471) q[0];
sx q[0];
rz(1.4824296) q[0];
rz(-0.99647635) q[1];
sx q[1];
rz(-1.6468208) q[1];
sx q[1];
rz(-1.5368808) q[1];
rz(2.773182) q[2];
sx q[2];
rz(-2.340292) q[2];
sx q[2];
rz(3.1030263) q[2];
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