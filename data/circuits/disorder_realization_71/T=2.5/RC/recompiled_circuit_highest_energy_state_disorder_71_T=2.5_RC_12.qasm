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
rz(0.55819297) q[0];
sx q[0];
rz(3.8605122) q[0];
sx q[0];
rz(10.100848) q[0];
rz(0.27322912) q[1];
sx q[1];
rz(-2.1844808) q[1];
sx q[1];
rz(0.91125429) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9118496) q[0];
sx q[0];
rz(-1.5485244) q[0];
sx q[0];
rz(0.020072083) q[0];
x q[1];
rz(1.5780294) q[2];
sx q[2];
rz(-2.1963495) q[2];
sx q[2];
rz(-0.48239732) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.377502) q[1];
sx q[1];
rz(-1.6030178) q[1];
sx q[1];
rz(2.568214) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0082672) q[3];
sx q[3];
rz(-0.76078868) q[3];
sx q[3];
rz(0.68851346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2216144) q[2];
sx q[2];
rz(-2.6875434) q[2];
sx q[2];
rz(-0.81217074) q[2];
rz(3.0657366) q[3];
sx q[3];
rz(-0.64801884) q[3];
sx q[3];
rz(-2.5465452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49038637) q[0];
sx q[0];
rz(-2.6631329) q[0];
sx q[0];
rz(1.867021) q[0];
rz(1.6365341) q[1];
sx q[1];
rz(-0.36205629) q[1];
sx q[1];
rz(1.9777745) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3904189) q[0];
sx q[0];
rz(-2.0295706) q[0];
sx q[0];
rz(2.075688) q[0];
x q[1];
rz(-1.1357434) q[2];
sx q[2];
rz(-2.0363303) q[2];
sx q[2];
rz(0.42143824) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5608985) q[1];
sx q[1];
rz(-0.95798641) q[1];
sx q[1];
rz(1.3152907) q[1];
rz(-pi) q[2];
rz(0.31812654) q[3];
sx q[3];
rz(-1.9788747) q[3];
sx q[3];
rz(2.7216743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.85094467) q[2];
sx q[2];
rz(-0.023357563) q[2];
sx q[2];
rz(-1.5811496) q[2];
rz(2.9110939) q[3];
sx q[3];
rz(-1.0483619) q[3];
sx q[3];
rz(-2.7010664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5480963) q[0];
sx q[0];
rz(-2.8289712) q[0];
sx q[0];
rz(0.12259677) q[0];
rz(-0.7971881) q[1];
sx q[1];
rz(-2.1295348) q[1];
sx q[1];
rz(3.0534993) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1377376) q[0];
sx q[0];
rz(-0.81043506) q[0];
sx q[0];
rz(1.9268981) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.083115301) q[2];
sx q[2];
rz(-1.3440455) q[2];
sx q[2];
rz(-0.47435681) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8574787) q[1];
sx q[1];
rz(-1.7076252) q[1];
sx q[1];
rz(-1.72422) q[1];
rz(-pi) q[2];
rz(0.53998472) q[3];
sx q[3];
rz(-1.267588) q[3];
sx q[3];
rz(-1.7005512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.88283551) q[2];
sx q[2];
rz(-1.5828524) q[2];
sx q[2];
rz(1.8176414) q[2];
rz(-2.6435408) q[3];
sx q[3];
rz(-0.96314722) q[3];
sx q[3];
rz(-0.74670416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28904706) q[0];
sx q[0];
rz(-2.9883224) q[0];
sx q[0];
rz(-0.14719851) q[0];
rz(2.4920801) q[1];
sx q[1];
rz(-1.5107061) q[1];
sx q[1];
rz(-1.6726327) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2788852) q[0];
sx q[0];
rz(-1.5636347) q[0];
sx q[0];
rz(1.347752) q[0];
rz(-0.18314731) q[2];
sx q[2];
rz(-1.9948655) q[2];
sx q[2];
rz(1.9522425) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8115639) q[1];
sx q[1];
rz(-2.4061476) q[1];
sx q[1];
rz(1.0388908) q[1];
rz(-pi) q[2];
rz(2.5178403) q[3];
sx q[3];
rz(-0.68188462) q[3];
sx q[3];
rz(-2.9870913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2088251) q[2];
sx q[2];
rz(-0.494445) q[2];
sx q[2];
rz(-0.24173582) q[2];
rz(2.406481) q[3];
sx q[3];
rz(-1.7416411) q[3];
sx q[3];
rz(2.476695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6786137) q[0];
sx q[0];
rz(-2.0942056) q[0];
sx q[0];
rz(0.12272923) q[0];
rz(-0.8853451) q[1];
sx q[1];
rz(-1.0095936) q[1];
sx q[1];
rz(-2.5811894) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3353676) q[0];
sx q[0];
rz(-2.0809552) q[0];
sx q[0];
rz(-0.2038184) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29717314) q[2];
sx q[2];
rz(-1.783833) q[2];
sx q[2];
rz(-1.7293255) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0593157) q[1];
sx q[1];
rz(-0.53566414) q[1];
sx q[1];
rz(1.0470864) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82920977) q[3];
sx q[3];
rz(-0.61625879) q[3];
sx q[3];
rz(2.5065638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.6178599) q[2];
sx q[2];
rz(-2.3475519) q[2];
sx q[2];
rz(0.20475556) q[2];
rz(0.93920416) q[3];
sx q[3];
rz(-1.771628) q[3];
sx q[3];
rz(0.088976629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8646249) q[0];
sx q[0];
rz(-3.0245916) q[0];
sx q[0];
rz(-2.4846039) q[0];
rz(-2.9981546) q[1];
sx q[1];
rz(-1.5564432) q[1];
sx q[1];
rz(2.4030446) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.448682) q[0];
sx q[0];
rz(-1.4553207) q[0];
sx q[0];
rz(-0.67489745) q[0];
x q[1];
rz(0.22566585) q[2];
sx q[2];
rz(-1.0899001) q[2];
sx q[2];
rz(-1.7250329) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2966531) q[1];
sx q[1];
rz(-0.9228188) q[1];
sx q[1];
rz(0.39877994) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.012771) q[3];
sx q[3];
rz(-2.2000562) q[3];
sx q[3];
rz(2.2037639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2584381) q[2];
sx q[2];
rz(-2.482735) q[2];
sx q[2];
rz(0.0030041791) q[2];
rz(-2.352412) q[3];
sx q[3];
rz(-2.7572258) q[3];
sx q[3];
rz(1.9743617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.060519144) q[0];
sx q[0];
rz(-2.7923212) q[0];
sx q[0];
rz(0.14807598) q[0];
rz(0.5332467) q[1];
sx q[1];
rz(-0.24594578) q[1];
sx q[1];
rz(1.6291133) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5083744) q[0];
sx q[0];
rz(-1.1092767) q[0];
sx q[0];
rz(1.2715696) q[0];
rz(-2.327269) q[2];
sx q[2];
rz(-1.0934798) q[2];
sx q[2];
rz(-2.5392591) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5606192) q[1];
sx q[1];
rz(-1.6915913) q[1];
sx q[1];
rz(1.9126152) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9307327) q[3];
sx q[3];
rz(-2.2818344) q[3];
sx q[3];
rz(0.27121997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.44275722) q[2];
sx q[2];
rz(-1.7008702) q[2];
sx q[2];
rz(-0.093078144) q[2];
rz(0.062151521) q[3];
sx q[3];
rz(-2.744894) q[3];
sx q[3];
rz(-1.2232346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1208948) q[0];
sx q[0];
rz(-2.5466205) q[0];
sx q[0];
rz(-2.4769532) q[0];
rz(-1.7918034) q[1];
sx q[1];
rz(-2.2638958) q[1];
sx q[1];
rz(-0.68914366) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4664554) q[0];
sx q[0];
rz(-0.27091089) q[0];
sx q[0];
rz(2.6938426) q[0];
x q[1];
rz(-2.1470239) q[2];
sx q[2];
rz(-2.3625018) q[2];
sx q[2];
rz(1.3887203) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7245111) q[1];
sx q[1];
rz(-1.2577836) q[1];
sx q[1];
rz(-1.723812) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6313305) q[3];
sx q[3];
rz(-1.7352967) q[3];
sx q[3];
rz(-1.9609083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0906715) q[2];
sx q[2];
rz(-2.4854269) q[2];
sx q[2];
rz(-1.3760759) q[2];
rz(-1.9164267) q[3];
sx q[3];
rz(-2.4296032) q[3];
sx q[3];
rz(0.041697748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10251481) q[0];
sx q[0];
rz(-3.097105) q[0];
sx q[0];
rz(-0.59120375) q[0];
rz(0.2419596) q[1];
sx q[1];
rz(-2.0457025) q[1];
sx q[1];
rz(-2.718149) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11277448) q[0];
sx q[0];
rz(-1.2679546) q[0];
sx q[0];
rz(-1.8870728) q[0];
x q[1];
rz(1.6744711) q[2];
sx q[2];
rz(-1.0846234) q[2];
sx q[2];
rz(2.1882265) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0656506) q[1];
sx q[1];
rz(-2.099875) q[1];
sx q[1];
rz(-0.34543646) q[1];
x q[2];
rz(-1.1657646) q[3];
sx q[3];
rz(-2.1486503) q[3];
sx q[3];
rz(-0.33194968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5847136) q[2];
sx q[2];
rz(-0.60606474) q[2];
sx q[2];
rz(-1.9752183) q[2];
rz(2.0139458) q[3];
sx q[3];
rz(-2.1731264) q[3];
sx q[3];
rz(-0.52496547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0135076) q[0];
sx q[0];
rz(-0.43376827) q[0];
sx q[0];
rz(0.67310131) q[0];
rz(-2.290944) q[1];
sx q[1];
rz(-1.3530082) q[1];
sx q[1];
rz(-2.8180715) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44064012) q[0];
sx q[0];
rz(-1.3214759) q[0];
sx q[0];
rz(-1.1855769) q[0];
rz(-pi) q[1];
rz(2.0704248) q[2];
sx q[2];
rz(-0.53348225) q[2];
sx q[2];
rz(-0.96772742) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9266745) q[1];
sx q[1];
rz(-1.2387382) q[1];
sx q[1];
rz(-2.3658793) q[1];
rz(-pi) q[2];
rz(-0.011908493) q[3];
sx q[3];
rz(-0.58387305) q[3];
sx q[3];
rz(-1.5990822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2181776) q[2];
sx q[2];
rz(-0.18109334) q[2];
sx q[2];
rz(3.0695445) q[2];
rz(-1.0142903) q[3];
sx q[3];
rz(-2.225596) q[3];
sx q[3];
rz(-0.54916507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2608248) q[0];
sx q[0];
rz(-1.4243955) q[0];
sx q[0];
rz(-2.0212174) q[0];
rz(-2.8063759) q[1];
sx q[1];
rz(-0.64246476) q[1];
sx q[1];
rz(-1.1186218) q[1];
rz(1.3041244) q[2];
sx q[2];
rz(-2.647807) q[2];
sx q[2];
rz(-0.37995445) q[2];
rz(0.14306457) q[3];
sx q[3];
rz(-1.5456556) q[3];
sx q[3];
rz(1.7329334) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
