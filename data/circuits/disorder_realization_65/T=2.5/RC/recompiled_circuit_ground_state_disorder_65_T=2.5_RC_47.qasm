OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.559691) q[0];
sx q[0];
rz(-0.087495916) q[0];
sx q[0];
rz(-3.0425332) q[0];
rz(-0.29095185) q[1];
sx q[1];
rz(-1.3909611) q[1];
sx q[1];
rz(-1.6265534) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5387548) q[0];
sx q[0];
rz(-1.8010209) q[0];
sx q[0];
rz(-1.099596) q[0];
rz(-pi) q[1];
rz(-1.0617822) q[2];
sx q[2];
rz(-1.1689593) q[2];
sx q[2];
rz(-0.33515829) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5196621) q[1];
sx q[1];
rz(-2.3708713) q[1];
sx q[1];
rz(-1.3407441) q[1];
x q[2];
rz(-2.4060449) q[3];
sx q[3];
rz(-2.7233019) q[3];
sx q[3];
rz(2.4052704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0809975) q[2];
sx q[2];
rz(-3.1298349) q[2];
sx q[2];
rz(0.11059977) q[2];
rz(3.1317173) q[3];
sx q[3];
rz(-0.014009352) q[3];
sx q[3];
rz(-0.88147718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.9843543) q[0];
sx q[0];
rz(-0.088033661) q[0];
sx q[0];
rz(-1.3798168) q[0];
rz(-3.1294322) q[1];
sx q[1];
rz(-2.1575243) q[1];
sx q[1];
rz(1.5252569) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3846958) q[0];
sx q[0];
rz(-1.257557) q[0];
sx q[0];
rz(1.5114853) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41864606) q[2];
sx q[2];
rz(-1.3867298) q[2];
sx q[2];
rz(-2.7238562) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.42022195) q[1];
sx q[1];
rz(-2.9530912) q[1];
sx q[1];
rz(-2.7367715) q[1];
x q[2];
rz(-0.44604519) q[3];
sx q[3];
rz(-1.6823263) q[3];
sx q[3];
rz(-2.1457637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81185594) q[2];
sx q[2];
rz(-0.024710329) q[2];
sx q[2];
rz(-0.032026637) q[2];
rz(0.61102593) q[3];
sx q[3];
rz(-1.150584) q[3];
sx q[3];
rz(-1.376763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5585612) q[0];
sx q[0];
rz(-0.91201454) q[0];
sx q[0];
rz(-1.6770021) q[0];
rz(1.5088082) q[1];
sx q[1];
rz(-1.2631515) q[1];
sx q[1];
rz(-0.060922932) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3896854) q[0];
sx q[0];
rz(-0.42493978) q[0];
sx q[0];
rz(-2.6838046) q[0];
rz(-pi) q[1];
rz(-3.0369122) q[2];
sx q[2];
rz(-1.6828244) q[2];
sx q[2];
rz(-0.45272967) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0197501) q[1];
sx q[1];
rz(-1.3977244) q[1];
sx q[1];
rz(-2.9815841) q[1];
rz(-pi) q[2];
rz(-0.64843602) q[3];
sx q[3];
rz(-0.47659527) q[3];
sx q[3];
rz(-2.4364382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9691201) q[2];
sx q[2];
rz(-3.1388333) q[2];
sx q[2];
rz(2.3976682) q[2];
rz(0.38532358) q[3];
sx q[3];
rz(-0.0016366882) q[3];
sx q[3];
rz(2.2374432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0124403) q[0];
sx q[0];
rz(-0.04239447) q[0];
sx q[0];
rz(2.2989035) q[0];
rz(-2.1878302) q[1];
sx q[1];
rz(-1.6396435) q[1];
sx q[1];
rz(-1.5488254) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5935338) q[0];
sx q[0];
rz(-1.6879544) q[0];
sx q[0];
rz(-1.5911932) q[0];
rz(-pi) q[1];
rz(-1.9325402) q[2];
sx q[2];
rz(-0.42297034) q[2];
sx q[2];
rz(-1.7646546) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6083513) q[1];
sx q[1];
rz(-2.5035739) q[1];
sx q[1];
rz(1.649943) q[1];
rz(-pi) q[2];
rz(-1.4709267) q[3];
sx q[3];
rz(-0.79260561) q[3];
sx q[3];
rz(2.2865393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.039975) q[2];
sx q[2];
rz(-0.016923252) q[2];
sx q[2];
rz(-2.7162111) q[2];
rz(1.3500805) q[3];
sx q[3];
rz(-1.5964419) q[3];
sx q[3];
rz(-2.8369956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72163433) q[0];
sx q[0];
rz(-0.0036792734) q[0];
sx q[0];
rz(-2.4256308) q[0];
rz(-1.6250027) q[1];
sx q[1];
rz(-0.45418987) q[1];
sx q[1];
rz(2.9862278) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10532285) q[0];
sx q[0];
rz(-1.5548517) q[0];
sx q[0];
rz(1.6769958) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4808319) q[2];
sx q[2];
rz(-1.9461921) q[2];
sx q[2];
rz(-1.9824972) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1137312) q[1];
sx q[1];
rz(-1.8370973) q[1];
sx q[1];
rz(-1.4425502) q[1];
rz(-pi) q[2];
rz(2.0955438) q[3];
sx q[3];
rz(-2.397846) q[3];
sx q[3];
rz(0.50713563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4510437) q[2];
sx q[2];
rz(-1.1210219) q[2];
sx q[2];
rz(1.6070018) q[2];
rz(1.0891886) q[3];
sx q[3];
rz(-0.023622731) q[3];
sx q[3];
rz(1.7781809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71517336) q[0];
sx q[0];
rz(-0.024113163) q[0];
sx q[0];
rz(0.56060785) q[0];
rz(-0.86588612) q[1];
sx q[1];
rz(-3.1325603) q[1];
sx q[1];
rz(-2.3057356) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3535806) q[0];
sx q[0];
rz(-0.29639527) q[0];
sx q[0];
rz(0.061876492) q[0];
rz(-pi) q[1];
rz(-1.4908815) q[2];
sx q[2];
rz(-0.40849388) q[2];
sx q[2];
rz(-1.6871883) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.33619704) q[1];
sx q[1];
rz(-1.8289036) q[1];
sx q[1];
rz(-1.6805524) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8172328) q[3];
sx q[3];
rz(-1.8772238) q[3];
sx q[3];
rz(0.6014733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5638289) q[2];
sx q[2];
rz(-1.8356297) q[2];
sx q[2];
rz(-1.5527661) q[2];
rz(-2.0896437) q[3];
sx q[3];
rz(-0.51783872) q[3];
sx q[3];
rz(0.53798211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2688667) q[0];
sx q[0];
rz(-2.5476542) q[0];
sx q[0];
rz(1.3237704) q[0];
rz(-2.2920604) q[1];
sx q[1];
rz(-3.140026) q[1];
sx q[1];
rz(-0.82226396) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36997312) q[0];
sx q[0];
rz(-1.7499588) q[0];
sx q[0];
rz(-0.097296299) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5157813) q[2];
sx q[2];
rz(-1.9085247) q[2];
sx q[2];
rz(1.6134451) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8755799) q[1];
sx q[1];
rz(-3.071738) q[1];
sx q[1];
rz(1.477398) q[1];
rz(-pi) q[2];
rz(-2.5093171) q[3];
sx q[3];
rz(-0.32514363) q[3];
sx q[3];
rz(-2.7141311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5795341) q[2];
sx q[2];
rz(-1.9066533) q[2];
sx q[2];
rz(-0.84260064) q[2];
rz(0.54251999) q[3];
sx q[3];
rz(-3.123057) q[3];
sx q[3];
rz(0.59244573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0845959) q[0];
sx q[0];
rz(-1.5942986) q[0];
sx q[0];
rz(-1.0513167) q[0];
rz(-1.3732789) q[1];
sx q[1];
rz(-0.028742464) q[1];
sx q[1];
rz(-1.2398667) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1662906) q[0];
sx q[0];
rz(-0.24595114) q[0];
sx q[0];
rz(2.9800132) q[0];
rz(-pi) q[1];
rz(-1.6592024) q[2];
sx q[2];
rz(-2.7475221) q[2];
sx q[2];
rz(1.566377) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6426601) q[1];
sx q[1];
rz(-0.98915425) q[1];
sx q[1];
rz(-10/(13*pi)) q[1];
rz(-1.8375592) q[3];
sx q[3];
rz(-1.3893616) q[3];
sx q[3];
rz(2.415588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5956868) q[2];
sx q[2];
rz(-0.02820708) q[2];
sx q[2];
rz(2.9941881) q[2];
rz(-1.5386511) q[3];
sx q[3];
rz(-1.7037062) q[3];
sx q[3];
rz(-2.9586207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1786757) q[0];
sx q[0];
rz(-0.33483949) q[0];
sx q[0];
rz(-1.8033002) q[0];
rz(1.6637404) q[1];
sx q[1];
rz(-2.0070952) q[1];
sx q[1];
rz(-0.17325625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6918744) q[0];
sx q[0];
rz(-1.0074179) q[0];
sx q[0];
rz(2.0823576) q[0];
rz(-0.40399) q[2];
sx q[2];
rz(-1.0183909) q[2];
sx q[2];
rz(-0.45726038) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.063364) q[1];
sx q[1];
rz(-1.8174108) q[1];
sx q[1];
rz(1.4117086) q[1];
rz(-pi) q[2];
rz(1.162649) q[3];
sx q[3];
rz(-0.014685304) q[3];
sx q[3];
rz(2.5774308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.133539) q[2];
sx q[2];
rz(-0.011913813) q[2];
sx q[2];
rz(2.1214205) q[2];
rz(3.0607439) q[3];
sx q[3];
rz(-3.1404218) q[3];
sx q[3];
rz(-1.1297869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6218277) q[0];
sx q[0];
rz(-0.28525678) q[0];
sx q[0];
rz(-1.6602302) q[0];
rz(0.18352428) q[1];
sx q[1];
rz(-2.8158975) q[1];
sx q[1];
rz(3.0537925) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7526898) q[0];
sx q[0];
rz(-1.4702893) q[0];
sx q[0];
rz(2.697562) q[0];
rz(2.8864809) q[2];
sx q[2];
rz(-2.8255551) q[2];
sx q[2];
rz(0.46816805) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1147881) q[1];
sx q[1];
rz(-2.8685796) q[1];
sx q[1];
rz(-2.6778912) q[1];
rz(-2.4724768) q[3];
sx q[3];
rz(-1.3312311) q[3];
sx q[3];
rz(-1.4544044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8537019) q[2];
sx q[2];
rz(-0.012306865) q[2];
sx q[2];
rz(-2.2575209) q[2];
rz(-2.4331802) q[3];
sx q[3];
rz(-0.00021472308) q[3];
sx q[3];
rz(0.29402548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52083279) q[0];
sx q[0];
rz(-1.5694869) q[0];
sx q[0];
rz(-1.6318305) q[0];
rz(0.028342551) q[1];
sx q[1];
rz(-0.62348532) q[1];
sx q[1];
rz(-3.0709406) q[1];
rz(2.4653409) q[2];
sx q[2];
rz(-0.55174232) q[2];
sx q[2];
rz(-0.55773013) q[2];
rz(-2.4130751) q[3];
sx q[3];
rz(-1.3195724) q[3];
sx q[3];
rz(2.6825678) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
