OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1146381) q[0];
sx q[0];
rz(-1.4517598) q[0];
sx q[0];
rz(-0.64557689) q[0];
rz(-2.7627856) q[1];
sx q[1];
rz(-1.768755) q[1];
sx q[1];
rz(-1.6436613) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49255532) q[0];
sx q[0];
rz(-1.7151924) q[0];
sx q[0];
rz(1.9306246) q[0];
rz(-pi) q[1];
rz(-1.475392) q[2];
sx q[2];
rz(-0.63268748) q[2];
sx q[2];
rz(-0.22104095) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.363417) q[1];
sx q[1];
rz(-0.84901224) q[1];
sx q[1];
rz(1.8587106) q[1];
x q[2];
rz(0.50819355) q[3];
sx q[3];
rz(-0.71532202) q[3];
sx q[3];
rz(3.1014266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.92007414) q[2];
sx q[2];
rz(-1.6211082) q[2];
sx q[2];
rz(0.34040889) q[2];
rz(-0.83299625) q[3];
sx q[3];
rz(-0.70327988) q[3];
sx q[3];
rz(-1.3566646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4689363) q[0];
sx q[0];
rz(-0.39009538) q[0];
sx q[0];
rz(0.54498589) q[0];
rz(-0.90822059) q[1];
sx q[1];
rz(-0.32749367) q[1];
sx q[1];
rz(0.82495904) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23404113) q[0];
sx q[0];
rz(-1.7579161) q[0];
sx q[0];
rz(-1.9812816) q[0];
x q[1];
rz(-0.40290515) q[2];
sx q[2];
rz(-1.5803792) q[2];
sx q[2];
rz(-0.314089) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.35926871) q[1];
sx q[1];
rz(-1.6869079) q[1];
sx q[1];
rz(2.5776723) q[1];
rz(0.19552688) q[3];
sx q[3];
rz(-1.0128847) q[3];
sx q[3];
rz(2.5415681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58723441) q[2];
sx q[2];
rz(-0.52296573) q[2];
sx q[2];
rz(2.9651802) q[2];
rz(-3.0107064) q[3];
sx q[3];
rz(-1.1754879) q[3];
sx q[3];
rz(-3.1089354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.762887) q[0];
sx q[0];
rz(-2.5580907) q[0];
sx q[0];
rz(2.6718455) q[0];
rz(-1.6167971) q[1];
sx q[1];
rz(-2.6641615) q[1];
sx q[1];
rz(0.038539561) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24667106) q[0];
sx q[0];
rz(-1.5630432) q[0];
sx q[0];
rz(-0.00064108032) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4171717) q[2];
sx q[2];
rz(-1.3426174) q[2];
sx q[2];
rz(-1.4128026) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4290532) q[1];
sx q[1];
rz(-2.6479811) q[1];
sx q[1];
rz(-1.7008971) q[1];
rz(-0.79629691) q[3];
sx q[3];
rz(-0.7634123) q[3];
sx q[3];
rz(0.89087668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5588351) q[2];
sx q[2];
rz(-2.0520515) q[2];
sx q[2];
rz(-0.91252404) q[2];
rz(-1.8330666) q[3];
sx q[3];
rz(-2.1378744) q[3];
sx q[3];
rz(-1.7002038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38347605) q[0];
sx q[0];
rz(-1.8261199) q[0];
sx q[0];
rz(-0.34969774) q[0];
rz(-1.8967459) q[1];
sx q[1];
rz(-2.8379776) q[1];
sx q[1];
rz(-2.9464338) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0965658) q[0];
sx q[0];
rz(-1.4782527) q[0];
sx q[0];
rz(-1.6916313) q[0];
rz(-pi) q[1];
rz(0.85907016) q[2];
sx q[2];
rz(-2.4675183) q[2];
sx q[2];
rz(-2.7328343) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.86707592) q[1];
sx q[1];
rz(-1.80596) q[1];
sx q[1];
rz(2.0348674) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5024662) q[3];
sx q[3];
rz(-1.5577321) q[3];
sx q[3];
rz(2.5845205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9106456) q[2];
sx q[2];
rz(-2.6516984) q[2];
sx q[2];
rz(1.8738497) q[2];
rz(-2.0729444) q[3];
sx q[3];
rz(-2.1290776) q[3];
sx q[3];
rz(2.3012565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0347663) q[0];
sx q[0];
rz(-1.6893457) q[0];
sx q[0];
rz(3.0019794) q[0];
rz(-1.0768249) q[1];
sx q[1];
rz(-2.1007517) q[1];
sx q[1];
rz(-0.37809125) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17079167) q[0];
sx q[0];
rz(-2.1149201) q[0];
sx q[0];
rz(2.8118954) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3230882) q[2];
sx q[2];
rz(-0.92220014) q[2];
sx q[2];
rz(1.6550145) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0873449) q[1];
sx q[1];
rz(-1.9964881) q[1];
sx q[1];
rz(0.40360968) q[1];
rz(1.3729172) q[3];
sx q[3];
rz(-1.1117001) q[3];
sx q[3];
rz(-1.0517694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.87171626) q[2];
sx q[2];
rz(-1.8173952) q[2];
sx q[2];
rz(-0.88796973) q[2];
rz(2.1652083) q[3];
sx q[3];
rz(-1.7247) q[3];
sx q[3];
rz(2.1358657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.3865005) q[0];
sx q[0];
rz(-3.0639102) q[0];
sx q[0];
rz(0.37762541) q[0];
rz(0.32304421) q[1];
sx q[1];
rz(-2.4781365) q[1];
sx q[1];
rz(2.2264218) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24419063) q[0];
sx q[0];
rz(-1.4892329) q[0];
sx q[0];
rz(0.32342644) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85257951) q[2];
sx q[2];
rz(-0.70194178) q[2];
sx q[2];
rz(-0.41350565) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9050161) q[1];
sx q[1];
rz(-1.6394776) q[1];
sx q[1];
rz(-2.828039) q[1];
rz(1.4573426) q[3];
sx q[3];
rz(-2.2621584) q[3];
sx q[3];
rz(-1.3464348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6043828) q[2];
sx q[2];
rz(-2.9769124) q[2];
sx q[2];
rz(0.79052314) q[2];
rz(2.8516155) q[3];
sx q[3];
rz(-0.73263779) q[3];
sx q[3];
rz(-2.524232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73615605) q[0];
sx q[0];
rz(-2.0744531) q[0];
sx q[0];
rz(-0.36703584) q[0];
rz(1.2332747) q[1];
sx q[1];
rz(-1.1739302) q[1];
sx q[1];
rz(1.4253915) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51222425) q[0];
sx q[0];
rz(-2.6656796) q[0];
sx q[0];
rz(-0.41643629) q[0];
rz(-3.0076214) q[2];
sx q[2];
rz(-2.6920762) q[2];
sx q[2];
rz(-1.9884895) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.62896171) q[1];
sx q[1];
rz(-2.1753864) q[1];
sx q[1];
rz(1.1058034) q[1];
x q[2];
rz(-0.013399259) q[3];
sx q[3];
rz(-2.074476) q[3];
sx q[3];
rz(-2.3434533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9795064) q[2];
sx q[2];
rz(-1.6493713) q[2];
sx q[2];
rz(1.9308176) q[2];
rz(-0.0018421729) q[3];
sx q[3];
rz(-0.76549923) q[3];
sx q[3];
rz(-0.65729284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2974671) q[0];
sx q[0];
rz(-2.5572889) q[0];
sx q[0];
rz(-0.076518245) q[0];
rz(-2.466295) q[1];
sx q[1];
rz(-2.8476871) q[1];
sx q[1];
rz(-2.3892367) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9936179) q[0];
sx q[0];
rz(-2.5320842) q[0];
sx q[0];
rz(0.57172914) q[0];
rz(-pi) q[1];
x q[1];
rz(1.760375) q[2];
sx q[2];
rz(-2.484212) q[2];
sx q[2];
rz(2.0331403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.016537746) q[1];
sx q[1];
rz(-2.9978831) q[1];
sx q[1];
rz(2.2347666) q[1];
rz(2.5520677) q[3];
sx q[3];
rz(-1.3988964) q[3];
sx q[3];
rz(0.37026438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.9383119) q[2];
sx q[2];
rz(-1.9613772) q[2];
sx q[2];
rz(-3.1414462) q[2];
rz(-2.032062) q[3];
sx q[3];
rz(-0.90679449) q[3];
sx q[3];
rz(-0.29763597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9922441) q[0];
sx q[0];
rz(-2.902817) q[0];
sx q[0];
rz(1.0060271) q[0];
rz(-2.8736615) q[1];
sx q[1];
rz(-1.3403099) q[1];
sx q[1];
rz(1.8267652) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6515503) q[0];
sx q[0];
rz(-0.74768066) q[0];
sx q[0];
rz(2.7689395) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4376274) q[2];
sx q[2];
rz(-1.2097934) q[2];
sx q[2];
rz(1.5732869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4970376) q[1];
sx q[1];
rz(-0.92004787) q[1];
sx q[1];
rz(-0.52849309) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1341348) q[3];
sx q[3];
rz(-2.3141626) q[3];
sx q[3];
rz(0.5479365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.43977794) q[2];
sx q[2];
rz(-2.5952227) q[2];
sx q[2];
rz(0.24547274) q[2];
rz(-0.43073511) q[3];
sx q[3];
rz(-2.0499178) q[3];
sx q[3];
rz(0.47732863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7219287) q[0];
sx q[0];
rz(-2.1491282) q[0];
sx q[0];
rz(2.8549109) q[0];
rz(-0.87896705) q[1];
sx q[1];
rz(-1.5022087) q[1];
sx q[1];
rz(-2.749696) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6673198) q[0];
sx q[0];
rz(-0.55756888) q[0];
sx q[0];
rz(1.5865109) q[0];
rz(-pi) q[1];
rz(2.2266065) q[2];
sx q[2];
rz(-1.6784385) q[2];
sx q[2];
rz(0.3273302) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0185768) q[1];
sx q[1];
rz(-1.3850817) q[1];
sx q[1];
rz(0.8155483) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48505731) q[3];
sx q[3];
rz(-0.61614803) q[3];
sx q[3];
rz(-2.171606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.60951704) q[2];
sx q[2];
rz(-0.993002) q[2];
sx q[2];
rz(1.9840476) q[2];
rz(-2.5907497) q[3];
sx q[3];
rz(-1.5778056) q[3];
sx q[3];
rz(1.9633861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3863603) q[0];
sx q[0];
rz(-1.7978783) q[0];
sx q[0];
rz(-1.88301) q[0];
rz(1.339284) q[1];
sx q[1];
rz(-2.5279999) q[1];
sx q[1];
rz(-2.7816714) q[1];
rz(2.3618868) q[2];
sx q[2];
rz(-1.69366) q[2];
sx q[2];
rz(2.2488307) q[2];
rz(-1.0212392) q[3];
sx q[3];
rz(-1.869535) q[3];
sx q[3];
rz(0.46365999) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];