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
rz(-1.018723) q[0];
sx q[0];
rz(-0.85917226) q[0];
sx q[0];
rz(2.3265042) q[0];
rz(-1.1852784) q[1];
sx q[1];
rz(-1.4108682) q[1];
sx q[1];
rz(-2.0739532) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1045221) q[0];
sx q[0];
rz(-0.43209729) q[0];
sx q[0];
rz(-1.62754) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52402988) q[2];
sx q[2];
rz(-1.2175778) q[2];
sx q[2];
rz(-2.2729682) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5112665) q[1];
sx q[1];
rz(-2.8095803) q[1];
sx q[1];
rz(1.270833) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8772763) q[3];
sx q[3];
rz(-2.2393763) q[3];
sx q[3];
rz(0.17449915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4472569) q[2];
sx q[2];
rz(-2.3981636) q[2];
sx q[2];
rz(-1.2747964) q[2];
rz(0.46191195) q[3];
sx q[3];
rz(-0.67449823) q[3];
sx q[3];
rz(-2.0089202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.3502515) q[0];
sx q[0];
rz(-0.30838648) q[0];
sx q[0];
rz(1.2530918) q[0];
rz(-0.14532267) q[1];
sx q[1];
rz(-1.3959613) q[1];
sx q[1];
rz(-1.0911509) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6515132) q[0];
sx q[0];
rz(-1.3567748) q[0];
sx q[0];
rz(1.1351311) q[0];
rz(-pi) q[1];
rz(1.2107466) q[2];
sx q[2];
rz(-2.0031679) q[2];
sx q[2];
rz(0.52456174) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5877643) q[1];
sx q[1];
rz(-1.3979646) q[1];
sx q[1];
rz(-0.21765222) q[1];
rz(-0.66295538) q[3];
sx q[3];
rz(-2.5409449) q[3];
sx q[3];
rz(-1.8523703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6609409) q[2];
sx q[2];
rz(-1.3903684) q[2];
sx q[2];
rz(2.6118028) q[2];
rz(-2.3482813) q[3];
sx q[3];
rz(-1.5373693) q[3];
sx q[3];
rz(-0.62354273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(2.6563501) q[0];
sx q[0];
rz(-0.95004496) q[0];
sx q[0];
rz(0.52870885) q[0];
rz(-2.5953925) q[1];
sx q[1];
rz(-2.1827953) q[1];
sx q[1];
rz(-2.8012457) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3211283) q[0];
sx q[0];
rz(-1.2437151) q[0];
sx q[0];
rz(-1.8595247) q[0];
rz(2.9067847) q[2];
sx q[2];
rz(-1.7124512) q[2];
sx q[2];
rz(1.0740785) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.18225154) q[1];
sx q[1];
rz(-1.4762344) q[1];
sx q[1];
rz(1.5675117) q[1];
rz(-pi) q[2];
rz(-2.765708) q[3];
sx q[3];
rz(-2.2760227) q[3];
sx q[3];
rz(0.24209596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1722374) q[2];
sx q[2];
rz(-0.41219553) q[2];
sx q[2];
rz(2.7117512) q[2];
rz(2.3853081) q[3];
sx q[3];
rz(-2.9679306) q[3];
sx q[3];
rz(0.82591301) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7569358) q[0];
sx q[0];
rz(-1.5058368) q[0];
sx q[0];
rz(1.4935619) q[0];
rz(-0.76796302) q[1];
sx q[1];
rz(-2.6710644) q[1];
sx q[1];
rz(2.5951662) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7872629) q[0];
sx q[0];
rz(-0.096219115) q[0];
sx q[0];
rz(0.83258273) q[0];
rz(-2.3713263) q[2];
sx q[2];
rz(-2.4663743) q[2];
sx q[2];
rz(-2.9807621) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.029303251) q[1];
sx q[1];
rz(-1.7186856) q[1];
sx q[1];
rz(1.6690955) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91937842) q[3];
sx q[3];
rz(-2.1438144) q[3];
sx q[3];
rz(2.1185377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4297318) q[2];
sx q[2];
rz(-1.6542566) q[2];
sx q[2];
rz(0.30409733) q[2];
rz(-1.3652623) q[3];
sx q[3];
rz(-1.2208166) q[3];
sx q[3];
rz(-1.0767267) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6292608) q[0];
sx q[0];
rz(-0.4442232) q[0];
sx q[0];
rz(0.00057922676) q[0];
rz(-1.0821139) q[1];
sx q[1];
rz(-2.6172456) q[1];
sx q[1];
rz(-2.2023315) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0359373) q[0];
sx q[0];
rz(-0.93917423) q[0];
sx q[0];
rz(1.0211591) q[0];
x q[1];
rz(-1.0910735) q[2];
sx q[2];
rz(-0.66893286) q[2];
sx q[2];
rz(-0.26593966) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0013922) q[1];
sx q[1];
rz(-1.0085229) q[1];
sx q[1];
rz(2.8881025) q[1];
rz(2.9143798) q[3];
sx q[3];
rz(-1.2649396) q[3];
sx q[3];
rz(1.7744886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1028563) q[2];
sx q[2];
rz(-1.880371) q[2];
sx q[2];
rz(-0.2612513) q[2];
rz(0.86483613) q[3];
sx q[3];
rz(-1.7295001) q[3];
sx q[3];
rz(0.29204667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.29702) q[0];
sx q[0];
rz(-1.7140056) q[0];
sx q[0];
rz(0.20508668) q[0];
rz(-1.0467485) q[1];
sx q[1];
rz(-1.3958684) q[1];
sx q[1];
rz(-0.37839016) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8119733) q[0];
sx q[0];
rz(-1.1866335) q[0];
sx q[0];
rz(-0.20150082) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8257337) q[2];
sx q[2];
rz(-2.6489885) q[2];
sx q[2];
rz(1.619018) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7415315) q[1];
sx q[1];
rz(-2.0963256) q[1];
sx q[1];
rz(2.7353213) q[1];
rz(-pi) q[2];
rz(-0.50726733) q[3];
sx q[3];
rz(-2.3823822) q[3];
sx q[3];
rz(1.678148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4947074) q[2];
sx q[2];
rz(-1.9834221) q[2];
sx q[2];
rz(2.0873439) q[2];
rz(-2.4833637) q[3];
sx q[3];
rz(-0.64526486) q[3];
sx q[3];
rz(-2.6575507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9996416) q[0];
sx q[0];
rz(-1.9331837) q[0];
sx q[0];
rz(-2.5010338) q[0];
rz(-2.7236252) q[1];
sx q[1];
rz(-1.2203981) q[1];
sx q[1];
rz(0.80498615) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1520839) q[0];
sx q[0];
rz(-1.594127) q[0];
sx q[0];
rz(0.70744608) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4067602) q[2];
sx q[2];
rz(-0.8443588) q[2];
sx q[2];
rz(-3.022829) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5857081) q[1];
sx q[1];
rz(-1.5991365) q[1];
sx q[1];
rz(-3.1217087) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2748313) q[3];
sx q[3];
rz(-2.1419542) q[3];
sx q[3];
rz(-2.7829602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.54887041) q[2];
sx q[2];
rz(-2.1465116) q[2];
sx q[2];
rz(-0.76118809) q[2];
rz(-0.74782863) q[3];
sx q[3];
rz(-1.8239832) q[3];
sx q[3];
rz(-0.67659155) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51650301) q[0];
sx q[0];
rz(-0.28446063) q[0];
sx q[0];
rz(1.6647343) q[0];
rz(0.45627108) q[1];
sx q[1];
rz(-1.4197333) q[1];
sx q[1];
rz(-2.2705073) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3696156) q[0];
sx q[0];
rz(-0.61423683) q[0];
sx q[0];
rz(-3.1200039) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4159059) q[2];
sx q[2];
rz(-1.911507) q[2];
sx q[2];
rz(2.129385) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.217389) q[1];
sx q[1];
rz(-2.743209) q[1];
sx q[1];
rz(-2.5474362) q[1];
x q[2];
rz(-3.0814287) q[3];
sx q[3];
rz(-1.1451266) q[3];
sx q[3];
rz(-0.7983467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2386834) q[2];
sx q[2];
rz(-2.6148655) q[2];
sx q[2];
rz(0.61863679) q[2];
rz(-0.020261852) q[3];
sx q[3];
rz(-2.198115) q[3];
sx q[3];
rz(2.3616135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(3.0778377) q[0];
sx q[0];
rz(-1.5851333) q[0];
sx q[0];
rz(0.42386398) q[0];
rz(2.9546812) q[1];
sx q[1];
rz(-2.321545) q[1];
sx q[1];
rz(-1.6835469) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8114421) q[0];
sx q[0];
rz(-1.4528265) q[0];
sx q[0];
rz(-1.5850204) q[0];
rz(1.1122658) q[2];
sx q[2];
rz(-0.56156172) q[2];
sx q[2];
rz(0.70365814) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8193389) q[1];
sx q[1];
rz(-1.5194494) q[1];
sx q[1];
rz(3.0027185) q[1];
rz(-2.7856876) q[3];
sx q[3];
rz(-1.487046) q[3];
sx q[3];
rz(0.29415302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3905048) q[2];
sx q[2];
rz(-2.0743399) q[2];
sx q[2];
rz(-1.0941774) q[2];
rz(1.4607726) q[3];
sx q[3];
rz(-1.3760309) q[3];
sx q[3];
rz(-1.338753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(1.3222892) q[0];
sx q[0];
rz(-1.7604473) q[0];
sx q[0];
rz(1.5018916) q[0];
rz(2.8627401) q[1];
sx q[1];
rz(-2.1809705) q[1];
sx q[1];
rz(-1.0241114) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6449924) q[0];
sx q[0];
rz(-1.4787714) q[0];
sx q[0];
rz(-2.651398) q[0];
rz(-0.050692888) q[2];
sx q[2];
rz(-1.8452574) q[2];
sx q[2];
rz(-1.3433742) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2453277) q[1];
sx q[1];
rz(-0.30054856) q[1];
sx q[1];
rz(-0.73722571) q[1];
rz(-pi) q[2];
rz(0.82641469) q[3];
sx q[3];
rz(-1.3152988) q[3];
sx q[3];
rz(2.5662553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9378822) q[2];
sx q[2];
rz(-1.8546162) q[2];
sx q[2];
rz(2.6046275) q[2];
rz(-1.2282061) q[3];
sx q[3];
rz(-2.1824586) q[3];
sx q[3];
rz(0.29135191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1163597) q[0];
sx q[0];
rz(-1.3675084) q[0];
sx q[0];
rz(1.0687923) q[0];
rz(-0.7069201) q[1];
sx q[1];
rz(-1.128935) q[1];
sx q[1];
rz(3.1405906) q[1];
rz(-1.3798643) q[2];
sx q[2];
rz(-1.1715062) q[2];
sx q[2];
rz(1.7572335) q[2];
rz(-1.7892006) q[3];
sx q[3];
rz(-1.7159749) q[3];
sx q[3];
rz(0.80930474) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
