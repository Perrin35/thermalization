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
rz(5.4240131) q[0];
sx q[0];
rz(11.751282) q[0];
rz(-1.1852784) q[1];
sx q[1];
rz(-1.4108682) q[1];
sx q[1];
rz(-2.0739532) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6237804) q[0];
sx q[0];
rz(-1.5470439) q[0];
sx q[0];
rz(-2.0022814) q[0];
rz(-pi) q[1];
rz(1.9733623) q[2];
sx q[2];
rz(-2.0595008) q[2];
sx q[2];
rz(2.2421064) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8274535) q[1];
sx q[1];
rz(-1.2541391) q[1];
sx q[1];
rz(0.10152557) q[1];
x q[2];
rz(-0.88495636) q[3];
sx q[3];
rz(-1.3643294) q[3];
sx q[3];
rz(-1.5790758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.69433576) q[2];
sx q[2];
rz(-0.7434291) q[2];
sx q[2];
rz(-1.2747964) q[2];
rz(-0.46191195) q[3];
sx q[3];
rz(-2.4670944) q[3];
sx q[3];
rz(-2.0089202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-2.3502515) q[0];
sx q[0];
rz(-2.8332062) q[0];
sx q[0];
rz(1.2530918) q[0];
rz(-0.14532267) q[1];
sx q[1];
rz(-1.3959613) q[1];
sx q[1];
rz(2.0504418) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6515132) q[0];
sx q[0];
rz(-1.7848178) q[0];
sx q[0];
rz(2.0064615) q[0];
x q[1];
rz(0.45812313) q[2];
sx q[2];
rz(-1.8964185) q[2];
sx q[2];
rz(1.9389012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.021046358) q[1];
sx q[1];
rz(-1.3564377) q[1];
sx q[1];
rz(1.3938741) q[1];
x q[2];
rz(-1.1717779) q[3];
sx q[3];
rz(-1.1091091) q[3];
sx q[3];
rz(2.6103013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.48065177) q[2];
sx q[2];
rz(-1.7512243) q[2];
sx q[2];
rz(-2.6118028) q[2];
rz(2.3482813) q[3];
sx q[3];
rz(-1.5373693) q[3];
sx q[3];
rz(0.62354273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48524258) q[0];
sx q[0];
rz(-2.1915477) q[0];
sx q[0];
rz(-2.6128838) q[0];
rz(2.5953925) q[1];
sx q[1];
rz(-2.1827953) q[1];
sx q[1];
rz(-0.34034696) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8204644) q[0];
sx q[0];
rz(-1.2437151) q[0];
sx q[0];
rz(-1.2820679) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54989782) q[2];
sx q[2];
rz(-2.8680414) q[2];
sx q[2];
rz(-0.036444681) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9941147) q[1];
sx q[1];
rz(-0.094618751) q[1];
sx q[1];
rz(-0.034618207) q[1];
rz(2.765708) q[3];
sx q[3];
rz(-0.86556992) q[3];
sx q[3];
rz(-2.8994967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9693552) q[2];
sx q[2];
rz(-0.41219553) q[2];
sx q[2];
rz(0.42984143) q[2];
rz(-0.75628453) q[3];
sx q[3];
rz(-0.1736621) q[3];
sx q[3];
rz(2.3156796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38465685) q[0];
sx q[0];
rz(-1.5058368) q[0];
sx q[0];
rz(1.6480308) q[0];
rz(2.3736296) q[1];
sx q[1];
rz(-2.6710644) q[1];
sx q[1];
rz(-0.54642645) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7872629) q[0];
sx q[0];
rz(-3.0453735) q[0];
sx q[0];
rz(0.83258273) q[0];
rz(0.52164061) q[2];
sx q[2];
rz(-1.120479) q[2];
sx q[2];
rz(0.7618103) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5560234) q[1];
sx q[1];
rz(-1.6680191) q[1];
sx q[1];
rz(-2.9929964) q[1];
rz(2.3872822) q[3];
sx q[3];
rz(-0.83900634) q[3];
sx q[3];
rz(-3.0712002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7118608) q[2];
sx q[2];
rz(-1.6542566) q[2];
sx q[2];
rz(2.8374953) q[2];
rz(1.3652623) q[3];
sx q[3];
rz(-1.2208166) q[3];
sx q[3];
rz(1.0767267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51233184) q[0];
sx q[0];
rz(-0.4442232) q[0];
sx q[0];
rz(3.1410134) q[0];
rz(2.0594788) q[1];
sx q[1];
rz(-2.6172456) q[1];
sx q[1];
rz(0.93926114) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2596596) q[0];
sx q[0];
rz(-1.1355917) q[0];
sx q[0];
rz(2.4324904) q[0];
rz(0.3498431) q[2];
sx q[2];
rz(-2.1533384) q[2];
sx q[2];
rz(2.2900641) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.68785948) q[1];
sx q[1];
rz(-0.61111585) q[1];
sx q[1];
rz(-1.9495717) q[1];
rz(-pi) q[2];
rz(-0.22721283) q[3];
sx q[3];
rz(-1.876653) q[3];
sx q[3];
rz(-1.7744886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1028563) q[2];
sx q[2];
rz(-1.880371) q[2];
sx q[2];
rz(0.2612513) q[2];
rz(-0.86483613) q[3];
sx q[3];
rz(-1.7295001) q[3];
sx q[3];
rz(-0.29204667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.29702) q[0];
sx q[0];
rz(-1.427587) q[0];
sx q[0];
rz(-2.936506) q[0];
rz(2.0948441) q[1];
sx q[1];
rz(-1.3958684) q[1];
sx q[1];
rz(2.7632025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9721824) q[0];
sx q[0];
rz(-2.7101288) q[0];
sx q[0];
rz(1.1110825) q[0];
rz(0.47176265) q[2];
sx q[2];
rz(-1.4233575) q[2];
sx q[2];
rz(-0.23213895) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.69120317) q[1];
sx q[1];
rz(-0.65237633) q[1];
sx q[1];
rz(2.1689586) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.50726733) q[3];
sx q[3];
rz(-0.75921042) q[3];
sx q[3];
rz(1.4634446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64688524) q[2];
sx q[2];
rz(-1.9834221) q[2];
sx q[2];
rz(2.0873439) q[2];
rz(2.4833637) q[3];
sx q[3];
rz(-2.4963278) q[3];
sx q[3];
rz(0.48404199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14195104) q[0];
sx q[0];
rz(-1.208409) q[0];
sx q[0];
rz(0.64055881) q[0];
rz(0.41796747) q[1];
sx q[1];
rz(-1.2203981) q[1];
sx q[1];
rz(-2.3366065) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55401257) q[0];
sx q[0];
rz(-2.4338284) q[0];
sx q[0];
rz(-3.1057024) q[0];
rz(-pi) q[1];
rz(-1.1518794) q[2];
sx q[2];
rz(-0.81406128) q[2];
sx q[2];
rz(0.69413041) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9737712) q[1];
sx q[1];
rz(-0.034618363) q[1];
sx q[1];
rz(0.95914118) q[1];
rz(-pi) q[2];
rz(-2.715519) q[3];
sx q[3];
rz(-0.63562993) q[3];
sx q[3];
rz(0.15492188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.54887041) q[2];
sx q[2];
rz(-0.9950811) q[2];
sx q[2];
rz(-2.3804046) q[2];
rz(0.74782863) q[3];
sx q[3];
rz(-1.8239832) q[3];
sx q[3];
rz(-2.4650011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51650301) q[0];
sx q[0];
rz(-0.28446063) q[0];
sx q[0];
rz(1.6647343) q[0];
rz(-2.6853216) q[1];
sx q[1];
rz(-1.4197333) q[1];
sx q[1];
rz(0.87108535) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3696156) q[0];
sx q[0];
rz(-0.61423683) q[0];
sx q[0];
rz(0.02158879) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4159059) q[2];
sx q[2];
rz(-1.2300856) q[2];
sx q[2];
rz(-1.0122077) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5567315) q[1];
sx q[1];
rz(-1.24354) q[1];
sx q[1];
rz(-1.3393988) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.060163946) q[3];
sx q[3];
rz(-1.1451266) q[3];
sx q[3];
rz(0.7983467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.90290922) q[2];
sx q[2];
rz(-2.6148655) q[2];
sx q[2];
rz(-0.61863679) q[2];
rz(-3.1213308) q[3];
sx q[3];
rz(-0.94347763) q[3];
sx q[3];
rz(-0.77997911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.063754931) q[0];
sx q[0];
rz(-1.5851333) q[0];
sx q[0];
rz(0.42386398) q[0];
rz(-0.18691143) q[1];
sx q[1];
rz(-2.321545) q[1];
sx q[1];
rz(1.4580457) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9317212) q[0];
sx q[0];
rz(-3.0227724) q[0];
sx q[0];
rz(-0.11943905) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1122658) q[2];
sx q[2];
rz(-2.5800309) q[2];
sx q[2];
rz(-2.4379345) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.600454) q[1];
sx q[1];
rz(-0.14800528) q[1];
sx q[1];
rz(-2.7861094) q[1];
rz(-pi) q[2];
rz(-2.7856876) q[3];
sx q[3];
rz(-1.6545466) q[3];
sx q[3];
rz(2.8474396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.75108782) q[2];
sx q[2];
rz(-1.0672528) q[2];
sx q[2];
rz(2.0474153) q[2];
rz(-1.4607726) q[3];
sx q[3];
rz(-1.7655617) q[3];
sx q[3];
rz(1.8028397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8193034) q[0];
sx q[0];
rz(-1.7604473) q[0];
sx q[0];
rz(-1.639701) q[0];
rz(-2.8627401) q[1];
sx q[1];
rz(-0.96062213) q[1];
sx q[1];
rz(2.1174812) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1163961) q[0];
sx q[0];
rz(-1.0828583) q[0];
sx q[0];
rz(1.4665718) q[0];
rz(-pi) q[1];
rz(3.0908998) q[2];
sx q[2];
rz(-1.8452574) q[2];
sx q[2];
rz(1.7982184) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0054202) q[1];
sx q[1];
rz(-1.7917624) q[1];
sx q[1];
rz(-1.7762089) q[1];
rz(-pi) q[2];
rz(-1.9387705) q[3];
sx q[3];
rz(-2.362613) q[3];
sx q[3];
rz(-0.72768962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9378822) q[2];
sx q[2];
rz(-1.8546162) q[2];
sx q[2];
rz(0.53696519) q[2];
rz(-1.9133866) q[3];
sx q[3];
rz(-0.95913404) q[3];
sx q[3];
rz(0.29135191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(1.1163597) q[0];
sx q[0];
rz(-1.7740842) q[0];
sx q[0];
rz(-2.0728003) q[0];
rz(0.7069201) q[1];
sx q[1];
rz(-2.0126577) q[1];
sx q[1];
rz(-0.001002034) q[1];
rz(-1.3798643) q[2];
sx q[2];
rz(-1.1715062) q[2];
sx q[2];
rz(1.7572335) q[2];
rz(-2.9929326) q[3];
sx q[3];
rz(-1.7868663) q[3];
sx q[3];
rz(-0.79358906) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
