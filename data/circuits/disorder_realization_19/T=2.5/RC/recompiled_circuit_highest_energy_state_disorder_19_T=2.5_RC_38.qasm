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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6237804) q[0];
sx q[0];
rz(-1.5470439) q[0];
sx q[0];
rz(-2.0022814) q[0];
rz(0.52402988) q[2];
sx q[2];
rz(-1.2175778) q[2];
sx q[2];
rz(-2.2729682) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2249445) q[1];
sx q[1];
rz(-1.6672581) q[1];
sx q[1];
rz(1.8889844) q[1];
x q[2];
rz(1.8901915) q[3];
sx q[3];
rz(-2.4301964) q[3];
sx q[3];
rz(2.9044202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.69433576) q[2];
sx q[2];
rz(-2.3981636) q[2];
sx q[2];
rz(-1.2747964) q[2];
rz(-2.6796807) q[3];
sx q[3];
rz(-0.67449823) q[3];
sx q[3];
rz(1.1326724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3502515) q[0];
sx q[0];
rz(-0.30838648) q[0];
sx q[0];
rz(-1.2530918) q[0];
rz(-0.14532267) q[1];
sx q[1];
rz(-1.7456313) q[1];
sx q[1];
rz(-2.0504418) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.017830124) q[0];
sx q[0];
rz(-1.995867) q[0];
sx q[0];
rz(2.9062889) q[0];
rz(-0.65204377) q[2];
sx q[2];
rz(-0.55527675) q[2];
sx q[2];
rz(-0.20737831) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1205463) q[1];
sx q[1];
rz(-1.7851549) q[1];
sx q[1];
rz(-1.3938741) q[1];
rz(-pi) q[2];
rz(1.1717779) q[3];
sx q[3];
rz(-1.1091091) q[3];
sx q[3];
rz(-2.6103013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.48065177) q[2];
sx q[2];
rz(-1.7512243) q[2];
sx q[2];
rz(0.52978984) q[2];
rz(0.79331136) q[3];
sx q[3];
rz(-1.6042234) q[3];
sx q[3];
rz(0.62354273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.6563501) q[0];
sx q[0];
rz(-2.1915477) q[0];
sx q[0];
rz(0.52870885) q[0];
rz(0.54620019) q[1];
sx q[1];
rz(-0.9587973) q[1];
sx q[1];
rz(2.8012457) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.074269144) q[0];
sx q[0];
rz(-0.43282291) q[0];
sx q[0];
rz(2.4433663) q[0];
x q[1];
rz(2.5916948) q[2];
sx q[2];
rz(-2.8680414) q[2];
sx q[2];
rz(0.036444681) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.18225154) q[1];
sx q[1];
rz(-1.6653582) q[1];
sx q[1];
rz(-1.5675117) q[1];
rz(-pi) q[2];
x q[2];
rz(2.765708) q[3];
sx q[3];
rz(-0.86556992) q[3];
sx q[3];
rz(-2.8994967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1722374) q[2];
sx q[2];
rz(-2.7293971) q[2];
sx q[2];
rz(-0.42984143) q[2];
rz(-0.75628453) q[3];
sx q[3];
rz(-0.1736621) q[3];
sx q[3];
rz(-0.82591301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38465685) q[0];
sx q[0];
rz(-1.5058368) q[0];
sx q[0];
rz(-1.6480308) q[0];
rz(0.76796302) q[1];
sx q[1];
rz(-2.6710644) q[1];
sx q[1];
rz(0.54642645) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1892197) q[0];
sx q[0];
rz(-1.5060987) q[0];
sx q[0];
rz(1.6420664) q[0];
rz(-0.77026639) q[2];
sx q[2];
rz(-2.4663743) q[2];
sx q[2];
rz(-0.16083052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.029303251) q[1];
sx q[1];
rz(-1.422907) q[1];
sx q[1];
rz(-1.6690955) q[1];
rz(-pi) q[2];
rz(0.68164556) q[3];
sx q[3];
rz(-1.036231) q[3];
sx q[3];
rz(0.93972423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4297318) q[2];
sx q[2];
rz(-1.4873361) q[2];
sx q[2];
rz(2.8374953) q[2];
rz(-1.3652623) q[3];
sx q[3];
rz(-1.920776) q[3];
sx q[3];
rz(-2.064866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6292608) q[0];
sx q[0];
rz(-2.6973695) q[0];
sx q[0];
rz(0.00057922676) q[0];
rz(-2.0594788) q[1];
sx q[1];
rz(-0.52434701) q[1];
sx q[1];
rz(-2.2023315) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9094641) q[0];
sx q[0];
rz(-2.3298023) q[0];
sx q[0];
rz(-0.62007298) q[0];
rz(0.95920697) q[2];
sx q[2];
rz(-1.2805174) q[2];
sx q[2];
rz(-0.91735754) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5737557) q[1];
sx q[1];
rz(-1.3569965) q[1];
sx q[1];
rz(2.1478199) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.22721283) q[3];
sx q[3];
rz(-1.2649396) q[3];
sx q[3];
rz(-1.3671041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1028563) q[2];
sx q[2];
rz(-1.2612217) q[2];
sx q[2];
rz(-2.8803414) q[2];
rz(0.86483613) q[3];
sx q[3];
rz(-1.7295001) q[3];
sx q[3];
rz(-2.849546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84457266) q[0];
sx q[0];
rz(-1.427587) q[0];
sx q[0];
rz(2.936506) q[0];
rz(-1.0467485) q[1];
sx q[1];
rz(-1.7457242) q[1];
sx q[1];
rz(0.37839016) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32961938) q[0];
sx q[0];
rz(-1.1866335) q[0];
sx q[0];
rz(2.9400918) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47176265) q[2];
sx q[2];
rz(-1.7182351) q[2];
sx q[2];
rz(0.23213895) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.69120317) q[1];
sx q[1];
rz(-0.65237633) q[1];
sx q[1];
rz(-2.1689586) q[1];
rz(0.50726733) q[3];
sx q[3];
rz(-0.75921042) q[3];
sx q[3];
rz(1.678148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4947074) q[2];
sx q[2];
rz(-1.9834221) q[2];
sx q[2];
rz(1.0542487) q[2];
rz(0.65822893) q[3];
sx q[3];
rz(-2.4963278) q[3];
sx q[3];
rz(-0.48404199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9996416) q[0];
sx q[0];
rz(-1.9331837) q[0];
sx q[0];
rz(0.64055881) q[0];
rz(-0.41796747) q[1];
sx q[1];
rz(-1.9211946) q[1];
sx q[1];
rz(0.80498615) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5875801) q[0];
sx q[0];
rz(-2.4338284) q[0];
sx q[0];
rz(-3.1057024) q[0];
rz(-pi) q[1];
rz(0.80193582) q[2];
sx q[2];
rz(-1.2705497) q[2];
sx q[2];
rz(1.173347) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0143483) q[1];
sx q[1];
rz(-1.5906723) q[1];
sx q[1];
rz(-1.5991421) q[1];
rz(-pi) q[2];
rz(1.8667614) q[3];
sx q[3];
rz(-0.99963847) q[3];
sx q[3];
rz(0.35863245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5927222) q[2];
sx q[2];
rz(-2.1465116) q[2];
sx q[2];
rz(-0.76118809) q[2];
rz(2.393764) q[3];
sx q[3];
rz(-1.8239832) q[3];
sx q[3];
rz(-0.67659155) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(0.87108535) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3251299) q[0];
sx q[0];
rz(-1.5583546) q[0];
sx q[0];
rz(0.61412707) q[0];
rz(-pi) q[1];
rz(0.34452166) q[2];
sx q[2];
rz(-1.7167175) q[2];
sx q[2];
rz(2.5308756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9242036) q[1];
sx q[1];
rz(-2.743209) q[1];
sx q[1];
rz(0.59415643) q[1];
rz(1.1444451) q[3];
sx q[3];
rz(-1.5160069) q[3];
sx q[3];
rz(0.74758119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.90290922) q[2];
sx q[2];
rz(-0.52672714) q[2];
sx q[2];
rz(-0.61863679) q[2];
rz(3.1213308) q[3];
sx q[3];
rz(-2.198115) q[3];
sx q[3];
rz(-0.77997911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0778377) q[0];
sx q[0];
rz(-1.5851333) q[0];
sx q[0];
rz(-2.7177287) q[0];
rz(-2.9546812) q[1];
sx q[1];
rz(-0.82004768) q[1];
sx q[1];
rz(-1.6835469) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8114421) q[0];
sx q[0];
rz(-1.6887661) q[0];
sx q[0];
rz(1.5850204) q[0];
rz(-pi) q[1];
rz(-1.0571613) q[2];
sx q[2];
rz(-1.8087401) q[2];
sx q[2];
rz(-0.4713716) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9002237) q[1];
sx q[1];
rz(-1.4321064) q[1];
sx q[1];
rz(1.6226416) q[1];
rz(-0.23641674) q[3];
sx q[3];
rz(-2.7763753) q[3];
sx q[3];
rz(-1.055298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3905048) q[2];
sx q[2];
rz(-1.0672528) q[2];
sx q[2];
rz(1.0941774) q[2];
rz(-1.68082) q[3];
sx q[3];
rz(-1.3760309) q[3];
sx q[3];
rz(1.8028397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8193034) q[0];
sx q[0];
rz(-1.3811454) q[0];
sx q[0];
rz(-1.5018916) q[0];
rz(0.27885258) q[1];
sx q[1];
rz(-0.96062213) q[1];
sx q[1];
rz(2.1174812) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0251966) q[0];
sx q[0];
rz(-2.0587344) q[0];
sx q[0];
rz(-1.6750209) q[0];
x q[1];
rz(3.0908998) q[2];
sx q[2];
rz(-1.2963352) q[2];
sx q[2];
rz(1.3433742) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1361724) q[1];
sx q[1];
rz(-1.3498303) q[1];
sx q[1];
rz(1.7762089) q[1];
rz(2.315178) q[3];
sx q[3];
rz(-1.3152988) q[3];
sx q[3];
rz(-2.5662553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9378822) q[2];
sx q[2];
rz(-1.2869765) q[2];
sx q[2];
rz(-2.6046275) q[2];
rz(1.2282061) q[3];
sx q[3];
rz(-0.95913404) q[3];
sx q[3];
rz(-2.8502407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0252329) q[0];
sx q[0];
rz(-1.7740842) q[0];
sx q[0];
rz(-2.0728003) q[0];
rz(-2.4346726) q[1];
sx q[1];
rz(-2.0126577) q[1];
sx q[1];
rz(-0.001002034) q[1];
rz(2.7189485) q[2];
sx q[2];
rz(-0.44036897) q[2];
sx q[2];
rz(-1.8457495) q[2];
rz(1.3523921) q[3];
sx q[3];
rz(-1.7159749) q[3];
sx q[3];
rz(0.80930474) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
