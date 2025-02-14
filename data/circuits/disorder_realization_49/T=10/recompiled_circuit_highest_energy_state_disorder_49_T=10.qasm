OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.385741204023361) q[0];
sx q[0];
rz(2.15856114228303) q[0];
sx q[0];
rz(9.98213628529712) q[0];
rz(-1.92692565917969) q[1];
sx q[1];
rz(7.13752237160737) q[1];
sx q[1];
rz(8.48272762297794) q[1];
cx q[1],q[0];
rz(2.1250114440918) q[0];
sx q[0];
rz(4.83102956612641) q[0];
sx q[0];
rz(12.3560290098111) q[0];
rz(0.932816922664642) q[2];
sx q[2];
rz(2.06174102623994) q[2];
sx q[2];
rz(9.94874701499149) q[2];
cx q[2],q[1];
rz(0.571919322013855) q[1];
sx q[1];
rz(1.50773564179475) q[1];
sx q[1];
rz(7.46772894858524) q[1];
rz(-0.837607860565186) q[3];
sx q[3];
rz(0.75467840035493) q[3];
sx q[3];
rz(11.1346715450208) q[3];
cx q[3],q[2];
rz(0.726309239864349) q[2];
sx q[2];
rz(4.90905896027619) q[2];
sx q[2];
rz(4.51258513926669) q[2];
rz(0.819094240665436) q[3];
sx q[3];
rz(2.29471251566941) q[3];
sx q[3];
rz(6.17292878626987) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.467219889163971) q[0];
sx q[0];
rz(1.90855661233003) q[0];
sx q[0];
rz(9.63065328299209) q[0];
rz(1.824183344841) q[1];
sx q[1];
rz(4.42339733441407) q[1];
sx q[1];
rz(15.2389921903531) q[1];
cx q[1],q[0];
rz(-0.436761111021042) q[0];
sx q[0];
rz(5.36738077004487) q[0];
sx q[0];
rz(6.25373313426181) q[0];
rz(3.55466842651367) q[2];
sx q[2];
rz(3.64294985135133) q[2];
sx q[2];
rz(8.50249329804584) q[2];
cx q[2],q[1];
rz(-2.93920207023621) q[1];
sx q[1];
rz(1.10310450394685) q[1];
sx q[1];
rz(13.1968738794248) q[1];
rz(-2.32483243942261) q[3];
sx q[3];
rz(5.15972772439057) q[3];
sx q[3];
rz(9.32559697925254) q[3];
cx q[3],q[2];
rz(2.5466902256012) q[2];
sx q[2];
rz(3.97420099576051) q[2];
sx q[2];
rz(8.7971948146741) q[2];
rz(1.07074308395386) q[3];
sx q[3];
rz(2.63797334034974) q[3];
sx q[3];
rz(12.8951883077542) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-4.54643583297729) q[0];
sx q[0];
rz(0.810318144159861) q[0];
sx q[0];
rz(10.2132623553197) q[0];
rz(-3.71271061897278) q[1];
sx q[1];
rz(1.32891372044618) q[1];
sx q[1];
rz(14.0123095273893) q[1];
cx q[1],q[0];
rz(2.33192372322083) q[0];
sx q[0];
rz(4.44371417363221) q[0];
sx q[0];
rz(14.1234631299894) q[0];
rz(4.17434787750244) q[2];
sx q[2];
rz(4.74554422696168) q[2];
sx q[2];
rz(13.3336632013242) q[2];
cx q[2],q[1];
rz(2.26479768753052) q[1];
sx q[1];
rz(4.81221965153749) q[1];
sx q[1];
rz(10.697099184982) q[1];
rz(1.05500447750092) q[3];
sx q[3];
rz(3.78030827839906) q[3];
sx q[3];
rz(10.5456787109296) q[3];
cx q[3],q[2];
rz(-1.04204046726227) q[2];
sx q[2];
rz(3.58400911291177) q[2];
sx q[2];
rz(9.04567000865146) q[2];
rz(1.37426173686981) q[3];
sx q[3];
rz(2.17024240096147) q[3];
sx q[3];
rz(8.4410387635152) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.207972601056099) q[0];
sx q[0];
rz(3.89638248284394) q[0];
sx q[0];
rz(7.19561669825717) q[0];
rz(4.53539228439331) q[1];
sx q[1];
rz(2.63639655907685) q[1];
sx q[1];
rz(6.58551261424228) q[1];
cx q[1],q[0];
rz(0.497198224067688) q[0];
sx q[0];
rz(3.28206023772294) q[0];
sx q[0];
rz(14.578684782974) q[0];
rz(1.9665492773056) q[2];
sx q[2];
rz(4.03757599194581) q[2];
sx q[2];
rz(3.99014184474155) q[2];
cx q[2],q[1];
rz(1.52795493602753) q[1];
sx q[1];
rz(9.63849607308442) q[1];
sx q[1];
rz(6.48848030566379) q[1];
rz(-0.693481683731079) q[3];
sx q[3];
rz(5.27354350884492) q[3];
sx q[3];
rz(8.31880948542758) q[3];
cx q[3],q[2];
rz(3.62590789794922) q[2];
sx q[2];
rz(3.84327551920945) q[2];
sx q[2];
rz(13.19575045108) q[2];
rz(1.7221919298172) q[3];
sx q[3];
rz(5.00200930436189) q[3];
sx q[3];
rz(10.1332329869191) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0929542705416679) q[0];
sx q[0];
rz(2.10436704953248) q[0];
sx q[0];
rz(10.7441268920819) q[0];
rz(2.05357789993286) q[1];
sx q[1];
rz(4.41337064106996) q[1];
sx q[1];
rz(7.96005103587314) q[1];
cx q[1],q[0];
rz(2.33254027366638) q[0];
sx q[0];
rz(2.27024242480332) q[0];
sx q[0];
rz(11.0375599622647) q[0];
rz(-0.0324499495327473) q[2];
sx q[2];
rz(2.20321491559083) q[2];
sx q[2];
rz(5.6411878824155) q[2];
cx q[2],q[1];
rz(-0.326809972524643) q[1];
sx q[1];
rz(4.37404993374879) q[1];
sx q[1];
rz(9.62571675180599) q[1];
rz(1.89280772209167) q[3];
sx q[3];
rz(2.51090410550172) q[3];
sx q[3];
rz(5.07277009486362) q[3];
cx q[3],q[2];
rz(-0.692708373069763) q[2];
sx q[2];
rz(0.118818910914012) q[2];
sx q[2];
rz(9.18785083889171) q[2];
rz(-4.1221227645874) q[3];
sx q[3];
rz(4.85389211972291) q[3];
sx q[3];
rz(7.22712895869418) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.15810814499855) q[0];
sx q[0];
rz(6.17369285424287) q[0];
sx q[0];
rz(9.27961545287772) q[0];
rz(1.62040865421295) q[1];
sx q[1];
rz(3.9415718634897) q[1];
sx q[1];
rz(9.41754384421884) q[1];
cx q[1],q[0];
rz(2.50993919372559) q[0];
sx q[0];
rz(4.47678199608857) q[0];
sx q[0];
rz(10.1345051884572) q[0];
rz(-3.98629117012024) q[2];
sx q[2];
rz(1.13960042794282) q[2];
sx q[2];
rz(10.2635951399724) q[2];
cx q[2],q[1];
rz(-1.59110939502716) q[1];
sx q[1];
rz(1.80460110505159) q[1];
sx q[1];
rz(10.8889590263288) q[1];
rz(-3.11298775672913) q[3];
sx q[3];
rz(4.54633835156495) q[3];
sx q[3];
rz(12.3744003534238) q[3];
cx q[3],q[2];
rz(-2.24685573577881) q[2];
sx q[2];
rz(7.47247615655) q[2];
sx q[2];
rz(8.84724858998462) q[2];
rz(1.15248107910156) q[3];
sx q[3];
rz(3.7265677173906) q[3];
sx q[3];
rz(8.01287434100314) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.56855273246765) q[0];
sx q[0];
rz(3.33724473615224) q[0];
sx q[0];
rz(7.34506604670688) q[0];
rz(1.38782393932343) q[1];
sx q[1];
rz(5.05271449883515) q[1];
sx q[1];
rz(10.9229409456174) q[1];
cx q[1],q[0];
rz(-3.95543575286865) q[0];
sx q[0];
rz(5.83812323411042) q[0];
sx q[0];
rz(8.1662202835004) q[0];
rz(1.14462804794312) q[2];
sx q[2];
rz(2.45056012471253) q[2];
sx q[2];
rz(9.86689046620532) q[2];
cx q[2],q[1];
rz(0.643872618675232) q[1];
sx q[1];
rz(1.79596975644166) q[1];
sx q[1];
rz(8.85560975073978) q[1];
rz(-2.05925536155701) q[3];
sx q[3];
rz(1.85912135441835) q[3];
sx q[3];
rz(11.847483611099) q[3];
cx q[3],q[2];
rz(1.74287033081055) q[2];
sx q[2];
rz(3.01694320340688) q[2];
sx q[2];
rz(11.6523442029874) q[2];
rz(2.38776087760925) q[3];
sx q[3];
rz(5.05152383645112) q[3];
sx q[3];
rz(8.96830258368655) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.13304257392883) q[0];
sx q[0];
rz(3.46678370435769) q[0];
sx q[0];
rz(9.4142112865965) q[0];
rz(-1.0039302110672) q[1];
sx q[1];
rz(4.55804017384584) q[1];
sx q[1];
rz(9.38579095750257) q[1];
cx q[1],q[0];
rz(-1.73321354389191) q[0];
sx q[0];
rz(2.77624023159082) q[0];
sx q[0];
rz(13.0747089147489) q[0];
rz(4.39558601379395) q[2];
sx q[2];
rz(4.70329454739625) q[2];
sx q[2];
rz(13.1938214063565) q[2];
cx q[2],q[1];
rz(0.403636395931244) q[1];
sx q[1];
rz(4.33412555058534) q[1];
sx q[1];
rz(11.7179853677671) q[1];
rz(-1.45639073848724) q[3];
sx q[3];
rz(8.75850787957246) q[3];
sx q[3];
rz(10.9270597457807) q[3];
cx q[3],q[2];
rz(-1.57993340492249) q[2];
sx q[2];
rz(4.37174478371675) q[2];
sx q[2];
rz(9.32523460536405) q[2];
rz(1.29334115982056) q[3];
sx q[3];
rz(4.99793222744996) q[3];
sx q[3];
rz(6.65874168872043) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.47159099578857) q[0];
sx q[0];
rz(4.91577854950959) q[0];
sx q[0];
rz(10.7369609832685) q[0];
rz(2.61477637290955) q[1];
sx q[1];
rz(5.52939644654328) q[1];
sx q[1];
rz(8.58801457881137) q[1];
cx q[1],q[0];
rz(3.14716291427612) q[0];
sx q[0];
rz(4.97244551976258) q[0];
sx q[0];
rz(9.12045822142764) q[0];
rz(1.0238801240921) q[2];
sx q[2];
rz(1.02670565445954) q[2];
sx q[2];
rz(10.0436507224958) q[2];
cx q[2],q[1];
rz(4.59722709655762) q[1];
sx q[1];
rz(5.33822551568086) q[1];
sx q[1];
rz(12.8910739183347) q[1];
rz(-1.189457654953) q[3];
sx q[3];
rz(2.01028886635835) q[3];
sx q[3];
rz(8.85568491219684) q[3];
cx q[3],q[2];
rz(-4.51224899291992) q[2];
sx q[2];
rz(6.89664188225801) q[2];
sx q[2];
rz(8.21754703520938) q[2];
rz(1.31974446773529) q[3];
sx q[3];
rz(4.71815958817536) q[3];
sx q[3];
rz(10.2210044026296) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.094548262655735) q[0];
sx q[0];
rz(2.66675090988214) q[0];
sx q[0];
rz(11.8497180700223) q[0];
rz(-1.56397497653961) q[1];
sx q[1];
rz(7.58690753777558) q[1];
sx q[1];
rz(11.9490232229154) q[1];
cx q[1],q[0];
rz(0.0940531492233276) q[0];
sx q[0];
rz(2.27290240128572) q[0];
sx q[0];
rz(8.07856533526584) q[0];
rz(-5.99361085891724) q[2];
sx q[2];
rz(0.485578449564525) q[2];
sx q[2];
rz(8.82046434878513) q[2];
cx q[2],q[1];
rz(0.523862838745117) q[1];
sx q[1];
rz(3.95551392634446) q[1];
sx q[1];
rz(11.6182439088742) q[1];
rz(2.50122928619385) q[3];
sx q[3];
rz(1.79882195790345) q[3];
sx q[3];
rz(11.2033177375714) q[3];
cx q[3],q[2];
rz(7.61274480819702) q[2];
sx q[2];
rz(2.75414076645906) q[2];
sx q[2];
rz(2.66239020823642) q[2];
rz(1.51743483543396) q[3];
sx q[3];
rz(4.54616359074647) q[3];
sx q[3];
rz(10.9242317438047) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.299782931804657) q[0];
sx q[0];
rz(2.12328925927217) q[0];
sx q[0];
rz(12.8211893796842) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-6.12543296813965) q[1];
sx q[1];
rz(3.20256163750822) q[1];
sx q[1];
rz(8.50348368882343) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(1.07523429393768) q[2];
sx q[2];
rz(7.64482990105683) q[2];
sx q[2];
rz(12.6833355188291) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(0.696326553821564) q[3];
sx q[3];
rz(1.06547919114167) q[3];
sx q[3];
rz(6.05331704615756) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
