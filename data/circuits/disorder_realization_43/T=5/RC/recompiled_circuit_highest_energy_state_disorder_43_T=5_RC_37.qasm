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
rz(2.9146258) q[0];
sx q[0];
rz(-0.6114971) q[0];
sx q[0];
rz(-2.5831232) q[0];
rz(0.59977579) q[1];
sx q[1];
rz(1.37473) q[1];
sx q[1];
rz(9.3461499) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1100548) q[0];
sx q[0];
rz(-2.1221836) q[0];
sx q[0];
rz(-2.4999775) q[0];
rz(-pi) q[1];
rz(1.6958649) q[2];
sx q[2];
rz(-0.89587051) q[2];
sx q[2];
rz(-1.1488016) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.47657644) q[1];
sx q[1];
rz(-1.4890097) q[1];
sx q[1];
rz(-0.25340124) q[1];
rz(-pi) q[2];
rz(-1.8923678) q[3];
sx q[3];
rz(-1.9059637) q[3];
sx q[3];
rz(0.17091076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2416396) q[2];
sx q[2];
rz(-0.048154801) q[2];
sx q[2];
rz(-2.107035) q[2];
rz(0.11040802) q[3];
sx q[3];
rz(-1.6610828) q[3];
sx q[3];
rz(2.835527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13650525) q[0];
sx q[0];
rz(-1.6965447) q[0];
sx q[0];
rz(-1.9553631) q[0];
rz(1.2677445) q[1];
sx q[1];
rz(-2.6823951) q[1];
sx q[1];
rz(1.7523821) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1339552) q[0];
sx q[0];
rz(-1.8764155) q[0];
sx q[0];
rz(2.1307178) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.21135881) q[2];
sx q[2];
rz(-2.6508001) q[2];
sx q[2];
rz(-2.9826814) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1139377) q[1];
sx q[1];
rz(-1.9009556) q[1];
sx q[1];
rz(0.3056674) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0458353) q[3];
sx q[3];
rz(-2.1928722) q[3];
sx q[3];
rz(-0.84703895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.17239751) q[2];
sx q[2];
rz(-2.0266504) q[2];
sx q[2];
rz(1.7051075) q[2];
rz(-1.8819594) q[3];
sx q[3];
rz(-2.6862222) q[3];
sx q[3];
rz(3.1148541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7773975) q[0];
sx q[0];
rz(-1.7561678) q[0];
sx q[0];
rz(1.9245603) q[0];
rz(-2.9265535) q[1];
sx q[1];
rz(-1.2478849) q[1];
sx q[1];
rz(1.7020114) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3453044) q[0];
sx q[0];
rz(-0.6718967) q[0];
sx q[0];
rz(-0.20033269) q[0];
x q[1];
rz(-2.0820175) q[2];
sx q[2];
rz(-1.4892092) q[2];
sx q[2];
rz(0.64543426) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9653757) q[1];
sx q[1];
rz(-1.6976408) q[1];
sx q[1];
rz(-1.3398257) q[1];
rz(1.2974627) q[3];
sx q[3];
rz(-1.36051) q[3];
sx q[3];
rz(-3.1137636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.8000468) q[2];
sx q[2];
rz(-1.6365106) q[2];
sx q[2];
rz(-0.81614196) q[2];
rz(-3.1331565) q[3];
sx q[3];
rz(-2.4237027) q[3];
sx q[3];
rz(1.5661904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7208045) q[0];
sx q[0];
rz(-1.1710465) q[0];
sx q[0];
rz(-0.26914832) q[0];
rz(-1.3828269) q[1];
sx q[1];
rz(-0.56126422) q[1];
sx q[1];
rz(2.6699578) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29971805) q[0];
sx q[0];
rz(-0.73900676) q[0];
sx q[0];
rz(-1.3625381) q[0];
rz(-3.0582171) q[2];
sx q[2];
rz(-0.40503392) q[2];
sx q[2];
rz(-1.1277367) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.1822575) q[1];
sx q[1];
rz(-1.253009) q[1];
sx q[1];
rz(0.55303581) q[1];
rz(-2.7876623) q[3];
sx q[3];
rz(-2.2056863) q[3];
sx q[3];
rz(2.3707182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65583324) q[2];
sx q[2];
rz(-1.7344319) q[2];
sx q[2];
rz(-1.887623) q[2];
rz(-2.8433825) q[3];
sx q[3];
rz(-1.2757653) q[3];
sx q[3];
rz(-1.6037174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0373847) q[0];
sx q[0];
rz(-2.2315114) q[0];
sx q[0];
rz(2.3846159) q[0];
rz(1.0589927) q[1];
sx q[1];
rz(-1.6964361) q[1];
sx q[1];
rz(2.7991378) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9152074) q[0];
sx q[0];
rz(-1.3801551) q[0];
sx q[0];
rz(-1.1676428) q[0];
rz(-pi) q[1];
rz(-2.8169698) q[2];
sx q[2];
rz(-2.2461736) q[2];
sx q[2];
rz(-1.7243232) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.40552822) q[1];
sx q[1];
rz(-1.2031415) q[1];
sx q[1];
rz(-0.50761948) q[1];
rz(2.5352372) q[3];
sx q[3];
rz(-0.49288756) q[3];
sx q[3];
rz(-0.78745251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.69514889) q[2];
sx q[2];
rz(-1.6764287) q[2];
sx q[2];
rz(-2.5214419) q[2];
rz(-2.5946963) q[3];
sx q[3];
rz(-2.9345025) q[3];
sx q[3];
rz(1.0774405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3243489) q[0];
sx q[0];
rz(-2.9381848) q[0];
sx q[0];
rz(1.9203583) q[0];
rz(-1.1034032) q[1];
sx q[1];
rz(-1.9788479) q[1];
sx q[1];
rz(-0.81072909) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9098783) q[0];
sx q[0];
rz(-0.30261573) q[0];
sx q[0];
rz(-2.581654) q[0];
rz(-pi) q[1];
rz(-2.3757908) q[2];
sx q[2];
rz(-0.81379902) q[2];
sx q[2];
rz(1.9794996) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.4524881) q[1];
sx q[1];
rz(-0.14912046) q[1];
sx q[1];
rz(-1.945687) q[1];
rz(-1.9493413) q[3];
sx q[3];
rz(-1.8486128) q[3];
sx q[3];
rz(1.5072106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.80060426) q[2];
sx q[2];
rz(-1.9438513) q[2];
sx q[2];
rz(-1.1654589) q[2];
rz(-0.10945877) q[3];
sx q[3];
rz(-1.0215003) q[3];
sx q[3];
rz(-0.060733184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3133746) q[0];
sx q[0];
rz(-0.010638588) q[0];
sx q[0];
rz(-2.4995372) q[0];
rz(0.24318801) q[1];
sx q[1];
rz(-0.41047341) q[1];
sx q[1];
rz(-0.64116716) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.083602) q[0];
sx q[0];
rz(-0.80785368) q[0];
sx q[0];
rz(0.74961787) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9320609) q[2];
sx q[2];
rz(-0.85339499) q[2];
sx q[2];
rz(3.1255869) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4879681) q[1];
sx q[1];
rz(-1.89978) q[1];
sx q[1];
rz(2.3423561) q[1];
rz(-2.9092714) q[3];
sx q[3];
rz(-1.9594176) q[3];
sx q[3];
rz(0.53750932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0082561) q[2];
sx q[2];
rz(-1.1587605) q[2];
sx q[2];
rz(0.94375098) q[2];
rz(-0.64822316) q[3];
sx q[3];
rz(-2.4002878) q[3];
sx q[3];
rz(-0.6664204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.42565313) q[0];
sx q[0];
rz(-2.3305927) q[0];
sx q[0];
rz(-2.7124523) q[0];
rz(-0.40402135) q[1];
sx q[1];
rz(-0.41530135) q[1];
sx q[1];
rz(2.2228352) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5024273) q[0];
sx q[0];
rz(-0.60898655) q[0];
sx q[0];
rz(-0.43242411) q[0];
x q[1];
rz(1.8317675) q[2];
sx q[2];
rz(-1.4497641) q[2];
sx q[2];
rz(0.82838917) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.017581878) q[1];
sx q[1];
rz(-0.38330829) q[1];
sx q[1];
rz(-0.95335754) q[1];
rz(1.8848583) q[3];
sx q[3];
rz(-1.5005642) q[3];
sx q[3];
rz(-2.8633871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1227485) q[2];
sx q[2];
rz(-0.65905535) q[2];
sx q[2];
rz(1.3502236) q[2];
rz(-0.33251479) q[3];
sx q[3];
rz(-0.88034383) q[3];
sx q[3];
rz(-1.3688709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3933082) q[0];
sx q[0];
rz(-2.4916593) q[0];
sx q[0];
rz(0.45676029) q[0];
rz(-1.7204174) q[1];
sx q[1];
rz(-1.2799809) q[1];
sx q[1];
rz(-3.0866887) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0676014) q[0];
sx q[0];
rz(-1.5419863) q[0];
sx q[0];
rz(2.3037698) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0012399) q[2];
sx q[2];
rz(-1.9483742) q[2];
sx q[2];
rz(-2.0098639) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0662288) q[1];
sx q[1];
rz(-1.3540596) q[1];
sx q[1];
rz(1.073887) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9076806) q[3];
sx q[3];
rz(-2.0816605) q[3];
sx q[3];
rz(0.0072016933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.88825893) q[2];
sx q[2];
rz(-0.71037018) q[2];
sx q[2];
rz(2.9448275) q[2];
rz(1.0737859) q[3];
sx q[3];
rz(-1.1353227) q[3];
sx q[3];
rz(-0.73608583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9842904) q[0];
sx q[0];
rz(-0.8534011) q[0];
sx q[0];
rz(0.89944696) q[0];
rz(-2.8304214) q[1];
sx q[1];
rz(-1.2093465) q[1];
sx q[1];
rz(-1.7412294) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3547826) q[0];
sx q[0];
rz(-1.3616832) q[0];
sx q[0];
rz(-1.0037006) q[0];
x q[1];
rz(0.24263361) q[2];
sx q[2];
rz(-1.0972692) q[2];
sx q[2];
rz(0.99019105) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8099762) q[1];
sx q[1];
rz(-2.4016233) q[1];
sx q[1];
rz(2.5472067) q[1];
x q[2];
rz(-2.1651044) q[3];
sx q[3];
rz(-2.1955238) q[3];
sx q[3];
rz(-0.67052602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.060999) q[2];
sx q[2];
rz(-2.521535) q[2];
sx q[2];
rz(1.9798123) q[2];
rz(-1.0787841) q[3];
sx q[3];
rz(-0.99083841) q[3];
sx q[3];
rz(-2.259528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5494736) q[0];
sx q[0];
rz(-1.619512) q[0];
sx q[0];
rz(-1.6721538) q[0];
rz(-0.1334162) q[1];
sx q[1];
rz(-1.3795556) q[1];
sx q[1];
rz(2.1605927) q[1];
rz(3.1405814) q[2];
sx q[2];
rz(-1.6242003) q[2];
sx q[2];
rz(1.5846197) q[2];
rz(-0.12984325) q[3];
sx q[3];
rz(-2.6363439) q[3];
sx q[3];
rz(1.5324203) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
