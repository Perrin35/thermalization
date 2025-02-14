OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6405606) q[0];
sx q[0];
rz(-0.80067331) q[0];
sx q[0];
rz(2.9963357) q[0];
rz(2.2486806) q[1];
sx q[1];
rz(-0.4141663) q[1];
sx q[1];
rz(-1.1069586) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1205638) q[0];
sx q[0];
rz(-1.0146838) q[0];
sx q[0];
rz(-2.4377941) q[0];
rz(-pi) q[1];
rz(-2.8584004) q[2];
sx q[2];
rz(-1.8841266) q[2];
sx q[2];
rz(-1.6746132) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9913276) q[1];
sx q[1];
rz(-1.2854165) q[1];
sx q[1];
rz(2.1726514) q[1];
x q[2];
rz(-1.1916394) q[3];
sx q[3];
rz(-1.772545) q[3];
sx q[3];
rz(0.74046046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5408111) q[2];
sx q[2];
rz(-1.7822219) q[2];
sx q[2];
rz(0.092279807) q[2];
rz(-1.7021092) q[3];
sx q[3];
rz(-1.1441792) q[3];
sx q[3];
rz(-2.2030742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0094902078) q[0];
sx q[0];
rz(-1.9785896) q[0];
sx q[0];
rz(-2.5376885) q[0];
rz(1.6101135) q[1];
sx q[1];
rz(-1.8096626) q[1];
sx q[1];
rz(1.5276705) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.912589) q[0];
sx q[0];
rz(-2.4618076) q[0];
sx q[0];
rz(1.0566105) q[0];
rz(2.1555734) q[2];
sx q[2];
rz(-2.2246317) q[2];
sx q[2];
rz(-0.80160917) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.81948796) q[1];
sx q[1];
rz(-2.8675962) q[1];
sx q[1];
rz(0.042983965) q[1];
rz(3.134733) q[3];
sx q[3];
rz(-0.09819542) q[3];
sx q[3];
rz(-0.47658893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.70832002) q[2];
sx q[2];
rz(-1.4488139) q[2];
sx q[2];
rz(-2.0937008) q[2];
rz(-2.4723054) q[3];
sx q[3];
rz(-0.71988121) q[3];
sx q[3];
rz(-1.1650813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.79353756) q[0];
sx q[0];
rz(-2.3568643) q[0];
sx q[0];
rz(-1.0070356) q[0];
rz(-2.8496565) q[1];
sx q[1];
rz(-2.597229) q[1];
sx q[1];
rz(-0.67063355) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6601283) q[0];
sx q[0];
rz(-2.1892836) q[0];
sx q[0];
rz(-3.0342558) q[0];
rz(-pi) q[1];
rz(0.79659454) q[2];
sx q[2];
rz(-2.0907986) q[2];
sx q[2];
rz(2.7869625) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5340928) q[1];
sx q[1];
rz(-2.1779446) q[1];
sx q[1];
rz(1.1578039) q[1];
x q[2];
rz(1.8995883) q[3];
sx q[3];
rz(-2.7955849) q[3];
sx q[3];
rz(-0.11869988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2900419) q[2];
sx q[2];
rz(-0.85323492) q[2];
sx q[2];
rz(-2.2500989) q[2];
rz(1.1024891) q[3];
sx q[3];
rz(-2.1474371) q[3];
sx q[3];
rz(1.7360784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
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
rz(0.82146984) q[0];
sx q[0];
rz(-1.3865043) q[0];
sx q[0];
rz(1.0572877) q[0];
rz(-0.0013466324) q[1];
sx q[1];
rz(-0.82095447) q[1];
sx q[1];
rz(2.3473306) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47424537) q[0];
sx q[0];
rz(-0.66215958) q[0];
sx q[0];
rz(-0.029900877) q[0];
rz(-pi) q[1];
rz(-0.35425953) q[2];
sx q[2];
rz(-1.2984167) q[2];
sx q[2];
rz(-2.543022) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22143198) q[1];
sx q[1];
rz(-2.2756612) q[1];
sx q[1];
rz(-0.95734289) q[1];
rz(-pi) q[2];
rz(0.36093386) q[3];
sx q[3];
rz(-0.31692255) q[3];
sx q[3];
rz(-1.288687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1069676) q[2];
sx q[2];
rz(-0.62109533) q[2];
sx q[2];
rz(2.1194439) q[2];
rz(1.243783) q[3];
sx q[3];
rz(-2.1918178) q[3];
sx q[3];
rz(-1.3589842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46955243) q[0];
sx q[0];
rz(-1.7615027) q[0];
sx q[0];
rz(-2.688038) q[0];
rz(-2.1039311) q[1];
sx q[1];
rz(-1.1353759) q[1];
sx q[1];
rz(-0.79197788) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2642518) q[0];
sx q[0];
rz(-1.8487572) q[0];
sx q[0];
rz(0.045601337) q[0];
rz(-2.7466082) q[2];
sx q[2];
rz(-1.028966) q[2];
sx q[2];
rz(-1.3196368) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9336613) q[1];
sx q[1];
rz(-1.9092036) q[1];
sx q[1];
rz(-1.4470571) q[1];
rz(-pi) q[2];
rz(-3.1230917) q[3];
sx q[3];
rz(-1.7089239) q[3];
sx q[3];
rz(0.02397315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8301293) q[2];
sx q[2];
rz(-0.66581231) q[2];
sx q[2];
rz(-2.8361481) q[2];
rz(2.9597802) q[3];
sx q[3];
rz(-1.5287377) q[3];
sx q[3];
rz(-2.4055433) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5859454) q[0];
sx q[0];
rz(-0.82252994) q[0];
sx q[0];
rz(0.12538759) q[0];
rz(-1.5665945) q[1];
sx q[1];
rz(-1.6927203) q[1];
sx q[1];
rz(-3.1094508) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9956995) q[0];
sx q[0];
rz(-0.75135485) q[0];
sx q[0];
rz(0.32352792) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3022167) q[2];
sx q[2];
rz(-1.617086) q[2];
sx q[2];
rz(2.3037825) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7508538) q[1];
sx q[1];
rz(-1.097569) q[1];
sx q[1];
rz(-0.79428947) q[1];
rz(-pi) q[2];
x q[2];
rz(3.016196) q[3];
sx q[3];
rz(-0.74955696) q[3];
sx q[3];
rz(0.22334418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1273758) q[2];
sx q[2];
rz(-1.2140423) q[2];
sx q[2];
rz(1.1227013) q[2];
rz(-3.0681916) q[3];
sx q[3];
rz(-0.95728907) q[3];
sx q[3];
rz(0.61298031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4323394) q[0];
sx q[0];
rz(-1.4839577) q[0];
sx q[0];
rz(0.51026979) q[0];
rz(0.36422745) q[1];
sx q[1];
rz(-0.42114708) q[1];
sx q[1];
rz(1.6417004) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0696568) q[0];
sx q[0];
rz(-0.47635117) q[0];
sx q[0];
rz(-1.8029965) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0232863) q[2];
sx q[2];
rz(-1.7276754) q[2];
sx q[2];
rz(-1.9594693) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31452306) q[1];
sx q[1];
rz(-2.0644958) q[1];
sx q[1];
rz(-1.3160454) q[1];
x q[2];
rz(2.4046201) q[3];
sx q[3];
rz(-1.342257) q[3];
sx q[3];
rz(-1.4481973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4437272) q[2];
sx q[2];
rz(-1.5865734) q[2];
sx q[2];
rz(-0.86722428) q[2];
rz(-1.5813658) q[3];
sx q[3];
rz(-0.19872228) q[3];
sx q[3];
rz(-1.8038512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7241868) q[0];
sx q[0];
rz(-0.066411821) q[0];
sx q[0];
rz(0.41626406) q[0];
rz(1.1626214) q[1];
sx q[1];
rz(-1.9527718) q[1];
sx q[1];
rz(0.75688854) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.074919393) q[0];
sx q[0];
rz(-2.6180589) q[0];
sx q[0];
rz(-2.174211) q[0];
rz(2.1483118) q[2];
sx q[2];
rz(-1.9406291) q[2];
sx q[2];
rz(0.12060697) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6357506) q[1];
sx q[1];
rz(-1.381784) q[1];
sx q[1];
rz(2.0057949) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59764909) q[3];
sx q[3];
rz(-0.72859287) q[3];
sx q[3];
rz(0.33580175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68652144) q[2];
sx q[2];
rz(-1.7073809) q[2];
sx q[2];
rz(-1.7835468) q[2];
rz(0.81651917) q[3];
sx q[3];
rz(-1.9101382) q[3];
sx q[3];
rz(-2.9787298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74129504) q[0];
sx q[0];
rz(-1.6930027) q[0];
sx q[0];
rz(-1.4073538) q[0];
rz(-0.74527144) q[1];
sx q[1];
rz(-2.407357) q[1];
sx q[1];
rz(-1.5819246) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2360412) q[0];
sx q[0];
rz(-2.2187382) q[0];
sx q[0];
rz(0.34867649) q[0];
x q[1];
rz(0.64777957) q[2];
sx q[2];
rz(-1.5967007) q[2];
sx q[2];
rz(1.4454973) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.386454) q[1];
sx q[1];
rz(-1.6999131) q[1];
sx q[1];
rz(1.1243724) q[1];
rz(0.060272597) q[3];
sx q[3];
rz(-0.99592402) q[3];
sx q[3];
rz(2.4337492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.72426307) q[2];
sx q[2];
rz(-1.686692) q[2];
sx q[2];
rz(-2.0261185) q[2];
rz(2.8880902) q[3];
sx q[3];
rz(-1.2616254) q[3];
sx q[3];
rz(3.1326262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1215006) q[0];
sx q[0];
rz(-1.2637063) q[0];
sx q[0];
rz(2.3790835) q[0];
rz(-2.1992042) q[1];
sx q[1];
rz(-2.0768879) q[1];
sx q[1];
rz(1.9047838) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4651637) q[0];
sx q[0];
rz(-1.0753514) q[0];
sx q[0];
rz(-1.5009319) q[0];
rz(-0.89434172) q[2];
sx q[2];
rz(-0.80781762) q[2];
sx q[2];
rz(-1.7730912) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2357465) q[1];
sx q[1];
rz(-1.5457834) q[1];
sx q[1];
rz(-0.025084875) q[1];
rz(-pi) q[2];
rz(1.9754161) q[3];
sx q[3];
rz(-1.5050833) q[3];
sx q[3];
rz(0.031377553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2463871) q[2];
sx q[2];
rz(-3.0463986) q[2];
sx q[2];
rz(3.134356) q[2];
rz(2.9712408) q[3];
sx q[3];
rz(-1.6536313) q[3];
sx q[3];
rz(0.65792221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8728747) q[0];
sx q[0];
rz(-1.3997411) q[0];
sx q[0];
rz(1.1135143) q[0];
rz(0.089182236) q[1];
sx q[1];
rz(-1.5166278) q[1];
sx q[1];
rz(1.2022432) q[1];
rz(-0.15684814) q[2];
sx q[2];
rz(-1.7050171) q[2];
sx q[2];
rz(-0.27324054) q[2];
rz(1.3842267) q[3];
sx q[3];
rz(-1.866478) q[3];
sx q[3];
rz(2.5685892) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
