OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.50103203) q[0];
sx q[0];
rz(-2.3409193) q[0];
sx q[0];
rz(-2.9963357) q[0];
rz(2.2486806) q[1];
sx q[1];
rz(-0.4141663) q[1];
sx q[1];
rz(-1.1069586) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.021028886) q[0];
sx q[0];
rz(-1.0146838) q[0];
sx q[0];
rz(-0.70379852) q[0];
rz(-pi) q[1];
rz(-0.28319226) q[2];
sx q[2];
rz(-1.8841266) q[2];
sx q[2];
rz(1.6746132) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9913276) q[1];
sx q[1];
rz(-1.8561761) q[1];
sx q[1];
rz(-0.96894126) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.065973) q[3];
sx q[3];
rz(-2.7144066) q[3];
sx q[3];
rz(-1.2963201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5408111) q[2];
sx q[2];
rz(-1.3593707) q[2];
sx q[2];
rz(-3.0493128) q[2];
rz(1.4394834) q[3];
sx q[3];
rz(-1.9974134) q[3];
sx q[3];
rz(-0.93851844) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0094902078) q[0];
sx q[0];
rz(-1.9785896) q[0];
sx q[0];
rz(-2.5376885) q[0];
rz(-1.5314792) q[1];
sx q[1];
rz(-1.8096626) q[1];
sx q[1];
rz(-1.6139222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92787111) q[0];
sx q[0];
rz(-1.8851213) q[0];
sx q[0];
rz(0.95753352) q[0];
x q[1];
rz(2.3983922) q[2];
sx q[2];
rz(-1.1172406) q[2];
sx q[2];
rz(-1.1519866) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.81948796) q[1];
sx q[1];
rz(-0.27399644) q[1];
sx q[1];
rz(-0.042983965) q[1];
rz(3.134733) q[3];
sx q[3];
rz(-3.0433972) q[3];
sx q[3];
rz(-2.6650037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4332726) q[2];
sx q[2];
rz(-1.4488139) q[2];
sx q[2];
rz(1.0478919) q[2];
rz(-0.6692872) q[3];
sx q[3];
rz(-0.71988121) q[3];
sx q[3];
rz(1.1650813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79353756) q[0];
sx q[0];
rz(-2.3568643) q[0];
sx q[0];
rz(2.134557) q[0];
rz(-0.29193613) q[1];
sx q[1];
rz(-2.597229) q[1];
sx q[1];
rz(-2.4709591) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9898674) q[0];
sx q[0];
rz(-1.6581931) q[0];
sx q[0];
rz(2.1920127) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.466356) q[2];
sx q[2];
rz(-2.2227382) q[2];
sx q[2];
rz(-2.3777131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2082767) q[1];
sx q[1];
rz(-1.9066992) q[1];
sx q[1];
rz(-2.4927054) q[1];
rz(-pi) q[2];
rz(1.8995883) q[3];
sx q[3];
rz(-0.34600779) q[3];
sx q[3];
rz(-3.0228928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85155073) q[2];
sx q[2];
rz(-0.85323492) q[2];
sx q[2];
rz(-2.2500989) q[2];
rz(2.0391035) q[3];
sx q[3];
rz(-2.1474371) q[3];
sx q[3];
rz(1.4055143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82146984) q[0];
sx q[0];
rz(-1.3865043) q[0];
sx q[0];
rz(-2.084305) q[0];
rz(3.140246) q[1];
sx q[1];
rz(-0.82095447) q[1];
sx q[1];
rz(2.3473306) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7052536) q[0];
sx q[0];
rz(-2.2326075) q[0];
sx q[0];
rz(1.5474942) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2813454) q[2];
sx q[2];
rz(-1.9114541) q[2];
sx q[2];
rz(-0.87305005) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2202178) q[1];
sx q[1];
rz(-2.0247321) q[1];
sx q[1];
rz(0.80516071) q[1];
rz(-1.6861071) q[3];
sx q[3];
rz(-1.274935) q[3];
sx q[3];
rz(1.4747696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1069676) q[2];
sx q[2];
rz(-2.5204973) q[2];
sx q[2];
rz(-2.1194439) q[2];
rz(-1.8978097) q[3];
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
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6720402) q[0];
sx q[0];
rz(-1.7615027) q[0];
sx q[0];
rz(2.688038) q[0];
rz(-1.0376616) q[1];
sx q[1];
rz(-1.1353759) q[1];
sx q[1];
rz(-2.3496148) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29402367) q[0];
sx q[0];
rz(-1.6146462) q[0];
sx q[0];
rz(-1.8490318) q[0];
x q[1];
rz(2.7466082) q[2];
sx q[2];
rz(-2.1126267) q[2];
sx q[2];
rz(1.8219558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9336613) q[1];
sx q[1];
rz(-1.232389) q[1];
sx q[1];
rz(-1.6945356) q[1];
x q[2];
rz(-1.7030992) q[3];
sx q[3];
rz(-3.0022394) q[3];
sx q[3];
rz(-2.9840368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3114634) q[2];
sx q[2];
rz(-0.66581231) q[2];
sx q[2];
rz(2.8361481) q[2];
rz(2.9597802) q[3];
sx q[3];
rz(-1.5287377) q[3];
sx q[3];
rz(-2.4055433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5859454) q[0];
sx q[0];
rz(-0.82252994) q[0];
sx q[0];
rz(-0.12538759) q[0];
rz(1.5665945) q[1];
sx q[1];
rz(-1.4488723) q[1];
sx q[1];
rz(-3.1094508) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9956995) q[0];
sx q[0];
rz(-0.75135485) q[0];
sx q[0];
rz(-0.32352792) q[0];
rz(-pi) q[1];
rz(-2.3022167) q[2];
sx q[2];
rz(-1.5245066) q[2];
sx q[2];
rz(-0.83781017) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.39073889) q[1];
sx q[1];
rz(-1.097569) q[1];
sx q[1];
rz(0.79428947) q[1];
x q[2];
rz(-0.12539668) q[3];
sx q[3];
rz(-0.74955696) q[3];
sx q[3];
rz(-2.9182485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0142168) q[2];
sx q[2];
rz(-1.2140423) q[2];
sx q[2];
rz(-2.0188913) q[2];
rz(-3.0681916) q[3];
sx q[3];
rz(-0.95728907) q[3];
sx q[3];
rz(-2.5286123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4323394) q[0];
sx q[0];
rz(-1.657635) q[0];
sx q[0];
rz(2.6313229) q[0];
rz(-2.7773652) q[1];
sx q[1];
rz(-2.7204456) q[1];
sx q[1];
rz(1.4998923) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0696568) q[0];
sx q[0];
rz(-2.6652415) q[0];
sx q[0];
rz(-1.3385962) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11830637) q[2];
sx q[2];
rz(-1.7276754) q[2];
sx q[2];
rz(1.1821234) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.31452306) q[1];
sx q[1];
rz(-1.0770969) q[1];
sx q[1];
rz(1.8255472) q[1];
rz(-pi) q[2];
rz(1.2664421) q[3];
sx q[3];
rz(-0.8571763) q[3];
sx q[3];
rz(-0.32538381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4437272) q[2];
sx q[2];
rz(-1.5865734) q[2];
sx q[2];
rz(2.2743684) q[2];
rz(1.5813658) q[3];
sx q[3];
rz(-0.19872228) q[3];
sx q[3];
rz(-1.3377415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7241868) q[0];
sx q[0];
rz(-3.0751808) q[0];
sx q[0];
rz(0.41626406) q[0];
rz(-1.1626214) q[1];
sx q[1];
rz(-1.9527718) q[1];
sx q[1];
rz(2.3847041) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3945321) q[0];
sx q[0];
rz(-1.9950657) q[0];
sx q[0];
rz(-0.31655689) q[0];
rz(-pi) q[1];
rz(0.95332884) q[2];
sx q[2];
rz(-2.4673415) q[2];
sx q[2];
rz(1.9566388) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6923995) q[1];
sx q[1];
rz(-0.47187065) q[1];
sx q[1];
rz(1.9969246) q[1];
rz(0.63558319) q[3];
sx q[3];
rz(-1.9548237) q[3];
sx q[3];
rz(-1.4366729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4550712) q[2];
sx q[2];
rz(-1.7073809) q[2];
sx q[2];
rz(1.7835468) q[2];
rz(-0.81651917) q[3];
sx q[3];
rz(-1.9101382) q[3];
sx q[3];
rz(-0.16286287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.73423568) q[1];
sx q[1];
rz(-1.5596681) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0228221) q[0];
sx q[0];
rz(-1.2948991) q[0];
sx q[0];
rz(-2.2488382) q[0];
rz(-0.042912622) q[2];
sx q[2];
rz(-0.6482228) q[2];
sx q[2];
rz(2.9820778) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8956977) q[1];
sx q[1];
rz(-2.0132397) q[1];
sx q[1];
rz(-2.9986283) q[1];
rz(0.060272597) q[3];
sx q[3];
rz(-2.1456686) q[3];
sx q[3];
rz(0.70784345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4173296) q[2];
sx q[2];
rz(-1.4549007) q[2];
sx q[2];
rz(1.1154741) q[2];
rz(-2.8880902) q[3];
sx q[3];
rz(-1.8799672) q[3];
sx q[3];
rz(3.1326262) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1215006) q[0];
sx q[0];
rz(-1.2637063) q[0];
sx q[0];
rz(-2.3790835) q[0];
rz(2.1992042) q[1];
sx q[1];
rz(-2.0768879) q[1];
sx q[1];
rz(-1.9047838) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6764289) q[0];
sx q[0];
rz(-1.0753514) q[0];
sx q[0];
rz(-1.5009319) q[0];
rz(-pi) q[1];
rz(-2.2472509) q[2];
sx q[2];
rz(-2.333775) q[2];
sx q[2];
rz(1.3685014) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0227388) q[1];
sx q[1];
rz(-3.10617) q[1];
sx q[1];
rz(2.3574748) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0701261) q[3];
sx q[3];
rz(-1.1671007) q[3];
sx q[3];
rz(-1.57406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8952055) q[2];
sx q[2];
rz(-3.0463986) q[2];
sx q[2];
rz(3.134356) q[2];
rz(-2.9712408) q[3];
sx q[3];
rz(-1.6536313) q[3];
sx q[3];
rz(-0.65792221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
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
rz(2.4287672) q[2];
sx q[2];
rz(-2.935514) q[2];
sx q[2];
rz(1.9998177) q[2];
rz(-2.5946272) q[3];
sx q[3];
rz(-2.7934358) q[3];
sx q[3];
rz(-1.1478333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
