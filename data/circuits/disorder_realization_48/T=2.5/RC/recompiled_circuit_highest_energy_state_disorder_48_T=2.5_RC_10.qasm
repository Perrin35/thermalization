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
rz(2.7895522) q[0];
sx q[0];
rz(-2.3166603) q[0];
sx q[0];
rz(2.6068249) q[0];
rz(0.98400247) q[1];
sx q[1];
rz(-2.6360631) q[1];
sx q[1];
rz(-1.8796138) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.004744) q[0];
sx q[0];
rz(-1.8333577) q[0];
sx q[0];
rz(-0.62341086) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1970799) q[2];
sx q[2];
rz(-2.0227832) q[2];
sx q[2];
rz(0.96679692) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.76319219) q[1];
sx q[1];
rz(-2.1175368) q[1];
sx q[1];
rz(-1.631447) q[1];
rz(-2.7522911) q[3];
sx q[3];
rz(-1.7598146) q[3];
sx q[3];
rz(-2.1019328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8925605) q[2];
sx q[2];
rz(-0.51237115) q[2];
sx q[2];
rz(-1.6585635) q[2];
rz(2.856971) q[3];
sx q[3];
rz(-0.62323815) q[3];
sx q[3];
rz(2.0774138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2738709) q[0];
sx q[0];
rz(-0.72672788) q[0];
sx q[0];
rz(0.04100767) q[0];
rz(1.2076591) q[1];
sx q[1];
rz(-2.9137847) q[1];
sx q[1];
rz(1.6809195) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2940833) q[0];
sx q[0];
rz(-0.14862157) q[0];
sx q[0];
rz(-0.20047184) q[0];
rz(-pi) q[1];
rz(-2.5778077) q[2];
sx q[2];
rz(-2.2905458) q[2];
sx q[2];
rz(1.9518746) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.41439287) q[1];
sx q[1];
rz(-1.314394) q[1];
sx q[1];
rz(-0.039915737) q[1];
rz(-pi) q[2];
rz(-1.8051265) q[3];
sx q[3];
rz(-0.70584345) q[3];
sx q[3];
rz(0.88283759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4760806) q[2];
sx q[2];
rz(-1.6802695) q[2];
sx q[2];
rz(1.3903138) q[2];
rz(0.0811854) q[3];
sx q[3];
rz(-0.66892162) q[3];
sx q[3];
rz(-1.4359052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7471033) q[0];
sx q[0];
rz(-3.1293226) q[0];
sx q[0];
rz(2.5315206) q[0];
rz(1.3129781) q[1];
sx q[1];
rz(-0.6956296) q[1];
sx q[1];
rz(2.5909766) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0673772) q[0];
sx q[0];
rz(-2.2226637) q[0];
sx q[0];
rz(-1.1135191) q[0];
x q[1];
rz(-2.8189893) q[2];
sx q[2];
rz(-1.3122953) q[2];
sx q[2];
rz(-0.58835627) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0455835) q[1];
sx q[1];
rz(-2.7852045) q[1];
sx q[1];
rz(2.136904) q[1];
rz(-pi) q[2];
rz(-2.9987644) q[3];
sx q[3];
rz(-0.83072829) q[3];
sx q[3];
rz(2.9252441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6637471) q[2];
sx q[2];
rz(-2.9691594) q[2];
sx q[2];
rz(-2.0752068) q[2];
rz(-1.9829228) q[3];
sx q[3];
rz(-1.7585157) q[3];
sx q[3];
rz(-0.042044736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31024194) q[0];
sx q[0];
rz(-0.61315918) q[0];
sx q[0];
rz(-0.9170652) q[0];
rz(-0.65545583) q[1];
sx q[1];
rz(-2.2323445) q[1];
sx q[1];
rz(-2.7520032) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661267) q[0];
sx q[0];
rz(-1.4891169) q[0];
sx q[0];
rz(-1.7577175) q[0];
rz(-pi) q[1];
rz(0.55928237) q[2];
sx q[2];
rz(-0.96508316) q[2];
sx q[2];
rz(0.16597151) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37871088) q[1];
sx q[1];
rz(-0.40232752) q[1];
sx q[1];
rz(-1.4961924) q[1];
rz(-0.97352685) q[3];
sx q[3];
rz(-2.0359267) q[3];
sx q[3];
rz(-1.2965508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3098844) q[2];
sx q[2];
rz(-2.1990621) q[2];
sx q[2];
rz(1.3523098) q[2];
rz(-2.3934707) q[3];
sx q[3];
rz(-1.3040521) q[3];
sx q[3];
rz(-1.5266533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30369502) q[0];
sx q[0];
rz(-0.56981531) q[0];
sx q[0];
rz(-1.8001528) q[0];
rz(-2.1603284) q[1];
sx q[1];
rz(-1.860268) q[1];
sx q[1];
rz(0.70995465) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34559743) q[0];
sx q[0];
rz(-1.5516722) q[0];
sx q[0];
rz(2.8157793) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1469028) q[2];
sx q[2];
rz(-1.7719977) q[2];
sx q[2];
rz(-2.3735769) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4932351) q[1];
sx q[1];
rz(-1.7410189) q[1];
sx q[1];
rz(-1.7864208) q[1];
x q[2];
rz(0.59779915) q[3];
sx q[3];
rz(-1.6493268) q[3];
sx q[3];
rz(0.30818916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.95825163) q[2];
sx q[2];
rz(-0.7603344) q[2];
sx q[2];
rz(0.90671268) q[2];
rz(0.18236154) q[3];
sx q[3];
rz(-1.7014818) q[3];
sx q[3];
rz(-2.2831634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7957423) q[0];
sx q[0];
rz(-2.1868571) q[0];
sx q[0];
rz(2.9398651) q[0];
rz(-1.3564823) q[1];
sx q[1];
rz(-2.4701665) q[1];
sx q[1];
rz(-1.1511525) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042554458) q[0];
sx q[0];
rz(-1.7793613) q[0];
sx q[0];
rz(-2.0046356) q[0];
x q[1];
rz(-2.5348659) q[2];
sx q[2];
rz(-1.9579685) q[2];
sx q[2];
rz(2.3558056) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5073302) q[1];
sx q[1];
rz(-1.1686022) q[1];
sx q[1];
rz(-2.5172154) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7986951) q[3];
sx q[3];
rz(-1.0712475) q[3];
sx q[3];
rz(1.7021021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1543033) q[2];
sx q[2];
rz(-0.58649784) q[2];
sx q[2];
rz(2.7941373) q[2];
rz(2.3197428) q[3];
sx q[3];
rz(-1.776266) q[3];
sx q[3];
rz(1.6177026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3407985) q[0];
sx q[0];
rz(-0.64685416) q[0];
sx q[0];
rz(-2.0183753) q[0];
rz(2.7128291) q[1];
sx q[1];
rz(-2.2518497) q[1];
sx q[1];
rz(-0.14895983) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2620747) q[0];
sx q[0];
rz(-1.5929758) q[0];
sx q[0];
rz(2.0878077) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2960213) q[2];
sx q[2];
rz(-1.3775577) q[2];
sx q[2];
rz(3.0770965) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.80854446) q[1];
sx q[1];
rz(-0.77189779) q[1];
sx q[1];
rz(-0.21265642) q[1];
rz(-pi) q[2];
rz(-0.034986939) q[3];
sx q[3];
rz(-1.9763499) q[3];
sx q[3];
rz(1.7008894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1296156) q[2];
sx q[2];
rz(-0.35795438) q[2];
sx q[2];
rz(-0.32148662) q[2];
rz(-2.0251515) q[3];
sx q[3];
rz(-2.2827086) q[3];
sx q[3];
rz(-1.727625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9880992) q[0];
sx q[0];
rz(-1.7970002) q[0];
sx q[0];
rz(-0.022911428) q[0];
rz(1.5494391) q[1];
sx q[1];
rz(-2.0982845) q[1];
sx q[1];
rz(-1.2124088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7717944) q[0];
sx q[0];
rz(-1.3098469) q[0];
sx q[0];
rz(2.1290995) q[0];
rz(-pi) q[1];
rz(-2.3018478) q[2];
sx q[2];
rz(-1.296412) q[2];
sx q[2];
rz(1.7626761) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3500994) q[1];
sx q[1];
rz(-1.1731358) q[1];
sx q[1];
rz(1.2418163) q[1];
x q[2];
rz(-1.4597662) q[3];
sx q[3];
rz(-0.71655203) q[3];
sx q[3];
rz(0.46886231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2108078) q[2];
sx q[2];
rz(-2.942694) q[2];
sx q[2];
rz(-0.99736324) q[2];
rz(1.8831683) q[3];
sx q[3];
rz(-1.1910028) q[3];
sx q[3];
rz(-1.2303801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1579943) q[0];
sx q[0];
rz(-2.4857434) q[0];
sx q[0];
rz(-0.47437814) q[0];
rz(-0.55167088) q[1];
sx q[1];
rz(-2.3730979) q[1];
sx q[1];
rz(-0.65753585) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0829959) q[0];
sx q[0];
rz(-1.5471282) q[0];
sx q[0];
rz(-0.039796782) q[0];
x q[1];
rz(-2.9716597) q[2];
sx q[2];
rz(-1.9714154) q[2];
sx q[2];
rz(-1.1378433) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8213013) q[1];
sx q[1];
rz(-1.5870915) q[1];
sx q[1];
rz(1.9391796) q[1];
rz(-pi) q[2];
rz(2.1980638) q[3];
sx q[3];
rz(-0.17788685) q[3];
sx q[3];
rz(1.8454144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.23003301) q[2];
sx q[2];
rz(-0.40731373) q[2];
sx q[2];
rz(1.3113021) q[2];
rz(2.4274872) q[3];
sx q[3];
rz(-1.2823391) q[3];
sx q[3];
rz(0.56984058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0062200935) q[0];
sx q[0];
rz(-2.1923809) q[0];
sx q[0];
rz(2.5332992) q[0];
rz(2.3237806) q[1];
sx q[1];
rz(-1.1759956) q[1];
sx q[1];
rz(-0.83293319) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3522763) q[0];
sx q[0];
rz(-2.7171405) q[0];
sx q[0];
rz(-0.12571521) q[0];
x q[1];
rz(2.8274546) q[2];
sx q[2];
rz(-1.8153662) q[2];
sx q[2];
rz(-2.7121674) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6842839) q[1];
sx q[1];
rz(-0.59694081) q[1];
sx q[1];
rz(-2.9039911) q[1];
rz(-pi) q[2];
rz(1.6021483) q[3];
sx q[3];
rz(-2.1636181) q[3];
sx q[3];
rz(-1.5167936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1828764) q[2];
sx q[2];
rz(-2.7887838) q[2];
sx q[2];
rz(-0.58615169) q[2];
rz(-2.2385249) q[3];
sx q[3];
rz(-2.51913) q[3];
sx q[3];
rz(-3.0557475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9027973) q[0];
sx q[0];
rz(-1.1302523) q[0];
sx q[0];
rz(-0.94373066) q[0];
rz(1.0974274) q[1];
sx q[1];
rz(-0.25249093) q[1];
sx q[1];
rz(2.9846334) q[1];
rz(-0.13931724) q[2];
sx q[2];
rz(-1.2323772) q[2];
sx q[2];
rz(2.7215794) q[2];
rz(-0.87370609) q[3];
sx q[3];
rz(-1.1813963) q[3];
sx q[3];
rz(1.2753857) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
