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
rz(1.9563142) q[1];
sx q[1];
rz(4.5524608) q[1];
sx q[1];
rz(8.3571385) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51781228) q[0];
sx q[0];
rz(-1.5470439) q[0];
sx q[0];
rz(-1.1393113) q[0];
rz(2.6175628) q[2];
sx q[2];
rz(-1.9240148) q[2];
sx q[2];
rz(0.86862446) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8274535) q[1];
sx q[1];
rz(-1.2541391) q[1];
sx q[1];
rz(3.0400671) q[1];
rz(2.8772763) q[3];
sx q[3];
rz(-0.90221635) q[3];
sx q[3];
rz(2.9670935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4472569) q[2];
sx q[2];
rz(-0.7434291) q[2];
sx q[2];
rz(1.2747964) q[2];
rz(-0.46191195) q[3];
sx q[3];
rz(-2.4670944) q[3];
sx q[3];
rz(1.1326724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79134113) q[0];
sx q[0];
rz(-2.8332062) q[0];
sx q[0];
rz(1.2530918) q[0];
rz(-2.99627) q[1];
sx q[1];
rz(-1.3959613) q[1];
sx q[1];
rz(1.0911509) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6328207) q[0];
sx q[0];
rz(-0.48235407) q[0];
sx q[0];
rz(-1.095196) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45812313) q[2];
sx q[2];
rz(-1.2451742) q[2];
sx q[2];
rz(1.9389012) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4636865) q[1];
sx q[1];
rz(-0.27707252) q[1];
sx q[1];
rz(2.4617247) q[1];
rz(-pi) q[2];
rz(-1.9698148) q[3];
sx q[3];
rz(-1.1091091) q[3];
sx q[3];
rz(-2.6103013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6609409) q[2];
sx q[2];
rz(-1.7512243) q[2];
sx q[2];
rz(0.52978984) q[2];
rz(-2.3482813) q[3];
sx q[3];
rz(-1.6042234) q[3];
sx q[3];
rz(0.62354273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6563501) q[0];
sx q[0];
rz(-0.95004496) q[0];
sx q[0];
rz(2.6128838) q[0];
rz(-2.5953925) q[1];
sx q[1];
rz(-2.1827953) q[1];
sx q[1];
rz(0.34034696) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2961194) q[0];
sx q[0];
rz(-1.8438135) q[0];
sx q[0];
rz(-2.8014328) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9067847) q[2];
sx q[2];
rz(-1.4291414) q[2];
sx q[2];
rz(-2.0675142) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9941147) q[1];
sx q[1];
rz(-3.0469739) q[1];
sx q[1];
rz(-3.1069744) q[1];
x q[2];
rz(1.163655) q[3];
sx q[3];
rz(-2.357886) q[3];
sx q[3];
rz(-2.3526085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9693552) q[2];
sx q[2];
rz(-2.7293971) q[2];
sx q[2];
rz(2.7117512) q[2];
rz(2.3853081) q[3];
sx q[3];
rz(-0.1736621) q[3];
sx q[3];
rz(-0.82591301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38465685) q[0];
sx q[0];
rz(-1.6357559) q[0];
sx q[0];
rz(1.4935619) q[0];
rz(2.3736296) q[1];
sx q[1];
rz(-0.47052828) q[1];
sx q[1];
rz(0.54642645) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38619216) q[0];
sx q[0];
rz(-1.4996756) q[0];
sx q[0];
rz(-0.064861809) q[0];
rz(1.0621319) q[2];
sx q[2];
rz(-2.0360247) q[2];
sx q[2];
rz(-2.0874799) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.56049743) q[1];
sx q[1];
rz(-2.9642134) q[1];
sx q[1];
rz(-2.5590798) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.68164556) q[3];
sx q[3];
rz(-2.1053616) q[3];
sx q[3];
rz(-2.2018684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4297318) q[2];
sx q[2];
rz(-1.4873361) q[2];
sx q[2];
rz(2.8374953) q[2];
rz(1.7763304) q[3];
sx q[3];
rz(-1.2208166) q[3];
sx q[3];
rz(-1.0767267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6292608) q[0];
sx q[0];
rz(-2.6973695) q[0];
sx q[0];
rz(0.00057922676) q[0];
rz(1.0821139) q[1];
sx q[1];
rz(-2.6172456) q[1];
sx q[1];
rz(-0.93926114) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88193306) q[0];
sx q[0];
rz(-1.1355917) q[0];
sx q[0];
rz(-0.70910221) q[0];
rz(0.95920697) q[2];
sx q[2];
rz(-1.8610753) q[2];
sx q[2];
rz(0.91735754) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4537332) q[1];
sx q[1];
rz(-0.61111585) q[1];
sx q[1];
rz(1.9495717) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1904508) q[3];
sx q[3];
rz(-2.7627146) q[3];
sx q[3];
rz(2.4293116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1028563) q[2];
sx q[2];
rz(-1.2612217) q[2];
sx q[2];
rz(-0.2612513) q[2];
rz(0.86483613) q[3];
sx q[3];
rz(-1.7295001) q[3];
sx q[3];
rz(-2.849546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.29702) q[0];
sx q[0];
rz(-1.7140056) q[0];
sx q[0];
rz(-0.20508668) q[0];
rz(-2.0948441) q[1];
sx q[1];
rz(-1.7457242) q[1];
sx q[1];
rz(2.7632025) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8119733) q[0];
sx q[0];
rz(-1.1866335) q[0];
sx q[0];
rz(-2.9400918) q[0];
rz(-pi) q[1];
rz(1.7360052) q[2];
sx q[2];
rz(-1.1045611) q[2];
sx q[2];
rz(1.8777443) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4000611) q[1];
sx q[1];
rz(-2.0963256) q[1];
sx q[1];
rz(-0.40627131) q[1];
rz(-pi) q[2];
rz(-2.6343253) q[3];
sx q[3];
rz(-0.75921042) q[3];
sx q[3];
rz(-1.4634446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.64688524) q[2];
sx q[2];
rz(-1.1581706) q[2];
sx q[2];
rz(2.0873439) q[2];
rz(-2.4833637) q[3];
sx q[3];
rz(-2.4963278) q[3];
sx q[3];
rz(2.6575507) q[3];
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
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9996416) q[0];
sx q[0];
rz(-1.208409) q[0];
sx q[0];
rz(-2.5010338) q[0];
rz(0.41796747) q[1];
sx q[1];
rz(-1.9211946) q[1];
sx q[1];
rz(-0.80498615) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5403596) q[0];
sx q[0];
rz(-2.2780097) q[0];
sx q[0];
rz(1.6014895) q[0];
x q[1];
rz(-1.1518794) q[2];
sx q[2];
rz(-0.81406128) q[2];
sx q[2];
rz(-2.4474622) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1678214) q[1];
sx q[1];
rz(-0.034618363) q[1];
sx q[1];
rz(0.95914118) q[1];
rz(-pi) q[2];
x q[2];
rz(0.42607362) q[3];
sx q[3];
rz(-2.5059627) q[3];
sx q[3];
rz(-0.15492188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54887041) q[2];
sx q[2];
rz(-0.9950811) q[2];
sx q[2];
rz(0.76118809) q[2];
rz(2.393764) q[3];
sx q[3];
rz(-1.8239832) q[3];
sx q[3];
rz(2.4650011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51650301) q[0];
sx q[0];
rz(-2.857132) q[0];
sx q[0];
rz(-1.4768584) q[0];
rz(2.6853216) q[1];
sx q[1];
rz(-1.4197333) q[1];
sx q[1];
rz(2.2705073) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74556158) q[0];
sx q[0];
rz(-2.1848688) q[0];
sx q[0];
rz(-1.5860193) q[0];
x q[1];
rz(-2.7311677) q[2];
sx q[2];
rz(-2.7685809) q[2];
sx q[2];
rz(-2.5665064) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.217389) q[1];
sx q[1];
rz(-2.743209) q[1];
sx q[1];
rz(-2.5474362) q[1];
rz(1.9971476) q[3];
sx q[3];
rz(-1.6255857) q[3];
sx q[3];
rz(-2.3940115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90290922) q[2];
sx q[2];
rz(-0.52672714) q[2];
sx q[2];
rz(2.5229559) q[2];
rz(-0.020261852) q[3];
sx q[3];
rz(-0.94347763) q[3];
sx q[3];
rz(0.77997911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063754931) q[0];
sx q[0];
rz(-1.5851333) q[0];
sx q[0];
rz(2.7177287) q[0];
rz(-0.18691143) q[1];
sx q[1];
rz(-2.321545) q[1];
sx q[1];
rz(-1.6835469) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8114421) q[0];
sx q[0];
rz(-1.6887661) q[0];
sx q[0];
rz(1.5565722) q[0];
rz(2.8700022) q[2];
sx q[2];
rz(-2.0686065) q[2];
sx q[2];
rz(1.23162) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5411387) q[1];
sx q[1];
rz(-0.14800528) q[1];
sx q[1];
rz(2.7861094) q[1];
rz(2.9051759) q[3];
sx q[3];
rz(-2.7763753) q[3];
sx q[3];
rz(2.0862947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.75108782) q[2];
sx q[2];
rz(-1.0672528) q[2];
sx q[2];
rz(-2.0474153) q[2];
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
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3222892) q[0];
sx q[0];
rz(-1.3811454) q[0];
sx q[0];
rz(-1.639701) q[0];
rz(-2.8627401) q[1];
sx q[1];
rz(-2.1809705) q[1];
sx q[1];
rz(1.0241114) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0251966) q[0];
sx q[0];
rz(-2.0587344) q[0];
sx q[0];
rz(1.4665718) q[0];
rz(-1.845593) q[2];
sx q[2];
rz(-1.5220023) q[2];
sx q[2];
rz(-2.9279207) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1361724) q[1];
sx q[1];
rz(-1.3498303) q[1];
sx q[1];
rz(1.3653838) q[1];
x q[2];
rz(1.9387705) q[3];
sx q[3];
rz(-0.77897969) q[3];
sx q[3];
rz(-0.72768962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2037105) q[2];
sx q[2];
rz(-1.8546162) q[2];
sx q[2];
rz(2.6046275) q[2];
rz(-1.9133866) q[3];
sx q[3];
rz(-2.1824586) q[3];
sx q[3];
rz(-0.29135191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0252329) q[0];
sx q[0];
rz(-1.7740842) q[0];
sx q[0];
rz(-2.0728003) q[0];
rz(2.4346726) q[1];
sx q[1];
rz(-1.128935) q[1];
sx q[1];
rz(3.1405906) q[1];
rz(-2.7189485) q[2];
sx q[2];
rz(-2.7012237) q[2];
sx q[2];
rz(1.2958432) q[2];
rz(1.7892006) q[3];
sx q[3];
rz(-1.4256178) q[3];
sx q[3];
rz(-2.3322879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
