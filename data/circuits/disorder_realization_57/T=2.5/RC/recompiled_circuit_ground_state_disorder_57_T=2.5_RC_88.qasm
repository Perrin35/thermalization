OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8785414) q[0];
sx q[0];
rz(-0.6113373) q[0];
sx q[0];
rz(1.4933458) q[0];
rz(-2.2892294) q[1];
sx q[1];
rz(-1.747793) q[1];
sx q[1];
rz(1.2300904) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48759197) q[0];
sx q[0];
rz(-2.0474259) q[0];
sx q[0];
rz(2.0311902) q[0];
rz(-pi) q[1];
rz(0.70582055) q[2];
sx q[2];
rz(-0.56669826) q[2];
sx q[2];
rz(-1.3659878) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.38065698) q[1];
sx q[1];
rz(-2.0208997) q[1];
sx q[1];
rz(-2.7633219) q[1];
rz(-1.2870077) q[3];
sx q[3];
rz(-1.5012245) q[3];
sx q[3];
rz(-2.8049935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40296951) q[2];
sx q[2];
rz(-1.7454742) q[2];
sx q[2];
rz(-1.1633066) q[2];
rz(2.2033384) q[3];
sx q[3];
rz(-2.3353751) q[3];
sx q[3];
rz(1.4279385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52629483) q[0];
sx q[0];
rz(-1.9087003) q[0];
sx q[0];
rz(1.1616608) q[0];
rz(-1.7192526) q[1];
sx q[1];
rz(-1.6484345) q[1];
sx q[1];
rz(-0.086440451) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42936642) q[0];
sx q[0];
rz(-1.4131568) q[0];
sx q[0];
rz(-1.0382354) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7941571) q[2];
sx q[2];
rz(-2.0513976) q[2];
sx q[2];
rz(-1.5011806) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5878488) q[1];
sx q[1];
rz(-1.5431197) q[1];
sx q[1];
rz(-1.5553655) q[1];
rz(-pi) q[2];
rz(-1.2041041) q[3];
sx q[3];
rz(-0.18309284) q[3];
sx q[3];
rz(-0.88475534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.57261673) q[2];
sx q[2];
rz(-0.47875753) q[2];
sx q[2];
rz(3.0968481) q[2];
rz(-0.71988002) q[3];
sx q[3];
rz(-1.6241112) q[3];
sx q[3];
rz(1.6911136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2429263) q[0];
sx q[0];
rz(-1.1761605) q[0];
sx q[0];
rz(-1.9157008) q[0];
rz(0.9179081) q[1];
sx q[1];
rz(-1.2914912) q[1];
sx q[1];
rz(1.8570522) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63375134) q[0];
sx q[0];
rz(-2.0720665) q[0];
sx q[0];
rz(-3.0847163) q[0];
rz(-0.62165181) q[2];
sx q[2];
rz(-1.5643411) q[2];
sx q[2];
rz(0.14989756) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1333993) q[1];
sx q[1];
rz(-0.8772521) q[1];
sx q[1];
rz(2.4213525) q[1];
rz(0.62326473) q[3];
sx q[3];
rz(-1.7812961) q[3];
sx q[3];
rz(2.7170469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7104177) q[2];
sx q[2];
rz(-1.147889) q[2];
sx q[2];
rz(0.61719027) q[2];
rz(2.2271473) q[3];
sx q[3];
rz(-0.72943288) q[3];
sx q[3];
rz(0.35703737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0964088) q[0];
sx q[0];
rz(-0.04627385) q[0];
sx q[0];
rz(-2.1531877) q[0];
rz(1.4811966) q[1];
sx q[1];
rz(-0.55440569) q[1];
sx q[1];
rz(2.5344417) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0074127277) q[0];
sx q[0];
rz(-1.3152342) q[0];
sx q[0];
rz(-1.8525396) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1902624) q[2];
sx q[2];
rz(-2.1759927) q[2];
sx q[2];
rz(0.91120341) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9937351) q[1];
sx q[1];
rz(-2.5566935) q[1];
sx q[1];
rz(2.2123076) q[1];
x q[2];
rz(1.7512583) q[3];
sx q[3];
rz(-1.4752098) q[3];
sx q[3];
rz(1.6570027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0854757) q[2];
sx q[2];
rz(-1.7530707) q[2];
sx q[2];
rz(-0.2307387) q[2];
rz(-0.76656109) q[3];
sx q[3];
rz(-2.9975588) q[3];
sx q[3];
rz(-1.2994331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4615999) q[0];
sx q[0];
rz(-1.5138641) q[0];
sx q[0];
rz(-0.067721279) q[0];
rz(-1.8374775) q[1];
sx q[1];
rz(-1.291052) q[1];
sx q[1];
rz(-1.7052604) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7505813) q[0];
sx q[0];
rz(-2.2428747) q[0];
sx q[0];
rz(2.6425022) q[0];
rz(-pi) q[1];
rz(-3.0313644) q[2];
sx q[2];
rz(-2.0617635) q[2];
sx q[2];
rz(1.5256745) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1597212) q[1];
sx q[1];
rz(-1.2471732) q[1];
sx q[1];
rz(-2.0702328) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.93238) q[3];
sx q[3];
rz(-1.3632332) q[3];
sx q[3];
rz(-2.521477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0962254) q[2];
sx q[2];
rz(-0.86288351) q[2];
sx q[2];
rz(-0.53675845) q[2];
rz(1.9219575) q[3];
sx q[3];
rz(-2.2370971) q[3];
sx q[3];
rz(1.9513244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9959975) q[0];
sx q[0];
rz(-2.1639316) q[0];
sx q[0];
rz(1.4033432) q[0];
rz(-2.6373236) q[1];
sx q[1];
rz(-1.206617) q[1];
sx q[1];
rz(-0.23955841) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0127691) q[0];
sx q[0];
rz(-1.7278624) q[0];
sx q[0];
rz(-1.3879018) q[0];
x q[1];
rz(-0.94953612) q[2];
sx q[2];
rz(-0.64709787) q[2];
sx q[2];
rz(0.71352772) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.319377) q[1];
sx q[1];
rz(-1.483157) q[1];
sx q[1];
rz(2.4725663) q[1];
rz(-2.2277545) q[3];
sx q[3];
rz(-0.60937928) q[3];
sx q[3];
rz(-2.1003124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2704894) q[2];
sx q[2];
rz(-1.5364545) q[2];
sx q[2];
rz(0.25320539) q[2];
rz(2.1829055) q[3];
sx q[3];
rz(-0.91512338) q[3];
sx q[3];
rz(1.8668713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7764928) q[0];
sx q[0];
rz(-0.72117844) q[0];
sx q[0];
rz(1.7505296) q[0];
rz(0.59025383) q[1];
sx q[1];
rz(-1.9275459) q[1];
sx q[1];
rz(-0.50277695) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3093291) q[0];
sx q[0];
rz(-1.676275) q[0];
sx q[0];
rz(-0.01450392) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46765695) q[2];
sx q[2];
rz(-2.4291647) q[2];
sx q[2];
rz(-1.4234655) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31820276) q[1];
sx q[1];
rz(-0.776459) q[1];
sx q[1];
rz(1.9253325) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8842949) q[3];
sx q[3];
rz(-2.0924221) q[3];
sx q[3];
rz(2.9056463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.34660029) q[2];
sx q[2];
rz(-1.184) q[2];
sx q[2];
rz(-2.6896175) q[2];
rz(0.31451264) q[3];
sx q[3];
rz(-1.9032685) q[3];
sx q[3];
rz(-3.0934635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(1.5101584) q[0];
sx q[0];
rz(-0.87961125) q[0];
sx q[0];
rz(-1.5055702) q[0];
rz(0.10558852) q[1];
sx q[1];
rz(-1.0786723) q[1];
sx q[1];
rz(-0.67965913) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9847947) q[0];
sx q[0];
rz(-2.121006) q[0];
sx q[0];
rz(2.2722428) q[0];
rz(-pi) q[1];
rz(-0.40186581) q[2];
sx q[2];
rz(-0.99152126) q[2];
sx q[2];
rz(1.4039672) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.41808203) q[1];
sx q[1];
rz(-0.46996221) q[1];
sx q[1];
rz(1.5638007) q[1];
x q[2];
rz(-0.8828396) q[3];
sx q[3];
rz(-0.2885836) q[3];
sx q[3];
rz(-1.2478873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.52123657) q[2];
sx q[2];
rz(-0.53400365) q[2];
sx q[2];
rz(0.83068577) q[2];
rz(1.0384809) q[3];
sx q[3];
rz(-0.65282789) q[3];
sx q[3];
rz(-1.4000819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1164383) q[0];
sx q[0];
rz(-1.9898299) q[0];
sx q[0];
rz(-2.108216) q[0];
rz(-2.0694464) q[1];
sx q[1];
rz(-1.1786345) q[1];
sx q[1];
rz(-0.56195608) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4406703) q[0];
sx q[0];
rz(-1.7096247) q[0];
sx q[0];
rz(-2.2825861) q[0];
rz(-pi) q[1];
rz(2.2267954) q[2];
sx q[2];
rz(-1.2848228) q[2];
sx q[2];
rz(-2.9867619) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1984673) q[1];
sx q[1];
rz(-2.2930305) q[1];
sx q[1];
rz(0.46015443) q[1];
x q[2];
rz(-1.5132684) q[3];
sx q[3];
rz(-1.5794601) q[3];
sx q[3];
rz(-1.5271562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0434887) q[2];
sx q[2];
rz(-1.6489112) q[2];
sx q[2];
rz(1.7334422) q[2];
rz(3.0812541) q[3];
sx q[3];
rz(-1.3121759) q[3];
sx q[3];
rz(-1.2618802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.9982346) q[0];
sx q[0];
rz(-0.77333212) q[0];
sx q[0];
rz(1.7711357) q[0];
rz(2.1342318) q[1];
sx q[1];
rz(-1.7528563) q[1];
sx q[1];
rz(-3.0000684) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62011224) q[0];
sx q[0];
rz(-0.65516383) q[0];
sx q[0];
rz(-0.76458365) q[0];
rz(-pi) q[1];
rz(-0.74574377) q[2];
sx q[2];
rz(-0.088938449) q[2];
sx q[2];
rz(2.3740785) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87301577) q[1];
sx q[1];
rz(-2.0076224) q[1];
sx q[1];
rz(-2.718513) q[1];
rz(-pi) q[2];
rz(-1.0753638) q[3];
sx q[3];
rz(-1.5513707) q[3];
sx q[3];
rz(2.7904449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0274028) q[2];
sx q[2];
rz(-1.9582615) q[2];
sx q[2];
rz(-1.2664504) q[2];
rz(-0.13628422) q[3];
sx q[3];
rz(-0.91629052) q[3];
sx q[3];
rz(-0.19792476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.65205735) q[0];
sx q[0];
rz(-1.0132402) q[0];
sx q[0];
rz(-1.5860438) q[0];
rz(-0.96792211) q[1];
sx q[1];
rz(-0.5236917) q[1];
sx q[1];
rz(2.0655469) q[1];
rz(0.10836149) q[2];
sx q[2];
rz(-2.32794) q[2];
sx q[2];
rz(1.9447504) q[2];
rz(1.3068253) q[3];
sx q[3];
rz(-0.19656678) q[3];
sx q[3];
rz(-1.8241775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
