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
rz(2.7338602) q[0];
sx q[0];
rz(-1.2091577) q[0];
sx q[0];
rz(-1.6785167) q[0];
rz(0.2585803) q[1];
sx q[1];
rz(-1.9039896) q[1];
sx q[1];
rz(3.0078476) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8834849) q[0];
sx q[0];
rz(-2.9843669) q[0];
sx q[0];
rz(1.5327461) q[0];
rz(2.1315246) q[2];
sx q[2];
rz(-1.0655996) q[2];
sx q[2];
rz(1.2743735) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3773777) q[1];
sx q[1];
rz(-1.0395607) q[1];
sx q[1];
rz(-1.2197391) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7236563) q[3];
sx q[3];
rz(-1.5709849) q[3];
sx q[3];
rz(2.3589695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3818843) q[2];
sx q[2];
rz(-0.81177652) q[2];
sx q[2];
rz(-1.831057) q[2];
rz(1.7740907) q[3];
sx q[3];
rz(-2.01367) q[3];
sx q[3];
rz(-0.016157063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97717184) q[0];
sx q[0];
rz(-1.1559957) q[0];
sx q[0];
rz(-2.1261862) q[0];
rz(-0.39335355) q[1];
sx q[1];
rz(-1.4721847) q[1];
sx q[1];
rz(-0.70158395) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.836228) q[0];
sx q[0];
rz(-0.71963632) q[0];
sx q[0];
rz(-1.4426065) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5364473) q[2];
sx q[2];
rz(-0.30314244) q[2];
sx q[2];
rz(0.042138635) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9522515) q[1];
sx q[1];
rz(-2.8165112) q[1];
sx q[1];
rz(1.8081244) q[1];
rz(-pi) q[2];
rz(0.75922482) q[3];
sx q[3];
rz(-2.1442778) q[3];
sx q[3];
rz(0.92347565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1141438) q[2];
sx q[2];
rz(-2.1178718) q[2];
sx q[2];
rz(-1.1055498) q[2];
rz(-2.8607199) q[3];
sx q[3];
rz(-2.5307145) q[3];
sx q[3];
rz(-0.15225473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1059145) q[0];
sx q[0];
rz(-2.4690101) q[0];
sx q[0];
rz(1.0726844) q[0];
rz(-0.013280344) q[1];
sx q[1];
rz(-1.5417136) q[1];
sx q[1];
rz(-0.57659155) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.428513) q[0];
sx q[0];
rz(-1.2831076) q[0];
sx q[0];
rz(-1.6384533) q[0];
rz(-pi) q[1];
rz(1.4378087) q[2];
sx q[2];
rz(-1.4879359) q[2];
sx q[2];
rz(1.8216009) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2943337) q[1];
sx q[1];
rz(-1.0953961) q[1];
sx q[1];
rz(-0.21902276) q[1];
rz(-1.0131888) q[3];
sx q[3];
rz(-1.7806781) q[3];
sx q[3];
rz(-0.55279532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2836634) q[2];
sx q[2];
rz(-2.2495705) q[2];
sx q[2];
rz(-2.499495) q[2];
rz(-0.56504956) q[3];
sx q[3];
rz(-2.8665906) q[3];
sx q[3];
rz(-2.9716085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86244407) q[0];
sx q[0];
rz(-1.2553517) q[0];
sx q[0];
rz(-0.24236648) q[0];
rz(0.53388059) q[1];
sx q[1];
rz(-0.68478525) q[1];
sx q[1];
rz(1.6530316) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80099312) q[0];
sx q[0];
rz(-1.6417608) q[0];
sx q[0];
rz(-2.5453955) q[0];
rz(-1.8901509) q[2];
sx q[2];
rz(-1.9803626) q[2];
sx q[2];
rz(1.7920272) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0341314) q[1];
sx q[1];
rz(-2.3337769) q[1];
sx q[1];
rz(-0.40523578) q[1];
x q[2];
rz(-3.0128128) q[3];
sx q[3];
rz(-0.30211651) q[3];
sx q[3];
rz(-1.0766034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.36294595) q[2];
sx q[2];
rz(-1.2641509) q[2];
sx q[2];
rz(-1.6710336) q[2];
rz(0.10739022) q[3];
sx q[3];
rz(-1.8348179) q[3];
sx q[3];
rz(1.2764527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11112467) q[0];
sx q[0];
rz(-2.680439) q[0];
sx q[0];
rz(-1.9386468) q[0];
rz(0.89490926) q[1];
sx q[1];
rz(-1.9590961) q[1];
sx q[1];
rz(2.9295909) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0184191) q[0];
sx q[0];
rz(-0.97476649) q[0];
sx q[0];
rz(-0.53431781) q[0];
x q[1];
rz(1.8708399) q[2];
sx q[2];
rz(-1.0677538) q[2];
sx q[2];
rz(-1.6479657) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.13276671) q[1];
sx q[1];
rz(-1.0100875) q[1];
sx q[1];
rz(-2.6067023) q[1];
x q[2];
rz(-0.0096036951) q[3];
sx q[3];
rz(-2.1240747) q[3];
sx q[3];
rz(-1.6320848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8733069) q[2];
sx q[2];
rz(-2.3640859) q[2];
sx q[2];
rz(3.0262465) q[2];
rz(2.1558732) q[3];
sx q[3];
rz(-1.7306354) q[3];
sx q[3];
rz(2.7058097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91144052) q[0];
sx q[0];
rz(-1.4608915) q[0];
sx q[0];
rz(0.69754115) q[0];
rz(-0.41360924) q[1];
sx q[1];
rz(-1.4613084) q[1];
sx q[1];
rz(2.4381309) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089982219) q[0];
sx q[0];
rz(-1.4282266) q[0];
sx q[0];
rz(-0.27326126) q[0];
rz(-1.037961) q[2];
sx q[2];
rz(-0.67801266) q[2];
sx q[2];
rz(-2.968766) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78402987) q[1];
sx q[1];
rz(-2.5270487) q[1];
sx q[1];
rz(-0.67484208) q[1];
rz(-pi) q[2];
x q[2];
rz(2.260842) q[3];
sx q[3];
rz(-0.38057571) q[3];
sx q[3];
rz(0.56169034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0387663) q[2];
sx q[2];
rz(-0.26831728) q[2];
sx q[2];
rz(-1.9290257) q[2];
rz(-2.8596527) q[3];
sx q[3];
rz(-1.5423256) q[3];
sx q[3];
rz(-2.6350002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7956227) q[0];
sx q[0];
rz(-2.3891734) q[0];
sx q[0];
rz(0.84683007) q[0];
rz(0.94547358) q[1];
sx q[1];
rz(-1.5676326) q[1];
sx q[1];
rz(0.6792773) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0196456) q[0];
sx q[0];
rz(-1.4942193) q[0];
sx q[0];
rz(-2.6640253) q[0];
rz(-pi) q[1];
rz(2.5176454) q[2];
sx q[2];
rz(-1.3610164) q[2];
sx q[2];
rz(-0.35614511) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5577364) q[1];
sx q[1];
rz(-0.5498372) q[1];
sx q[1];
rz(1.932895) q[1];
rz(-pi) q[2];
rz(1.8550917) q[3];
sx q[3];
rz(-1.5522027) q[3];
sx q[3];
rz(0.0068243703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2519553) q[2];
sx q[2];
rz(-2.1869982) q[2];
sx q[2];
rz(-1.4944705) q[2];
rz(2.5877118) q[3];
sx q[3];
rz(-1.5364105) q[3];
sx q[3];
rz(-1.4166098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3225591) q[0];
sx q[0];
rz(-2.366134) q[0];
sx q[0];
rz(-1.092528) q[0];
rz(-2.20772) q[1];
sx q[1];
rz(-1.7364419) q[1];
sx q[1];
rz(0.78528231) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99150554) q[0];
sx q[0];
rz(-1.1508216) q[0];
sx q[0];
rz(-0.49862592) q[0];
rz(0.44194989) q[2];
sx q[2];
rz(-1.0681249) q[2];
sx q[2];
rz(1.3208226) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31639659) q[1];
sx q[1];
rz(-1.8421208) q[1];
sx q[1];
rz(-1.8261938) q[1];
x q[2];
rz(2.5232072) q[3];
sx q[3];
rz(-1.9342285) q[3];
sx q[3];
rz(-1.5314764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87172047) q[2];
sx q[2];
rz(-1.4464658) q[2];
sx q[2];
rz(0.63228697) q[2];
rz(0.91442433) q[3];
sx q[3];
rz(-1.2605896) q[3];
sx q[3];
rz(-2.9818592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3226427) q[0];
sx q[0];
rz(-0.85700789) q[0];
sx q[0];
rz(0.44347611) q[0];
rz(2.8387198) q[1];
sx q[1];
rz(-0.58285204) q[1];
sx q[1];
rz(1.3076967) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6752077) q[0];
sx q[0];
rz(-1.8716011) q[0];
sx q[0];
rz(1.3928901) q[0];
rz(-pi) q[1];
rz(0.50755395) q[2];
sx q[2];
rz(-1.2967464) q[2];
sx q[2];
rz(-0.7407032) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9759226) q[1];
sx q[1];
rz(-2.5745086) q[1];
sx q[1];
rz(2.0815064) q[1];
rz(1.100622) q[3];
sx q[3];
rz(-1.8368145) q[3];
sx q[3];
rz(-3.0373552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0206454) q[2];
sx q[2];
rz(-1.0409313) q[2];
sx q[2];
rz(0.61332235) q[2];
rz(-2.3239418) q[3];
sx q[3];
rz(-0.48912564) q[3];
sx q[3];
rz(-2.8221655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9870975) q[0];
sx q[0];
rz(-1.3767865) q[0];
sx q[0];
rz(2.7823271) q[0];
rz(2.477395) q[1];
sx q[1];
rz(-1.3281053) q[1];
sx q[1];
rz(-0.42116234) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0764685) q[0];
sx q[0];
rz(-1.9389922) q[0];
sx q[0];
rz(-2.9028331) q[0];
rz(-pi) q[1];
x q[1];
rz(1.498327) q[2];
sx q[2];
rz(-1.5654828) q[2];
sx q[2];
rz(1.1491691) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1637472) q[1];
sx q[1];
rz(-2.6085269) q[1];
sx q[1];
rz(-1.0372437) q[1];
rz(0.99160414) q[3];
sx q[3];
rz(-1.0458071) q[3];
sx q[3];
rz(0.92574173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8757214) q[2];
sx q[2];
rz(-3.044812) q[2];
sx q[2];
rz(1.9270012) q[2];
rz(2.752979) q[3];
sx q[3];
rz(-2.2190084) q[3];
sx q[3];
rz(1.0987561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.494396) q[0];
sx q[0];
rz(-1.5955441) q[0];
sx q[0];
rz(1.2936976) q[0];
rz(-2.9248059) q[1];
sx q[1];
rz(-2.4516791) q[1];
sx q[1];
rz(-2.6286415) q[1];
rz(-0.51914712) q[2];
sx q[2];
rz(-2.3943974) q[2];
sx q[2];
rz(1.129528) q[2];
rz(-0.46617266) q[3];
sx q[3];
rz(-1.331658) q[3];
sx q[3];
rz(-0.030203947) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
