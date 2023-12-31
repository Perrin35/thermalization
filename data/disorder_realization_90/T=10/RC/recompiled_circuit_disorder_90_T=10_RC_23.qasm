OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7251627) q[0];
sx q[0];
rz(-3.0017612) q[0];
sx q[0];
rz(2.5319985) q[0];
rz(-2.4729589) q[1];
sx q[1];
rz(-0.86548391) q[1];
sx q[1];
rz(-3.0545711) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10843033) q[0];
sx q[0];
rz(-0.96956367) q[0];
sx q[0];
rz(2.6866459) q[0];
x q[1];
rz(0.1331698) q[2];
sx q[2];
rz(-0.77169092) q[2];
sx q[2];
rz(2.0193677) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.39057186) q[1];
sx q[1];
rz(-2.6387072) q[1];
sx q[1];
rz(-1.2807756) q[1];
rz(-3.0356785) q[3];
sx q[3];
rz(-0.4142136) q[3];
sx q[3];
rz(-2.8416882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0027851) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(-0.37386093) q[2];
rz(0.3368245) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(0.22836223) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8841298) q[0];
sx q[0];
rz(-2.7682436) q[0];
sx q[0];
rz(-1.194838) q[0];
rz(-0.082611235) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(3.1412178) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99129358) q[0];
sx q[0];
rz(-0.55643686) q[0];
sx q[0];
rz(2.2668905) q[0];
rz(-pi) q[1];
rz(0.20034321) q[2];
sx q[2];
rz(-1.3140972) q[2];
sx q[2];
rz(2.8368907) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13818571) q[1];
sx q[1];
rz(-2.1669743) q[1];
sx q[1];
rz(1.8562993) q[1];
rz(-1.5156636) q[3];
sx q[3];
rz(-2.030636) q[3];
sx q[3];
rz(-0.22051792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3327545) q[2];
sx q[2];
rz(-2.9465582) q[2];
sx q[2];
rz(2.9648798) q[2];
rz(0.79408944) q[3];
sx q[3];
rz(-0.77787557) q[3];
sx q[3];
rz(-0.55364048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90213838) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(-2.7368271) q[0];
rz(1.8602712) q[1];
sx q[1];
rz(-2.5679913) q[1];
sx q[1];
rz(-1.8331029) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.20298) q[0];
sx q[0];
rz(-0.70883195) q[0];
sx q[0];
rz(2.3908486) q[0];
rz(-1.3345509) q[2];
sx q[2];
rz(-1.0457195) q[2];
sx q[2];
rz(2.6785148) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82061003) q[1];
sx q[1];
rz(-0.87249407) q[1];
sx q[1];
rz(-0.33336158) q[1];
rz(-1.4806467) q[3];
sx q[3];
rz(-1.1336859) q[3];
sx q[3];
rz(2.3428832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9222766) q[2];
sx q[2];
rz(-0.52307659) q[2];
sx q[2];
rz(-0.02040872) q[2];
rz(-2.0698047) q[3];
sx q[3];
rz(-1.1276378) q[3];
sx q[3];
rz(2.4782457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9999009) q[0];
sx q[0];
rz(-2.2321556) q[0];
sx q[0];
rz(-2.2256057) q[0];
rz(-0.46332106) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(1.0571009) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4823285) q[0];
sx q[0];
rz(-1.5963012) q[0];
sx q[0];
rz(2.7042424) q[0];
rz(-pi) q[1];
rz(-2.1787203) q[2];
sx q[2];
rz(-2.0553556) q[2];
sx q[2];
rz(-2.059666) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4790198) q[1];
sx q[1];
rz(-2.9111007) q[1];
sx q[1];
rz(1.6307994) q[1];
x q[2];
rz(-1.7980868) q[3];
sx q[3];
rz(-0.82633457) q[3];
sx q[3];
rz(2.8518761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2477734) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(-2.7988953) q[2];
rz(1.4438859) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3559568) q[0];
sx q[0];
rz(-0.56348339) q[0];
sx q[0];
rz(-0.086439565) q[0];
rz(-1.3899639) q[1];
sx q[1];
rz(-2.2098863) q[1];
sx q[1];
rz(-0.46868971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86779224) q[0];
sx q[0];
rz(-2.3228354) q[0];
sx q[0];
rz(2.1208982) q[0];
rz(-pi) q[1];
rz(2.8332963) q[2];
sx q[2];
rz(-1.7436276) q[2];
sx q[2];
rz(1.0002713) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3846181) q[1];
sx q[1];
rz(-2.9006835) q[1];
sx q[1];
rz(0.5730281) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.080805578) q[3];
sx q[3];
rz(-0.81597933) q[3];
sx q[3];
rz(2.9007343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9995352) q[2];
sx q[2];
rz(-2.2613328) q[2];
sx q[2];
rz(2.2772677) q[2];
rz(2.9243829) q[3];
sx q[3];
rz(-2.5253798) q[3];
sx q[3];
rz(-0.29278452) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7140759) q[0];
sx q[0];
rz(-1.7345411) q[0];
sx q[0];
rz(-0.19700225) q[0];
rz(1.7794094) q[1];
sx q[1];
rz(-1.4809337) q[1];
sx q[1];
rz(-0.33624712) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2339736) q[0];
sx q[0];
rz(-2.189878) q[0];
sx q[0];
rz(1.0236077) q[0];
rz(-pi) q[1];
rz(2.3669846) q[2];
sx q[2];
rz(-2.3850072) q[2];
sx q[2];
rz(0.87385439) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.044144883) q[1];
sx q[1];
rz(-2.1598585) q[1];
sx q[1];
rz(1.6277905) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4451609) q[3];
sx q[3];
rz(-1.7034354) q[3];
sx q[3];
rz(-1.2030676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6062935) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(2.2952648) q[2];
rz(-1.9019295) q[3];
sx q[3];
rz(-1.9097493) q[3];
sx q[3];
rz(3.1404176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3939312) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(-0.25892648) q[0];
rz(-1.3461643) q[1];
sx q[1];
rz(-1.7566453) q[1];
sx q[1];
rz(1.0940201) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0703697) q[0];
sx q[0];
rz(-0.56814146) q[0];
sx q[0];
rz(2.2413261) q[0];
x q[1];
rz(-1.5161683) q[2];
sx q[2];
rz(-2.9187181) q[2];
sx q[2];
rz(-0.54578997) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.5309696) q[1];
sx q[1];
rz(-1.5429284) q[1];
sx q[1];
rz(-1.665297) q[1];
x q[2];
rz(0.098071531) q[3];
sx q[3];
rz(-1.5489998) q[3];
sx q[3];
rz(-1.1636213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.1743494) q[2];
sx q[2];
rz(-1.7847585) q[2];
sx q[2];
rz(2.8209177) q[2];
rz(2.705412) q[3];
sx q[3];
rz(-2.6685721) q[3];
sx q[3];
rz(-0.44803739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86673474) q[0];
sx q[0];
rz(-2.6651356) q[0];
sx q[0];
rz(2.136769) q[0];
rz(0.59016219) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(-2.337713) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76669508) q[0];
sx q[0];
rz(-2.156213) q[0];
sx q[0];
rz(-1.5234408) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29381807) q[2];
sx q[2];
rz(-0.57379469) q[2];
sx q[2];
rz(2.3125355) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4623973) q[1];
sx q[1];
rz(-0.089226626) q[1];
sx q[1];
rz(1.6855082) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6611093) q[3];
sx q[3];
rz(-2.2323425) q[3];
sx q[3];
rz(2.6722398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.33264318) q[2];
sx q[2];
rz(-2.2968764) q[2];
sx q[2];
rz(1.0104898) q[2];
rz(2.5381952) q[3];
sx q[3];
rz(-1.6167275) q[3];
sx q[3];
rz(-0.55019125) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6036966) q[0];
sx q[0];
rz(-1.2768856) q[0];
sx q[0];
rz(-0.4075152) q[0];
rz(2.852476) q[1];
sx q[1];
rz(-1.1228077) q[1];
sx q[1];
rz(-2.3908652) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2931965) q[0];
sx q[0];
rz(-1.2788532) q[0];
sx q[0];
rz(1.2743203) q[0];
x q[1];
rz(0.27751343) q[2];
sx q[2];
rz(-0.29007426) q[2];
sx q[2];
rz(-1.8857423) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.99414413) q[1];
sx q[1];
rz(-1.3683043) q[1];
sx q[1];
rz(1.2575498) q[1];
rz(-0.30131807) q[3];
sx q[3];
rz(-1.5457866) q[3];
sx q[3];
rz(-0.1015639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0687381) q[2];
sx q[2];
rz(-2.4170503) q[2];
sx q[2];
rz(1.1661952) q[2];
rz(-1.7769622) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(-0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5831379) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(-1.7512084) q[0];
rz(-0.33069262) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(-1.4601382) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64797348) q[0];
sx q[0];
rz(-2.4612392) q[0];
sx q[0];
rz(2.2749659) q[0];
rz(1.5286469) q[2];
sx q[2];
rz(-2.6372006) q[2];
sx q[2];
rz(1.1500037) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9863661) q[1];
sx q[1];
rz(-0.44253293) q[1];
sx q[1];
rz(0.44154422) q[1];
x q[2];
rz(0.86942418) q[3];
sx q[3];
rz(-0.93360177) q[3];
sx q[3];
rz(2.3616392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7918487) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(1.4617408) q[2];
rz(1.1200303) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(-2.6859443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5158952) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(-1.3810146) q[1];
sx q[1];
rz(-1.8534503) q[1];
sx q[1];
rz(1.9402515) q[1];
rz(-1.21576) q[2];
sx q[2];
rz(-0.81510614) q[2];
sx q[2];
rz(0.19744273) q[2];
rz(1.1827042) q[3];
sx q[3];
rz(-2.2065065) q[3];
sx q[3];
rz(1.828215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
