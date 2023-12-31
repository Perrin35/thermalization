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
rz(0.13983146) q[0];
sx q[0];
rz(10.034372) q[0];
rz(-2.4729589) q[1];
sx q[1];
rz(-0.86548391) q[1];
sx q[1];
rz(-3.0545711) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5370109) q[0];
sx q[0];
rz(-0.7365948) q[0];
sx q[0];
rz(-2.1405311) q[0];
rz(-pi) q[1];
rz(1.6992703) q[2];
sx q[2];
rz(-2.3339084) q[2];
sx q[2];
rz(1.8345923) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.062292369) q[1];
sx q[1];
rz(-2.0508517) q[1];
sx q[1];
rz(0.15602195) q[1];
rz(-pi) q[2];
rz(2.7294455) q[3];
sx q[3];
rz(-1.6133568) q[3];
sx q[3];
rz(1.7736848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0027851) q[2];
sx q[2];
rz(-1.7613208) q[2];
sx q[2];
rz(0.37386093) q[2];
rz(-2.8047681) q[3];
sx q[3];
rz(-1.5461494) q[3];
sx q[3];
rz(-0.22836223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8841298) q[0];
sx q[0];
rz(-2.7682436) q[0];
sx q[0];
rz(-1.9467547) q[0];
rz(3.0589814) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(-0.00037489051) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9449687) q[0];
sx q[0];
rz(-1.225291) q[0];
sx q[0];
rz(-2.0161122) q[0];
rz(-1.3090918) q[2];
sx q[2];
rz(-1.3771025) q[2];
sx q[2];
rz(-1.2145834) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2692711) q[1];
sx q[1];
rz(-1.806013) q[1];
sx q[1];
rz(-2.5260731) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0307816) q[3];
sx q[3];
rz(-2.6786945) q[3];
sx q[3];
rz(0.09679951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80883819) q[2];
sx q[2];
rz(-2.9465582) q[2];
sx q[2];
rz(-2.9648798) q[2];
rz(2.3475032) q[3];
sx q[3];
rz(-0.77787557) q[3];
sx q[3];
rz(-2.5879522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90213838) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(-2.7368271) q[0];
rz(-1.8602712) q[1];
sx q[1];
rz(-2.5679913) q[1];
sx q[1];
rz(-1.3084897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8264016) q[0];
sx q[0];
rz(-2.0668525) q[0];
sx q[0];
rz(2.1000923) q[0];
rz(-pi) q[1];
rz(1.3345509) q[2];
sx q[2];
rz(-2.0958732) q[2];
sx q[2];
rz(-0.46307785) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8150427) q[1];
sx q[1];
rz(-0.76154852) q[1];
sx q[1];
rz(1.1990859) q[1];
x q[2];
rz(-2.951252) q[3];
sx q[3];
rz(-0.44572178) q[3];
sx q[3];
rz(1.0090855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.219316) q[2];
sx q[2];
rz(-0.52307659) q[2];
sx q[2];
rz(-0.02040872) q[2];
rz(-1.071788) q[3];
sx q[3];
rz(-1.1276378) q[3];
sx q[3];
rz(-2.4782457) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14169176) q[0];
sx q[0];
rz(-2.2321556) q[0];
sx q[0];
rz(2.2256057) q[0];
rz(-0.46332106) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(1.0571009) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9986344) q[0];
sx q[0];
rz(-0.43804533) q[0];
sx q[0];
rz(3.0814339) q[0];
x q[1];
rz(2.1787203) q[2];
sx q[2];
rz(-1.086237) q[2];
sx q[2];
rz(-2.059666) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6009439) q[1];
sx q[1];
rz(-1.3407267) q[1];
sx q[1];
rz(-0.01407108) q[1];
x q[2];
rz(0.23987694) q[3];
sx q[3];
rz(-0.7719709) q[3];
sx q[3];
rz(-3.1023657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89381924) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(0.34269732) q[2];
rz(1.4438859) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(0.7152344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(2.6729029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86779224) q[0];
sx q[0];
rz(-2.3228354) q[0];
sx q[0];
rz(-1.0206945) q[0];
rz(-pi) q[1];
rz(0.52207559) q[2];
sx q[2];
rz(-0.35208382) q[2];
sx q[2];
rz(3.0662231) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3846181) q[1];
sx q[1];
rz(-0.24090919) q[1];
sx q[1];
rz(-0.5730281) q[1];
x q[2];
rz(-0.080805578) q[3];
sx q[3];
rz(-2.3256133) q[3];
sx q[3];
rz(-2.9007343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1420574) q[2];
sx q[2];
rz(-0.88025981) q[2];
sx q[2];
rz(-0.86432499) q[2];
rz(-2.9243829) q[3];
sx q[3];
rz(-2.5253798) q[3];
sx q[3];
rz(-2.8488081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42751673) q[0];
sx q[0];
rz(-1.4070516) q[0];
sx q[0];
rz(2.9445904) q[0];
rz(-1.7794094) q[1];
sx q[1];
rz(-1.660659) q[1];
sx q[1];
rz(2.8053455) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67664458) q[0];
sx q[0];
rz(-2.0083545) q[0];
sx q[0];
rz(-2.4462571) q[0];
x q[1];
rz(-0.98724987) q[2];
sx q[2];
rz(-1.0580214) q[2];
sx q[2];
rz(-3.0836881) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5583401) q[1];
sx q[1];
rz(-1.5234158) q[1];
sx q[1];
rz(-2.5517795) q[1];
rz(-pi) q[2];
rz(-2.4451609) q[3];
sx q[3];
rz(-1.7034354) q[3];
sx q[3];
rz(-1.9385251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6062935) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(0.84632787) q[2];
rz(1.9019295) q[3];
sx q[3];
rz(-1.9097493) q[3];
sx q[3];
rz(0.0011750778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3939312) q[0];
sx q[0];
rz(-1.0378391) q[0];
sx q[0];
rz(2.8826662) q[0];
rz(-1.7954284) q[1];
sx q[1];
rz(-1.7566453) q[1];
sx q[1];
rz(-1.0940201) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8261995) q[0];
sx q[0];
rz(-1.1356192) q[0];
sx q[0];
rz(0.3776334) q[0];
rz(1.3482434) q[2];
sx q[2];
rz(-1.5828653) q[2];
sx q[2];
rz(-1.0782858) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6106231) q[1];
sx q[1];
rz(-1.5429284) q[1];
sx q[1];
rz(-1.4762957) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.098071531) q[3];
sx q[3];
rz(-1.5489998) q[3];
sx q[3];
rz(1.1636213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9672433) q[2];
sx q[2];
rz(-1.3568342) q[2];
sx q[2];
rz(-2.8209177) q[2];
rz(-2.705412) q[3];
sx q[3];
rz(-2.6685721) q[3];
sx q[3];
rz(-2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2748579) q[0];
sx q[0];
rz(-0.47645706) q[0];
sx q[0];
rz(2.136769) q[0];
rz(-0.59016219) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(-0.80387962) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3636712) q[0];
sx q[0];
rz(-1.6102618) q[0];
sx q[0];
rz(0.58593336) q[0];
rz(-pi) q[1];
rz(-2.5875823) q[2];
sx q[2];
rz(-1.7286574) q[2];
sx q[2];
rz(-0.49288921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56402928) q[1];
sx q[1];
rz(-1.4821577) q[1];
sx q[1];
rz(0.010239756) q[1];
rz(-pi) q[2];
rz(-3.026268) q[3];
sx q[3];
rz(-0.66676312) q[3];
sx q[3];
rz(-2.525884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8089495) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(1.0104898) q[2];
rz(-2.5381952) q[3];
sx q[3];
rz(-1.5248652) q[3];
sx q[3];
rz(-0.55019125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6036966) q[0];
sx q[0];
rz(-1.2768856) q[0];
sx q[0];
rz(-2.7340775) q[0];
rz(2.852476) q[1];
sx q[1];
rz(-2.018785) q[1];
sx q[1];
rz(2.3908652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4780316) q[0];
sx q[0];
rz(-2.7285517) q[0];
sx q[0];
rz(-0.77126276) q[0];
x q[1];
rz(-2.8640792) q[2];
sx q[2];
rz(-2.8515184) q[2];
sx q[2];
rz(-1.2558503) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6299906) q[1];
sx q[1];
rz(-1.2641608) q[1];
sx q[1];
rz(-2.929045) q[1];
rz(-pi) q[2];
x q[2];
rz(0.084089355) q[3];
sx q[3];
rz(-2.8392699) q[3];
sx q[3];
rz(-1.7526527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0687381) q[2];
sx q[2];
rz(-2.4170503) q[2];
sx q[2];
rz(1.1661952) q[2];
rz(-1.7769622) q[3];
sx q[3];
rz(-0.77562538) q[3];
sx q[3];
rz(0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5584548) q[0];
sx q[0];
rz(-0.82644176) q[0];
sx q[0];
rz(1.3903842) q[0];
rz(0.33069262) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(-1.6814544) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64797348) q[0];
sx q[0];
rz(-0.68035347) q[0];
sx q[0];
rz(2.2749659) q[0];
x q[1];
rz(1.6129458) q[2];
sx q[2];
rz(-2.6372006) q[2];
sx q[2];
rz(-1.1500037) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8192484) q[1];
sx q[1];
rz(-1.3867612) q[1];
sx q[1];
rz(-0.40477246) q[1];
x q[2];
rz(2.2721685) q[3];
sx q[3];
rz(-0.93360177) q[3];
sx q[3];
rz(-2.3616392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34974393) q[2];
sx q[2];
rz(-1.3246374) q[2];
sx q[2];
rz(-1.6798518) q[2];
rz(-1.1200303) q[3];
sx q[3];
rz(-2.5199065) q[3];
sx q[3];
rz(-2.6859443) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5158952) q[0];
sx q[0];
rz(-1.5213756) q[0];
sx q[0];
rz(-1.9949927) q[0];
rz(1.3810146) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(0.78787055) q[2];
sx q[2];
rz(-1.8265767) q[2];
sx q[2];
rz(1.5192601) q[2];
rz(1.9588884) q[3];
sx q[3];
rz(-0.93508616) q[3];
sx q[3];
rz(-1.3133776) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
