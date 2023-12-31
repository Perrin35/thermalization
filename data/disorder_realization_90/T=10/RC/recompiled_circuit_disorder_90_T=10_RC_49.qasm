OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4164299) q[0];
sx q[0];
rz(-0.13983146) q[0];
sx q[0];
rz(-2.5319985) q[0];
rz(0.66863376) q[1];
sx q[1];
rz(-2.2761087) q[1];
sx q[1];
rz(3.0545711) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5370109) q[0];
sx q[0];
rz(-0.7365948) q[0];
sx q[0];
rz(2.1405311) q[0];
rz(2.3743462) q[2];
sx q[2];
rz(-1.6635206) q[2];
sx q[2];
rz(0.35284943) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.39057186) q[1];
sx q[1];
rz(-2.6387072) q[1];
sx q[1];
rz(-1.2807756) q[1];
rz(-pi) q[2];
rz(2.7294455) q[3];
sx q[3];
rz(-1.6133568) q[3];
sx q[3];
rz(1.7736848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0027851) q[2];
sx q[2];
rz(-1.3802718) q[2];
sx q[2];
rz(-0.37386093) q[2];
rz(-2.8047681) q[3];
sx q[3];
rz(-1.5954433) q[3];
sx q[3];
rz(-2.9132304) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2574629) q[0];
sx q[0];
rz(-2.7682436) q[0];
sx q[0];
rz(1.9467547) q[0];
rz(0.082611235) q[1];
sx q[1];
rz(-1.9742842) q[1];
sx q[1];
rz(0.00037489051) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1502991) q[0];
sx q[0];
rz(-2.5851558) q[0];
sx q[0];
rz(2.2668905) q[0];
rz(-pi) q[1];
rz(-0.20034321) q[2];
sx q[2];
rz(-1.8274954) q[2];
sx q[2];
rz(2.8368907) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5217168) q[1];
sx q[1];
rz(-2.4881425) q[1];
sx q[1];
rz(0.39342777) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.625929) q[3];
sx q[3];
rz(-2.030636) q[3];
sx q[3];
rz(0.22051792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80883819) q[2];
sx q[2];
rz(-2.9465582) q[2];
sx q[2];
rz(2.9648798) q[2];
rz(2.3475032) q[3];
sx q[3];
rz(-0.77787557) q[3];
sx q[3];
rz(-2.5879522) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90213838) q[0];
sx q[0];
rz(-2.4376526) q[0];
sx q[0];
rz(-0.40476558) q[0];
rz(1.2813214) q[1];
sx q[1];
rz(-0.57360137) q[1];
sx q[1];
rz(1.3084897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3151911) q[0];
sx q[0];
rz(-1.0747402) q[0];
sx q[0];
rz(1.0415003) q[0];
x q[1];
rz(1.8070418) q[2];
sx q[2];
rz(-2.0958732) q[2];
sx q[2];
rz(-2.6785148) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8150427) q[1];
sx q[1];
rz(-0.76154852) q[1];
sx q[1];
rz(1.9425068) q[1];
x q[2];
rz(-0.19034068) q[3];
sx q[3];
rz(-2.6958709) q[3];
sx q[3];
rz(1.0090855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9222766) q[2];
sx q[2];
rz(-2.6185161) q[2];
sx q[2];
rz(3.1211839) q[2];
rz(2.0698047) q[3];
sx q[3];
rz(-2.0139549) q[3];
sx q[3];
rz(-0.66334692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14169176) q[0];
sx q[0];
rz(-2.2321556) q[0];
sx q[0];
rz(2.2256057) q[0];
rz(0.46332106) q[1];
sx q[1];
rz(-2.0756192) q[1];
sx q[1];
rz(2.0844918) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4823285) q[0];
sx q[0];
rz(-1.5452914) q[0];
sx q[0];
rz(0.43735023) q[0];
x q[1];
rz(-0.57025036) q[2];
sx q[2];
rz(-2.1006561) q[2];
sx q[2];
rz(0.8023163) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.108236) q[1];
sx q[1];
rz(-1.557096) q[1];
sx q[1];
rz(1.800888) q[1];
x q[2];
rz(-0.23987694) q[3];
sx q[3];
rz(-2.3696218) q[3];
sx q[3];
rz(0.039226942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.89381924) q[2];
sx q[2];
rz(-2.5017068) q[2];
sx q[2];
rz(2.7988953) q[2];
rz(-1.6977067) q[3];
sx q[3];
rz(-2.2659437) q[3];
sx q[3];
rz(0.7152344) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3559568) q[0];
sx q[0];
rz(-2.5781093) q[0];
sx q[0];
rz(0.086439565) q[0];
rz(-1.3899639) q[1];
sx q[1];
rz(-2.2098863) q[1];
sx q[1];
rz(2.6729029) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2738004) q[0];
sx q[0];
rz(-2.3228354) q[0];
sx q[0];
rz(1.0206945) q[0];
rz(-pi) q[1];
rz(-2.6195171) q[2];
sx q[2];
rz(-2.7895088) q[2];
sx q[2];
rz(0.075369518) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9710755) q[1];
sx q[1];
rz(-1.3689539) q[1];
sx q[1];
rz(1.438373) q[1];
rz(3.0607871) q[3];
sx q[3];
rz(-2.3256133) q[3];
sx q[3];
rz(0.24085837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9995352) q[2];
sx q[2];
rz(-0.88025981) q[2];
sx q[2];
rz(2.2772677) q[2];
rz(-2.9243829) q[3];
sx q[3];
rz(-2.5253798) q[3];
sx q[3];
rz(0.29278452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7140759) q[0];
sx q[0];
rz(-1.4070516) q[0];
sx q[0];
rz(0.19700225) q[0];
rz(-1.3621832) q[1];
sx q[1];
rz(-1.660659) q[1];
sx q[1];
rz(-2.8053455) q[1];
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
rz(-2.3669846) q[2];
sx q[2];
rz(-2.3850072) q[2];
sx q[2];
rz(-0.87385439) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0834004) q[1];
sx q[1];
rz(-0.5914878) q[1];
sx q[1];
rz(3.0565492) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9365262) q[3];
sx q[3];
rz(-2.4347217) q[3];
sx q[3];
rz(2.6170078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.53529915) q[2];
sx q[2];
rz(-1.4761816) q[2];
sx q[2];
rz(2.2952648) q[2];
rz(-1.2396631) q[3];
sx q[3];
rz(-1.2318434) q[3];
sx q[3];
rz(-0.0011750778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3939312) q[0];
sx q[0];
rz(-2.1037536) q[0];
sx q[0];
rz(-2.8826662) q[0];
rz(1.7954284) q[1];
sx q[1];
rz(-1.3849473) q[1];
sx q[1];
rz(2.0475725) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31539311) q[0];
sx q[0];
rz(-1.1356192) q[0];
sx q[0];
rz(-0.3776334) q[0];
rz(1.3482434) q[2];
sx q[2];
rz(-1.5828653) q[2];
sx q[2];
rz(-1.0782858) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8158405) q[1];
sx q[1];
rz(-0.098512352) q[1];
sx q[1];
rz(-1.2835531) q[1];
rz(-1.592698) q[3];
sx q[3];
rz(-1.6688445) q[3];
sx q[3];
rz(2.736562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.1743494) q[2];
sx q[2];
rz(-1.7847585) q[2];
sx q[2];
rz(-2.8209177) q[2];
rz(0.43618068) q[3];
sx q[3];
rz(-0.47302055) q[3];
sx q[3];
rz(2.6935553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86673474) q[0];
sx q[0];
rz(-2.6651356) q[0];
sx q[0];
rz(-2.136769) q[0];
rz(-2.5514305) q[1];
sx q[1];
rz(-2.2361123) q[1];
sx q[1];
rz(-2.337713) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2893387) q[0];
sx q[0];
rz(-2.5544871) q[0];
sx q[0];
rz(0.071285204) q[0];
x q[1];
rz(-0.29381807) q[2];
sx q[2];
rz(-2.567798) q[2];
sx q[2];
rz(2.3125355) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4623973) q[1];
sx q[1];
rz(-0.089226626) q[1];
sx q[1];
rz(-1.4560844) q[1];
rz(-pi) q[2];
rz(-3.026268) q[3];
sx q[3];
rz(-0.66676312) q[3];
sx q[3];
rz(-2.525884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.33264318) q[2];
sx q[2];
rz(-0.84471622) q[2];
sx q[2];
rz(2.1311029) q[2];
rz(-2.5381952) q[3];
sx q[3];
rz(-1.5248652) q[3];
sx q[3];
rz(2.5914014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6036966) q[0];
sx q[0];
rz(-1.2768856) q[0];
sx q[0];
rz(2.7340775) q[0];
rz(-2.852476) q[1];
sx q[1];
rz(-2.018785) q[1];
sx q[1];
rz(-2.3908652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5068881) q[0];
sx q[0];
rz(-1.2872211) q[0];
sx q[0];
rz(0.30446913) q[0];
rz(1.4892011) q[2];
sx q[2];
rz(-1.8494693) q[2];
sx q[2];
rz(-0.96688731) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51160204) q[1];
sx q[1];
rz(-1.8774319) q[1];
sx q[1];
rz(-0.2125477) q[1];
rz(-pi) q[2];
rz(-0.30131807) q[3];
sx q[3];
rz(-1.5958061) q[3];
sx q[3];
rz(-3.0400288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0728545) q[2];
sx q[2];
rz(-2.4170503) q[2];
sx q[2];
rz(-1.9753974) q[2];
rz(-1.3646305) q[3];
sx q[3];
rz(-2.3659673) q[3];
sx q[3];
rz(0.021818074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5831379) q[0];
sx q[0];
rz(-2.3151509) q[0];
sx q[0];
rz(-1.7512084) q[0];
rz(2.8109) q[1];
sx q[1];
rz(-2.3762517) q[1];
sx q[1];
rz(1.6814544) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5064142) q[0];
sx q[0];
rz(-1.9902475) q[0];
sx q[0];
rz(-1.0181396) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6129458) q[2];
sx q[2];
rz(-2.6372006) q[2];
sx q[2];
rz(-1.1500037) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8148988) q[1];
sx q[1];
rz(-1.9683451) q[1];
sx q[1];
rz(-1.3709929) q[1];
rz(-pi) q[2];
rz(-2.4246033) q[3];
sx q[3];
rz(-2.2319712) q[3];
sx q[3];
rz(-0.17718525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7918487) q[2];
sx q[2];
rz(-1.8169553) q[2];
sx q[2];
rz(1.4617408) q[2];
rz(2.0215624) q[3];
sx q[3];
rz(-0.62168613) q[3];
sx q[3];
rz(-0.45564836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6256975) q[0];
sx q[0];
rz(-1.6202171) q[0];
sx q[0];
rz(1.1465999) q[0];
rz(-1.760578) q[1];
sx q[1];
rz(-1.2881423) q[1];
sx q[1];
rz(-1.2013411) q[1];
rz(-2.3537221) q[2];
sx q[2];
rz(-1.8265767) q[2];
sx q[2];
rz(1.5192601) q[2];
rz(-2.667726) q[3];
sx q[3];
rz(-2.4110473) q[3];
sx q[3];
rz(-1.9163781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
