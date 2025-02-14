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
rz(-2.3186853) q[0];
sx q[0];
rz(3.5003852) q[0];
sx q[0];
rz(8.556463) q[0];
rz(1.7262285) q[1];
sx q[1];
rz(-2.0077029) q[1];
sx q[1];
rz(0.99376065) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9159587) q[0];
sx q[0];
rz(-1.4580112) q[0];
sx q[0];
rz(-0.40333545) q[0];
rz(-pi) q[1];
rz(0.5014855) q[2];
sx q[2];
rz(-1.9619313) q[2];
sx q[2];
rz(-0.30539612) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.704761) q[1];
sx q[1];
rz(-2.0109573) q[1];
sx q[1];
rz(-2.6635567) q[1];
rz(1.1339784) q[3];
sx q[3];
rz(-1.7941495) q[3];
sx q[3];
rz(-0.88296605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.22780861) q[2];
sx q[2];
rz(-1.8081534) q[2];
sx q[2];
rz(2.559973) q[2];
rz(-2.4070814) q[3];
sx q[3];
rz(-1.651265) q[3];
sx q[3];
rz(-0.41828004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2468579) q[0];
sx q[0];
rz(-1.7623836) q[0];
sx q[0];
rz(0.49767622) q[0];
rz(1.0408164) q[1];
sx q[1];
rz(-0.40319315) q[1];
sx q[1];
rz(1.9042447) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5631249) q[0];
sx q[0];
rz(-1.5284561) q[0];
sx q[0];
rz(-2.4189831) q[0];
rz(-pi) q[1];
rz(0.18816973) q[2];
sx q[2];
rz(-1.9261179) q[2];
sx q[2];
rz(2.4057092) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7844441) q[1];
sx q[1];
rz(-1.7116575) q[1];
sx q[1];
rz(2.0037093) q[1];
rz(-pi) q[2];
rz(1.7877616) q[3];
sx q[3];
rz(-1.3111909) q[3];
sx q[3];
rz(1.5322154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.70025468) q[2];
sx q[2];
rz(-2.4105218) q[2];
sx q[2];
rz(-0.31044427) q[2];
rz(1.9849518) q[3];
sx q[3];
rz(-2.0168596) q[3];
sx q[3];
rz(-1.1667075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1336841) q[0];
sx q[0];
rz(-0.5846566) q[0];
sx q[0];
rz(-2.7957918) q[0];
rz(-0.24468228) q[1];
sx q[1];
rz(-0.9318277) q[1];
sx q[1];
rz(1.4366879) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0416243) q[0];
sx q[0];
rz(-1.9910553) q[0];
sx q[0];
rz(0.34456518) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57049306) q[2];
sx q[2];
rz(-1.2353503) q[2];
sx q[2];
rz(0.48656175) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86307058) q[1];
sx q[1];
rz(-2.3211042) q[1];
sx q[1];
rz(-2.8500975) q[1];
rz(1.9448024) q[3];
sx q[3];
rz(-2.0197649) q[3];
sx q[3];
rz(1.6351007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1338542) q[2];
sx q[2];
rz(-2.014092) q[2];
sx q[2];
rz(-0.63068843) q[2];
rz(-0.33356365) q[3];
sx q[3];
rz(-1.0015229) q[3];
sx q[3];
rz(1.001531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0665322) q[0];
sx q[0];
rz(-2.1933031) q[0];
sx q[0];
rz(-1.7370268) q[0];
rz(-2.7724077) q[1];
sx q[1];
rz(-1.7351979) q[1];
sx q[1];
rz(0.085748347) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2959901) q[0];
sx q[0];
rz(-1.1628502) q[0];
sx q[0];
rz(-1.244581) q[0];
x q[1];
rz(-0.46041885) q[2];
sx q[2];
rz(-1.9398089) q[2];
sx q[2];
rz(-2.2031914) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4032674) q[1];
sx q[1];
rz(-1.8083296) q[1];
sx q[1];
rz(1.0942142) q[1];
rz(-1.2530009) q[3];
sx q[3];
rz(-1.1927342) q[3];
sx q[3];
rz(1.8658569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.839445) q[2];
sx q[2];
rz(-1.2371233) q[2];
sx q[2];
rz(-2.6117924) q[2];
rz(0.038330404) q[3];
sx q[3];
rz(-0.72961346) q[3];
sx q[3];
rz(1.5420325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2235276) q[0];
sx q[0];
rz(-0.46850884) q[0];
sx q[0];
rz(-1.202762) q[0];
rz(-0.27944061) q[1];
sx q[1];
rz(-2.1165106) q[1];
sx q[1];
rz(-0.7739982) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27085486) q[0];
sx q[0];
rz(-1.8841198) q[0];
sx q[0];
rz(-0.82672755) q[0];
rz(-pi) q[1];
rz(1.4882795) q[2];
sx q[2];
rz(-1.3386209) q[2];
sx q[2];
rz(-1.2007932) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8355636) q[1];
sx q[1];
rz(-1.4317498) q[1];
sx q[1];
rz(3.1285888) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0353885) q[3];
sx q[3];
rz(-1.5683953) q[3];
sx q[3];
rz(0.46214275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1191795) q[2];
sx q[2];
rz(-3.000562) q[2];
sx q[2];
rz(-2.5775583) q[2];
rz(0.65308475) q[3];
sx q[3];
rz(-2.0061195) q[3];
sx q[3];
rz(-0.45026067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7804724) q[0];
sx q[0];
rz(-0.55279624) q[0];
sx q[0];
rz(-0.64315382) q[0];
rz(1.9505352) q[1];
sx q[1];
rz(-1.4570313) q[1];
sx q[1];
rz(-2.4868884) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7793286) q[0];
sx q[0];
rz(-1.6342499) q[0];
sx q[0];
rz(-0.16880798) q[0];
x q[1];
rz(1.5327318) q[2];
sx q[2];
rz(-2.2543779) q[2];
sx q[2];
rz(0.75929196) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5703598) q[1];
sx q[1];
rz(-2.1132937) q[1];
sx q[1];
rz(0.076537655) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7406171) q[3];
sx q[3];
rz(-0.39547503) q[3];
sx q[3];
rz(-0.91646376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5074629) q[2];
sx q[2];
rz(-2.7644988) q[2];
sx q[2];
rz(1.6667574) q[2];
rz(-2.6569488) q[3];
sx q[3];
rz(-2.1438997) q[3];
sx q[3];
rz(-1.2887597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5019048) q[0];
sx q[0];
rz(-2.8785093) q[0];
sx q[0];
rz(-2.4581773) q[0];
rz(-3.0112093) q[1];
sx q[1];
rz(-1.5510635) q[1];
sx q[1];
rz(0.032616671) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17856971) q[0];
sx q[0];
rz(-1.7462329) q[0];
sx q[0];
rz(-1.001872) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8993699) q[2];
sx q[2];
rz(-0.70021473) q[2];
sx q[2];
rz(-1.2716573) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4144145) q[1];
sx q[1];
rz(-1.7373996) q[1];
sx q[1];
rz(1.3018621) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8835589) q[3];
sx q[3];
rz(-0.88538187) q[3];
sx q[3];
rz(0.097921927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42328295) q[2];
sx q[2];
rz(-1.1085359) q[2];
sx q[2];
rz(0.56524593) q[2];
rz(1.7806753) q[3];
sx q[3];
rz(-2.8748685) q[3];
sx q[3];
rz(2.3274073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3703506) q[0];
sx q[0];
rz(-3.0841565) q[0];
sx q[0];
rz(-2.8420319) q[0];
rz(1.4467422) q[1];
sx q[1];
rz(-0.87211496) q[1];
sx q[1];
rz(0.95132336) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46630105) q[0];
sx q[0];
rz(-1.6032752) q[0];
sx q[0];
rz(-1.7759274) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0092333992) q[2];
sx q[2];
rz(-1.7823185) q[2];
sx q[2];
rz(2.4516425) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0971157) q[1];
sx q[1];
rz(-1.591178) q[1];
sx q[1];
rz(-2.7173032) q[1];
x q[2];
rz(-0.087835066) q[3];
sx q[3];
rz(-1.8316275) q[3];
sx q[3];
rz(1.6139327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1424554) q[2];
sx q[2];
rz(-0.74644011) q[2];
sx q[2];
rz(0.26774055) q[2];
rz(-1.0033876) q[3];
sx q[3];
rz(-0.95971862) q[3];
sx q[3];
rz(0.80823922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35850152) q[0];
sx q[0];
rz(-2.6434904) q[0];
sx q[0];
rz(2.1844693) q[0];
rz(2.762291) q[1];
sx q[1];
rz(-0.4117659) q[1];
sx q[1];
rz(-2.9764825) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9921761) q[0];
sx q[0];
rz(-0.78127978) q[0];
sx q[0];
rz(-2.6447456) q[0];
rz(-pi) q[1];
rz(-1.6696641) q[2];
sx q[2];
rz(-2.6672088) q[2];
sx q[2];
rz(-2.6883467) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2735426) q[1];
sx q[1];
rz(-0.99813491) q[1];
sx q[1];
rz(-0.21577253) q[1];
rz(-pi) q[2];
rz(-2.9028724) q[3];
sx q[3];
rz(-2.1607499) q[3];
sx q[3];
rz(1.1600398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2740606) q[2];
sx q[2];
rz(-1.910285) q[2];
sx q[2];
rz(2.3629698) q[2];
rz(-2.8042931) q[3];
sx q[3];
rz(-1.0236579) q[3];
sx q[3];
rz(1.5844828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.857665) q[0];
sx q[0];
rz(-1.090467) q[0];
sx q[0];
rz(2.0528059) q[0];
rz(2.7133443) q[1];
sx q[1];
rz(-2.1183522) q[1];
sx q[1];
rz(-0.64839378) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89598362) q[0];
sx q[0];
rz(-1.147384) q[0];
sx q[0];
rz(-1.8603252) q[0];
rz(-pi) q[1];
rz(1.9087547) q[2];
sx q[2];
rz(-1.612886) q[2];
sx q[2];
rz(-1.011285) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.32247816) q[1];
sx q[1];
rz(-1.6129061) q[1];
sx q[1];
rz(-2.4933715) q[1];
rz(-pi) q[2];
rz(-1.6850059) q[3];
sx q[3];
rz(-1.0100967) q[3];
sx q[3];
rz(1.4637092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.26583656) q[2];
sx q[2];
rz(-1.0886322) q[2];
sx q[2];
rz(-0.17975532) q[2];
rz(1.1579375) q[3];
sx q[3];
rz(-1.4561184) q[3];
sx q[3];
rz(1.2402844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4614048) q[0];
sx q[0];
rz(-2.9869933) q[0];
sx q[0];
rz(-2.3416478) q[0];
rz(1.9752621) q[1];
sx q[1];
rz(-0.98465289) q[1];
sx q[1];
rz(-2.226895) q[1];
rz(0.24576743) q[2];
sx q[2];
rz(-1.2817597) q[2];
sx q[2];
rz(1.3195932) q[2];
rz(2.8020482) q[3];
sx q[3];
rz(-1.3673269) q[3];
sx q[3];
rz(0.23585503) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
