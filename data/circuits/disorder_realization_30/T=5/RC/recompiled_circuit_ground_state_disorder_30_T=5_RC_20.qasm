OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.339191) q[0];
sx q[0];
rz(5.1270687) q[0];
sx q[0];
rz(13.083885) q[0];
rz(1.327688) q[1];
sx q[1];
rz(7.1751243) q[1];
sx q[1];
rz(8.8827477) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21236496) q[0];
sx q[0];
rz(-2.3327418) q[0];
sx q[0];
rz(-1.8083284) q[0];
x q[1];
rz(-2.3752604) q[2];
sx q[2];
rz(-1.7558985) q[2];
sx q[2];
rz(1.2183777) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0412647) q[1];
sx q[1];
rz(-0.81721837) q[1];
sx q[1];
rz(-0.27268091) q[1];
rz(-pi) q[2];
rz(2.9085701) q[3];
sx q[3];
rz(-1.7400422) q[3];
sx q[3];
rz(-0.47692933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1338256) q[2];
sx q[2];
rz(-2.0962891) q[2];
sx q[2];
rz(-0.93057752) q[2];
rz(-3.0104356) q[3];
sx q[3];
rz(-0.37808642) q[3];
sx q[3];
rz(1.650943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5432878) q[0];
sx q[0];
rz(-0.85619339) q[0];
sx q[0];
rz(1.8970066) q[0];
rz(-0.92567956) q[1];
sx q[1];
rz(-1.2558179) q[1];
sx q[1];
rz(1.4941039) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0605436) q[0];
sx q[0];
rz(-0.88969066) q[0];
sx q[0];
rz(-0.42411719) q[0];
rz(-0.67333198) q[2];
sx q[2];
rz(-1.3048561) q[2];
sx q[2];
rz(0.8348726) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0706961) q[1];
sx q[1];
rz(-0.58377171) q[1];
sx q[1];
rz(2.4811005) q[1];
x q[2];
rz(3.1140987) q[3];
sx q[3];
rz(-2.8373233) q[3];
sx q[3];
rz(-0.19434838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99574789) q[2];
sx q[2];
rz(-2.6971942) q[2];
sx q[2];
rz(0.90163976) q[2];
rz(1.3580648) q[3];
sx q[3];
rz(-1.3172904) q[3];
sx q[3];
rz(-2.5942514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3667592) q[0];
sx q[0];
rz(-1.1792553) q[0];
sx q[0];
rz(-1.1460079) q[0];
rz(0.92578069) q[1];
sx q[1];
rz(-1.0113167) q[1];
sx q[1];
rz(2.944223) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5615047) q[0];
sx q[0];
rz(-2.5338123) q[0];
sx q[0];
rz(1.32919) q[0];
x q[1];
rz(2.8497265) q[2];
sx q[2];
rz(-1.2780683) q[2];
sx q[2];
rz(1.1605934) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6655953) q[1];
sx q[1];
rz(-1.531812) q[1];
sx q[1];
rz(-1.798756) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8160964) q[3];
sx q[3];
rz(-0.56762689) q[3];
sx q[3];
rz(-2.5095255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7388514) q[2];
sx q[2];
rz(-2.1612942) q[2];
sx q[2];
rz(-0.19035467) q[2];
rz(-0.061577408) q[3];
sx q[3];
rz(-0.96934167) q[3];
sx q[3];
rz(2.0891345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9848118) q[0];
sx q[0];
rz(-1.7233912) q[0];
sx q[0];
rz(2.5132827) q[0];
rz(-1.5208987) q[1];
sx q[1];
rz(-1.9338806) q[1];
sx q[1];
rz(0.77888387) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.223004) q[0];
sx q[0];
rz(-2.1445334) q[0];
sx q[0];
rz(2.6911435) q[0];
rz(-pi) q[1];
rz(2.4905048) q[2];
sx q[2];
rz(-2.6991803) q[2];
sx q[2];
rz(0.96786066) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0592093) q[1];
sx q[1];
rz(-1.6184855) q[1];
sx q[1];
rz(1.5573182) q[1];
rz(-pi) q[2];
rz(0.811399) q[3];
sx q[3];
rz(-0.86199844) q[3];
sx q[3];
rz(-0.78042316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4817619) q[2];
sx q[2];
rz(-1.0783106) q[2];
sx q[2];
rz(-4.3241186e-05) q[2];
rz(1.0846042) q[3];
sx q[3];
rz(-1.1559887) q[3];
sx q[3];
rz(-2.675975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67149177) q[0];
sx q[0];
rz(-1.4441613) q[0];
sx q[0];
rz(-1.8430365) q[0];
rz(0.57688722) q[1];
sx q[1];
rz(-1.8678317) q[1];
sx q[1];
rz(0.44688046) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80432207) q[0];
sx q[0];
rz(-1.7511586) q[0];
sx q[0];
rz(2.4867663) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8088008) q[2];
sx q[2];
rz(-2.8341475) q[2];
sx q[2];
rz(0.36556903) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2858823) q[1];
sx q[1];
rz(-1.9706763) q[1];
sx q[1];
rz(0.22004308) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6303667) q[3];
sx q[3];
rz(-1.3856264) q[3];
sx q[3];
rz(-0.51051729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.4345066) q[2];
sx q[2];
rz(-2.862317) q[2];
sx q[2];
rz(0.2198098) q[2];
rz(1.6541727) q[3];
sx q[3];
rz(-1.4498962) q[3];
sx q[3];
rz(2.4421104) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.952482) q[0];
sx q[0];
rz(-0.43287745) q[0];
sx q[0];
rz(-2.2440946) q[0];
rz(-1.3709566) q[1];
sx q[1];
rz(-1.0400583) q[1];
sx q[1];
rz(0.8955566) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8093256) q[0];
sx q[0];
rz(-2.4054962) q[0];
sx q[0];
rz(2.9311137) q[0];
rz(-pi) q[1];
x q[1];
rz(0.61195749) q[2];
sx q[2];
rz(-2.0921807) q[2];
sx q[2];
rz(0.88980161) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.59437737) q[1];
sx q[1];
rz(-1.9037694) q[1];
sx q[1];
rz(-0.39526387) q[1];
x q[2];
rz(-0.79781161) q[3];
sx q[3];
rz(-0.63029248) q[3];
sx q[3];
rz(0.52784656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.36876496) q[2];
sx q[2];
rz(-2.0389098) q[2];
sx q[2];
rz(-2.9844798) q[2];
rz(-2.469192) q[3];
sx q[3];
rz(-1.7417358) q[3];
sx q[3];
rz(0.11252832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.86033487) q[0];
sx q[0];
rz(-0.66547886) q[0];
sx q[0];
rz(1.459664) q[0];
rz(-2.7809987) q[1];
sx q[1];
rz(-1.4692042) q[1];
sx q[1];
rz(-1.2999387) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0361745) q[0];
sx q[0];
rz(-1.4486533) q[0];
sx q[0];
rz(2.7226177) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1088013) q[2];
sx q[2];
rz(-2.1401775) q[2];
sx q[2];
rz(1.5608567) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9323401) q[1];
sx q[1];
rz(-0.31267088) q[1];
sx q[1];
rz(0.31950848) q[1];
rz(-pi) q[2];
rz(-2.594527) q[3];
sx q[3];
rz(-0.69098847) q[3];
sx q[3];
rz(-0.83741659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2099057) q[2];
sx q[2];
rz(-1.0301215) q[2];
sx q[2];
rz(-0.21200655) q[2];
rz(2.0906406) q[3];
sx q[3];
rz(-1.7651599) q[3];
sx q[3];
rz(-3.0534548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(1.5365005) q[0];
sx q[0];
rz(-0.78992805) q[0];
sx q[0];
rz(2.9686046) q[0];
rz(2.0269003) q[1];
sx q[1];
rz(-2.098691) q[1];
sx q[1];
rz(3.0317422) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82634559) q[0];
sx q[0];
rz(-1.3777121) q[0];
sx q[0];
rz(-0.56435926) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82288701) q[2];
sx q[2];
rz(-1.7639065) q[2];
sx q[2];
rz(2.2597974) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5067271) q[1];
sx q[1];
rz(-2.8322729) q[1];
sx q[1];
rz(0.97815467) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19211003) q[3];
sx q[3];
rz(-0.56280901) q[3];
sx q[3];
rz(0.61432381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4032119) q[2];
sx q[2];
rz(-1.2028368) q[2];
sx q[2];
rz(-0.66486764) q[2];
rz(0.46857771) q[3];
sx q[3];
rz(-1.5237619) q[3];
sx q[3];
rz(2.328228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.453603) q[0];
sx q[0];
rz(-0.60992321) q[0];
sx q[0];
rz(-2.4000121) q[0];
rz(1.1874416) q[1];
sx q[1];
rz(-1.5267742) q[1];
sx q[1];
rz(2.5724519) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3305107) q[0];
sx q[0];
rz(-1.2449322) q[0];
sx q[0];
rz(-1.261607) q[0];
rz(-pi) q[1];
rz(1.6422317) q[2];
sx q[2];
rz(-2.0234081) q[2];
sx q[2];
rz(0.49734383) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4148767) q[1];
sx q[1];
rz(-1.1099713) q[1];
sx q[1];
rz(2.0445092) q[1];
x q[2];
rz(2.3027116) q[3];
sx q[3];
rz(-1.0195707) q[3];
sx q[3];
rz(-1.4428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6471214) q[2];
sx q[2];
rz(-2.1529866) q[2];
sx q[2];
rz(-0.76623255) q[2];
rz(1.9689485) q[3];
sx q[3];
rz(-1.5394883) q[3];
sx q[3];
rz(-2.9197781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2072993) q[0];
sx q[0];
rz(-0.3485637) q[0];
sx q[0];
rz(0.19185129) q[0];
rz(-0.81709298) q[1];
sx q[1];
rz(-2.8849738) q[1];
sx q[1];
rz(-2.4339035) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1908042) q[0];
sx q[0];
rz(-0.30256264) q[0];
sx q[0];
rz(-0.62015073) q[0];
rz(-2.8079671) q[2];
sx q[2];
rz(-1.0051491) q[2];
sx q[2];
rz(0.27496342) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0470815) q[1];
sx q[1];
rz(-0.99267497) q[1];
sx q[1];
rz(-2.6972527) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1698506) q[3];
sx q[3];
rz(-1.069456) q[3];
sx q[3];
rz(-2.3399692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2231458) q[2];
sx q[2];
rz(-0.94474363) q[2];
sx q[2];
rz(-2.4165912) q[2];
rz(-0.21507344) q[3];
sx q[3];
rz(-0.30078617) q[3];
sx q[3];
rz(0.71896368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9216777) q[0];
sx q[0];
rz(-0.81659962) q[0];
sx q[0];
rz(0.30246977) q[0];
rz(1.1014145) q[1];
sx q[1];
rz(-1.9093724) q[1];
sx q[1];
rz(-2.2473635) q[1];
rz(-2.106059) q[2];
sx q[2];
rz(-1.6720094) q[2];
sx q[2];
rz(-1.592247) q[2];
rz(3.0871747) q[3];
sx q[3];
rz(-2.6734753) q[3];
sx q[3];
rz(-2.5772167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
