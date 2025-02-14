OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.1640373) q[0];
sx q[0];
rz(-1.9174175) q[0];
sx q[0];
rz(-1.2807711) q[0];
rz(-0.7723074) q[1];
sx q[1];
rz(-0.63900715) q[1];
sx q[1];
rz(-2.8746936) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7909662) q[0];
sx q[0];
rz(-1.4746951) q[0];
sx q[0];
rz(-3.0477357) q[0];
rz(1.2615777) q[2];
sx q[2];
rz(-2.4291933) q[2];
sx q[2];
rz(1.0605896) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7952629) q[1];
sx q[1];
rz(-0.85546934) q[1];
sx q[1];
rz(-0.99154179) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4679174) q[3];
sx q[3];
rz(-2.0112546) q[3];
sx q[3];
rz(2.6552424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12307564) q[2];
sx q[2];
rz(-2.1776431) q[2];
sx q[2];
rz(0.71259552) q[2];
rz(-0.16768843) q[3];
sx q[3];
rz(-0.84961397) q[3];
sx q[3];
rz(0.4445506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0406822) q[0];
sx q[0];
rz(-1.3586783) q[0];
sx q[0];
rz(-1.4605301) q[0];
rz(-1.679861) q[1];
sx q[1];
rz(-1.0667543) q[1];
sx q[1];
rz(2.0251958) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24943251) q[0];
sx q[0];
rz(-1.6139702) q[0];
sx q[0];
rz(0.024726111) q[0];
x q[1];
rz(-1.8113715) q[2];
sx q[2];
rz(-0.5690388) q[2];
sx q[2];
rz(2.9441058) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.18552407) q[1];
sx q[1];
rz(-1.9466234) q[1];
sx q[1];
rz(1.8517428) q[1];
x q[2];
rz(1.7061911) q[3];
sx q[3];
rz(-1.97768) q[3];
sx q[3];
rz(-3.1171006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1560893) q[2];
sx q[2];
rz(-2.6458793) q[2];
sx q[2];
rz(2.8893341) q[2];
rz(1.0559399) q[3];
sx q[3];
rz(-2.2838433) q[3];
sx q[3];
rz(2.1411538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58241874) q[0];
sx q[0];
rz(-0.90455872) q[0];
sx q[0];
rz(-2.3140123) q[0];
rz(-2.6170392) q[1];
sx q[1];
rz(-1.2453715) q[1];
sx q[1];
rz(2.1771199) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8480769) q[0];
sx q[0];
rz(-1.6381253) q[0];
sx q[0];
rz(-2.1490554) q[0];
rz(-pi) q[1];
rz(1.7894423) q[2];
sx q[2];
rz(-2.7579569) q[2];
sx q[2];
rz(-0.84886692) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.8207404) q[1];
sx q[1];
rz(-0.94106442) q[1];
sx q[1];
rz(1.6321502) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12196864) q[3];
sx q[3];
rz(-0.44443529) q[3];
sx q[3];
rz(2.0214863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7437637) q[2];
sx q[2];
rz(-2.9167794) q[2];
sx q[2];
rz(1.1986097) q[2];
rz(1.4926636) q[3];
sx q[3];
rz(-2.1677833) q[3];
sx q[3];
rz(-2.5627513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.32597932) q[0];
sx q[0];
rz(-0.15376832) q[0];
sx q[0];
rz(-1.2157259) q[0];
rz(2.1513596) q[1];
sx q[1];
rz(-2.4001887) q[1];
sx q[1];
rz(-0.56891099) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.617482) q[0];
sx q[0];
rz(-2.4883399) q[0];
sx q[0];
rz(-1.0934483) q[0];
rz(-pi) q[1];
rz(1.9401476) q[2];
sx q[2];
rz(-1.8180038) q[2];
sx q[2];
rz(1.2170685) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7987631) q[1];
sx q[1];
rz(-2.0137278) q[1];
sx q[1];
rz(-2.8196067) q[1];
x q[2];
rz(-0.91671555) q[3];
sx q[3];
rz(-0.93610969) q[3];
sx q[3];
rz(-0.095528729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8625921) q[2];
sx q[2];
rz(-2.180884) q[2];
sx q[2];
rz(2.847239) q[2];
rz(-0.66458464) q[3];
sx q[3];
rz(-0.75988257) q[3];
sx q[3];
rz(1.6871281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6890474) q[0];
sx q[0];
rz(-1.713546) q[0];
sx q[0];
rz(-2.8277165) q[0];
rz(2.2398056) q[1];
sx q[1];
rz(-1.3548464) q[1];
sx q[1];
rz(1.4656167) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2814502) q[0];
sx q[0];
rz(-1.6897795) q[0];
sx q[0];
rz(3.1156179) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1526665) q[2];
sx q[2];
rz(-0.2731495) q[2];
sx q[2];
rz(1.404431) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.213433) q[1];
sx q[1];
rz(-1.4959236) q[1];
sx q[1];
rz(2.7533745) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0907822) q[3];
sx q[3];
rz(-0.63119315) q[3];
sx q[3];
rz(2.7630591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.19447154) q[2];
sx q[2];
rz(-2.3919545) q[2];
sx q[2];
rz(1.0844024) q[2];
rz(-2.8179152) q[3];
sx q[3];
rz(-0.71666986) q[3];
sx q[3];
rz(-2.9173541) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1502458) q[0];
sx q[0];
rz(-0.24557376) q[0];
sx q[0];
rz(-0.41159758) q[0];
rz(-2.4161074) q[1];
sx q[1];
rz(-1.2440224) q[1];
sx q[1];
rz(1.4010319) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0822015) q[0];
sx q[0];
rz(-2.2103346) q[0];
sx q[0];
rz(1.8138422) q[0];
rz(3.0040222) q[2];
sx q[2];
rz(-0.74512945) q[2];
sx q[2];
rz(-0.18411769) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7159375) q[1];
sx q[1];
rz(-1.9129941) q[1];
sx q[1];
rz(-1.9075431) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9855491) q[3];
sx q[3];
rz(-2.6855009) q[3];
sx q[3];
rz(1.2233011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1650042) q[2];
sx q[2];
rz(-0.62825957) q[2];
sx q[2];
rz(1.4264301) q[2];
rz(-0.12380883) q[3];
sx q[3];
rz(-2.0721469) q[3];
sx q[3];
rz(2.6482705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8641149) q[0];
sx q[0];
rz(-1.9954229) q[0];
sx q[0];
rz(1.3364828) q[0];
rz(-2.2026964) q[1];
sx q[1];
rz(-1.122033) q[1];
sx q[1];
rz(2.8575361) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2920609) q[0];
sx q[0];
rz(-2.1725328) q[0];
sx q[0];
rz(-1.1547778) q[0];
rz(-1.797126) q[2];
sx q[2];
rz(-1.411419) q[2];
sx q[2];
rz(-1.2633737) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3160243) q[1];
sx q[1];
rz(-3.1085759) q[1];
sx q[1];
rz(-1.2570639) q[1];
x q[2];
rz(1.7098996) q[3];
sx q[3];
rz(-1.749732) q[3];
sx q[3];
rz(0.034497189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.17860086) q[2];
sx q[2];
rz(-2.4602349) q[2];
sx q[2];
rz(-0.94092384) q[2];
rz(3.1407147) q[3];
sx q[3];
rz(-1.9160756) q[3];
sx q[3];
rz(3.0555449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0451999) q[0];
sx q[0];
rz(-2.2470076) q[0];
sx q[0];
rz(-2.9659502) q[0];
rz(-1.9215709) q[1];
sx q[1];
rz(-2.2717387) q[1];
sx q[1];
rz(-0.87179914) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5822495) q[0];
sx q[0];
rz(-1.7212369) q[0];
sx q[0];
rz(1.3036672) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7824509) q[2];
sx q[2];
rz(-1.3030714) q[2];
sx q[2];
rz(-2.8504782) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3094191) q[1];
sx q[1];
rz(-1.7298352) q[1];
sx q[1];
rz(-0.60950233) q[1];
rz(-2.5848735) q[3];
sx q[3];
rz(-1.5079323) q[3];
sx q[3];
rz(-1.5134461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4153727) q[2];
sx q[2];
rz(-1.8612334) q[2];
sx q[2];
rz(-1.0225164) q[2];
rz(0.38698777) q[3];
sx q[3];
rz(-1.3849247) q[3];
sx q[3];
rz(-2.3302087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.33685327) q[0];
sx q[0];
rz(-1.9751208) q[0];
sx q[0];
rz(2.8785896) q[0];
rz(-0.31245843) q[1];
sx q[1];
rz(-1.0954906) q[1];
sx q[1];
rz(1.1824664) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9179494) q[0];
sx q[0];
rz(-1.765476) q[0];
sx q[0];
rz(-2.0854092) q[0];
rz(-pi) q[1];
rz(2.5850641) q[2];
sx q[2];
rz(-0.91224837) q[2];
sx q[2];
rz(1.6375033) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.39315614) q[1];
sx q[1];
rz(-1.9416205) q[1];
sx q[1];
rz(-1.7840339) q[1];
rz(-pi) q[2];
rz(0.5236756) q[3];
sx q[3];
rz(-2.0770235) q[3];
sx q[3];
rz(0.62577137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.075228127) q[2];
sx q[2];
rz(-1.790739) q[2];
sx q[2];
rz(-0.24825516) q[2];
rz(-0.17812854) q[3];
sx q[3];
rz(-1.0823366) q[3];
sx q[3];
rz(2.3586912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42661509) q[0];
sx q[0];
rz(-0.39418945) q[0];
sx q[0];
rz(-2.07975) q[0];
rz(-1.8408403) q[1];
sx q[1];
rz(-1.6164833) q[1];
sx q[1];
rz(2.010407) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12515629) q[0];
sx q[0];
rz(-1.1068692) q[0];
sx q[0];
rz(-0.50986503) q[0];
rz(0.71604095) q[2];
sx q[2];
rz(-2.5114905) q[2];
sx q[2];
rz(1.014647) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6415875) q[1];
sx q[1];
rz(-1.4785071) q[1];
sx q[1];
rz(-0.39938836) q[1];
rz(-2.5401073) q[3];
sx q[3];
rz(-1.3772528) q[3];
sx q[3];
rz(2.7880965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0608369) q[2];
sx q[2];
rz(-2.8751825) q[2];
sx q[2];
rz(-2.5401435) q[2];
rz(2.1001749) q[3];
sx q[3];
rz(-1.6626549) q[3];
sx q[3];
rz(-1.7634332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1489442) q[0];
sx q[0];
rz(-0.61275488) q[0];
sx q[0];
rz(-2.3339094) q[0];
rz(0.07769892) q[1];
sx q[1];
rz(-0.54010375) q[1];
sx q[1];
rz(1.7355951) q[1];
rz(-0.40602691) q[2];
sx q[2];
rz(-2.1320504) q[2];
sx q[2];
rz(-0.8036094) q[2];
rz(0.79176767) q[3];
sx q[3];
rz(-2.0590012) q[3];
sx q[3];
rz(1.7760359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
