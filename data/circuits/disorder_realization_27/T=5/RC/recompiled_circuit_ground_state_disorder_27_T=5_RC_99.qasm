OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7967427) q[0];
sx q[0];
rz(-2.8673708) q[0];
sx q[0];
rz(2.5728777) q[0];
rz(1.2110127) q[1];
sx q[1];
rz(-2.14415) q[1];
sx q[1];
rz(-0.2675736) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28963213) q[0];
sx q[0];
rz(-1.3697764) q[0];
sx q[0];
rz(2.5876849) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2368343) q[2];
sx q[2];
rz(-1.4638454) q[2];
sx q[2];
rz(3.1283875) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2530816) q[1];
sx q[1];
rz(-3.1329229) q[1];
sx q[1];
rz(1.289283) q[1];
x q[2];
rz(1.9851793) q[3];
sx q[3];
rz(-0.15358812) q[3];
sx q[3];
rz(1.6382662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.25207511) q[2];
sx q[2];
rz(-1.0675665) q[2];
sx q[2];
rz(2.0102823) q[2];
rz(-0.45025292) q[3];
sx q[3];
rz(-0.69142747) q[3];
sx q[3];
rz(0.94436193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4985519) q[0];
sx q[0];
rz(-2.2096071) q[0];
sx q[0];
rz(1.407628) q[0];
rz(-2.6990926) q[1];
sx q[1];
rz(-1.4235556) q[1];
sx q[1];
rz(0.59534591) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.567826) q[0];
sx q[0];
rz(-1.8086595) q[0];
sx q[0];
rz(-0.25517558) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48974719) q[2];
sx q[2];
rz(-2.1648266) q[2];
sx q[2];
rz(-0.56383946) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7253984) q[1];
sx q[1];
rz(-0.42896118) q[1];
sx q[1];
rz(-3.129175) q[1];
x q[2];
rz(0.61025783) q[3];
sx q[3];
rz(-1.9034981) q[3];
sx q[3];
rz(-1.4833029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.88196102) q[2];
sx q[2];
rz(-0.61820784) q[2];
sx q[2];
rz(0.76914966) q[2];
rz(-1.361557) q[3];
sx q[3];
rz(-0.96746126) q[3];
sx q[3];
rz(2.5462525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24460569) q[0];
sx q[0];
rz(-0.91490442) q[0];
sx q[0];
rz(-0.33682522) q[0];
rz(1.4312076) q[1];
sx q[1];
rz(-0.84588784) q[1];
sx q[1];
rz(3.0444042) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1203994) q[0];
sx q[0];
rz(-2.057907) q[0];
sx q[0];
rz(-0.22986408) q[0];
x q[1];
rz(-2.7601542) q[2];
sx q[2];
rz(-1.6614011) q[2];
sx q[2];
rz(-2.7862797) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5640024) q[1];
sx q[1];
rz(-1.2256943) q[1];
sx q[1];
rz(-1.4240828) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52690701) q[3];
sx q[3];
rz(-2.1873173) q[3];
sx q[3];
rz(-1.34672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1691957) q[2];
sx q[2];
rz(-1.7652067) q[2];
sx q[2];
rz(2.3248559) q[2];
rz(-0.9225325) q[3];
sx q[3];
rz(-0.43729344) q[3];
sx q[3];
rz(1.9492662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.8836477) q[0];
sx q[0];
rz(-1.8035996) q[0];
sx q[0];
rz(2.2747967) q[0];
rz(-1.2358933) q[1];
sx q[1];
rz(-2.0265323) q[1];
sx q[1];
rz(2.8996276) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4075027) q[0];
sx q[0];
rz(-2.3781812) q[0];
sx q[0];
rz(-1.9620738) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4093269) q[2];
sx q[2];
rz(-1.4982371) q[2];
sx q[2];
rz(-0.41997318) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.45631187) q[1];
sx q[1];
rz(-1.1667098) q[1];
sx q[1];
rz(-1.4220974) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5210152) q[3];
sx q[3];
rz(-1.0743273) q[3];
sx q[3];
rz(0.97126006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8786826) q[2];
sx q[2];
rz(-1.5107369) q[2];
sx q[2];
rz(2.6061457) q[2];
rz(2.6483436) q[3];
sx q[3];
rz(-2.2507164) q[3];
sx q[3];
rz(2.6935327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.066147476) q[0];
sx q[0];
rz(-2.9904521) q[0];
sx q[0];
rz(-2.2976663) q[0];
rz(-1.6663724) q[1];
sx q[1];
rz(-1.4862783) q[1];
sx q[1];
rz(-2.7484238) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.524784) q[0];
sx q[0];
rz(-0.83833414) q[0];
sx q[0];
rz(-2.7969267) q[0];
x q[1];
rz(2.8507502) q[2];
sx q[2];
rz(-0.65724361) q[2];
sx q[2];
rz(2.2387981) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1816493) q[1];
sx q[1];
rz(-1.8487329) q[1];
sx q[1];
rz(-0.47680579) q[1];
rz(0.78759463) q[3];
sx q[3];
rz(-1.2114085) q[3];
sx q[3];
rz(-1.6000634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7273442) q[2];
sx q[2];
rz(-2.9372637) q[2];
sx q[2];
rz(2.7434529) q[2];
rz(2.5967755) q[3];
sx q[3];
rz(-2.3530493) q[3];
sx q[3];
rz(-1.8383693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9141465) q[0];
sx q[0];
rz(-2.1307724) q[0];
sx q[0];
rz(0.8557125) q[0];
rz(2.4489467) q[1];
sx q[1];
rz(-0.99450642) q[1];
sx q[1];
rz(2.8725502) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0318258) q[0];
sx q[0];
rz(-2.9485011) q[0];
sx q[0];
rz(-1.4433799) q[0];
rz(-pi) q[1];
rz(1.5085717) q[2];
sx q[2];
rz(-0.50743689) q[2];
sx q[2];
rz(-2.1560046) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0768041) q[1];
sx q[1];
rz(-0.85569438) q[1];
sx q[1];
rz(-2.8984757) q[1];
x q[2];
rz(-1.3790491) q[3];
sx q[3];
rz(-1.7320398) q[3];
sx q[3];
rz(1.92056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0272224) q[2];
sx q[2];
rz(-1.8222858) q[2];
sx q[2];
rz(-0.11030062) q[2];
rz(-0.86841622) q[3];
sx q[3];
rz(-1.9649558) q[3];
sx q[3];
rz(-1.7051914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7480943) q[0];
sx q[0];
rz(-0.49999923) q[0];
sx q[0];
rz(-2.2802343) q[0];
rz(-1.512108) q[1];
sx q[1];
rz(-1.6810828) q[1];
sx q[1];
rz(-0.82383627) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.829946) q[0];
sx q[0];
rz(-2.5426546) q[0];
sx q[0];
rz(1.7700559) q[0];
rz(-pi) q[1];
rz(-2.3917213) q[2];
sx q[2];
rz(-1.2629384) q[2];
sx q[2];
rz(-1.4726382) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8331404) q[1];
sx q[1];
rz(-0.96953934) q[1];
sx q[1];
rz(-1.8300301) q[1];
rz(2.4895913) q[3];
sx q[3];
rz(-1.3175512) q[3];
sx q[3];
rz(0.30919231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5065826) q[2];
sx q[2];
rz(-2.625605) q[2];
sx q[2];
rz(0.79279509) q[2];
rz(3.1363764) q[3];
sx q[3];
rz(-0.79013932) q[3];
sx q[3];
rz(-1.691157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1714627) q[0];
sx q[0];
rz(-1.1928394) q[0];
sx q[0];
rz(-0.18950732) q[0];
rz(0.7729404) q[1];
sx q[1];
rz(-0.49630061) q[1];
sx q[1];
rz(-2.5453087) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4204105) q[0];
sx q[0];
rz(-0.89119688) q[0];
sx q[0];
rz(-1.7880938) q[0];
rz(-pi) q[1];
rz(-2.2484915) q[2];
sx q[2];
rz(-0.89701954) q[2];
sx q[2];
rz(2.386415) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.90702) q[1];
sx q[1];
rz(-1.4661745) q[1];
sx q[1];
rz(1.6638882) q[1];
x q[2];
rz(-2.1607481) q[3];
sx q[3];
rz(-2.337237) q[3];
sx q[3];
rz(0.70589069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4160055) q[2];
sx q[2];
rz(-3.1276939) q[2];
sx q[2];
rz(1.2131946) q[2];
rz(0.98617918) q[3];
sx q[3];
rz(-1.7257907) q[3];
sx q[3];
rz(1.2590316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1213433) q[0];
sx q[0];
rz(-1.5859402) q[0];
sx q[0];
rz(-2.0315309) q[0];
rz(-0.32866651) q[1];
sx q[1];
rz(-1.5549436) q[1];
sx q[1];
rz(-1.2967671) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9220306) q[0];
sx q[0];
rz(-0.91518213) q[0];
sx q[0];
rz(-3.0232885) q[0];
rz(-pi) q[1];
rz(1.1743714) q[2];
sx q[2];
rz(-1.6280481) q[2];
sx q[2];
rz(-0.020584189) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3039141) q[1];
sx q[1];
rz(-1.078036) q[1];
sx q[1];
rz(1.4999092) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27196692) q[3];
sx q[3];
rz(-0.68373954) q[3];
sx q[3];
rz(0.34599397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0686331) q[2];
sx q[2];
rz(-0.98560846) q[2];
sx q[2];
rz(-0.10406058) q[2];
rz(-1.1348628) q[3];
sx q[3];
rz(-1.3645423) q[3];
sx q[3];
rz(0.22741905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65570152) q[0];
sx q[0];
rz(-0.98485297) q[0];
sx q[0];
rz(-2.4110598) q[0];
rz(2.5841374) q[1];
sx q[1];
rz(-1.1703706) q[1];
sx q[1];
rz(2.7117859) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24048478) q[0];
sx q[0];
rz(-0.76178023) q[0];
sx q[0];
rz(-0.00073379993) q[0];
rz(-2.7676959) q[2];
sx q[2];
rz(-1.3329525) q[2];
sx q[2];
rz(1.4542945) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3278663) q[1];
sx q[1];
rz(-1.9414492) q[1];
sx q[1];
rz(-0.085538105) q[1];
x q[2];
rz(1.2005931) q[3];
sx q[3];
rz(-1.4552081) q[3];
sx q[3];
rz(-0.59578958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80959117) q[2];
sx q[2];
rz(-1.0217228) q[2];
sx q[2];
rz(0.29279718) q[2];
rz(3.0012567) q[3];
sx q[3];
rz(-2.1122746) q[3];
sx q[3];
rz(-0.50104195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3705227) q[0];
sx q[0];
rz(-1.6765544) q[0];
sx q[0];
rz(0.21677207) q[0];
rz(2.4222005) q[1];
sx q[1];
rz(-1.8665301) q[1];
sx q[1];
rz(-2.9449609) q[1];
rz(-1.0631845) q[2];
sx q[2];
rz(-1.7287935) q[2];
sx q[2];
rz(-2.4274735) q[2];
rz(1.2483531) q[3];
sx q[3];
rz(-1.4436246) q[3];
sx q[3];
rz(1.1818813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
