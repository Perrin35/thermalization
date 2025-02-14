OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34484997) q[0];
sx q[0];
rz(-0.27422187) q[0];
sx q[0];
rz(-2.5728777) q[0];
rz(1.2110127) q[1];
sx q[1];
rz(-2.14415) q[1];
sx q[1];
rz(2.8740191) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8519605) q[0];
sx q[0];
rz(-1.3697764) q[0];
sx q[0];
rz(0.55390771) q[0];
x q[1];
rz(1.9047584) q[2];
sx q[2];
rz(-1.4638454) q[2];
sx q[2];
rz(-0.013205139) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5423746) q[1];
sx q[1];
rz(-1.5683878) q[1];
sx q[1];
rz(1.5791248) q[1];
rz(-1.4300284) q[3];
sx q[3];
rz(-1.632431) q[3];
sx q[3];
rz(0.34256645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8895175) q[2];
sx q[2];
rz(-1.0675665) q[2];
sx q[2];
rz(2.0102823) q[2];
rz(-2.6913397) q[3];
sx q[3];
rz(-2.4501652) q[3];
sx q[3];
rz(-2.1972307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4985519) q[0];
sx q[0];
rz(-2.2096071) q[0];
sx q[0];
rz(-1.7339647) q[0];
rz(-2.6990926) q[1];
sx q[1];
rz(-1.4235556) q[1];
sx q[1];
rz(0.59534591) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93564088) q[0];
sx q[0];
rz(-1.3229587) q[0];
sx q[0];
rz(1.3252844) q[0];
rz(2.2240665) q[2];
sx q[2];
rz(-1.9712312) q[2];
sx q[2];
rz(-1.296907) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.40253958) q[1];
sx q[1];
rz(-1.1418704) q[1];
sx q[1];
rz(1.5651171) q[1];
x q[2];
rz(-1.1717623) q[3];
sx q[3];
rz(-2.1431987) q[3];
sx q[3];
rz(-0.31203285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2596316) q[2];
sx q[2];
rz(-0.61820784) q[2];
sx q[2];
rz(0.76914966) q[2];
rz(-1.7800356) q[3];
sx q[3];
rz(-0.96746126) q[3];
sx q[3];
rz(-2.5462525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24460569) q[0];
sx q[0];
rz(-2.2266882) q[0];
sx q[0];
rz(-2.8047674) q[0];
rz(1.4312076) q[1];
sx q[1];
rz(-0.84588784) q[1];
sx q[1];
rz(3.0444042) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.021193209) q[0];
sx q[0];
rz(-1.0836856) q[0];
sx q[0];
rz(-0.22986408) q[0];
rz(-pi) q[1];
rz(-0.23938208) q[2];
sx q[2];
rz(-2.7500543) q[2];
sx q[2];
rz(1.4373612) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5640024) q[1];
sx q[1];
rz(-1.2256943) q[1];
sx q[1];
rz(-1.4240828) q[1];
x q[2];
rz(2.2575374) q[3];
sx q[3];
rz(-1.1480322) q[3];
sx q[3];
rz(-2.5930269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.97239697) q[2];
sx q[2];
rz(-1.376386) q[2];
sx q[2];
rz(2.3248559) q[2];
rz(-0.9225325) q[3];
sx q[3];
rz(-2.7042992) q[3];
sx q[3];
rz(-1.9492662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8836477) q[0];
sx q[0];
rz(-1.3379931) q[0];
sx q[0];
rz(2.2747967) q[0];
rz(-1.2358933) q[1];
sx q[1];
rz(-2.0265323) q[1];
sx q[1];
rz(-0.24196504) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6886238) q[0];
sx q[0];
rz(-1.8376266) q[0];
sx q[0];
rz(0.84653207) q[0];
x q[1];
rz(0.073512065) q[2];
sx q[2];
rz(-1.7318372) q[2];
sx q[2];
rz(1.1390151) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6852808) q[1];
sx q[1];
rz(-1.1667098) q[1];
sx q[1];
rz(-1.4220974) q[1];
rz(-3.0499942) q[3];
sx q[3];
rz(-0.49875101) q[3];
sx q[3];
rz(2.2745511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8786826) q[2];
sx q[2];
rz(-1.5107369) q[2];
sx q[2];
rz(-0.53544694) q[2];
rz(0.49324909) q[3];
sx q[3];
rz(-0.89087629) q[3];
sx q[3];
rz(2.6935327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0754452) q[0];
sx q[0];
rz(-2.9904521) q[0];
sx q[0];
rz(-2.2976663) q[0];
rz(-1.4752202) q[1];
sx q[1];
rz(-1.6553144) q[1];
sx q[1];
rz(0.39316887) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9520023) q[0];
sx q[0];
rz(-1.31685) q[0];
sx q[0];
rz(-2.3334731) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8507502) q[2];
sx q[2];
rz(-0.65724361) q[2];
sx q[2];
rz(-2.2387981) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6715309) q[1];
sx q[1];
rz(-1.1137149) q[1];
sx q[1];
rz(1.8815329) q[1];
rz(-pi) q[2];
x q[2];
rz(2.353998) q[3];
sx q[3];
rz(-1.9301842) q[3];
sx q[3];
rz(-1.6000634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4142485) q[2];
sx q[2];
rz(-2.9372637) q[2];
sx q[2];
rz(-0.3981398) q[2];
rz(-0.54481715) q[3];
sx q[3];
rz(-0.78854338) q[3];
sx q[3];
rz(1.8383693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9141465) q[0];
sx q[0];
rz(-2.1307724) q[0];
sx q[0];
rz(-0.8557125) q[0];
rz(-0.69264597) q[1];
sx q[1];
rz(-0.99450642) q[1];
sx q[1];
rz(-0.2690424) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5554898) q[0];
sx q[0];
rz(-1.5951831) q[0];
sx q[0];
rz(1.3792319) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5085717) q[2];
sx q[2];
rz(-2.6341558) q[2];
sx q[2];
rz(-0.98558805) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.70316305) q[1];
sx q[1];
rz(-0.7483349) q[1];
sx q[1];
rz(1.3003883) q[1];
rz(-0.16420047) q[3];
sx q[3];
rz(-1.3815666) q[3];
sx q[3];
rz(2.7606719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0272224) q[2];
sx q[2];
rz(-1.8222858) q[2];
sx q[2];
rz(-3.031292) q[2];
rz(0.86841622) q[3];
sx q[3];
rz(-1.9649558) q[3];
sx q[3];
rz(1.7051914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39349839) q[0];
sx q[0];
rz(-2.6415934) q[0];
sx q[0];
rz(-2.2802343) q[0];
rz(-1.6294847) q[1];
sx q[1];
rz(-1.4605099) q[1];
sx q[1];
rz(-0.82383627) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0938823) q[0];
sx q[0];
rz(-1.4589696) q[0];
sx q[0];
rz(0.98112962) q[0];
x q[1];
rz(-2.3917213) q[2];
sx q[2];
rz(-1.2629384) q[2];
sx q[2];
rz(-1.4726382) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8331404) q[1];
sx q[1];
rz(-0.96953934) q[1];
sx q[1];
rz(-1.3115626) q[1];
rz(0.40315513) q[3];
sx q[3];
rz(-0.69270999) q[3];
sx q[3];
rz(1.5787391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6350101) q[2];
sx q[2];
rz(-0.51598769) q[2];
sx q[2];
rz(0.79279509) q[2];
rz(0.005216287) q[3];
sx q[3];
rz(-0.79013932) q[3];
sx q[3];
rz(1.691157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.97013) q[0];
sx q[0];
rz(-1.9487533) q[0];
sx q[0];
rz(0.18950732) q[0];
rz(2.3686523) q[1];
sx q[1];
rz(-2.645292) q[1];
sx q[1];
rz(0.596284) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7211821) q[0];
sx q[0];
rz(-0.89119688) q[0];
sx q[0];
rz(-1.3534989) q[0];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3459732) q[1];
sx q[1];
rz(-1.6633777) q[1];
sx q[1];
rz(0.10507344) q[1];
rz(-pi) q[2];
rz(-0.52395405) q[3];
sx q[3];
rz(-0.92902459) q[3];
sx q[3];
rz(-0.061836035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.72558713) q[2];
sx q[2];
rz(-0.013898762) q[2];
sx q[2];
rz(1.2131946) q[2];
rz(-2.1554135) q[3];
sx q[3];
rz(-1.7257907) q[3];
sx q[3];
rz(-1.8825611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1213433) q[0];
sx q[0];
rz(-1.5859402) q[0];
sx q[0];
rz(1.1100618) q[0];
rz(-2.8129261) q[1];
sx q[1];
rz(-1.5549436) q[1];
sx q[1];
rz(1.2967671) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4121108) q[0];
sx q[0];
rz(-0.66464948) q[0];
sx q[0];
rz(-1.4185216) q[0];
x q[1];
rz(-1.1743714) q[2];
sx q[2];
rz(-1.6280481) q[2];
sx q[2];
rz(-3.1210085) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.98666603) q[1];
sx q[1];
rz(-0.49741751) q[1];
sx q[1];
rz(3.0104396) q[1];
x q[2];
rz(1.3553008) q[3];
sx q[3];
rz(-0.91665506) q[3];
sx q[3];
rz(-0.69132016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0729596) q[2];
sx q[2];
rz(-2.1559842) q[2];
sx q[2];
rz(-3.0375321) q[2];
rz(1.1348628) q[3];
sx q[3];
rz(-1.7770504) q[3];
sx q[3];
rz(-2.9141736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.4858911) q[0];
sx q[0];
rz(-0.98485297) q[0];
sx q[0];
rz(2.4110598) q[0];
rz(-2.5841374) q[1];
sx q[1];
rz(-1.1703706) q[1];
sx q[1];
rz(0.4298068) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24048478) q[0];
sx q[0];
rz(-2.3798124) q[0];
sx q[0];
rz(0.00073379993) q[0];
x q[1];
rz(2.5556106) q[2];
sx q[2];
rz(-0.44011099) q[2];
sx q[2];
rz(2.7173079) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.046172) q[1];
sx q[1];
rz(-2.7616427) q[1];
sx q[1];
rz(-1.7871961) q[1];
x q[2];
rz(1.9409995) q[3];
sx q[3];
rz(-1.6863846) q[3];
sx q[3];
rz(2.5458031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3320015) q[2];
sx q[2];
rz(-1.0217228) q[2];
sx q[2];
rz(2.8487955) q[2];
rz(-0.14033595) q[3];
sx q[3];
rz(-1.0293181) q[3];
sx q[3];
rz(0.50104195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.77107) q[0];
sx q[0];
rz(-1.4650383) q[0];
sx q[0];
rz(-2.9248206) q[0];
rz(-2.4222005) q[1];
sx q[1];
rz(-1.2750625) q[1];
sx q[1];
rz(0.19663179) q[1];
rz(2.9612598) q[2];
sx q[2];
rz(-2.0714932) q[2];
sx q[2];
rz(-0.76938236) q[2];
rz(1.9543129) q[3];
sx q[3];
rz(-2.7957932) q[3];
sx q[3];
rz(2.3898706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
