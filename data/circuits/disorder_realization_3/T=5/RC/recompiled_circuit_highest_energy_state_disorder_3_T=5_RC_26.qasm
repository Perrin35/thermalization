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
rz(0.02778223) q[0];
sx q[0];
rz(-1.2454998) q[0];
sx q[0];
rz(-2.0552638) q[0];
rz(2.0431986) q[1];
sx q[1];
rz(-1.9280704) q[1];
sx q[1];
rz(-1.1462125) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1725464) q[0];
sx q[0];
rz(-1.2451225) q[0];
sx q[0];
rz(2.9708751) q[0];
x q[1];
rz(0.59046794) q[2];
sx q[2];
rz(-1.2256978) q[2];
sx q[2];
rz(-0.1869299) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0806345) q[1];
sx q[1];
rz(-0.74392156) q[1];
sx q[1];
rz(-2.3638636) q[1];
x q[2];
rz(-2.0729154) q[3];
sx q[3];
rz(-1.3716591) q[3];
sx q[3];
rz(2.5728383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22831297) q[2];
sx q[2];
rz(-1.5659119) q[2];
sx q[2];
rz(2.7794465) q[2];
rz(-0.95237887) q[3];
sx q[3];
rz(-0.55914545) q[3];
sx q[3];
rz(2.4366116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97149616) q[0];
sx q[0];
rz(-0.2146475) q[0];
sx q[0];
rz(-1.5317408) q[0];
rz(1.2486628) q[1];
sx q[1];
rz(-1.5208533) q[1];
sx q[1];
rz(-0.85266399) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.06518) q[0];
sx q[0];
rz(-1.914304) q[0];
sx q[0];
rz(-1.0344013) q[0];
rz(-pi) q[1];
x q[1];
rz(0.12262362) q[2];
sx q[2];
rz(-1.2957591) q[2];
sx q[2];
rz(-3.1240535) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.98140796) q[1];
sx q[1];
rz(-1.369919) q[1];
sx q[1];
rz(2.1698879) q[1];
rz(-1.0107105) q[3];
sx q[3];
rz(-1.9017856) q[3];
sx q[3];
rz(0.50332848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2109005) q[2];
sx q[2];
rz(-2.1810668) q[2];
sx q[2];
rz(1.9453913) q[2];
rz(-0.026737468) q[3];
sx q[3];
rz(-2.6452711) q[3];
sx q[3];
rz(-1.4727717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6579599) q[0];
sx q[0];
rz(-0.25068972) q[0];
sx q[0];
rz(-2.2657917) q[0];
rz(2.7818413) q[1];
sx q[1];
rz(-1.3641554) q[1];
sx q[1];
rz(2.13805) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26311603) q[0];
sx q[0];
rz(-0.56385332) q[0];
sx q[0];
rz(3.0082358) q[0];
x q[1];
rz(2.7458756) q[2];
sx q[2];
rz(-0.96683244) q[2];
sx q[2];
rz(0.25091083) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0599614) q[1];
sx q[1];
rz(-1.4278894) q[1];
sx q[1];
rz(-0.68024723) q[1];
rz(-0.64522532) q[3];
sx q[3];
rz(-0.78515437) q[3];
sx q[3];
rz(-0.19715362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.935219) q[2];
sx q[2];
rz(-2.424365) q[2];
sx q[2];
rz(2.3865872) q[2];
rz(0.098091789) q[3];
sx q[3];
rz(-0.57049975) q[3];
sx q[3];
rz(-0.3977972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0574684) q[0];
sx q[0];
rz(-1.5115154) q[0];
sx q[0];
rz(0.64495069) q[0];
rz(-1.0954674) q[1];
sx q[1];
rz(-2.7332833) q[1];
sx q[1];
rz(-2.0300949) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31331691) q[0];
sx q[0];
rz(-3.0544825) q[0];
sx q[0];
rz(-0.75096186) q[0];
rz(-pi) q[1];
rz(0.45504163) q[2];
sx q[2];
rz(-0.18559449) q[2];
sx q[2];
rz(1.1374823) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9605691) q[1];
sx q[1];
rz(-1.7180859) q[1];
sx q[1];
rz(0.21206124) q[1];
rz(-2.4006149) q[3];
sx q[3];
rz(-1.40687) q[3];
sx q[3];
rz(-1.0865097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.633454) q[2];
sx q[2];
rz(-0.78712574) q[2];
sx q[2];
rz(2.656929) q[2];
rz(2.700452) q[3];
sx q[3];
rz(-1.7287247) q[3];
sx q[3];
rz(-1.887656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6345374) q[0];
sx q[0];
rz(-0.84369722) q[0];
sx q[0];
rz(-1.7150568) q[0];
rz(0.85975319) q[1];
sx q[1];
rz(-0.25329241) q[1];
sx q[1];
rz(-2.3111129) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62522307) q[0];
sx q[0];
rz(-2.255356) q[0];
sx q[0];
rz(-0.76098085) q[0];
rz(0.93797063) q[2];
sx q[2];
rz(-1.5245078) q[2];
sx q[2];
rz(-2.7673495) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66248686) q[1];
sx q[1];
rz(-1.431152) q[1];
sx q[1];
rz(0.73151155) q[1];
rz(-pi) q[2];
rz(0.40404218) q[3];
sx q[3];
rz(-0.38876226) q[3];
sx q[3];
rz(-0.59170149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6972203) q[2];
sx q[2];
rz(-2.1264919) q[2];
sx q[2];
rz(-2.3946136) q[2];
rz(-0.30484453) q[3];
sx q[3];
rz(-1.3957142) q[3];
sx q[3];
rz(1.2386809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9656669) q[0];
sx q[0];
rz(-1.8400064) q[0];
sx q[0];
rz(2.9081705) q[0];
rz(2.6790791) q[1];
sx q[1];
rz(-1.1902483) q[1];
sx q[1];
rz(-0.40406427) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3041026) q[0];
sx q[0];
rz(-1.7902052) q[0];
sx q[0];
rz(0.93774937) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42521221) q[2];
sx q[2];
rz(-1.325144) q[2];
sx q[2];
rz(-2.6157891) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0928468) q[1];
sx q[1];
rz(-1.6786075) q[1];
sx q[1];
rz(2.5786933) q[1];
x q[2];
rz(2.8430953) q[3];
sx q[3];
rz(-0.65234557) q[3];
sx q[3];
rz(2.9708095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.81249753) q[2];
sx q[2];
rz(-1.3500682) q[2];
sx q[2];
rz(-2.4883032) q[2];
rz(1.4981883) q[3];
sx q[3];
rz(-2.523962) q[3];
sx q[3];
rz(0.60630715) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6752328) q[0];
sx q[0];
rz(-1.2132069) q[0];
sx q[0];
rz(2.8208222) q[0];
rz(-0.65387154) q[1];
sx q[1];
rz(-2.8120698) q[1];
sx q[1];
rz(3.0810862) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4757547) q[0];
sx q[0];
rz(-1.8403957) q[0];
sx q[0];
rz(2.9682387) q[0];
rz(-pi) q[1];
rz(-0.94580067) q[2];
sx q[2];
rz(-2.201718) q[2];
sx q[2];
rz(-2.6827902) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6937032) q[1];
sx q[1];
rz(-2.1809019) q[1];
sx q[1];
rz(2.542716) q[1];
x q[2];
rz(-2.7869446) q[3];
sx q[3];
rz(-1.9546247) q[3];
sx q[3];
rz(-0.5634748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3939646) q[2];
sx q[2];
rz(-1.4614033) q[2];
sx q[2];
rz(-0.72236577) q[2];
rz(-0.28313053) q[3];
sx q[3];
rz(-0.88909283) q[3];
sx q[3];
rz(2.6319671) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79001456) q[0];
sx q[0];
rz(-0.55792648) q[0];
sx q[0];
rz(2.956692) q[0];
rz(2.6376873) q[1];
sx q[1];
rz(-1.7920707) q[1];
sx q[1];
rz(2.6769743) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7085717) q[0];
sx q[0];
rz(-1.4056217) q[0];
sx q[0];
rz(-0.1689706) q[0];
rz(-pi) q[1];
rz(0.41733317) q[2];
sx q[2];
rz(-1.482903) q[2];
sx q[2];
rz(2.7869239) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.020430001) q[1];
sx q[1];
rz(-1.9489701) q[1];
sx q[1];
rz(0.44831033) q[1];
rz(-pi) q[2];
rz(1.3724324) q[3];
sx q[3];
rz(-2.2806205) q[3];
sx q[3];
rz(1.7539303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.23286143) q[2];
sx q[2];
rz(-1.6249012) q[2];
sx q[2];
rz(-1.0535343) q[2];
rz(0.72707027) q[3];
sx q[3];
rz(-1.2289685) q[3];
sx q[3];
rz(2.5963636) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19186774) q[0];
sx q[0];
rz(-1.9504915) q[0];
sx q[0];
rz(-0.18549347) q[0];
rz(-0.56974757) q[1];
sx q[1];
rz(-2.2173917) q[1];
sx q[1];
rz(-1.1572908) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7737024) q[0];
sx q[0];
rz(-1.4285402) q[0];
sx q[0];
rz(0.41928798) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5026538) q[2];
sx q[2];
rz(-0.63387094) q[2];
sx q[2];
rz(-0.14574742) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.19543513) q[1];
sx q[1];
rz(-1.7845673) q[1];
sx q[1];
rz(2.5740088) q[1];
x q[2];
rz(1.1284105) q[3];
sx q[3];
rz(-0.81656045) q[3];
sx q[3];
rz(0.26536322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.40413228) q[2];
sx q[2];
rz(-2.9170333) q[2];
sx q[2];
rz(1.7411211) q[2];
rz(2.776966) q[3];
sx q[3];
rz(-0.76321634) q[3];
sx q[3];
rz(-2.3188685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6309361) q[0];
sx q[0];
rz(-0.83844227) q[0];
sx q[0];
rz(-3.1070218) q[0];
rz(-2.8055387) q[1];
sx q[1];
rz(-0.91468179) q[1];
sx q[1];
rz(-0.58759442) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1711271) q[0];
sx q[0];
rz(-1.1538527) q[0];
sx q[0];
rz(-1.72669) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73070902) q[2];
sx q[2];
rz(-1.8160422) q[2];
sx q[2];
rz(0.8697378) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.529189) q[1];
sx q[1];
rz(-1.8029787) q[1];
sx q[1];
rz(-3.0983689) q[1];
x q[2];
rz(0.67983277) q[3];
sx q[3];
rz(-1.9380016) q[3];
sx q[3];
rz(-0.96585551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8079638) q[2];
sx q[2];
rz(-0.89666349) q[2];
sx q[2];
rz(-2.4133852) q[2];
rz(0.41094574) q[3];
sx q[3];
rz(-0.22308895) q[3];
sx q[3];
rz(-1.2724916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67644607) q[0];
sx q[0];
rz(-1.5560173) q[0];
sx q[0];
rz(-1.7898855) q[0];
rz(-0.43988718) q[1];
sx q[1];
rz(-2.7688347) q[1];
sx q[1];
rz(-0.29040029) q[1];
rz(0.34392233) q[2];
sx q[2];
rz(-1.8544329) q[2];
sx q[2];
rz(0.0022619958) q[2];
rz(-2.7308913) q[3];
sx q[3];
rz(-3.0252785) q[3];
sx q[3];
rz(-2.5276285) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
