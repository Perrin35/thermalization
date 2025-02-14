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
rz(-3.1138104) q[0];
sx q[0];
rz(-1.8960928) q[0];
sx q[0];
rz(-1.0863289) q[0];
rz(2.0431986) q[1];
sx q[1];
rz(1.2135222) q[1];
sx q[1];
rz(10.57099) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1725464) q[0];
sx q[0];
rz(-1.8964701) q[0];
sx q[0];
rz(0.17071755) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.97922) q[2];
sx q[2];
rz(-2.1222489) q[2];
sx q[2];
rz(-1.1609032) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0806345) q[1];
sx q[1];
rz(-0.74392156) q[1];
sx q[1];
rz(0.77772909) q[1];
x q[2];
rz(-0.2262874) q[3];
sx q[3];
rz(-2.0620966) q[3];
sx q[3];
rz(-2.0313583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9132797) q[2];
sx q[2];
rz(-1.5756807) q[2];
sx q[2];
rz(-0.3621462) q[2];
rz(-0.95237887) q[3];
sx q[3];
rz(-0.55914545) q[3];
sx q[3];
rz(2.4366116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1700965) q[0];
sx q[0];
rz(-0.2146475) q[0];
sx q[0];
rz(-1.5317408) q[0];
rz(1.8929298) q[1];
sx q[1];
rz(-1.6207393) q[1];
sx q[1];
rz(2.2889287) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.296761) q[0];
sx q[0];
rz(-2.0728025) q[0];
sx q[0];
rz(2.7472578) q[0];
rz(-pi) q[1];
rz(1.1617848) q[2];
sx q[2];
rz(-0.30050052) q[2];
sx q[2];
rz(2.7331293) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1601847) q[1];
sx q[1];
rz(-1.369919) q[1];
sx q[1];
rz(-2.1698879) q[1];
rz(-pi) q[2];
rz(-2.1308822) q[3];
sx q[3];
rz(-1.9017856) q[3];
sx q[3];
rz(-0.50332848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2109005) q[2];
sx q[2];
rz(-0.96052581) q[2];
sx q[2];
rz(1.9453913) q[2];
rz(-3.1148552) q[3];
sx q[3];
rz(-0.49632159) q[3];
sx q[3];
rz(-1.4727717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4836327) q[0];
sx q[0];
rz(-2.8909029) q[0];
sx q[0];
rz(0.87580097) q[0];
rz(-0.35975131) q[1];
sx q[1];
rz(-1.3641554) q[1];
sx q[1];
rz(-1.0035427) q[1];
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
rz(-0.92873145) q[2];
sx q[2];
rz(-1.2479397) q[2];
sx q[2];
rz(2.0546469) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0599614) q[1];
sx q[1];
rz(-1.7137032) q[1];
sx q[1];
rz(0.68024723) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67387132) q[3];
sx q[3];
rz(-2.0099075) q[3];
sx q[3];
rz(1.278745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.935219) q[2];
sx q[2];
rz(-2.424365) q[2];
sx q[2];
rz(-0.75500542) q[2];
rz(-0.098091789) q[3];
sx q[3];
rz(-0.57049975) q[3];
sx q[3];
rz(-2.7437955) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0841242) q[0];
sx q[0];
rz(-1.6300772) q[0];
sx q[0];
rz(2.496642) q[0];
rz(1.0954674) q[1];
sx q[1];
rz(-2.7332833) q[1];
sx q[1];
rz(-1.1114978) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1350461) q[0];
sx q[0];
rz(-1.6301951) q[0];
sx q[0];
rz(3.0778372) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45504163) q[2];
sx q[2];
rz(-0.18559449) q[2];
sx q[2];
rz(-1.1374823) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7202338) q[1];
sx q[1];
rz(-1.780527) q[1];
sx q[1];
rz(-1.420182) q[1];
rz(-pi) q[2];
rz(-1.7913412) q[3];
sx q[3];
rz(-0.84201563) q[3];
sx q[3];
rz(0.33607863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.633454) q[2];
sx q[2];
rz(-0.78712574) q[2];
sx q[2];
rz(0.48466361) q[2];
rz(2.700452) q[3];
sx q[3];
rz(-1.7287247) q[3];
sx q[3];
rz(-1.887656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6345374) q[0];
sx q[0];
rz(-0.84369722) q[0];
sx q[0];
rz(1.7150568) q[0];
rz(-0.85975319) q[1];
sx q[1];
rz(-0.25329241) q[1];
sx q[1];
rz(-0.8304798) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7380421) q[0];
sx q[0];
rz(-1.0071686) q[0];
sx q[0];
rz(0.72569816) q[0];
rz(-pi) q[1];
rz(1.4926339) q[2];
sx q[2];
rz(-0.63428464) q[2];
sx q[2];
rz(1.8820349) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4791058) q[1];
sx q[1];
rz(-1.431152) q[1];
sx q[1];
rz(-0.73151155) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4111341) q[3];
sx q[3];
rz(-1.9267907) q[3];
sx q[3];
rz(-2.1170664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6972203) q[2];
sx q[2];
rz(-2.1264919) q[2];
sx q[2];
rz(-2.3946136) q[2];
rz(-2.8367481) q[3];
sx q[3];
rz(-1.3957142) q[3];
sx q[3];
rz(-1.2386809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1759258) q[0];
sx q[0];
rz(-1.8400064) q[0];
sx q[0];
rz(-2.9081705) q[0];
rz(-2.6790791) q[1];
sx q[1];
rz(-1.9513444) q[1];
sx q[1];
rz(-0.40406427) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8916805) q[0];
sx q[0];
rz(-0.95522987) q[0];
sx q[0];
rz(-2.8717442) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5954757) q[2];
sx q[2];
rz(-0.48729333) q[2];
sx q[2];
rz(-0.55216002) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3531472) q[1];
sx q[1];
rz(-2.569558) q[1];
sx q[1];
rz(-0.20010179) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5109291) q[3];
sx q[3];
rz(-1.3913097) q[3];
sx q[3];
rz(-1.6398182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.81249753) q[2];
sx q[2];
rz(-1.3500682) q[2];
sx q[2];
rz(-0.65328944) q[2];
rz(1.6434044) q[3];
sx q[3];
rz(-0.61763063) q[3];
sx q[3];
rz(0.60630715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6752328) q[0];
sx q[0];
rz(-1.2132069) q[0];
sx q[0];
rz(2.8208222) q[0];
rz(-2.4877211) q[1];
sx q[1];
rz(-2.8120698) q[1];
sx q[1];
rz(-3.0810862) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.084235926) q[0];
sx q[0];
rz(-0.31937803) q[0];
sx q[0];
rz(-1.0127629) q[0];
rz(-2.195792) q[2];
sx q[2];
rz(-0.93987465) q[2];
sx q[2];
rz(0.45880246) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.49566073) q[1];
sx q[1];
rz(-2.0510608) q[1];
sx q[1];
rz(0.86841327) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7869446) q[3];
sx q[3];
rz(-1.186968) q[3];
sx q[3];
rz(0.5634748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3939646) q[2];
sx q[2];
rz(-1.4614033) q[2];
sx q[2];
rz(-2.4192269) q[2];
rz(-2.8584621) q[3];
sx q[3];
rz(-2.2524998) q[3];
sx q[3];
rz(2.6319671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79001456) q[0];
sx q[0];
rz(-0.55792648) q[0];
sx q[0];
rz(-0.18490069) q[0];
rz(0.50390538) q[1];
sx q[1];
rz(-1.7920707) q[1];
sx q[1];
rz(0.46461836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9757742) q[0];
sx q[0];
rz(-1.4041472) q[0];
sx q[0];
rz(-1.7383132) q[0];
rz(-2.9275168) q[2];
sx q[2];
rz(-0.42595562) q[2];
sx q[2];
rz(1.4115806) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7261207) q[1];
sx q[1];
rz(-1.1562043) q[1];
sx q[1];
rz(1.1555671) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7691602) q[3];
sx q[3];
rz(-2.2806205) q[3];
sx q[3];
rz(-1.7539303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.23286143) q[2];
sx q[2];
rz(-1.5166914) q[2];
sx q[2];
rz(1.0535343) q[2];
rz(2.4145224) q[3];
sx q[3];
rz(-1.2289685) q[3];
sx q[3];
rz(0.54522902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.9497249) q[0];
sx q[0];
rz(-1.1911012) q[0];
sx q[0];
rz(2.9560992) q[0];
rz(-0.56974757) q[1];
sx q[1];
rz(-2.2173917) q[1];
sx q[1];
rz(-1.1572908) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51087166) q[0];
sx q[0];
rz(-2.7001885) q[0];
sx q[0];
rz(2.8033103) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5026538) q[2];
sx q[2];
rz(-0.63387094) q[2];
sx q[2];
rz(-2.9958452) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.696582) q[1];
sx q[1];
rz(-0.60234501) q[1];
sx q[1];
rz(-0.38378832) q[1];
rz(-pi) q[2];
rz(0.80482153) q[3];
sx q[3];
rz(-1.8880883) q[3];
sx q[3];
rz(-2.1498093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.40413228) q[2];
sx q[2];
rz(-2.9170333) q[2];
sx q[2];
rz(-1.7411211) q[2];
rz(0.36462668) q[3];
sx q[3];
rz(-2.3783763) q[3];
sx q[3];
rz(0.82272416) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6309361) q[0];
sx q[0];
rz(-2.3031504) q[0];
sx q[0];
rz(-3.1070218) q[0];
rz(-2.8055387) q[1];
sx q[1];
rz(-0.91468179) q[1];
sx q[1];
rz(2.5539982) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97046554) q[0];
sx q[0];
rz(-1.1538527) q[0];
sx q[0];
rz(1.4149026) q[0];
rz(2.7828092) q[2];
sx q[2];
rz(-0.76350313) q[2];
sx q[2];
rz(0.43649188) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.9683444) q[1];
sx q[1];
rz(-1.6128596) q[1];
sx q[1];
rz(-1.8031881) q[1];
rz(-1.1114656) q[3];
sx q[3];
rz(-0.94375347) q[3];
sx q[3];
rz(-2.2541943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8079638) q[2];
sx q[2];
rz(-0.89666349) q[2];
sx q[2];
rz(-2.4133852) q[2];
rz(-2.7306469) q[3];
sx q[3];
rz(-2.9185037) q[3];
sx q[3];
rz(-1.8691011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4651466) q[0];
sx q[0];
rz(-1.5855753) q[0];
sx q[0];
rz(1.3517071) q[0];
rz(0.43988718) q[1];
sx q[1];
rz(-0.37275795) q[1];
sx q[1];
rz(2.8511924) q[1];
rz(2.4287379) q[2];
sx q[2];
rz(-0.44217449) q[2];
sx q[2];
rz(0.90978734) q[2];
rz(-1.5241809) q[3];
sx q[3];
rz(-1.4641932) q[3];
sx q[3];
rz(-2.1144397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
