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
rz(-0.40773243) q[0];
sx q[0];
rz(4.3507504) q[0];
sx q[0];
rz(11.103295) q[0];
rz(-2.8830124) q[1];
sx q[1];
rz(-1.2376031) q[1];
sx q[1];
rz(0.13374506) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8834849) q[0];
sx q[0];
rz(-0.15722577) q[0];
sx q[0];
rz(1.6088465) q[0];
x q[1];
rz(2.1315246) q[2];
sx q[2];
rz(-2.0759931) q[2];
sx q[2];
rz(1.8672191) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7642149) q[1];
sx q[1];
rz(-1.0395607) q[1];
sx q[1];
rz(-1.9218535) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4179363) q[3];
sx q[3];
rz(-1.5709849) q[3];
sx q[3];
rz(2.3589695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3818843) q[2];
sx q[2];
rz(-0.81177652) q[2];
sx q[2];
rz(-1.3105357) q[2];
rz(1.7740907) q[3];
sx q[3];
rz(-1.1279227) q[3];
sx q[3];
rz(0.016157063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1644208) q[0];
sx q[0];
rz(-1.9855969) q[0];
sx q[0];
rz(1.0154065) q[0];
rz(2.7482391) q[1];
sx q[1];
rz(-1.669408) q[1];
sx q[1];
rz(0.70158395) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9727952) q[0];
sx q[0];
rz(-1.4864362) q[0];
sx q[0];
rz(-2.286351) q[0];
rz(-pi) q[1];
x q[1];
rz(0.010741269) q[2];
sx q[2];
rz(-1.2678384) q[2];
sx q[2];
rz(0.078127351) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.9522515) q[1];
sx q[1];
rz(-0.32508141) q[1];
sx q[1];
rz(1.8081244) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75922482) q[3];
sx q[3];
rz(-2.1442778) q[3];
sx q[3];
rz(-0.92347565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.027448805) q[2];
sx q[2];
rz(-2.1178718) q[2];
sx q[2];
rz(-1.1055498) q[2];
rz(-0.28087273) q[3];
sx q[3];
rz(-2.5307145) q[3];
sx q[3];
rz(0.15225473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0356782) q[0];
sx q[0];
rz(-2.4690101) q[0];
sx q[0];
rz(-2.0689082) q[0];
rz(0.013280344) q[1];
sx q[1];
rz(-1.5417136) q[1];
sx q[1];
rz(0.57659155) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.428513) q[0];
sx q[0];
rz(-1.2831076) q[0];
sx q[0];
rz(1.6384533) q[0];
rz(-3.0579975) q[2];
sx q[2];
rz(-1.7033249) q[2];
sx q[2];
rz(2.8797163) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.84725897) q[1];
sx q[1];
rz(-1.0953961) q[1];
sx q[1];
rz(-0.21902276) q[1];
rz(1.1880831) q[3];
sx q[3];
rz(-0.59186223) q[3];
sx q[3];
rz(0.6955717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2836634) q[2];
sx q[2];
rz(-2.2495705) q[2];
sx q[2];
rz(2.499495) q[2];
rz(-0.56504956) q[3];
sx q[3];
rz(-0.27500209) q[3];
sx q[3];
rz(-0.16998418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86244407) q[0];
sx q[0];
rz(-1.886241) q[0];
sx q[0];
rz(-2.8992262) q[0];
rz(0.53388059) q[1];
sx q[1];
rz(-0.68478525) q[1];
sx q[1];
rz(1.6530316) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72172644) q[0];
sx q[0];
rz(-0.97630608) q[0];
sx q[0];
rz(1.4851024) q[0];
x q[1];
rz(-0.42885355) q[2];
sx q[2];
rz(-1.2786713) q[2];
sx q[2];
rz(-3.0512864) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3900534) q[1];
sx q[1];
rz(-1.8597457) q[1];
sx q[1];
rz(-2.3759936) q[1];
rz(1.6107992) q[3];
sx q[3];
rz(-1.2712595) q[3];
sx q[3];
rz(-1.9301723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.36294595) q[2];
sx q[2];
rz(-1.8774418) q[2];
sx q[2];
rz(1.6710336) q[2];
rz(3.0342024) q[3];
sx q[3];
rz(-1.8348179) q[3];
sx q[3];
rz(-1.2764527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11112467) q[0];
sx q[0];
rz(-0.46115369) q[0];
sx q[0];
rz(-1.2029458) q[0];
rz(-0.89490926) q[1];
sx q[1];
rz(-1.1824965) q[1];
sx q[1];
rz(2.9295909) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.014678) q[0];
sx q[0];
rz(-2.0058332) q[0];
sx q[0];
rz(-0.90330173) q[0];
x q[1];
rz(-1.8708399) q[2];
sx q[2];
rz(-1.0677538) q[2];
sx q[2];
rz(1.6479657) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3983237) q[1];
sx q[1];
rz(-1.124427) q[1];
sx q[1];
rz(-0.94029398) q[1];
rz(-pi) q[2];
rz(-3.131989) q[3];
sx q[3];
rz(-2.1240747) q[3];
sx q[3];
rz(-1.5095078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8733069) q[2];
sx q[2];
rz(-0.77750677) q[2];
sx q[2];
rz(3.0262465) q[2];
rz(-2.1558732) q[3];
sx q[3];
rz(-1.4109572) q[3];
sx q[3];
rz(2.7058097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2301521) q[0];
sx q[0];
rz(-1.4608915) q[0];
sx q[0];
rz(-2.4440515) q[0];
rz(0.41360924) q[1];
sx q[1];
rz(-1.6802843) q[1];
sx q[1];
rz(-0.70346171) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1915775) q[0];
sx q[0];
rz(-0.30739014) q[0];
sx q[0];
rz(2.6527575) q[0];
rz(-pi) q[1];
rz(-2.7532534) q[2];
sx q[2];
rz(-2.1415798) q[2];
sx q[2];
rz(2.6663189) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5589964) q[1];
sx q[1];
rz(-1.1038053) q[1];
sx q[1];
rz(-1.9860616) q[1];
x q[2];
rz(2.260842) q[3];
sx q[3];
rz(-2.7610169) q[3];
sx q[3];
rz(2.5799023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.10282639) q[2];
sx q[2];
rz(-2.8732754) q[2];
sx q[2];
rz(-1.9290257) q[2];
rz(2.8596527) q[3];
sx q[3];
rz(-1.599267) q[3];
sx q[3];
rz(-2.6350002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7956227) q[0];
sx q[0];
rz(-0.75241929) q[0];
sx q[0];
rz(2.2947626) q[0];
rz(2.1961191) q[1];
sx q[1];
rz(-1.5739601) q[1];
sx q[1];
rz(0.6792773) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69791171) q[0];
sx q[0];
rz(-0.4831995) q[0];
sx q[0];
rz(-0.16541055) q[0];
rz(-2.7921259) q[2];
sx q[2];
rz(-2.4878056) q[2];
sx q[2];
rz(0.9330627) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8161387) q[1];
sx q[1];
rz(-1.3846163) q[1];
sx q[1];
rz(-1.0503616) q[1];
rz(-pi) q[2];
rz(1.6369989) q[3];
sx q[3];
rz(-2.8567064) q[3];
sx q[3];
rz(-1.6411622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.88963738) q[2];
sx q[2];
rz(-0.95459443) q[2];
sx q[2];
rz(-1.4944705) q[2];
rz(-2.5877118) q[3];
sx q[3];
rz(-1.5364105) q[3];
sx q[3];
rz(1.4166098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3225591) q[0];
sx q[0];
rz(-0.77545866) q[0];
sx q[0];
rz(1.092528) q[0];
rz(-2.20772) q[1];
sx q[1];
rz(-1.7364419) q[1];
sx q[1];
rz(-2.3563103) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3438246) q[0];
sx q[0];
rz(-1.1189119) q[0];
sx q[0];
rz(-2.0411819) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6996428) q[2];
sx q[2];
rz(-2.0734678) q[2];
sx q[2];
rz(1.3208226) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3242625) q[1];
sx q[1];
rz(-1.3249389) q[1];
sx q[1];
rz(-2.861633) q[1];
rz(-pi) q[2];
rz(-2.5232072) q[3];
sx q[3];
rz(-1.9342285) q[3];
sx q[3];
rz(-1.6101162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.87172047) q[2];
sx q[2];
rz(-1.4464658) q[2];
sx q[2];
rz(-2.5093057) q[2];
rz(-0.91442433) q[3];
sx q[3];
rz(-1.881003) q[3];
sx q[3];
rz(-2.9818592) q[3];
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
rz(0.81894994) q[0];
sx q[0];
rz(-0.85700789) q[0];
sx q[0];
rz(0.44347611) q[0];
rz(0.30287287) q[1];
sx q[1];
rz(-2.5587406) q[1];
sx q[1];
rz(1.3076967) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0904059) q[0];
sx q[0];
rz(-1.4009579) q[0];
sx q[0];
rz(-0.30533653) q[0];
x q[1];
rz(2.6340387) q[2];
sx q[2];
rz(-1.8448463) q[2];
sx q[2];
rz(-0.7407032) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84660599) q[1];
sx q[1];
rz(-1.8364803) q[1];
sx q[1];
rz(1.0636661) q[1];
x q[2];
rz(0.29662432) q[3];
sx q[3];
rz(-1.1184177) q[3];
sx q[3];
rz(1.5993702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0206454) q[2];
sx q[2];
rz(-2.1006613) q[2];
sx q[2];
rz(-0.61332235) q[2];
rz(-0.81765085) q[3];
sx q[3];
rz(-0.48912564) q[3];
sx q[3];
rz(2.8221655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.9870975) q[0];
sx q[0];
rz(-1.7648062) q[0];
sx q[0];
rz(-0.35926551) q[0];
rz(-2.477395) q[1];
sx q[1];
rz(-1.3281053) q[1];
sx q[1];
rz(0.42116234) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065124113) q[0];
sx q[0];
rz(-1.9389922) q[0];
sx q[0];
rz(2.9028331) q[0];
rz(-pi) q[1];
rz(1.498327) q[2];
sx q[2];
rz(-1.5761099) q[2];
sx q[2];
rz(-1.1491691) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5790184) q[1];
sx q[1];
rz(-2.0236602) q[1];
sx q[1];
rz(-2.8500506) q[1];
rz(-pi) q[2];
rz(-2.3844658) q[3];
sx q[3];
rz(-0.76092488) q[3];
sx q[3];
rz(-1.2989375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8757214) q[2];
sx q[2];
rz(-0.096780626) q[2];
sx q[2];
rz(1.9270012) q[2];
rz(-0.38861361) q[3];
sx q[3];
rz(-0.92258421) q[3];
sx q[3];
rz(2.0428366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.494396) q[0];
sx q[0];
rz(-1.5460486) q[0];
sx q[0];
rz(-1.8478951) q[0];
rz(-0.21678674) q[1];
sx q[1];
rz(-0.68991359) q[1];
sx q[1];
rz(0.51295113) q[1];
rz(0.51914712) q[2];
sx q[2];
rz(-0.74719528) q[2];
sx q[2];
rz(-2.0120646) q[2];
rz(-1.3043591) q[3];
sx q[3];
rz(-1.1188917) q[3];
sx q[3];
rz(1.6592142) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
