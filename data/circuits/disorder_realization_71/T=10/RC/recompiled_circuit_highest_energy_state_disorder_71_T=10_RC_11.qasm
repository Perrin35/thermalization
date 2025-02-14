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
rz(-0.072862536) q[0];
sx q[0];
rz(-2.1003567) q[0];
sx q[0];
rz(0.6676724) q[0];
rz(3.0296037) q[1];
sx q[1];
rz(-1.4854687) q[1];
sx q[1];
rz(2.8356584) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65016215) q[0];
sx q[0];
rz(-1.9823649) q[0];
sx q[0];
rz(1.6603907) q[0];
rz(-pi) q[1];
rz(2.5391766) q[2];
sx q[2];
rz(-2.4309046) q[2];
sx q[2];
rz(-1.7809694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.48253912) q[1];
sx q[1];
rz(-1.6371563) q[1];
sx q[1];
rz(2.769196) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.080360595) q[3];
sx q[3];
rz(-0.70801641) q[3];
sx q[3];
rz(0.77786672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.84451199) q[2];
sx q[2];
rz(-1.7252555) q[2];
sx q[2];
rz(0.97999209) q[2];
rz(-0.32198191) q[3];
sx q[3];
rz(-1.201509) q[3];
sx q[3];
rz(0.15596786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.62671536) q[0];
sx q[0];
rz(-1.6139655) q[0];
sx q[0];
rz(-2.2514586) q[0];
rz(-2.0251677) q[1];
sx q[1];
rz(-2.2869488) q[1];
sx q[1];
rz(1.1847624) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1370727) q[0];
sx q[0];
rz(-1.3810087) q[0];
sx q[0];
rz(-0.50317554) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1211195) q[2];
sx q[2];
rz(-1.2492078) q[2];
sx q[2];
rz(-2.2905397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55966942) q[1];
sx q[1];
rz(-1.6861177) q[1];
sx q[1];
rz(1.6827464) q[1];
x q[2];
rz(1.5542278) q[3];
sx q[3];
rz(-1.145716) q[3];
sx q[3];
rz(2.1147902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4778135) q[2];
sx q[2];
rz(-2.9728643) q[2];
sx q[2];
rz(-1.4230049) q[2];
rz(2.6723828) q[3];
sx q[3];
rz(-1.4956632) q[3];
sx q[3];
rz(0.34876987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94217268) q[0];
sx q[0];
rz(-1.9439789) q[0];
sx q[0];
rz(1.0823826) q[0];
rz(-2.9669145) q[1];
sx q[1];
rz(-1.382261) q[1];
sx q[1];
rz(-0.47294012) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2286005) q[0];
sx q[0];
rz(-0.45683041) q[0];
sx q[0];
rz(-1.940371) q[0];
rz(-1.9835112) q[2];
sx q[2];
rz(-0.7447401) q[2];
sx q[2];
rz(-1.2117529) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.14607695) q[1];
sx q[1];
rz(-0.56997609) q[1];
sx q[1];
rz(-1.6126627) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6430506) q[3];
sx q[3];
rz(-0.99851552) q[3];
sx q[3];
rz(-1.3839433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6579154) q[2];
sx q[2];
rz(-1.6562853) q[2];
sx q[2];
rz(1.8243054) q[2];
rz(-0.45051908) q[3];
sx q[3];
rz(-0.91166383) q[3];
sx q[3];
rz(-1.3594422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4860185) q[0];
sx q[0];
rz(-1.065932) q[0];
sx q[0];
rz(-2.4701212) q[0];
rz(-2.2010522) q[1];
sx q[1];
rz(-2.7214919) q[1];
sx q[1];
rz(0.95983517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3726138) q[0];
sx q[0];
rz(-0.44071925) q[0];
sx q[0];
rz(2.0924241) q[0];
rz(-pi) q[1];
rz(0.66633983) q[2];
sx q[2];
rz(-1.7494836) q[2];
sx q[2];
rz(1.7354244) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9588622) q[1];
sx q[1];
rz(-2.7736543) q[1];
sx q[1];
rz(-2.5065305) q[1];
x q[2];
rz(-2.6996451) q[3];
sx q[3];
rz(-1.4728896) q[3];
sx q[3];
rz(2.4349495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.81022108) q[2];
sx q[2];
rz(-1.081531) q[2];
sx q[2];
rz(0.01037154) q[2];
rz(-1.2822019) q[3];
sx q[3];
rz(-2.5344262) q[3];
sx q[3];
rz(-0.31266323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66987592) q[0];
sx q[0];
rz(-0.80729055) q[0];
sx q[0];
rz(-2.5655991) q[0];
rz(2.5895789) q[1];
sx q[1];
rz(-1.3060528) q[1];
sx q[1];
rz(0.920151) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3475889) q[0];
sx q[0];
rz(-0.92692843) q[0];
sx q[0];
rz(0.76066239) q[0];
rz(-pi) q[1];
rz(1.5414562) q[2];
sx q[2];
rz(-1.27904) q[2];
sx q[2];
rz(-2.2911366) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0219974) q[1];
sx q[1];
rz(-1.9174077) q[1];
sx q[1];
rz(0.24878169) q[1];
rz(0.12432762) q[3];
sx q[3];
rz(-1.6724068) q[3];
sx q[3];
rz(-1.6989482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0803904) q[2];
sx q[2];
rz(-1.7549425) q[2];
sx q[2];
rz(0.45026711) q[2];
rz(-0.14084147) q[3];
sx q[3];
rz(-1.2248657) q[3];
sx q[3];
rz(2.5950477) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8394354) q[0];
sx q[0];
rz(-1.7385229) q[0];
sx q[0];
rz(-2.9490525) q[0];
rz(2.4977066) q[1];
sx q[1];
rz(-1.5778678) q[1];
sx q[1];
rz(-0.083316915) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2177888) q[0];
sx q[0];
rz(-1.4155651) q[0];
sx q[0];
rz(-0.96331994) q[0];
rz(-pi) q[1];
rz(1.6896137) q[2];
sx q[2];
rz(-1.4897222) q[2];
sx q[2];
rz(0.18866779) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0472368) q[1];
sx q[1];
rz(-1.8810913) q[1];
sx q[1];
rz(0.17961802) q[1];
x q[2];
rz(-2.7658773) q[3];
sx q[3];
rz(-1.8109115) q[3];
sx q[3];
rz(-2.1841336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7328428) q[2];
sx q[2];
rz(-0.3877483) q[2];
sx q[2];
rz(-2.5852618) q[2];
rz(2.2522816) q[3];
sx q[3];
rz(-1.5143062) q[3];
sx q[3];
rz(0.48243943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0675875) q[0];
sx q[0];
rz(-2.6370625) q[0];
sx q[0];
rz(2.0489847) q[0];
rz(-0.69674528) q[1];
sx q[1];
rz(-0.69843355) q[1];
sx q[1];
rz(0.651407) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0677867) q[0];
sx q[0];
rz(-1.4069562) q[0];
sx q[0];
rz(-0.6730754) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90187855) q[2];
sx q[2];
rz(-1.9142391) q[2];
sx q[2];
rz(0.65367389) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9735865) q[1];
sx q[1];
rz(-2.1128707) q[1];
sx q[1];
rz(0.31967052) q[1];
x q[2];
rz(-1.8508554) q[3];
sx q[3];
rz(-0.055563888) q[3];
sx q[3];
rz(1.1058713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.50218987) q[2];
sx q[2];
rz(-2.4263589) q[2];
sx q[2];
rz(2.5779842) q[2];
rz(3.0230076) q[3];
sx q[3];
rz(-1.010681) q[3];
sx q[3];
rz(2.3732869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80657715) q[0];
sx q[0];
rz(-1.7105569) q[0];
sx q[0];
rz(-0.12120506) q[0];
rz(-2.99446) q[1];
sx q[1];
rz(-1.1558665) q[1];
sx q[1];
rz(2.3187231) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2911573) q[0];
sx q[0];
rz(-1.4802093) q[0];
sx q[0];
rz(-0.25048243) q[0];
rz(3.0953832) q[2];
sx q[2];
rz(-2.5760057) q[2];
sx q[2];
rz(-2.2639745) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40892021) q[1];
sx q[1];
rz(-0.80574811) q[1];
sx q[1];
rz(-0.21269704) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3579509) q[3];
sx q[3];
rz(-1.6506288) q[3];
sx q[3];
rz(-0.24536192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0269028) q[2];
sx q[2];
rz(-0.59640408) q[2];
sx q[2];
rz(0.58843311) q[2];
rz(0.19649291) q[3];
sx q[3];
rz(-1.4339707) q[3];
sx q[3];
rz(2.4832723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7115985) q[0];
sx q[0];
rz(-0.64084941) q[0];
sx q[0];
rz(0.95630056) q[0];
rz(-0.20005964) q[1];
sx q[1];
rz(-1.450489) q[1];
sx q[1];
rz(-0.53093451) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5945608) q[0];
sx q[0];
rz(-2.6480125) q[0];
sx q[0];
rz(1.7112687) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42163452) q[2];
sx q[2];
rz(-2.4608825) q[2];
sx q[2];
rz(1.9599078) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0772452) q[1];
sx q[1];
rz(-0.56659154) q[1];
sx q[1];
rz(-2.8342325) q[1];
rz(2.8538029) q[3];
sx q[3];
rz(-2.6734201) q[3];
sx q[3];
rz(2.5361907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7134573) q[2];
sx q[2];
rz(-2.0811446) q[2];
sx q[2];
rz(-1.1327845) q[2];
rz(-0.97635859) q[3];
sx q[3];
rz(-2.4553757) q[3];
sx q[3];
rz(1.2767701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3499632) q[0];
sx q[0];
rz(-0.11194734) q[0];
sx q[0];
rz(1.5524701) q[0];
rz(0.38996977) q[1];
sx q[1];
rz(-1.5682805) q[1];
sx q[1];
rz(-0.76036298) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0977265) q[0];
sx q[0];
rz(-2.0368493) q[0];
sx q[0];
rz(-1.4938484) q[0];
x q[1];
rz(-2.6796723) q[2];
sx q[2];
rz(-1.5963781) q[2];
sx q[2];
rz(-1.2939351) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2901569) q[1];
sx q[1];
rz(-2.7692215) q[1];
sx q[1];
rz(-1.464616) q[1];
rz(0.30837183) q[3];
sx q[3];
rz(-1.0532612) q[3];
sx q[3];
rz(2.8116159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4752263) q[2];
sx q[2];
rz(-0.97576371) q[2];
sx q[2];
rz(-2.547612) q[2];
rz(-0.62371671) q[3];
sx q[3];
rz(-1.2902322) q[3];
sx q[3];
rz(-2.311603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8754616) q[0];
sx q[0];
rz(-0.13127357) q[0];
sx q[0];
rz(1.776478) q[0];
rz(0.020137067) q[1];
sx q[1];
rz(-0.29300856) q[1];
sx q[1];
rz(-1.551052) q[1];
rz(1.9556373) q[2];
sx q[2];
rz(-1.5784752) q[2];
sx q[2];
rz(-1.2316647) q[2];
rz(-0.094811335) q[3];
sx q[3];
rz(-0.1764134) q[3];
sx q[3];
rz(-1.3865948) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
