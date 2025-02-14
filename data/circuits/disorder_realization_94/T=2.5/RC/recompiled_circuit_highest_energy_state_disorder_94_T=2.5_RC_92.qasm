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
rz(-1.0903519) q[0];
sx q[0];
rz(-0.89972377) q[0];
sx q[0];
rz(-1.1296912) q[0];
rz(0.55454412) q[1];
sx q[1];
rz(-0.50538844) q[1];
sx q[1];
rz(2.5880421) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81766291) q[0];
sx q[0];
rz(-2.637902) q[0];
sx q[0];
rz(-1.5970741) q[0];
rz(-0.0099382691) q[2];
sx q[2];
rz(-1.1608089) q[2];
sx q[2];
rz(0.25914295) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5419936) q[1];
sx q[1];
rz(-2.2570199) q[1];
sx q[1];
rz(2.0938796) q[1];
x q[2];
rz(-3.0985988) q[3];
sx q[3];
rz(-1.3519796) q[3];
sx q[3];
rz(0.18312632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0022137) q[2];
sx q[2];
rz(-1.6422108) q[2];
sx q[2];
rz(-1.1688983) q[2];
rz(-2.4161731) q[3];
sx q[3];
rz(-1.9099216) q[3];
sx q[3];
rz(1.5861082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91189522) q[0];
sx q[0];
rz(-0.40922368) q[0];
sx q[0];
rz(-1.0503861) q[0];
rz(-1.7627675) q[1];
sx q[1];
rz(-1.7345112) q[1];
sx q[1];
rz(-2.1398267) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79369583) q[0];
sx q[0];
rz(-0.2783723) q[0];
sx q[0];
rz(-1.9547561) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1304847) q[2];
sx q[2];
rz(-2.610865) q[2];
sx q[2];
rz(2.2953334) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6607501) q[1];
sx q[1];
rz(-0.4933115) q[1];
sx q[1];
rz(-0.50220614) q[1];
x q[2];
rz(-2.7185387) q[3];
sx q[3];
rz(-1.3025369) q[3];
sx q[3];
rz(0.34706193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.072814552) q[2];
sx q[2];
rz(-0.38663703) q[2];
sx q[2];
rz(-1.6437423) q[2];
rz(0.98613286) q[3];
sx q[3];
rz(-0.96753263) q[3];
sx q[3];
rz(-2.6223474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3471974) q[0];
sx q[0];
rz(-2.6674542) q[0];
sx q[0];
rz(0.073632181) q[0];
rz(-0.73078784) q[1];
sx q[1];
rz(-1.7981139) q[1];
sx q[1];
rz(-2.1582019) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2692719) q[0];
sx q[0];
rz(-1.2896363) q[0];
sx q[0];
rz(-2.9806656) q[0];
rz(-0.47731303) q[2];
sx q[2];
rz(-0.92318857) q[2];
sx q[2];
rz(-1.2106193) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3427998) q[1];
sx q[1];
rz(-1.4451318) q[1];
sx q[1];
rz(1.8428414) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0737209) q[3];
sx q[3];
rz(-1.5024324) q[3];
sx q[3];
rz(-1.9828486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.694515) q[2];
sx q[2];
rz(-2.647001) q[2];
sx q[2];
rz(0.05973235) q[2];
rz(0.51074243) q[3];
sx q[3];
rz(-0.35664883) q[3];
sx q[3];
rz(-1.7662778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8259976) q[0];
sx q[0];
rz(-0.70846486) q[0];
sx q[0];
rz(0.84061709) q[0];
rz(0.21545848) q[1];
sx q[1];
rz(-2.0550199) q[1];
sx q[1];
rz(0.1159018) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0784779) q[0];
sx q[0];
rz(-1.3097094) q[0];
sx q[0];
rz(-1.9223708) q[0];
rz(1.7809107) q[2];
sx q[2];
rz(-1.8571915) q[2];
sx q[2];
rz(2.1330733) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.85899587) q[1];
sx q[1];
rz(-2.1934185) q[1];
sx q[1];
rz(2.7864561) q[1];
x q[2];
rz(0.23502879) q[3];
sx q[3];
rz(-0.42065603) q[3];
sx q[3];
rz(3.0893081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0243715) q[2];
sx q[2];
rz(-1.8481959) q[2];
sx q[2];
rz(-2.2852066) q[2];
rz(-0.74445009) q[3];
sx q[3];
rz(-0.88862935) q[3];
sx q[3];
rz(-0.3652679) q[3];
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
rz(2.982548) q[0];
sx q[0];
rz(-2.5578975) q[0];
sx q[0];
rz(0.99345508) q[0];
rz(2.0221201) q[1];
sx q[1];
rz(-1.1886202) q[1];
sx q[1];
rz(3.1192034) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3684621) q[0];
sx q[0];
rz(-1.4903632) q[0];
sx q[0];
rz(-3.0162433) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3359433) q[2];
sx q[2];
rz(-0.73084007) q[2];
sx q[2];
rz(-1.8455055) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.34848475) q[1];
sx q[1];
rz(-0.99103084) q[1];
sx q[1];
rz(2.0391885) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0898354) q[3];
sx q[3];
rz(-0.93624712) q[3];
sx q[3];
rz(0.1635199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.28919724) q[2];
sx q[2];
rz(-2.9183233) q[2];
sx q[2];
rz(-1.7247464) q[2];
rz(-2.1220186) q[3];
sx q[3];
rz(-2.151078) q[3];
sx q[3];
rz(2.8250601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1932909) q[0];
sx q[0];
rz(-0.010007771) q[0];
sx q[0];
rz(-0.77578068) q[0];
rz(-1.4099482) q[1];
sx q[1];
rz(-1.8312788) q[1];
sx q[1];
rz(1.5595248) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821302) q[0];
sx q[0];
rz(-3.0728615) q[0];
sx q[0];
rz(2.0267846) q[0];
rz(-2.1072793) q[2];
sx q[2];
rz(-0.63648495) q[2];
sx q[2];
rz(-0.45490593) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.75943817) q[1];
sx q[1];
rz(-1.4389137) q[1];
sx q[1];
rz(-1.1005172) q[1];
rz(2.6234145) q[3];
sx q[3];
rz(-1.2478531) q[3];
sx q[3];
rz(-2.2707828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5887624) q[2];
sx q[2];
rz(-1.7071807) q[2];
sx q[2];
rz(0.32858783) q[2];
rz(-1.5594679) q[3];
sx q[3];
rz(-2.4170473) q[3];
sx q[3];
rz(2.4174378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62411031) q[0];
sx q[0];
rz(-2.3243853) q[0];
sx q[0];
rz(-2.4716614) q[0];
rz(-1.6239369) q[1];
sx q[1];
rz(-1.0847963) q[1];
sx q[1];
rz(-2.0673015) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23465189) q[0];
sx q[0];
rz(-2.2379506) q[0];
sx q[0];
rz(-1.2950598) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45410486) q[2];
sx q[2];
rz(-2.7412446) q[2];
sx q[2];
rz(0.12794447) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5902293) q[1];
sx q[1];
rz(-1.4783715) q[1];
sx q[1];
rz(1.2552099) q[1];
rz(-pi) q[2];
x q[2];
rz(1.303745) q[3];
sx q[3];
rz(-0.69483583) q[3];
sx q[3];
rz(-1.1649433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7488148) q[2];
sx q[2];
rz(-0.98069507) q[2];
sx q[2];
rz(-1.7559715) q[2];
rz(-2.4393926) q[3];
sx q[3];
rz(-0.32679138) q[3];
sx q[3];
rz(2.0328111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7053407) q[0];
sx q[0];
rz(-1.5441283) q[0];
sx q[0];
rz(-0.77823773) q[0];
rz(1.8477919) q[1];
sx q[1];
rz(-1.8840645) q[1];
sx q[1];
rz(0.055880849) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5326753) q[0];
sx q[0];
rz(-1.4977542) q[0];
sx q[0];
rz(-0.29204412) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7588024) q[2];
sx q[2];
rz(-1.3488608) q[2];
sx q[2];
rz(1.7772016) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1712492) q[1];
sx q[1];
rz(-1.9546485) q[1];
sx q[1];
rz(-1.7136736) q[1];
rz(-pi) q[2];
rz(-1.0033402) q[3];
sx q[3];
rz(-2.3262089) q[3];
sx q[3];
rz(1.585683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0539661) q[2];
sx q[2];
rz(-2.1199333) q[2];
sx q[2];
rz(0.6915687) q[2];
rz(0.21928731) q[3];
sx q[3];
rz(-1.6577474) q[3];
sx q[3];
rz(1.9482435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-1.2174862) q[0];
sx q[0];
rz(-0.036343887) q[0];
sx q[0];
rz(-2.0935667) q[0];
rz(2.3955087) q[1];
sx q[1];
rz(-1.289295) q[1];
sx q[1];
rz(0.41864905) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11514535) q[0];
sx q[0];
rz(-1.5799684) q[0];
sx q[0];
rz(-1.5546726) q[0];
rz(-pi) q[1];
rz(1.8739088) q[2];
sx q[2];
rz(-2.4753548) q[2];
sx q[2];
rz(-2.4002838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7884085) q[1];
sx q[1];
rz(-1.4338081) q[1];
sx q[1];
rz(-1.4333588) q[1];
rz(0.86781921) q[3];
sx q[3];
rz(-1.6074925) q[3];
sx q[3];
rz(-2.2040895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6333542) q[2];
sx q[2];
rz(-0.60428047) q[2];
sx q[2];
rz(-2.5227127) q[2];
rz(-0.45281705) q[3];
sx q[3];
rz(-1.489095) q[3];
sx q[3];
rz(2.648491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.598269) q[0];
sx q[0];
rz(-3.0040574) q[0];
sx q[0];
rz(1.2602873) q[0];
rz(0.29809412) q[1];
sx q[1];
rz(-0.9577848) q[1];
sx q[1];
rz(1.6424929) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74841046) q[0];
sx q[0];
rz(-1.5481189) q[0];
sx q[0];
rz(-2.3298954) q[0];
x q[1];
rz(-0.28147667) q[2];
sx q[2];
rz(-1.4325085) q[2];
sx q[2];
rz(-2.2654576) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3819653) q[1];
sx q[1];
rz(-1.5873199) q[1];
sx q[1];
rz(-1.3445059) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9177336) q[3];
sx q[3];
rz(-0.41439842) q[3];
sx q[3];
rz(0.5887659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5160211) q[2];
sx q[2];
rz(-1.8578153) q[2];
sx q[2];
rz(0.18216356) q[2];
rz(-2.9158909) q[3];
sx q[3];
rz(-1.00939) q[3];
sx q[3];
rz(1.3753447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7094649) q[0];
sx q[0];
rz(-1.3668677) q[0];
sx q[0];
rz(0.9912542) q[0];
rz(-1.2200914) q[1];
sx q[1];
rz(-0.5443926) q[1];
sx q[1];
rz(2.3701445) q[1];
rz(0.091339672) q[2];
sx q[2];
rz(-1.7484574) q[2];
sx q[2];
rz(-1.4469528) q[2];
rz(-2.3867802) q[3];
sx q[3];
rz(-2.1960498) q[3];
sx q[3];
rz(-1.4084569) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
