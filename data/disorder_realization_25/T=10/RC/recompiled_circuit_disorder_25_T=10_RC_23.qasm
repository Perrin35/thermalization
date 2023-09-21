OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(-2.893653) q[0];
sx q[0];
rz(0.82013446) q[0];
rz(-3.1007383) q[1];
sx q[1];
rz(-0.78894579) q[1];
sx q[1];
rz(-0.087892428) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1782921) q[0];
sx q[0];
rz(-1.1305729) q[0];
sx q[0];
rz(0.51142366) q[0];
rz(-pi) q[1];
rz(0.90859969) q[2];
sx q[2];
rz(-2.5052921) q[2];
sx q[2];
rz(1.5696021) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.81401285) q[1];
sx q[1];
rz(-2.6373133) q[1];
sx q[1];
rz(0.0045279702) q[1];
rz(0.0058562972) q[3];
sx q[3];
rz(-1.3121038) q[3];
sx q[3];
rz(-0.70219016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1408046) q[2];
sx q[2];
rz(-1.6892097) q[2];
sx q[2];
rz(-0.56935707) q[2];
rz(-1.6128444) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(1.3585842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5434718) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(2.5879481) q[0];
rz(1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(1.233261) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20576661) q[0];
sx q[0];
rz(-2.2201949) q[0];
sx q[0];
rz(-2.4539024) q[0];
rz(2.4670062) q[2];
sx q[2];
rz(-1.470675) q[2];
sx q[2];
rz(-1.1652201) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.43014363) q[1];
sx q[1];
rz(-1.5866382) q[1];
sx q[1];
rz(1.364691) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29970701) q[3];
sx q[3];
rz(-0.68972677) q[3];
sx q[3];
rz(-2.0852065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1374986) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(-3.1393576) q[2];
rz(-2.3114752) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(-1.8921651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24901351) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(0.85154831) q[0];
rz(0.7710723) q[1];
sx q[1];
rz(-1.4665736) q[1];
sx q[1];
rz(0.071203701) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035633798) q[0];
sx q[0];
rz(-1.7674315) q[0];
sx q[0];
rz(-1.3282302) q[0];
rz(-pi) q[1];
rz(-1.115828) q[2];
sx q[2];
rz(-1.8075382) q[2];
sx q[2];
rz(-0.67982212) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5350266) q[1];
sx q[1];
rz(-2.2177794) q[1];
sx q[1];
rz(1.7391298) q[1];
rz(-pi) q[2];
rz(2.2749388) q[3];
sx q[3];
rz(-2.1979077) q[3];
sx q[3];
rz(1.31124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3123793) q[2];
sx q[2];
rz(-2.2177314) q[2];
sx q[2];
rz(-2.6181347) q[2];
rz(-2.8097025) q[3];
sx q[3];
rz(-3.0811946) q[3];
sx q[3];
rz(-2.6383242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(1.7317384) q[0];
sx q[0];
rz(-2.8926352) q[0];
sx q[0];
rz(0.82558924) q[0];
rz(-1.9748953) q[1];
sx q[1];
rz(-0.86667934) q[1];
sx q[1];
rz(1.4473787) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3576413) q[0];
sx q[0];
rz(-2.8581627) q[0];
sx q[0];
rz(2.1901603) q[0];
rz(1.3924613) q[2];
sx q[2];
rz(-1.5797838) q[2];
sx q[2];
rz(2.044878) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4059658) q[1];
sx q[1];
rz(-1.4117129) q[1];
sx q[1];
rz(-0.87079485) q[1];
rz(1.3284239) q[3];
sx q[3];
rz(-1.216785) q[3];
sx q[3];
rz(0.44204516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6772785) q[2];
sx q[2];
rz(-1.7769287) q[2];
sx q[2];
rz(2.9966667) q[2];
rz(-2.1302917) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(-0.01468006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99575627) q[0];
sx q[0];
rz(-3.0881112) q[0];
sx q[0];
rz(-0.68403912) q[0];
rz(1.1902635) q[1];
sx q[1];
rz(-2.4763156) q[1];
sx q[1];
rz(-1.2129983) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4764403) q[0];
sx q[0];
rz(-1.1993009) q[0];
sx q[0];
rz(-1.8508338) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91765399) q[2];
sx q[2];
rz(-1.0461763) q[2];
sx q[2];
rz(1.4337002) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7291521) q[1];
sx q[1];
rz(-1.3153331) q[1];
sx q[1];
rz(1.4490119) q[1];
x q[2];
rz(2.1608859) q[3];
sx q[3];
rz(-2.8299035) q[3];
sx q[3];
rz(-0.89524549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6325536) q[2];
sx q[2];
rz(-1.2712487) q[2];
sx q[2];
rz(2.5115013) q[2];
rz(2.5693494) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(-2.2563289) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0260139) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(0.28636006) q[0];
rz(0.59965602) q[1];
sx q[1];
rz(-2.0136132) q[1];
sx q[1];
rz(-2.5568331) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6385348) q[0];
sx q[0];
rz(-0.98031822) q[0];
sx q[0];
rz(-1.1919341) q[0];
rz(-0.037126304) q[2];
sx q[2];
rz(-2.526844) q[2];
sx q[2];
rz(-2.6908713) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4510182) q[1];
sx q[1];
rz(-1.8585397) q[1];
sx q[1];
rz(0.87177999) q[1];
x q[2];
rz(1.5201969) q[3];
sx q[3];
rz(-1.4046602) q[3];
sx q[3];
rz(-1.43464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.16453234) q[2];
sx q[2];
rz(-2.2571199) q[2];
sx q[2];
rz(1.5131081) q[2];
rz(2.1785054) q[3];
sx q[3];
rz(-1.77308) q[3];
sx q[3];
rz(-0.62455463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0528089) q[0];
sx q[0];
rz(-1.4316906) q[0];
sx q[0];
rz(-2.5928296) q[0];
rz(0.051963003) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(2.5351977) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0604027) q[0];
sx q[0];
rz(-2.5572544) q[0];
sx q[0];
rz(1.1330182) q[0];
x q[1];
rz(-2.7658471) q[2];
sx q[2];
rz(-0.94259113) q[2];
sx q[2];
rz(1.9919765) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4013306) q[1];
sx q[1];
rz(-1.6832502) q[1];
sx q[1];
rz(2.3349891) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2783979) q[3];
sx q[3];
rz(-2.6226225) q[3];
sx q[3];
rz(2.2591345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2018044) q[2];
sx q[2];
rz(-0.99210343) q[2];
sx q[2];
rz(-0.50841224) q[2];
rz(-1.1172179) q[3];
sx q[3];
rz(-2.9682187) q[3];
sx q[3];
rz(-1.4846876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(2.8023119) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(-0.75396496) q[0];
rz(2.1915961) q[1];
sx q[1];
rz(-1.4818622) q[1];
sx q[1];
rz(2.9139013) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84275093) q[0];
sx q[0];
rz(-1.1584873) q[0];
sx q[0];
rz(-0.077667872) q[0];
rz(-pi) q[1];
rz(-0.67543244) q[2];
sx q[2];
rz(-2.5917412) q[2];
sx q[2];
rz(0.39381248) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88340064) q[1];
sx q[1];
rz(-1.0419783) q[1];
sx q[1];
rz(-1.8316815) q[1];
x q[2];
rz(-0.55898198) q[3];
sx q[3];
rz(-1.535927) q[3];
sx q[3];
rz(-2.1502286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.078538744) q[2];
sx q[2];
rz(-0.20919122) q[2];
sx q[2];
rz(-2.3030346) q[2];
rz(2.7770384) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(-0.84028876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7446328) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(-2.3378085) q[0];
rz(-1.0546168) q[1];
sx q[1];
rz(-2.5024253) q[1];
sx q[1];
rz(1.1484336) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.046612) q[0];
sx q[0];
rz(-1.5946832) q[0];
sx q[0];
rz(1.4987962) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36275136) q[2];
sx q[2];
rz(-1.0748378) q[2];
sx q[2];
rz(2.9438058) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2003186) q[1];
sx q[1];
rz(-0.66220821) q[1];
sx q[1];
rz(0.83076417) q[1];
x q[2];
rz(-1.6059444) q[3];
sx q[3];
rz(-1.2190378) q[3];
sx q[3];
rz(-0.37946057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.049872963) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(2.3633374) q[2];
rz(-2.9459279) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(1.0218609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8675999) q[0];
sx q[0];
rz(-0.99853981) q[0];
sx q[0];
rz(2.8883949) q[0];
rz(-0.7397488) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(-2.443312) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65933278) q[0];
sx q[0];
rz(-1.2685304) q[0];
sx q[0];
rz(-1.6460101) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6449498) q[2];
sx q[2];
rz(-1.5955131) q[2];
sx q[2];
rz(2.2029684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.89430289) q[1];
sx q[1];
rz(-1.786788) q[1];
sx q[1];
rz(-1.1029907) q[1];
x q[2];
rz(2.9000557) q[3];
sx q[3];
rz(-1.0165443) q[3];
sx q[3];
rz(-2.9709771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.11761052) q[2];
sx q[2];
rz(-2.2587946) q[2];
sx q[2];
rz(0.84021604) q[2];
rz(2.501287) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(-1.9742112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-1.8828076) q[0];
sx q[0];
rz(-0.85492815) q[0];
sx q[0];
rz(1.7229765) q[0];
rz(-3.1124658) q[1];
sx q[1];
rz(-0.11321414) q[1];
sx q[1];
rz(-1.3197457) q[1];
rz(-2.6559033) q[2];
sx q[2];
rz(-0.66076856) q[2];
sx q[2];
rz(3.0271157) q[2];
rz(-1.7520262) q[3];
sx q[3];
rz(-2.443234) q[3];
sx q[3];
rz(0.20886226) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
