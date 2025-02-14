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
rz(2.3866374) q[0];
sx q[0];
rz(-0.75140262) q[0];
sx q[0];
rz(-1.0487392) q[0];
rz(-2.7222848) q[1];
sx q[1];
rz(4.5276548) q[1];
sx q[1];
rz(6.164896) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2631419) q[0];
sx q[0];
rz(-2.1241444) q[0];
sx q[0];
rz(0.31991495) q[0];
rz(-0.30960942) q[2];
sx q[2];
rz(-0.053507858) q[2];
sx q[2];
rz(0.87602304) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.58913) q[1];
sx q[1];
rz(-1.5877866) q[1];
sx q[1];
rz(0.029689992) q[1];
rz(-pi) q[2];
rz(-2.8382073) q[3];
sx q[3];
rz(-1.4974563) q[3];
sx q[3];
rz(2.908542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2179541) q[2];
sx q[2];
rz(-3.1380234) q[2];
sx q[2];
rz(2.9812319) q[2];
rz(1.977836) q[3];
sx q[3];
rz(-2.0573503) q[3];
sx q[3];
rz(-2.3273996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9775951) q[0];
sx q[0];
rz(-1.6544592) q[0];
sx q[0];
rz(2.1556222) q[0];
rz(-1.5522955) q[1];
sx q[1];
rz(-0.28737107) q[1];
sx q[1];
rz(-1.5812965) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4673104) q[0];
sx q[0];
rz(-0.63100609) q[0];
sx q[0];
rz(2.254807) q[0];
x q[1];
rz(-2.5471117) q[2];
sx q[2];
rz(-1.5882601) q[2];
sx q[2];
rz(-1.524615) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1007019) q[1];
sx q[1];
rz(-1.9305848) q[1];
sx q[1];
rz(-1.3967819) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12388568) q[3];
sx q[3];
rz(-2.200211) q[3];
sx q[3];
rz(0.38485011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.61957773) q[2];
sx q[2];
rz(-1.3238944) q[2];
sx q[2];
rz(1.0349234) q[2];
rz(-0.50152913) q[3];
sx q[3];
rz(-3.0540255) q[3];
sx q[3];
rz(2.1653304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4185249) q[0];
sx q[0];
rz(-1.9880966) q[0];
sx q[0];
rz(1.9368517) q[0];
rz(-1.0459666) q[1];
sx q[1];
rz(-3.0515262) q[1];
sx q[1];
rz(0.14030309) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36049309) q[0];
sx q[0];
rz(-1.3083959) q[0];
sx q[0];
rz(-2.187832) q[0];
rz(-pi) q[1];
rz(0.84425462) q[2];
sx q[2];
rz(-1.8945969) q[2];
sx q[2];
rz(-2.5716788) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.54934422) q[1];
sx q[1];
rz(-0.98824938) q[1];
sx q[1];
rz(-1.2986533) q[1];
rz(-pi) q[2];
rz(3.0012673) q[3];
sx q[3];
rz(-1.232339) q[3];
sx q[3];
rz(2.942977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3673765) q[2];
sx q[2];
rz(-0.99952951) q[2];
sx q[2];
rz(-0.2224758) q[2];
rz(-3.0761278) q[3];
sx q[3];
rz(-1.2832578) q[3];
sx q[3];
rz(2.2953667) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9533933) q[0];
sx q[0];
rz(-3.1173752) q[0];
sx q[0];
rz(-0.54704332) q[0];
rz(2.8833) q[1];
sx q[1];
rz(-3.1196085) q[1];
sx q[1];
rz(0.34150728) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19778457) q[0];
sx q[0];
rz(-3.1349234) q[0];
sx q[0];
rz(-1.3640553) q[0];
x q[1];
rz(-1.9893622) q[2];
sx q[2];
rz(-1.2593049) q[2];
sx q[2];
rz(0.77889393) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.472802) q[1];
sx q[1];
rz(-2.3040132) q[1];
sx q[1];
rz(1.1146783) q[1];
rz(2.6686752) q[3];
sx q[3];
rz(-1.1120136) q[3];
sx q[3];
rz(-0.44675752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.42939886) q[2];
sx q[2];
rz(-1.2960351) q[2];
sx q[2];
rz(-0.94924259) q[2];
rz(-2.3927169) q[3];
sx q[3];
rz(-1.2810992) q[3];
sx q[3];
rz(0.062189814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
rz(0.6147989) q[0];
sx q[0];
rz(-0.034788046) q[0];
sx q[0];
rz(1.590796) q[0];
rz(-1.7855478) q[1];
sx q[1];
rz(-0.0043914774) q[1];
sx q[1];
rz(-3.0785676) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7558799) q[0];
sx q[0];
rz(-1.5405419) q[0];
sx q[0];
rz(-3.0755733) q[0];
rz(2.3895415) q[2];
sx q[2];
rz(-1.1630485) q[2];
sx q[2];
rz(2.4501981) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.121513) q[1];
sx q[1];
rz(-0.20993349) q[1];
sx q[1];
rz(-1.6990183) q[1];
x q[2];
rz(-2.6659417) q[3];
sx q[3];
rz(-0.84053549) q[3];
sx q[3];
rz(-2.2851839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.65681347) q[2];
sx q[2];
rz(-1.8317089) q[2];
sx q[2];
rz(0.64378929) q[2];
rz(0.8748318) q[3];
sx q[3];
rz(-2.8372786) q[3];
sx q[3];
rz(0.8005825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0463882) q[0];
sx q[0];
rz(-0.055618532) q[0];
sx q[0];
rz(-2.7860506) q[0];
rz(-0.19861673) q[1];
sx q[1];
rz(-3.1348517) q[1];
sx q[1];
rz(-0.14828646) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3902238) q[0];
sx q[0];
rz(-1.7538188) q[0];
sx q[0];
rz(1.5733294) q[0];
x q[1];
rz(0.4112377) q[2];
sx q[2];
rz(-1.6401575) q[2];
sx q[2];
rz(-2.3313521) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0394888) q[1];
sx q[1];
rz(-1.436427) q[1];
sx q[1];
rz(2.9062382) q[1];
x q[2];
rz(-0.96252302) q[3];
sx q[3];
rz(-1.4279162) q[3];
sx q[3];
rz(-2.8382728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6716914) q[2];
sx q[2];
rz(-0.24158676) q[2];
sx q[2];
rz(0.12413231) q[2];
rz(-0.5747253) q[3];
sx q[3];
rz(-2.997213) q[3];
sx q[3];
rz(-0.1709443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2588876) q[0];
sx q[0];
rz(-0.12508617) q[0];
sx q[0];
rz(-2.4001154) q[0];
rz(2.8575836) q[1];
sx q[1];
rz(-3.1378855) q[1];
sx q[1];
rz(0.31518087) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18703546) q[0];
sx q[0];
rz(-1.5439568) q[0];
sx q[0];
rz(1.5052273) q[0];
x q[1];
rz(-1.1229321) q[2];
sx q[2];
rz(-1.7583946) q[2];
sx q[2];
rz(-0.59567398) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6446723) q[1];
sx q[1];
rz(-0.61750353) q[1];
sx q[1];
rz(-0.27568494) q[1];
rz(1.4676845) q[3];
sx q[3];
rz(-0.22480837) q[3];
sx q[3];
rz(-1.6112427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.86889851) q[2];
sx q[2];
rz(-1.0947451) q[2];
sx q[2];
rz(0.72186738) q[2];
rz(-2.7591211) q[3];
sx q[3];
rz(-1.9964652) q[3];
sx q[3];
rz(1.0684048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5751936) q[0];
sx q[0];
rz(-3.1167751) q[0];
sx q[0];
rz(-1.5665293) q[0];
rz(0.20340915) q[1];
sx q[1];
rz(-1.8433488) q[1];
sx q[1];
rz(-0.64483109) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6015907) q[0];
sx q[0];
rz(-2.9338957) q[0];
sx q[0];
rz(1.6112441) q[0];
x q[1];
rz(3.105509) q[2];
sx q[2];
rz(-0.60257327) q[2];
sx q[2];
rz(0.35843231) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.37721021) q[1];
sx q[1];
rz(-1.6444211) q[1];
sx q[1];
rz(-1.0318569) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4025492) q[3];
sx q[3];
rz(-0.62646455) q[3];
sx q[3];
rz(1.6221969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5975534) q[2];
sx q[2];
rz(-2.7880703) q[2];
sx q[2];
rz(-0.37975797) q[2];
rz(2.0960268) q[3];
sx q[3];
rz(-1.2328204) q[3];
sx q[3];
rz(1.1988962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3658635) q[0];
sx q[0];
rz(-3.1079223) q[0];
sx q[0];
rz(-1.7849543) q[0];
rz(2.7011073) q[1];
sx q[1];
rz(-1.0904652) q[1];
sx q[1];
rz(2.4408565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25862274) q[0];
sx q[0];
rz(-2.077266) q[0];
sx q[0];
rz(-2.3372786) q[0];
x q[1];
rz(-1.8092625) q[2];
sx q[2];
rz(-0.70435134) q[2];
sx q[2];
rz(0.65504247) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8776013) q[1];
sx q[1];
rz(-0.7670247) q[1];
sx q[1];
rz(-0.33543889) q[1];
rz(-pi) q[2];
rz(1.5666991) q[3];
sx q[3];
rz(-1.1381686) q[3];
sx q[3];
rz(3.0790902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35357722) q[2];
sx q[2];
rz(-2.7699296) q[2];
sx q[2];
rz(-1.8549982) q[2];
rz(-2.6364117) q[3];
sx q[3];
rz(-0.44613871) q[3];
sx q[3];
rz(-1.9193468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5196359) q[0];
sx q[0];
rz(-0.049787909) q[0];
sx q[0];
rz(-1.5995481) q[0];
rz(-0.75625769) q[1];
sx q[1];
rz(-0.007096346) q[1];
sx q[1];
rz(-2.8047628) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23754691) q[0];
sx q[0];
rz(-1.2531452) q[0];
sx q[0];
rz(1.5685515) q[0];
rz(2.9354503) q[2];
sx q[2];
rz(-1.2954324) q[2];
sx q[2];
rz(-2.5554267) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.017180786) q[1];
sx q[1];
rz(-1.4596816) q[1];
sx q[1];
rz(0.085787048) q[1];
rz(-pi) q[2];
x q[2];
rz(2.182951) q[3];
sx q[3];
rz(-2.4869182) q[3];
sx q[3];
rz(2.502113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51356641) q[2];
sx q[2];
rz(-0.94945532) q[2];
sx q[2];
rz(2.8619518) q[2];
rz(-0.63129342) q[3];
sx q[3];
rz(-0.93835962) q[3];
sx q[3];
rz(2.5173371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8414128) q[0];
sx q[0];
rz(-1.5501839) q[0];
sx q[0];
rz(-1.3612904) q[0];
rz(2.367876) q[1];
sx q[1];
rz(-2.5061889) q[1];
sx q[1];
rz(-2.9255964) q[1];
rz(0.86376376) q[2];
sx q[2];
rz(-2.2939199) q[2];
sx q[2];
rz(-1.6242956) q[2];
rz(-1.5480883) q[3];
sx q[3];
rz(-1.1946214) q[3];
sx q[3];
rz(-0.0047636845) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
