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
rz(1.9396012) q[0];
sx q[0];
rz(-0.48299462) q[0];
sx q[0];
rz(-1.510409) q[0];
rz(-0.13934879) q[1];
sx q[1];
rz(5.723602) q[1];
sx q[1];
rz(8.6948123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48929991) q[0];
sx q[0];
rz(-2.518939) q[0];
sx q[0];
rz(1.3176509) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59487409) q[2];
sx q[2];
rz(-1.7584137) q[2];
sx q[2];
rz(1.7930195) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3986721) q[1];
sx q[1];
rz(-1.6026596) q[1];
sx q[1];
rz(0.58018654) q[1];
x q[2];
rz(2.9128051) q[3];
sx q[3];
rz(-1.2902033) q[3];
sx q[3];
rz(0.56786637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0240747) q[2];
sx q[2];
rz(-2.8705609) q[2];
sx q[2];
rz(2.3226341) q[2];
rz(-3.135318) q[3];
sx q[3];
rz(-1.2333906) q[3];
sx q[3];
rz(2.1591469) q[3];
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
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821424) q[0];
sx q[0];
rz(-1.5511976) q[0];
sx q[0];
rz(-0.5994125) q[0];
rz(-0.86743152) q[1];
sx q[1];
rz(-2.1116833) q[1];
sx q[1];
rz(0.74554602) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2363385) q[0];
sx q[0];
rz(-2.8239692) q[0];
sx q[0];
rz(2.7228217) q[0];
rz(2.8023984) q[2];
sx q[2];
rz(-0.789398) q[2];
sx q[2];
rz(-1.2123002) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6709969) q[1];
sx q[1];
rz(-1.7299011) q[1];
sx q[1];
rz(-2.5039423) q[1];
rz(0.87941951) q[3];
sx q[3];
rz(-1.9718277) q[3];
sx q[3];
rz(1.5123617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6856689) q[2];
sx q[2];
rz(-2.8440639) q[2];
sx q[2];
rz(-1.7342742) q[2];
rz(-0.84434858) q[3];
sx q[3];
rz(-0.83365369) q[3];
sx q[3];
rz(-1.0367397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5582964) q[0];
sx q[0];
rz(-1.1953657) q[0];
sx q[0];
rz(1.012828) q[0];
rz(-0.60802513) q[1];
sx q[1];
rz(-1.5826179) q[1];
sx q[1];
rz(-1.8720522) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0758666) q[0];
sx q[0];
rz(-1.1320496) q[0];
sx q[0];
rz(2.7164396) q[0];
x q[1];
rz(-2.24176) q[2];
sx q[2];
rz(-1.6353893) q[2];
sx q[2];
rz(-3.10499) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0996272) q[1];
sx q[1];
rz(-1.7668952) q[1];
sx q[1];
rz(-2.4134497) q[1];
x q[2];
rz(-0.57165159) q[3];
sx q[3];
rz(-2.5749675) q[3];
sx q[3];
rz(-1.0039312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6262007) q[2];
sx q[2];
rz(-2.0833368) q[2];
sx q[2];
rz(1.2916279) q[2];
rz(-3.1332704) q[3];
sx q[3];
rz(-0.95299995) q[3];
sx q[3];
rz(0.21361175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0030647) q[0];
sx q[0];
rz(-2.5862638) q[0];
sx q[0];
rz(1.3767161) q[0];
rz(1.5273013) q[1];
sx q[1];
rz(-1.2034028) q[1];
sx q[1];
rz(-2.6208904) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6277058) q[0];
sx q[0];
rz(-1.9235652) q[0];
sx q[0];
rz(1.1350495) q[0];
rz(-pi) q[1];
rz(2.7815656) q[2];
sx q[2];
rz(-1.0485149) q[2];
sx q[2];
rz(1.6619267) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6286875) q[1];
sx q[1];
rz(-2.4081552) q[1];
sx q[1];
rz(-2.4705486) q[1];
x q[2];
rz(2.2124452) q[3];
sx q[3];
rz(-0.7693253) q[3];
sx q[3];
rz(-1.3206467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.5248096) q[2];
sx q[2];
rz(-2.3526134) q[2];
sx q[2];
rz(1.3516124) q[2];
rz(-2.1221519) q[3];
sx q[3];
rz(-0.56988684) q[3];
sx q[3];
rz(-1.1714237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5919507) q[0];
sx q[0];
rz(-1.6787981) q[0];
sx q[0];
rz(3.0426262) q[0];
rz(0.70676604) q[1];
sx q[1];
rz(-0.87044972) q[1];
sx q[1];
rz(-2.6720572) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3833475) q[0];
sx q[0];
rz(-1.7086141) q[0];
sx q[0];
rz(3.1175749) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37920538) q[2];
sx q[2];
rz(-1.0673041) q[2];
sx q[2];
rz(0.45631726) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8626735) q[1];
sx q[1];
rz(-1.1871183) q[1];
sx q[1];
rz(2.5088599) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.585586) q[3];
sx q[3];
rz(-1.4463498) q[3];
sx q[3];
rz(-0.57226244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2609451) q[2];
sx q[2];
rz(-1.3346883) q[2];
sx q[2];
rz(-1.3060695) q[2];
rz(-0.4246873) q[3];
sx q[3];
rz(-1.3771907) q[3];
sx q[3];
rz(0.88596058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(3.0223618) q[0];
sx q[0];
rz(-0.62218085) q[0];
sx q[0];
rz(-1.8192044) q[0];
rz(-1.127683) q[1];
sx q[1];
rz(-1.5195945) q[1];
sx q[1];
rz(1.3816396) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7776426) q[0];
sx q[0];
rz(-0.76350437) q[0];
sx q[0];
rz(1.4457153) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53333111) q[2];
sx q[2];
rz(-1.0600277) q[2];
sx q[2];
rz(-3.0549218) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7834251) q[1];
sx q[1];
rz(-1.7562215) q[1];
sx q[1];
rz(-0.76134759) q[1];
rz(-pi) q[2];
rz(0.23412496) q[3];
sx q[3];
rz(-2.8322599) q[3];
sx q[3];
rz(-0.36946378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3809001) q[2];
sx q[2];
rz(-1.5194632) q[2];
sx q[2];
rz(-2.6603928) q[2];
rz(0.5091269) q[3];
sx q[3];
rz(-0.23734084) q[3];
sx q[3];
rz(0.38207644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5803489) q[0];
sx q[0];
rz(-0.48155293) q[0];
sx q[0];
rz(0.089381889) q[0];
rz(-2.0629758) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(-0.99348974) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.908113) q[0];
sx q[0];
rz(-2.3589239) q[0];
sx q[0];
rz(-0.054305768) q[0];
x q[1];
rz(2.2682796) q[2];
sx q[2];
rz(-0.9160348) q[2];
sx q[2];
rz(-2.8755434) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.74718432) q[1];
sx q[1];
rz(-1.4239053) q[1];
sx q[1];
rz(1.9373075) q[1];
x q[2];
rz(2.5105623) q[3];
sx q[3];
rz(-2.0027805) q[3];
sx q[3];
rz(-0.61842266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3420458) q[2];
sx q[2];
rz(-0.61508238) q[2];
sx q[2];
rz(2.7395524) q[2];
rz(-1.0953995) q[3];
sx q[3];
rz(-1.2145372) q[3];
sx q[3];
rz(2.5636165) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47356975) q[0];
sx q[0];
rz(-2.3325925) q[0];
sx q[0];
rz(-1.3336257) q[0];
rz(-0.85583055) q[1];
sx q[1];
rz(-1.7355093) q[1];
sx q[1];
rz(-0.91748253) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3429392) q[0];
sx q[0];
rz(-2.9173033) q[0];
sx q[0];
rz(-2.0464315) q[0];
rz(-3.0227674) q[2];
sx q[2];
rz(-0.95044327) q[2];
sx q[2];
rz(0.22681776) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.92261295) q[1];
sx q[1];
rz(-2.726788) q[1];
sx q[1];
rz(0.99402438) q[1];
rz(-pi) q[2];
rz(2.0794271) q[3];
sx q[3];
rz(-0.98782238) q[3];
sx q[3];
rz(2.3140918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5069919) q[2];
sx q[2];
rz(-1.736182) q[2];
sx q[2];
rz(1.290192) q[2];
rz(-2.2318132) q[3];
sx q[3];
rz(-1.6981643) q[3];
sx q[3];
rz(-3.1006052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23583394) q[0];
sx q[0];
rz(-1.9941149) q[0];
sx q[0];
rz(2.8644417) q[0];
rz(1.2241036) q[1];
sx q[1];
rz(-1.6033019) q[1];
sx q[1];
rz(1.7701497) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7200206) q[0];
sx q[0];
rz(-2.4457481) q[0];
sx q[0];
rz(-0.30465845) q[0];
rz(-pi) q[1];
rz(-2.0061699) q[2];
sx q[2];
rz(-2.0757542) q[2];
sx q[2];
rz(1.4471042) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4491357) q[1];
sx q[1];
rz(-0.68329158) q[1];
sx q[1];
rz(-1.9016674) q[1];
x q[2];
rz(1.3731433) q[3];
sx q[3];
rz(-2.4098793) q[3];
sx q[3];
rz(-0.55303516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93854967) q[2];
sx q[2];
rz(-0.86157346) q[2];
sx q[2];
rz(-0.68823632) q[2];
rz(-0.31442434) q[3];
sx q[3];
rz(-2.7721072) q[3];
sx q[3];
rz(-1.8800053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37453434) q[0];
sx q[0];
rz(-0.86807591) q[0];
sx q[0];
rz(0.57149291) q[0];
rz(-0.68069619) q[1];
sx q[1];
rz(-1.9711767) q[1];
sx q[1];
rz(-0.64819711) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2637973) q[0];
sx q[0];
rz(-1.5476942) q[0];
sx q[0];
rz(0.013851555) q[0];
rz(-pi) q[1];
rz(0.86096455) q[2];
sx q[2];
rz(-0.87536821) q[2];
sx q[2];
rz(-0.26317138) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.56420621) q[1];
sx q[1];
rz(-1.5930452) q[1];
sx q[1];
rz(0.69907) q[1];
rz(-2.6246214) q[3];
sx q[3];
rz(-1.3992157) q[3];
sx q[3];
rz(-1.3763381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8615243) q[2];
sx q[2];
rz(-1.0736977) q[2];
sx q[2];
rz(-1.6746707) q[2];
rz(2.9649949) q[3];
sx q[3];
rz(-0.64591518) q[3];
sx q[3];
rz(-1.8471898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27547729) q[0];
sx q[0];
rz(-1.3963516) q[0];
sx q[0];
rz(1.8657952) q[0];
rz(-0.38446174) q[1];
sx q[1];
rz(-1.6356331) q[1];
sx q[1];
rz(2.5148139) q[1];
rz(-2.7914417) q[2];
sx q[2];
rz(-1.2731009) q[2];
sx q[2];
rz(0.44558744) q[2];
rz(-2.1351142) q[3];
sx q[3];
rz(-2.0340393) q[3];
sx q[3];
rz(-1.5599193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
