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
rz(3.0022439) q[1];
sx q[1];
rz(-2.5820093) q[1];
sx q[1];
rz(-2.411627) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18099526) q[0];
sx q[0];
rz(-2.1707524) q[0];
sx q[0];
rz(2.9636895) q[0];
x q[1];
rz(-1.7961305) q[2];
sx q[2];
rz(-2.1538434) q[2];
sx q[2];
rz(3.0449113) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0180359) q[1];
sx q[1];
rz(-0.58096051) q[1];
sx q[1];
rz(3.0835129) q[1];
rz(-pi) q[2];
rz(0.9040622) q[3];
sx q[3];
rz(-0.36012563) q[3];
sx q[3];
rz(-1.2670497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0240747) q[2];
sx q[2];
rz(-0.2710318) q[2];
sx q[2];
rz(2.3226341) q[2];
rz(3.135318) q[3];
sx q[3];
rz(-1.9082021) q[3];
sx q[3];
rz(-0.98244572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5821424) q[0];
sx q[0];
rz(-1.5903951) q[0];
sx q[0];
rz(2.5421802) q[0];
rz(-2.2741611) q[1];
sx q[1];
rz(-1.0299094) q[1];
sx q[1];
rz(-2.3960466) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2363385) q[0];
sx q[0];
rz(-2.8239692) q[0];
sx q[0];
rz(2.7228217) q[0];
rz(-0.7600766) q[2];
sx q[2];
rz(-1.3323297) q[2];
sx q[2];
rz(-2.5395405) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.47059575) q[1];
sx q[1];
rz(-1.4116916) q[1];
sx q[1];
rz(-0.63765031) q[1];
rz(-pi) q[2];
rz(2.638444) q[3];
sx q[3];
rz(-0.9434349) q[3];
sx q[3];
rz(-0.37093758) q[3];
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
rz(2.2972441) q[3];
sx q[3];
rz(-0.83365369) q[3];
sx q[3];
rz(2.1048529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5582964) q[0];
sx q[0];
rz(-1.946227) q[0];
sx q[0];
rz(-1.012828) q[0];
rz(2.5335675) q[1];
sx q[1];
rz(-1.5826179) q[1];
sx q[1];
rz(1.2695405) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4465049) q[0];
sx q[0];
rz(-1.9534612) q[0];
sx q[0];
rz(2.0464567) q[0];
x q[1];
rz(-0.082398675) q[2];
sx q[2];
rz(-2.2401056) q[2];
sx q[2];
rz(-1.5562039) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.041965466) q[1];
sx q[1];
rz(-1.3746975) q[1];
sx q[1];
rz(-2.4134497) q[1];
rz(-pi) q[2];
rz(-0.57165159) q[3];
sx q[3];
rz(-0.56662512) q[3];
sx q[3];
rz(-2.1376615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.515392) q[2];
sx q[2];
rz(-1.0582558) q[2];
sx q[2];
rz(-1.8499648) q[2];
rz(-0.0083222566) q[3];
sx q[3];
rz(-0.95299995) q[3];
sx q[3];
rz(-0.21361175) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0030647) q[0];
sx q[0];
rz(-2.5862638) q[0];
sx q[0];
rz(1.3767161) q[0];
rz(-1.6142913) q[1];
sx q[1];
rz(-1.9381899) q[1];
sx q[1];
rz(-0.52070224) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5138869) q[0];
sx q[0];
rz(-1.9235652) q[0];
sx q[0];
rz(-1.1350495) q[0];
rz(-2.7815656) q[2];
sx q[2];
rz(-2.0930778) q[2];
sx q[2];
rz(-1.479666) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5508073) q[1];
sx q[1];
rz(-1.1414764) q[1];
sx q[1];
rz(-2.5270259) q[1];
x q[2];
rz(2.6163231) q[3];
sx q[3];
rz(-2.161918) q[3];
sx q[3];
rz(2.1256465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.5248096) q[2];
sx q[2];
rz(-0.7889792) q[2];
sx q[2];
rz(-1.7899803) q[2];
rz(2.1221519) q[3];
sx q[3];
rz(-0.56988684) q[3];
sx q[3];
rz(1.1714237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.549642) q[0];
sx q[0];
rz(-1.6787981) q[0];
sx q[0];
rz(3.0426262) q[0];
rz(-2.4348266) q[1];
sx q[1];
rz(-0.87044972) q[1];
sx q[1];
rz(0.4695355) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5851368) q[0];
sx q[0];
rz(-3.0017108) q[0];
sx q[0];
rz(1.742247) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0355613) q[2];
sx q[2];
rz(-1.9010086) q[2];
sx q[2];
rz(2.2170628) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3217056) q[1];
sx q[1];
rz(-0.72607909) q[1];
sx q[1];
rz(-0.5989845) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55600663) q[3];
sx q[3];
rz(-1.4463498) q[3];
sx q[3];
rz(2.5693302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2609451) q[2];
sx q[2];
rz(-1.3346883) q[2];
sx q[2];
rz(-1.3060695) q[2];
rz(-0.4246873) q[3];
sx q[3];
rz(-1.7644019) q[3];
sx q[3];
rz(2.2556321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11923085) q[0];
sx q[0];
rz(-2.5194118) q[0];
sx q[0];
rz(-1.8192044) q[0];
rz(-1.127683) q[1];
sx q[1];
rz(-1.6219982) q[1];
sx q[1];
rz(-1.3816396) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19162543) q[0];
sx q[0];
rz(-0.81474308) q[0];
sx q[0];
rz(-0.11884584) q[0];
rz(-pi) q[1];
rz(-0.83397978) q[2];
sx q[2];
rz(-0.72089855) q[2];
sx q[2];
rz(-2.3490459) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.35816757) q[1];
sx q[1];
rz(-1.3853711) q[1];
sx q[1];
rz(-0.76134759) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30140169) q[3];
sx q[3];
rz(-1.5001138) q[3];
sx q[3];
rz(1.4247198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7606925) q[2];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5612438) q[0];
sx q[0];
rz(-2.6600397) q[0];
sx q[0];
rz(-3.0522108) q[0];
rz(2.0629758) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(-2.1481029) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3099965) q[0];
sx q[0];
rz(-0.78959268) q[0];
sx q[0];
rz(1.6247276) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69665945) q[2];
sx q[2];
rz(-2.2242745) q[2];
sx q[2];
rz(-0.6763263) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.76748874) q[1];
sx q[1];
rz(-1.933177) q[1];
sx q[1];
rz(-0.15717536) q[1];
x q[2];
rz(-2.5105623) q[3];
sx q[3];
rz(-2.0027805) q[3];
sx q[3];
rz(0.61842266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3420458) q[2];
sx q[2];
rz(-2.5265103) q[2];
sx q[2];
rz(2.7395524) q[2];
rz(-2.0461931) q[3];
sx q[3];
rz(-1.9270555) q[3];
sx q[3];
rz(2.5636165) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47356975) q[0];
sx q[0];
rz(-2.3325925) q[0];
sx q[0];
rz(1.8079669) q[0];
rz(0.85583055) q[1];
sx q[1];
rz(-1.4060833) q[1];
sx q[1];
rz(2.2241101) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3429392) q[0];
sx q[0];
rz(-0.22428939) q[0];
sx q[0];
rz(2.0464315) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4063666) q[2];
sx q[2];
rz(-2.5114369) q[2];
sx q[2];
rz(-0.024261628) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0304347) q[1];
sx q[1];
rz(-1.7923755) q[1];
sx q[1];
rz(1.2171924) q[1];
rz(-1.0621655) q[3];
sx q[3];
rz(-2.1537703) q[3];
sx q[3];
rz(0.82750083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5069919) q[2];
sx q[2];
rz(-1.736182) q[2];
sx q[2];
rz(-1.8514006) q[2];
rz(-0.90977943) q[3];
sx q[3];
rz(-1.4434283) q[3];
sx q[3];
rz(0.040987404) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9057587) q[0];
sx q[0];
rz(-1.1474778) q[0];
sx q[0];
rz(-2.8644417) q[0];
rz(1.9174891) q[1];
sx q[1];
rz(-1.6033019) q[1];
sx q[1];
rz(-1.7701497) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4215721) q[0];
sx q[0];
rz(-0.69584457) q[0];
sx q[0];
rz(-2.8369342) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4898289) q[2];
sx q[2];
rz(-2.4874176) q[2];
sx q[2];
rz(-0.68133611) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7595948) q[1];
sx q[1];
rz(-1.7773668) q[1];
sx q[1];
rz(-2.2269511) q[1];
rz(-pi) q[2];
x q[2];
rz(0.17454608) q[3];
sx q[3];
rz(-0.8564328) q[3];
sx q[3];
rz(-2.3256231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.93854967) q[2];
sx q[2];
rz(-0.86157346) q[2];
sx q[2];
rz(-0.68823632) q[2];
rz(0.31442434) q[3];
sx q[3];
rz(-0.36948547) q[3];
sx q[3];
rz(-1.8800053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7670583) q[0];
sx q[0];
rz(-0.86807591) q[0];
sx q[0];
rz(0.57149291) q[0];
rz(2.4608965) q[1];
sx q[1];
rz(-1.1704159) q[1];
sx q[1];
rz(0.64819711) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72361904) q[0];
sx q[0];
rz(-0.026935808) q[0];
sx q[0];
rz(2.1108146) q[0];
x q[1];
rz(-0.86096455) q[2];
sx q[2];
rz(-0.87536821) q[2];
sx q[2];
rz(-2.8784213) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.025291) q[1];
sx q[1];
rz(-2.2696583) q[1];
sx q[1];
rz(1.5998597) q[1];
rz(-1.7675507) q[3];
sx q[3];
rz(-2.0794387) q[3];
sx q[3];
rz(3.0439049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8615243) q[2];
sx q[2];
rz(-2.0678949) q[2];
sx q[2];
rz(-1.466922) q[2];
rz(-2.9649949) q[3];
sx q[3];
rz(-0.64591518) q[3];
sx q[3];
rz(-1.2944029) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27547729) q[0];
sx q[0];
rz(-1.3963516) q[0];
sx q[0];
rz(1.8657952) q[0];
rz(2.7571309) q[1];
sx q[1];
rz(-1.6356331) q[1];
sx q[1];
rz(2.5148139) q[1];
rz(-2.4118829) q[2];
sx q[2];
rz(-0.45558246) q[2];
sx q[2];
rz(2.6930553) q[2];
rz(0.53388673) q[3];
sx q[3];
rz(-2.0697099) q[3];
sx q[3];
rz(2.8768215) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
