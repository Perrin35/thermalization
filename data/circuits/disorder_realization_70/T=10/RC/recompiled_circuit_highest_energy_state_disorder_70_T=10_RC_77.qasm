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
rz(1.7477859) q[0];
sx q[0];
rz(-2.2414247) q[0];
sx q[0];
rz(1.409344) q[0];
rz(0.70977587) q[1];
sx q[1];
rz(3.4748454) q[1];
sx q[1];
rz(9.8069053) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47412518) q[0];
sx q[0];
rz(-0.84185266) q[0];
sx q[0];
rz(-1.390662) q[0];
rz(-pi) q[1];
rz(-3.0067031) q[2];
sx q[2];
rz(-1.3185878) q[2];
sx q[2];
rz(1.8602961) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.00013) q[1];
sx q[1];
rz(-1.9381818) q[1];
sx q[1];
rz(0.39333435) q[1];
x q[2];
rz(0.87859359) q[3];
sx q[3];
rz(-1.6631563) q[3];
sx q[3];
rz(-1.1042995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0673151) q[2];
sx q[2];
rz(-1.5019608) q[2];
sx q[2];
rz(3.066646) q[2];
rz(1.2429169) q[3];
sx q[3];
rz(-2.8800745) q[3];
sx q[3];
rz(0.3956795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13260929) q[0];
sx q[0];
rz(-2.7217396) q[0];
sx q[0];
rz(-0.62927759) q[0];
rz(-2.3043326) q[1];
sx q[1];
rz(-2.3910797) q[1];
sx q[1];
rz(0.57964051) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0934187) q[0];
sx q[0];
rz(-0.77147986) q[0];
sx q[0];
rz(-0.29399612) q[0];
rz(-pi) q[1];
x q[1];
rz(0.57168269) q[2];
sx q[2];
rz(-2.6808156) q[2];
sx q[2];
rz(-0.95651115) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.77387735) q[1];
sx q[1];
rz(-1.9019777) q[1];
sx q[1];
rz(-2.7298688) q[1];
x q[2];
rz(2.9778538) q[3];
sx q[3];
rz(-1.5511912) q[3];
sx q[3];
rz(-2.1480008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.89977449) q[2];
sx q[2];
rz(-2.9782229) q[2];
sx q[2];
rz(2.3972798) q[2];
rz(-0.36738473) q[3];
sx q[3];
rz(-1.2920047) q[3];
sx q[3];
rz(1.415409) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.431417) q[0];
sx q[0];
rz(-2.1911868) q[0];
sx q[0];
rz(1.0071734) q[0];
rz(-0.011064359) q[1];
sx q[1];
rz(-2.8304351) q[1];
sx q[1];
rz(-0.99753582) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1384093) q[0];
sx q[0];
rz(-1.54502) q[0];
sx q[0];
rz(1.4332214) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11409161) q[2];
sx q[2];
rz(-0.81159752) q[2];
sx q[2];
rz(1.6676192) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.026873206) q[1];
sx q[1];
rz(-0.51843868) q[1];
sx q[1];
rz(0.015109574) q[1];
rz(-0.51181958) q[3];
sx q[3];
rz(-1.5580172) q[3];
sx q[3];
rz(2.6731773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.02115383) q[2];
sx q[2];
rz(-2.8595371) q[2];
sx q[2];
rz(2.2998478) q[2];
rz(-2.4833931) q[3];
sx q[3];
rz(-2.2704312) q[3];
sx q[3];
rz(3.1200718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4206674) q[0];
sx q[0];
rz(-0.1249211) q[0];
sx q[0];
rz(0.7290054) q[0];
rz(-0.76338243) q[1];
sx q[1];
rz(-2.5114676) q[1];
sx q[1];
rz(0.36300945) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7621146) q[0];
sx q[0];
rz(-1.861137) q[0];
sx q[0];
rz(-1.3069673) q[0];
rz(-pi) q[1];
rz(0.55018376) q[2];
sx q[2];
rz(-2.2386058) q[2];
sx q[2];
rz(-0.097152348) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22717366) q[1];
sx q[1];
rz(-2.8110831) q[1];
sx q[1];
rz(2.6918489) q[1];
x q[2];
rz(0.12673817) q[3];
sx q[3];
rz(-1.5180382) q[3];
sx q[3];
rz(2.3873752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9558941) q[2];
sx q[2];
rz(-0.82315767) q[2];
sx q[2];
rz(1.4996747) q[2];
rz(0.58756346) q[3];
sx q[3];
rz(-2.1524119) q[3];
sx q[3];
rz(0.51830083) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3041621) q[0];
sx q[0];
rz(-2.4346133) q[0];
sx q[0];
rz(-0.28513232) q[0];
rz(0.25310165) q[1];
sx q[1];
rz(-2.1123835) q[1];
sx q[1];
rz(-1.0458127) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2946729) q[0];
sx q[0];
rz(-0.6631279) q[0];
sx q[0];
rz(-0.66594932) q[0];
rz(2.3836993) q[2];
sx q[2];
rz(-0.57595384) q[2];
sx q[2];
rz(2.8754701) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6599413) q[1];
sx q[1];
rz(-1.5143947) q[1];
sx q[1];
rz(2.7475471) q[1];
x q[2];
rz(2.8941514) q[3];
sx q[3];
rz(-0.69290042) q[3];
sx q[3];
rz(-0.019236658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73695856) q[2];
sx q[2];
rz(-1.0598695) q[2];
sx q[2];
rz(0.57445478) q[2];
rz(-2.5308841) q[3];
sx q[3];
rz(-2.6210531) q[3];
sx q[3];
rz(-1.0569388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8059998) q[0];
sx q[0];
rz(-1.3845504) q[0];
sx q[0];
rz(0.35032508) q[0];
rz(-2.8490745) q[1];
sx q[1];
rz(-3.0151093) q[1];
sx q[1];
rz(-0.77936053) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.821601) q[0];
sx q[0];
rz(-1.4519094) q[0];
sx q[0];
rz(-0.058811896) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33213385) q[2];
sx q[2];
rz(-2.6382425) q[2];
sx q[2];
rz(2.533874) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36679572) q[1];
sx q[1];
rz(-1.8908943) q[1];
sx q[1];
rz(1.3572378) q[1];
rz(-pi) q[2];
rz(2.8661714) q[3];
sx q[3];
rz(-0.53715992) q[3];
sx q[3];
rz(-2.8583683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.23955762) q[2];
sx q[2];
rz(-1.6263447) q[2];
sx q[2];
rz(2.1774192) q[2];
rz(-2.9686109) q[3];
sx q[3];
rz(-2.1292584) q[3];
sx q[3];
rz(3.0103736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59703374) q[0];
sx q[0];
rz(-1.1976765) q[0];
sx q[0];
rz(-0.069393754) q[0];
rz(1.3736877) q[1];
sx q[1];
rz(-1.2488139) q[1];
sx q[1];
rz(-0.51838851) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4189506) q[0];
sx q[0];
rz(-1.7403474) q[0];
sx q[0];
rz(-0.8954312) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47641192) q[2];
sx q[2];
rz(-2.0134543) q[2];
sx q[2];
rz(-1.2141643) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.035210755) q[1];
sx q[1];
rz(-1.7089881) q[1];
sx q[1];
rz(1.5968678) q[1];
rz(-1.5126918) q[3];
sx q[3];
rz(-1.4062506) q[3];
sx q[3];
rz(-1.8376415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.18846506) q[2];
sx q[2];
rz(-1.943482) q[2];
sx q[2];
rz(0.24448621) q[2];
rz(-1.9237579) q[3];
sx q[3];
rz(-3.0471424) q[3];
sx q[3];
rz(0.34734669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0511047) q[0];
sx q[0];
rz(-0.29260391) q[0];
sx q[0];
rz(-0.69945139) q[0];
rz(0.72169101) q[1];
sx q[1];
rz(-1.7803918) q[1];
sx q[1];
rz(-1.1915421) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22868294) q[0];
sx q[0];
rz(-1.3587202) q[0];
sx q[0];
rz(1.4430844) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9279153) q[2];
sx q[2];
rz(-1.4136633) q[2];
sx q[2];
rz(-0.84136277) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4556132) q[1];
sx q[1];
rz(-2.1081104) q[1];
sx q[1];
rz(-0.098453589) q[1];
x q[2];
rz(-0.44765391) q[3];
sx q[3];
rz(-1.8869683) q[3];
sx q[3];
rz(-0.74861909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2193853) q[2];
sx q[2];
rz(-2.3228513) q[2];
sx q[2];
rz(1.4156263) q[2];
rz(0.43736941) q[3];
sx q[3];
rz(-0.24961095) q[3];
sx q[3];
rz(2.6387446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1943844) q[0];
sx q[0];
rz(-1.2057065) q[0];
sx q[0];
rz(-2.6956287) q[0];
rz(-1.3392316) q[1];
sx q[1];
rz(-2.4868592) q[1];
sx q[1];
rz(1.9816678) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8561607) q[0];
sx q[0];
rz(-0.92575476) q[0];
sx q[0];
rz(1.1354394) q[0];
rz(3.0662698) q[2];
sx q[2];
rz(-1.3100071) q[2];
sx q[2];
rz(1.0950077) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36148188) q[1];
sx q[1];
rz(-2.105753) q[1];
sx q[1];
rz(1.8423716) q[1];
x q[2];
rz(-2.2307617) q[3];
sx q[3];
rz(-1.0928705) q[3];
sx q[3];
rz(-2.7895989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7356073) q[2];
sx q[2];
rz(-2.1129825) q[2];
sx q[2];
rz(-0.73060161) q[2];
rz(0.90100151) q[3];
sx q[3];
rz(-0.42947072) q[3];
sx q[3];
rz(-0.034339529) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5071097) q[0];
sx q[0];
rz(-2.6431838) q[0];
sx q[0];
rz(0.46257567) q[0];
rz(-2.0787461) q[1];
sx q[1];
rz(-1.7741508) q[1];
sx q[1];
rz(-3.0752693) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0339113) q[0];
sx q[0];
rz(-2.4723791) q[0];
sx q[0];
rz(2.6341556) q[0];
rz(-pi) q[1];
x q[1];
rz(0.1758258) q[2];
sx q[2];
rz(-1.9011902) q[2];
sx q[2];
rz(1.8595472) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5196701) q[1];
sx q[1];
rz(-1.038045) q[1];
sx q[1];
rz(-2.4581535) q[1];
rz(-pi) q[2];
rz(0.22610449) q[3];
sx q[3];
rz(-1.1563099) q[3];
sx q[3];
rz(-1.7630486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1041717) q[2];
sx q[2];
rz(-0.53356844) q[2];
sx q[2];
rz(1.9083692) q[2];
rz(0.45796606) q[3];
sx q[3];
rz(-0.27126867) q[3];
sx q[3];
rz(2.7428194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0290699) q[0];
sx q[0];
rz(-1.4334913) q[0];
sx q[0];
rz(1.7423472) q[0];
rz(0.15432547) q[1];
sx q[1];
rz(-1.7185153) q[1];
sx q[1];
rz(-1.1850866) q[1];
rz(-0.90171705) q[2];
sx q[2];
rz(-2.6606885) q[2];
sx q[2];
rz(-0.056405141) q[2];
rz(0.48404398) q[3];
sx q[3];
rz(-2.0917907) q[3];
sx q[3];
rz(1.2863822) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
