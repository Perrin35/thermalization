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
rz(-0.67796081) q[0];
sx q[0];
rz(-0.27422658) q[0];
sx q[0];
rz(1.8507313) q[0];
rz(-3.3554606) q[1];
sx q[1];
rz(5.9670347) q[1];
sx q[1];
rz(11.066758) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3363269) q[0];
sx q[0];
rz(-1.299843) q[0];
sx q[0];
rz(-1.3265557) q[0];
rz(-pi) q[1];
rz(-2.847509) q[2];
sx q[2];
rz(-2.3399647) q[2];
sx q[2];
rz(-1.557136) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1823481) q[1];
sx q[1];
rz(-0.79234353) q[1];
sx q[1];
rz(0.46087973) q[1];
rz(0.47874449) q[3];
sx q[3];
rz(-0.9968206) q[3];
sx q[3];
rz(-2.6285302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9700254) q[2];
sx q[2];
rz(-2.1000523) q[2];
sx q[2];
rz(-2.2402666) q[2];
rz(-2.6775635) q[3];
sx q[3];
rz(-1.8601067) q[3];
sx q[3];
rz(-3.071781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0483911) q[0];
sx q[0];
rz(-0.036660107) q[0];
sx q[0];
rz(-1.2777591) q[0];
rz(0.094206421) q[1];
sx q[1];
rz(-0.46368805) q[1];
sx q[1];
rz(-1.6069848) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8600991) q[0];
sx q[0];
rz(-1.5380385) q[0];
sx q[0];
rz(-0.97086914) q[0];
x q[1];
rz(1.7530551) q[2];
sx q[2];
rz(-1.6267952) q[2];
sx q[2];
rz(-0.1268498) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0329735) q[1];
sx q[1];
rz(-2.0619644) q[1];
sx q[1];
rz(-1.2742548) q[1];
x q[2];
rz(0.8022763) q[3];
sx q[3];
rz(-2.5893988) q[3];
sx q[3];
rz(1.7084165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4216807) q[2];
sx q[2];
rz(-1.9399425) q[2];
sx q[2];
rz(-2.9812532) q[2];
rz(2.4028589) q[3];
sx q[3];
rz(-0.84638798) q[3];
sx q[3];
rz(-1.7140478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5147603) q[0];
sx q[0];
rz(-1.7034096) q[0];
sx q[0];
rz(-2.6318188) q[0];
rz(-1.7104507) q[1];
sx q[1];
rz(-1.1809843) q[1];
sx q[1];
rz(-2.3950155) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7187481) q[0];
sx q[0];
rz(-1.5507409) q[0];
sx q[0];
rz(1.0182747) q[0];
x q[1];
rz(0.66303894) q[2];
sx q[2];
rz(-2.2929077) q[2];
sx q[2];
rz(2.7977365) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4164734) q[1];
sx q[1];
rz(-1.2095569) q[1];
sx q[1];
rz(-0.84049417) q[1];
rz(-pi) q[2];
rz(0.38807773) q[3];
sx q[3];
rz(-2.2142234) q[3];
sx q[3];
rz(0.79659792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9868682) q[2];
sx q[2];
rz(-2.8747323) q[2];
sx q[2];
rz(1.223684) q[2];
rz(2.0679421) q[3];
sx q[3];
rz(-1.4370388) q[3];
sx q[3];
rz(-2.2435718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8680442) q[0];
sx q[0];
rz(-0.88669625) q[0];
sx q[0];
rz(-0.37937382) q[0];
rz(-0.1768449) q[1];
sx q[1];
rz(-1.6601446) q[1];
sx q[1];
rz(-2.3462229) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5569755) q[0];
sx q[0];
rz(-1.4780227) q[0];
sx q[0];
rz(1.2014821) q[0];
rz(-pi) q[1];
rz(1.8537117) q[2];
sx q[2];
rz(-2.3672315) q[2];
sx q[2];
rz(1.1224358) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8898192) q[1];
sx q[1];
rz(-1.3135513) q[1];
sx q[1];
rz(2.1104269) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0159608) q[3];
sx q[3];
rz(-1.5141308) q[3];
sx q[3];
rz(1.980933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.4312326) q[2];
sx q[2];
rz(-1.7962339) q[2];
sx q[2];
rz(-1.7809407) q[2];
rz(2.1121173) q[3];
sx q[3];
rz(-1.6607213) q[3];
sx q[3];
rz(1.3220538) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4904356) q[0];
sx q[0];
rz(-1.54162) q[0];
sx q[0];
rz(-1.5929476) q[0];
rz(0.40052888) q[1];
sx q[1];
rz(-1.488204) q[1];
sx q[1];
rz(-1.4097479) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1574788) q[0];
sx q[0];
rz(-1.3712008) q[0];
sx q[0];
rz(2.4019813) q[0];
rz(0.39938853) q[2];
sx q[2];
rz(-0.3872954) q[2];
sx q[2];
rz(-0.072810955) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9292007) q[1];
sx q[1];
rz(-0.98941313) q[1];
sx q[1];
rz(2.6049032) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4758167) q[3];
sx q[3];
rz(-2.5174369) q[3];
sx q[3];
rz(-0.74172663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3132402) q[2];
sx q[2];
rz(-1.7007549) q[2];
sx q[2];
rz(-0.43133119) q[2];
rz(-3.0865772) q[3];
sx q[3];
rz(-2.8300245) q[3];
sx q[3];
rz(-0.55606786) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3738275) q[0];
sx q[0];
rz(-2.8887833) q[0];
sx q[0];
rz(2.1164236) q[0];
rz(-0.45285666) q[1];
sx q[1];
rz(-0.57944524) q[1];
sx q[1];
rz(-2.2377009) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9867803) q[0];
sx q[0];
rz(-0.93789369) q[0];
sx q[0];
rz(2.7557719) q[0];
rz(-2.507693) q[2];
sx q[2];
rz(-1.0154775) q[2];
sx q[2];
rz(-2.7125396) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3749668) q[1];
sx q[1];
rz(-1.8410676) q[1];
sx q[1];
rz(-0.75281669) q[1];
rz(-pi) q[2];
rz(-2.970231) q[3];
sx q[3];
rz(-1.4079861) q[3];
sx q[3];
rz(0.57933116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30770939) q[2];
sx q[2];
rz(-1.4557975) q[2];
sx q[2];
rz(-2.9252388) q[2];
rz(0.89573914) q[3];
sx q[3];
rz(-2.4560865) q[3];
sx q[3];
rz(-1.640813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1322587) q[0];
sx q[0];
rz(-1.1067156) q[0];
sx q[0];
rz(2.4427781) q[0];
rz(1.9578594) q[1];
sx q[1];
rz(-1.4212757) q[1];
sx q[1];
rz(-0.85174495) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4119371) q[0];
sx q[0];
rz(-1.0348399) q[0];
sx q[0];
rz(-0.044117731) q[0];
rz(-pi) q[1];
x q[1];
rz(0.039489517) q[2];
sx q[2];
rz(-2.563884) q[2];
sx q[2];
rz(-2.5414987) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.70515436) q[1];
sx q[1];
rz(-0.90696224) q[1];
sx q[1];
rz(1.1052119) q[1];
rz(-0.64166358) q[3];
sx q[3];
rz(-1.9681276) q[3];
sx q[3];
rz(-1.7375377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8375887) q[2];
sx q[2];
rz(-1.8303266) q[2];
sx q[2];
rz(-2.6336929) q[2];
rz(-1.0287644) q[3];
sx q[3];
rz(-0.84380904) q[3];
sx q[3];
rz(0.60104162) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69865882) q[0];
sx q[0];
rz(-1.826094) q[0];
sx q[0];
rz(-3.1319295) q[0];
rz(-1.5254321) q[1];
sx q[1];
rz(-1.3776255) q[1];
sx q[1];
rz(-2.756871) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5949109) q[0];
sx q[0];
rz(-1.4230886) q[0];
sx q[0];
rz(1.3315931) q[0];
rz(-pi) q[1];
rz(-1.184395) q[2];
sx q[2];
rz(-1.3829872) q[2];
sx q[2];
rz(-1.1593429) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2710451) q[1];
sx q[1];
rz(-1.7779852) q[1];
sx q[1];
rz(-1.2882339) q[1];
x q[2];
rz(0.77120292) q[3];
sx q[3];
rz(-1.1386385) q[3];
sx q[3];
rz(1.766253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.559451) q[2];
sx q[2];
rz(-0.95953512) q[2];
sx q[2];
rz(-3.1257296) q[2];
rz(1.9150241) q[3];
sx q[3];
rz(-0.80569402) q[3];
sx q[3];
rz(-2.5198643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3847619) q[0];
sx q[0];
rz(-0.66487304) q[0];
sx q[0];
rz(-1.0913947) q[0];
rz(2.3233844) q[1];
sx q[1];
rz(-1.810377) q[1];
sx q[1];
rz(2.4726726) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92907897) q[0];
sx q[0];
rz(-0.90796472) q[0];
sx q[0];
rz(0.79848358) q[0];
rz(2.0537655) q[2];
sx q[2];
rz(-2.3685799) q[2];
sx q[2];
rz(0.20394606) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3260169) q[1];
sx q[1];
rz(-1.4998124) q[1];
sx q[1];
rz(-2.1287588) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1534641) q[3];
sx q[3];
rz(-1.8496397) q[3];
sx q[3];
rz(-3.0089965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6341256) q[2];
sx q[2];
rz(-0.52906817) q[2];
sx q[2];
rz(0.96366209) q[2];
rz(-1.4108747) q[3];
sx q[3];
rz(-1.4286634) q[3];
sx q[3];
rz(1.7760407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1249579) q[0];
sx q[0];
rz(-0.5683012) q[0];
sx q[0];
rz(-1.6868663) q[0];
rz(1.1085054) q[1];
sx q[1];
rz(-1.5364372) q[1];
sx q[1];
rz(0.32807168) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0963225) q[0];
sx q[0];
rz(-0.12841405) q[0];
sx q[0];
rz(-1.1044109) q[0];
x q[1];
rz(-0.180822) q[2];
sx q[2];
rz(-1.6852813) q[2];
sx q[2];
rz(-0.78071981) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0340424) q[1];
sx q[1];
rz(-0.75354105) q[1];
sx q[1];
rz(0.57489245) q[1];
rz(-2.1222018) q[3];
sx q[3];
rz(-0.20272045) q[3];
sx q[3];
rz(-2.3644476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6803117) q[2];
sx q[2];
rz(-1.2993456) q[2];
sx q[2];
rz(-0.4168365) q[2];
rz(-2.4428115) q[3];
sx q[3];
rz(-0.67658934) q[3];
sx q[3];
rz(-0.39066395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1998491) q[0];
sx q[0];
rz(-1.9244292) q[0];
sx q[0];
rz(1.695965) q[0];
rz(-2.4327714) q[1];
sx q[1];
rz(-0.91239057) q[1];
sx q[1];
rz(-0.053587996) q[1];
rz(0.93408549) q[2];
sx q[2];
rz(-1.7333422) q[2];
sx q[2];
rz(-2.588803) q[2];
rz(2.3117456) q[3];
sx q[3];
rz(-1.1329069) q[3];
sx q[3];
rz(-2.2923242) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
