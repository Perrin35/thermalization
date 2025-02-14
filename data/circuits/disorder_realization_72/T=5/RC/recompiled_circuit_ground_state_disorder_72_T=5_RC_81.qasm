OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.44719625) q[0];
sx q[0];
rz(-1.5153272) q[0];
sx q[0];
rz(0.84994999) q[0];
rz(1.3810459) q[1];
sx q[1];
rz(-0.84264207) q[1];
sx q[1];
rz(-3.091264) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2272826) q[0];
sx q[0];
rz(-1.7512055) q[0];
sx q[0];
rz(-1.4854027) q[0];
rz(-pi) q[1];
rz(-0.16762208) q[2];
sx q[2];
rz(-0.95778886) q[2];
sx q[2];
rz(-2.5000664) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.23230339) q[1];
sx q[1];
rz(-2.9348845) q[1];
sx q[1];
rz(2.7951393) q[1];
rz(-pi) q[2];
rz(-2.5011236) q[3];
sx q[3];
rz(-1.5239118) q[3];
sx q[3];
rz(-1.3848927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3927268) q[2];
sx q[2];
rz(-1.8957081) q[2];
sx q[2];
rz(2.5457814) q[2];
rz(0.81123224) q[3];
sx q[3];
rz(-1.2657284) q[3];
sx q[3];
rz(-0.59734145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8314464) q[0];
sx q[0];
rz(-0.50066384) q[0];
sx q[0];
rz(1.3432107) q[0];
rz(1.7559747) q[1];
sx q[1];
rz(-1.132553) q[1];
sx q[1];
rz(-2.0766506) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8174521) q[0];
sx q[0];
rz(-1.5498501) q[0];
sx q[0];
rz(-0.073168427) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3716613) q[2];
sx q[2];
rz(-2.0624447) q[2];
sx q[2];
rz(2.2624598) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8044195) q[1];
sx q[1];
rz(-1.7111756) q[1];
sx q[1];
rz(-1.5765203) q[1];
x q[2];
rz(-0.029964793) q[3];
sx q[3];
rz(-1.6707722) q[3];
sx q[3];
rz(-2.4520017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6407577) q[2];
sx q[2];
rz(-2.1186327) q[2];
sx q[2];
rz(-0.41066059) q[2];
rz(-1.143035) q[3];
sx q[3];
rz(-0.7468907) q[3];
sx q[3];
rz(0.66543287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2086585) q[0];
sx q[0];
rz(-1.7993131) q[0];
sx q[0];
rz(0.68516532) q[0];
rz(-2.4655828) q[1];
sx q[1];
rz(-1.6513377) q[1];
sx q[1];
rz(-0.73838678) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91440141) q[0];
sx q[0];
rz(-1.8912328) q[0];
sx q[0];
rz(-2.6242424) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1072363) q[2];
sx q[2];
rz(-0.7329251) q[2];
sx q[2];
rz(2.422621) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1967839) q[1];
sx q[1];
rz(-2.9843708) q[1];
sx q[1];
rz(1.2320089) q[1];
rz(-pi) q[2];
rz(0.97279588) q[3];
sx q[3];
rz(-2.1913652) q[3];
sx q[3];
rz(2.6120294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43397778) q[2];
sx q[2];
rz(-1.6699764) q[2];
sx q[2];
rz(2.1493256) q[2];
rz(-0.17539242) q[3];
sx q[3];
rz(-2.6447191) q[3];
sx q[3];
rz(-2.8309256) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6076412) q[0];
sx q[0];
rz(-1.1247922) q[0];
sx q[0];
rz(-2.9006309) q[0];
rz(2.2843649) q[1];
sx q[1];
rz(-2.3576184) q[1];
sx q[1];
rz(-1.8704174) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4963112) q[0];
sx q[0];
rz(-1.4110312) q[0];
sx q[0];
rz(-1.1259354) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3996467) q[2];
sx q[2];
rz(-2.4782054) q[2];
sx q[2];
rz(-2.9398769) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1041455) q[1];
sx q[1];
rz(-1.1318143) q[1];
sx q[1];
rz(-2.6789915) q[1];
rz(1.9828834) q[3];
sx q[3];
rz(-0.85972393) q[3];
sx q[3];
rz(-2.3726316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.50095373) q[2];
sx q[2];
rz(-0.95879889) q[2];
sx q[2];
rz(-1.0183838) q[2];
rz(-1.1650677) q[3];
sx q[3];
rz(-0.94912052) q[3];
sx q[3];
rz(-0.44786662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.1027706) q[0];
sx q[0];
rz(-0.7907246) q[0];
sx q[0];
rz(0.08654174) q[0];
rz(1.6697007) q[1];
sx q[1];
rz(-2.4457928) q[1];
sx q[1];
rz(-2.5873628) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.897096) q[0];
sx q[0];
rz(-0.88567552) q[0];
sx q[0];
rz(1.8573631) q[0];
rz(-pi) q[1];
rz(2.9678095) q[2];
sx q[2];
rz(-2.4638468) q[2];
sx q[2];
rz(0.89231561) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1405331) q[1];
sx q[1];
rz(-2.3591318) q[1];
sx q[1];
rz(-2.4665574) q[1];
x q[2];
rz(1.0029129) q[3];
sx q[3];
rz(-2.6315303) q[3];
sx q[3];
rz(-2.3167852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.294813) q[2];
sx q[2];
rz(-2.7816935) q[2];
sx q[2];
rz(-0.93152085) q[2];
rz(-0.31342634) q[3];
sx q[3];
rz(-0.64160186) q[3];
sx q[3];
rz(2.6581367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1035006) q[0];
sx q[0];
rz(-0.8041389) q[0];
sx q[0];
rz(-1.7622129) q[0];
rz(-2.3992959) q[1];
sx q[1];
rz(-1.392044) q[1];
sx q[1];
rz(-1.0662063) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48165769) q[0];
sx q[0];
rz(-1.602442) q[0];
sx q[0];
rz(-3.0040347) q[0];
rz(-pi) q[1];
rz(1.3308238) q[2];
sx q[2];
rz(-0.53388816) q[2];
sx q[2];
rz(1.1660341) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4263947) q[1];
sx q[1];
rz(-2.4860588) q[1];
sx q[1];
rz(-1.1525201) q[1];
rz(-1.4857852) q[3];
sx q[3];
rz(-1.8411311) q[3];
sx q[3];
rz(-2.6694989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1833056) q[2];
sx q[2];
rz(-0.97465193) q[2];
sx q[2];
rz(-2.8948114) q[2];
rz(2.0731481) q[3];
sx q[3];
rz(-1.1698134) q[3];
sx q[3];
rz(2.0717513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5088155) q[0];
sx q[0];
rz(-2.1822073) q[0];
sx q[0];
rz(-0.36163133) q[0];
rz(1.5879141) q[1];
sx q[1];
rz(-1.5648774) q[1];
sx q[1];
rz(-1.9379001) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3643575) q[0];
sx q[0];
rz(-2.2493752) q[0];
sx q[0];
rz(-0.70058544) q[0];
rz(1.3839528) q[2];
sx q[2];
rz(-1.0770385) q[2];
sx q[2];
rz(-0.46390033) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7261914) q[1];
sx q[1];
rz(-2.5091689) q[1];
sx q[1];
rz(1.9167621) q[1];
x q[2];
rz(0.64048119) q[3];
sx q[3];
rz(-2.6121217) q[3];
sx q[3];
rz(2.01628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2162073) q[2];
sx q[2];
rz(-1.6758827) q[2];
sx q[2];
rz(2.1290131) q[2];
rz(-0.25977627) q[3];
sx q[3];
rz(-0.24661073) q[3];
sx q[3];
rz(2.7369505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15553661) q[0];
sx q[0];
rz(-1.2825092) q[0];
sx q[0];
rz(2.3543661) q[0];
rz(2.2975445) q[1];
sx q[1];
rz(-0.35875741) q[1];
sx q[1];
rz(-1.2740096) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9106261) q[0];
sx q[0];
rz(-1.9073124) q[0];
sx q[0];
rz(1.6331768) q[0];
x q[1];
rz(-0.34113555) q[2];
sx q[2];
rz(-2.3722509) q[2];
sx q[2];
rz(0.064521964) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40416557) q[1];
sx q[1];
rz(-2.3336172) q[1];
sx q[1];
rz(-0.91973181) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9926821) q[3];
sx q[3];
rz(-2.732548) q[3];
sx q[3];
rz(2.8623476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0969703) q[2];
sx q[2];
rz(-1.5644093) q[2];
sx q[2];
rz(0.38193646) q[2];
rz(1.1218128) q[3];
sx q[3];
rz(-0.33676454) q[3];
sx q[3];
rz(-3.0776265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0762416) q[0];
sx q[0];
rz(-0.92472804) q[0];
sx q[0];
rz(1.4724154) q[0];
rz(-2.8352101) q[1];
sx q[1];
rz(-0.52426052) q[1];
sx q[1];
rz(0.27149567) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9005405) q[0];
sx q[0];
rz(-2.0203622) q[0];
sx q[0];
rz(0.44067128) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36868544) q[2];
sx q[2];
rz(-2.2767608) q[2];
sx q[2];
rz(0.60312229) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2047386) q[1];
sx q[1];
rz(-0.7619155) q[1];
sx q[1];
rz(-0.73128766) q[1];
rz(-pi) q[2];
rz(-2.697345) q[3];
sx q[3];
rz(-1.0716728) q[3];
sx q[3];
rz(0.80330144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2740115) q[2];
sx q[2];
rz(-1.1271971) q[2];
sx q[2];
rz(2.4809044) q[2];
rz(0.88179669) q[3];
sx q[3];
rz(-0.56395689) q[3];
sx q[3];
rz(0.86625117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26447403) q[0];
sx q[0];
rz(-1.5333804) q[0];
sx q[0];
rz(1.3382753) q[0];
rz(-2.0954466) q[1];
sx q[1];
rz(-1.0271065) q[1];
sx q[1];
rz(-2.963692) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4433817) q[0];
sx q[0];
rz(-2.215235) q[0];
sx q[0];
rz(-0.81581913) q[0];
x q[1];
rz(-2.1523802) q[2];
sx q[2];
rz(-1.4228369) q[2];
sx q[2];
rz(1.6616247) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2819351) q[1];
sx q[1];
rz(-1.9023696) q[1];
sx q[1];
rz(-0.34032199) q[1];
rz(-0.33938395) q[3];
sx q[3];
rz(-2.748162) q[3];
sx q[3];
rz(0.14190488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.88343128) q[2];
sx q[2];
rz(-2.0475755) q[2];
sx q[2];
rz(2.2056313) q[2];
rz(-0.59518138) q[3];
sx q[3];
rz(-1.429957) q[3];
sx q[3];
rz(-0.87966758) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4277988) q[0];
sx q[0];
rz(-1.6811163) q[0];
sx q[0];
rz(1.2291193) q[0];
rz(-1.8145369) q[1];
sx q[1];
rz(-1.6358903) q[1];
sx q[1];
rz(-1.9930175) q[1];
rz(1.3934515) q[2];
sx q[2];
rz(-1.4589571) q[2];
sx q[2];
rz(-2.2082885) q[2];
rz(2.8758373) q[3];
sx q[3];
rz(-1.8674217) q[3];
sx q[3];
rz(-2.7056497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
