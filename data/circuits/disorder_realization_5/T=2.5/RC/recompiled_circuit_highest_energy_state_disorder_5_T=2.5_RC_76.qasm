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
rz(-0.37762168) q[0];
sx q[0];
rz(-2.7130337) q[0];
sx q[0];
rz(-0.23634401) q[0];
rz(-1.6726681) q[1];
sx q[1];
rz(-1.5863215) q[1];
sx q[1];
rz(-0.16170391) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2387213) q[0];
sx q[0];
rz(-1.4618902) q[0];
sx q[0];
rz(-0.98734537) q[0];
x q[1];
rz(-1.5125844) q[2];
sx q[2];
rz(-1.7663284) q[2];
sx q[2];
rz(-0.19042507) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.573805) q[1];
sx q[1];
rz(-1.5638788) q[1];
sx q[1];
rz(1.5439537) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6246965) q[3];
sx q[3];
rz(-0.42498838) q[3];
sx q[3];
rz(1.4864511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.24254313) q[2];
sx q[2];
rz(-2.264302) q[2];
sx q[2];
rz(-0.38929942) q[2];
rz(0.23046514) q[3];
sx q[3];
rz(-0.018229818) q[3];
sx q[3];
rz(-0.61483312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56735754) q[0];
sx q[0];
rz(-2.1901972) q[0];
sx q[0];
rz(1.6472598) q[0];
rz(1.5556473) q[1];
sx q[1];
rz(-2.9289398) q[1];
sx q[1];
rz(-1.1344604) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26274219) q[0];
sx q[0];
rz(-1.8244184) q[0];
sx q[0];
rz(2.4154759) q[0];
rz(-pi) q[1];
rz(-2.4058543) q[2];
sx q[2];
rz(-1.904907) q[2];
sx q[2];
rz(0.10057893) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6210584) q[1];
sx q[1];
rz(-2.1466594) q[1];
sx q[1];
rz(0.49382468) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5364356) q[3];
sx q[3];
rz(-0.89563771) q[3];
sx q[3];
rz(-1.1881811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0906585) q[2];
sx q[2];
rz(-2.2218573) q[2];
sx q[2];
rz(-1.301379) q[2];
rz(-2.0672412) q[3];
sx q[3];
rz(-0.31000546) q[3];
sx q[3];
rz(1.5460792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.018547) q[0];
sx q[0];
rz(-2.8392241) q[0];
sx q[0];
rz(2.5240335) q[0];
rz(-2.0410208) q[1];
sx q[1];
rz(-0.01958422) q[1];
sx q[1];
rz(-0.40357959) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27173938) q[0];
sx q[0];
rz(-1.6855441) q[0];
sx q[0];
rz(0.010557584) q[0];
rz(-pi) q[1];
rz(-2.9681272) q[2];
sx q[2];
rz(-2.6092165) q[2];
sx q[2];
rz(2.8796822) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.25572919) q[1];
sx q[1];
rz(-1.7049864) q[1];
sx q[1];
rz(1.7291452) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37637122) q[3];
sx q[3];
rz(-2.5507002) q[3];
sx q[3];
rz(1.8666238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6214211) q[2];
sx q[2];
rz(-1.6941864) q[2];
sx q[2];
rz(-0.78110313) q[2];
rz(2.805294) q[3];
sx q[3];
rz(-1.4399485) q[3];
sx q[3];
rz(2.4380016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6853365) q[0];
sx q[0];
rz(-2.6254613) q[0];
sx q[0];
rz(1.9235032) q[0];
rz(1.3457899) q[1];
sx q[1];
rz(-3.1333874) q[1];
sx q[1];
rz(1.9075314) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4505554) q[0];
sx q[0];
rz(-0.19007401) q[0];
sx q[0];
rz(0.051817373) q[0];
rz(-2.7379964) q[2];
sx q[2];
rz(-0.6535614) q[2];
sx q[2];
rz(2.5338478) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.36653301) q[1];
sx q[1];
rz(-1.9294191) q[1];
sx q[1];
rz(0.26945646) q[1];
rz(-pi) q[2];
rz(0.59438057) q[3];
sx q[3];
rz(-0.012033741) q[3];
sx q[3];
rz(-1.8569225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9464843) q[2];
sx q[2];
rz(-2.7184964) q[2];
sx q[2];
rz(1.1484324) q[2];
rz(0.75442433) q[3];
sx q[3];
rz(-1.1634588) q[3];
sx q[3];
rz(-2.8783126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(0.36963439) q[0];
sx q[0];
rz(-0.088005528) q[0];
sx q[0];
rz(-2.5391286) q[0];
rz(0.98214904) q[1];
sx q[1];
rz(-3.1352477) q[1];
sx q[1];
rz(-2.6314661) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9038123) q[0];
sx q[0];
rz(-2.8717715) q[0];
sx q[0];
rz(2.3348007) q[0];
x q[1];
rz(-3.0544364) q[2];
sx q[2];
rz(-1.4498324) q[2];
sx q[2];
rz(-0.5145413) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9822471) q[1];
sx q[1];
rz(-0.67894113) q[1];
sx q[1];
rz(0.62980501) q[1];
x q[2];
rz(-0.00058393936) q[3];
sx q[3];
rz(-2.3272132) q[3];
sx q[3];
rz(1.7581122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0977352) q[2];
sx q[2];
rz(-1.1172224) q[2];
sx q[2];
rz(-1.0350234) q[2];
rz(-1.6739316) q[3];
sx q[3];
rz(-0.36164713) q[3];
sx q[3];
rz(-2.5586832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8408836) q[0];
sx q[0];
rz(-2.1568334) q[0];
sx q[0];
rz(-1.0915407) q[0];
rz(2.7409399) q[1];
sx q[1];
rz(-3.1392097) q[1];
sx q[1];
rz(-0.56652743) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0877675) q[0];
sx q[0];
rz(-1.4364409) q[0];
sx q[0];
rz(-2.3162486) q[0];
rz(2.285901) q[2];
sx q[2];
rz(-2.1314179) q[2];
sx q[2];
rz(-2.4770346) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8174584) q[1];
sx q[1];
rz(-2.5310881) q[1];
sx q[1];
rz(-0.39892674) q[1];
rz(-0.57265307) q[3];
sx q[3];
rz(-1.476247) q[3];
sx q[3];
rz(-1.270164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1953676) q[2];
sx q[2];
rz(-0.75104284) q[2];
sx q[2];
rz(-1.6132149) q[2];
rz(1.001312) q[3];
sx q[3];
rz(-0.87037218) q[3];
sx q[3];
rz(2.9010469) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85750759) q[0];
sx q[0];
rz(-0.79428285) q[0];
sx q[0];
rz(1.1450144) q[0];
rz(1.5223632) q[1];
sx q[1];
rz(-3.122819) q[1];
sx q[1];
rz(-1.1633263) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49646851) q[0];
sx q[0];
rz(-1.2815406) q[0];
sx q[0];
rz(0.048286322) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2448213) q[2];
sx q[2];
rz(-1.687703) q[2];
sx q[2];
rz(-0.74806556) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.12194828) q[1];
sx q[1];
rz(-0.62133978) q[1];
sx q[1];
rz(1.7985733) q[1];
rz(-pi) q[2];
rz(1.4520423) q[3];
sx q[3];
rz(-0.94014478) q[3];
sx q[3];
rz(-1.5662367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6254977) q[2];
sx q[2];
rz(-2.1558546) q[2];
sx q[2];
rz(0.37334785) q[2];
rz(0.18800023) q[3];
sx q[3];
rz(-1.5550273) q[3];
sx q[3];
rz(-1.9787623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2242551) q[0];
sx q[0];
rz(-0.52977109) q[0];
sx q[0];
rz(-0.24714558) q[0];
rz(2.9180134) q[1];
sx q[1];
rz(-0.0037184628) q[1];
sx q[1];
rz(1.6252958) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.926439) q[0];
sx q[0];
rz(-1.5064539) q[0];
sx q[0];
rz(-1.6872726) q[0];
rz(-pi) q[1];
rz(1.9886279) q[2];
sx q[2];
rz(-2.2925809) q[2];
sx q[2];
rz(1.6642451) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5471955) q[1];
sx q[1];
rz(-1.5548348) q[1];
sx q[1];
rz(2.7475824) q[1];
x q[2];
rz(-0.070930158) q[3];
sx q[3];
rz(-2.5462357) q[3];
sx q[3];
rz(-1.8728674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.94459263) q[2];
sx q[2];
rz(-0.71114117) q[2];
sx q[2];
rz(-1.5468583) q[2];
rz(-2.366015) q[3];
sx q[3];
rz(-1.8577134) q[3];
sx q[3];
rz(-1.4625134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32770661) q[0];
sx q[0];
rz(-2.1729108) q[0];
sx q[0];
rz(1.1073329) q[0];
rz(2.2164717) q[1];
sx q[1];
rz(-0.0019625891) q[1];
sx q[1];
rz(2.3855239) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30413142) q[0];
sx q[0];
rz(-0.63429773) q[0];
sx q[0];
rz(1.0360121) q[0];
rz(-pi) q[1];
rz(-1.0399516) q[2];
sx q[2];
rz(-1.565298) q[2];
sx q[2];
rz(-1.0571684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9097164) q[1];
sx q[1];
rz(-1.8617587) q[1];
sx q[1];
rz(-2.5582486) q[1];
rz(-2.954079) q[3];
sx q[3];
rz(-0.57996677) q[3];
sx q[3];
rz(2.7051276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6206616) q[2];
sx q[2];
rz(-2.2563917) q[2];
sx q[2];
rz(-1.0010285) q[2];
rz(1.3641317) q[3];
sx q[3];
rz(-0.93327156) q[3];
sx q[3];
rz(2.6076243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.119568) q[0];
sx q[0];
rz(-1.3725932) q[0];
sx q[0];
rz(0.46510988) q[0];
rz(-1.3348835) q[1];
sx q[1];
rz(-0.3723793) q[1];
sx q[1];
rz(-1.5660657) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50349277) q[0];
sx q[0];
rz(-1.3412807) q[0];
sx q[0];
rz(-3.0969308) q[0];
x q[1];
rz(2.8303873) q[2];
sx q[2];
rz(-1.0619783) q[2];
sx q[2];
rz(2.6151163) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.59908341) q[1];
sx q[1];
rz(-1.570382) q[1];
sx q[1];
rz(3.1387657) q[1];
rz(-pi) q[2];
rz(0.09327831) q[3];
sx q[3];
rz(-0.3808379) q[3];
sx q[3];
rz(0.24254984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.89432013) q[2];
sx q[2];
rz(-0.044581052) q[2];
sx q[2];
rz(1.0856005) q[2];
rz(1.9191437) q[3];
sx q[3];
rz(-0.51210755) q[3];
sx q[3];
rz(1.8634169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4492252) q[0];
sx q[0];
rz(-1.4341555) q[0];
sx q[0];
rz(1.7644802) q[0];
rz(-1.5583246) q[1];
sx q[1];
rz(-2.227034) q[1];
sx q[1];
rz(-2.9169678) q[1];
rz(-0.10075154) q[2];
sx q[2];
rz(-1.5755079) q[2];
sx q[2];
rz(1.8625349) q[2];
rz(-1.1034154) q[3];
sx q[3];
rz(-2.6513908) q[3];
sx q[3];
rz(-1.8309616) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
