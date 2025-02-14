OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.14804949) q[0];
sx q[0];
rz(-0.43819675) q[0];
sx q[0];
rz(10.155805) q[0];
rz(-2.464715) q[1];
sx q[1];
rz(-0.69898611) q[1];
sx q[1];
rz(-2.5875523) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.06023887) q[0];
sx q[0];
rz(-1.0805655) q[0];
sx q[0];
rz(-1.0361363) q[0];
rz(1.0921539) q[2];
sx q[2];
rz(-1.3502099) q[2];
sx q[2];
rz(1.5696862) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89008625) q[1];
sx q[1];
rz(-2.4533669) q[1];
sx q[1];
rz(1.5176956) q[1];
rz(-pi) q[2];
rz(-0.4783614) q[3];
sx q[3];
rz(-0.94650062) q[3];
sx q[3];
rz(-2.9251298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76250166) q[2];
sx q[2];
rz(-1.5662301) q[2];
sx q[2];
rz(-2.9762034) q[2];
rz(-0.23073828) q[3];
sx q[3];
rz(-2.8985891) q[3];
sx q[3];
rz(2.2802584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47925258) q[0];
sx q[0];
rz(-0.32821822) q[0];
sx q[0];
rz(1.7829371) q[0];
rz(-1.6183629) q[1];
sx q[1];
rz(-0.75444573) q[1];
sx q[1];
rz(-2.3014136) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1217499) q[0];
sx q[0];
rz(-0.54980924) q[0];
sx q[0];
rz(2.0998831) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39583037) q[2];
sx q[2];
rz(-2.2326222) q[2];
sx q[2];
rz(-1.0368376) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.11586861) q[1];
sx q[1];
rz(-0.77889538) q[1];
sx q[1];
rz(-0.42515691) q[1];
rz(-pi) q[2];
rz(1.4235557) q[3];
sx q[3];
rz(-2.3309532) q[3];
sx q[3];
rz(-2.7642706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.75300616) q[2];
sx q[2];
rz(-1.9084088) q[2];
sx q[2];
rz(2.2071655) q[2];
rz(1.2855444) q[3];
sx q[3];
rz(-0.779874) q[3];
sx q[3];
rz(0.88750315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0282106) q[0];
sx q[0];
rz(-1.0209571) q[0];
sx q[0];
rz(1.5895948) q[0];
rz(-1.7734843) q[1];
sx q[1];
rz(-1.869447) q[1];
sx q[1];
rz(2.0210463) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89884886) q[0];
sx q[0];
rz(-1.1087497) q[0];
sx q[0];
rz(0.36926134) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2539961) q[2];
sx q[2];
rz(-1.996737) q[2];
sx q[2];
rz(-0.48522247) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2255133) q[1];
sx q[1];
rz(-2.597993) q[1];
sx q[1];
rz(-2.9617642) q[1];
rz(1.8561826) q[3];
sx q[3];
rz(-2.719398) q[3];
sx q[3];
rz(-0.58870047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0997194) q[2];
sx q[2];
rz(-2.9283044) q[2];
sx q[2];
rz(-0.29207692) q[2];
rz(1.9000351) q[3];
sx q[3];
rz(-1.7009267) q[3];
sx q[3];
rz(0.45274538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.1555772) q[0];
sx q[0];
rz(-0.41322511) q[0];
sx q[0];
rz(-2.7098932) q[0];
rz(-0.15402928) q[1];
sx q[1];
rz(-1.5700211) q[1];
sx q[1];
rz(2.0808751) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4765324) q[0];
sx q[0];
rz(-1.8899584) q[0];
sx q[0];
rz(1.9299797) q[0];
rz(-pi) q[1];
rz(-1.8155104) q[2];
sx q[2];
rz(-1.7438466) q[2];
sx q[2];
rz(-1.3428549) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.65217291) q[1];
sx q[1];
rz(-1.8967581) q[1];
sx q[1];
rz(-2.9054545) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4134365) q[3];
sx q[3];
rz(-1.3829136) q[3];
sx q[3];
rz(1.431501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3920307) q[2];
sx q[2];
rz(-1.531484) q[2];
sx q[2];
rz(0.57382601) q[2];
rz(3.0841893) q[3];
sx q[3];
rz(-1.4797689) q[3];
sx q[3];
rz(-2.3315232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5350128) q[0];
sx q[0];
rz(-2.8856394) q[0];
sx q[0];
rz(-2.8888597) q[0];
rz(-0.17929721) q[1];
sx q[1];
rz(-1.7959692) q[1];
sx q[1];
rz(-1.7326694) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94570075) q[0];
sx q[0];
rz(-1.333814) q[0];
sx q[0];
rz(-3.0653565) q[0];
rz(-1.1834572) q[2];
sx q[2];
rz(-1.172386) q[2];
sx q[2];
rz(0.18926316) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9086011) q[1];
sx q[1];
rz(-1.4484693) q[1];
sx q[1];
rz(-0.78559383) q[1];
rz(-1.9305655) q[3];
sx q[3];
rz(-1.0749987) q[3];
sx q[3];
rz(-1.825765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1285642) q[2];
sx q[2];
rz(-1.3821673) q[2];
sx q[2];
rz(-1.2010835) q[2];
rz(2.3342093) q[3];
sx q[3];
rz(-1.3599334) q[3];
sx q[3];
rz(-1.132563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0627237) q[0];
sx q[0];
rz(-1.6428592) q[0];
sx q[0];
rz(2.3003182) q[0];
rz(1.9367283) q[1];
sx q[1];
rz(-1.3179904) q[1];
sx q[1];
rz(-1.0296317) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1618274) q[0];
sx q[0];
rz(-2.0812391) q[0];
sx q[0];
rz(0.011988601) q[0];
rz(-pi) q[1];
rz(-0.51635833) q[2];
sx q[2];
rz(-1.8022924) q[2];
sx q[2];
rz(-1.4524492) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8358546) q[1];
sx q[1];
rz(-1.5785913) q[1];
sx q[1];
rz(-2.9628997) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5978686) q[3];
sx q[3];
rz(-0.35055509) q[3];
sx q[3];
rz(-1.9026827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.131375) q[2];
sx q[2];
rz(-2.3848332) q[2];
sx q[2];
rz(-0.49794623) q[2];
rz(0.52305269) q[3];
sx q[3];
rz(-1.4191351) q[3];
sx q[3];
rz(0.12652346) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2987357) q[0];
sx q[0];
rz(-0.13164483) q[0];
sx q[0];
rz(2.329622) q[0];
rz(0.89093351) q[1];
sx q[1];
rz(-1.9782601) q[1];
sx q[1];
rz(-2.1422211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97702128) q[0];
sx q[0];
rz(-2.0334822) q[0];
sx q[0];
rz(-2.0624731) q[0];
rz(-pi) q[1];
rz(-1.6738724) q[2];
sx q[2];
rz(-1.4749881) q[2];
sx q[2];
rz(2.9745548) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.68518406) q[1];
sx q[1];
rz(-2.051609) q[1];
sx q[1];
rz(-2.494704) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8511397) q[3];
sx q[3];
rz(-0.92602611) q[3];
sx q[3];
rz(1.5034663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.823395) q[2];
sx q[2];
rz(-0.91788936) q[2];
sx q[2];
rz(1.2981752) q[2];
rz(-0.038481742) q[3];
sx q[3];
rz(-2.0352071) q[3];
sx q[3];
rz(-1.6160256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0549952) q[0];
sx q[0];
rz(-1.1308068) q[0];
sx q[0];
rz(0.31563345) q[0];
rz(1.9075958) q[1];
sx q[1];
rz(-0.71593586) q[1];
sx q[1];
rz(1.2589781) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6427073) q[0];
sx q[0];
rz(-2.2872637) q[0];
sx q[0];
rz(-1.9509208) q[0];
rz(-1.744049) q[2];
sx q[2];
rz(-0.8129186) q[2];
sx q[2];
rz(3.0305741) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96871829) q[1];
sx q[1];
rz(-0.58191381) q[1];
sx q[1];
rz(2.3496778) q[1];
rz(-3.0397814) q[3];
sx q[3];
rz(-1.57668) q[3];
sx q[3];
rz(2.7525097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.2235609) q[2];
sx q[2];
rz(-2.3980902) q[2];
sx q[2];
rz(3.1136759) q[2];
rz(2.3935086) q[3];
sx q[3];
rz(-1.7223765) q[3];
sx q[3];
rz(-2.9885651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7276723) q[0];
sx q[0];
rz(-2.1509009) q[0];
sx q[0];
rz(-2.682611) q[0];
rz(-1.1363632) q[1];
sx q[1];
rz(-1.4308948) q[1];
sx q[1];
rz(2.9232025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1735575) q[0];
sx q[0];
rz(-1.5474678) q[0];
sx q[0];
rz(-1.5541881) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2480466) q[2];
sx q[2];
rz(-2.7392929) q[2];
sx q[2];
rz(0.93076555) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7228411) q[1];
sx q[1];
rz(-1.4237849) q[1];
sx q[1];
rz(-0.49572368) q[1];
rz(-pi) q[2];
rz(1.8828527) q[3];
sx q[3];
rz(-0.40570212) q[3];
sx q[3];
rz(-2.4108374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8454664) q[2];
sx q[2];
rz(-2.8870236) q[2];
sx q[2];
rz(2.0938342) q[2];
rz(2.3962077) q[3];
sx q[3];
rz(-1.5785917) q[3];
sx q[3];
rz(1.6489702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02136852) q[0];
sx q[0];
rz(-2.5701018) q[0];
sx q[0];
rz(-2.7987203) q[0];
rz(-0.19035569) q[1];
sx q[1];
rz(-2.1326667) q[1];
sx q[1];
rz(-0.51317936) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29268713) q[0];
sx q[0];
rz(-0.93528895) q[0];
sx q[0];
rz(-2.9440895) q[0];
x q[1];
rz(-0.97810271) q[2];
sx q[2];
rz(-2.5524271) q[2];
sx q[2];
rz(1.1617253) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.045195508) q[1];
sx q[1];
rz(-1.1151114) q[1];
sx q[1];
rz(1.1229188) q[1];
rz(0.14182183) q[3];
sx q[3];
rz(-1.1537227) q[3];
sx q[3];
rz(0.50167044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43789431) q[2];
sx q[2];
rz(-0.63592211) q[2];
sx q[2];
rz(-2.5305914) q[2];
rz(2.8827187) q[3];
sx q[3];
rz(-1.2114108) q[3];
sx q[3];
rz(-1.3795616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7543058) q[0];
sx q[0];
rz(-1.5504693) q[0];
sx q[0];
rz(-1.5691527) q[0];
rz(0.98987956) q[1];
sx q[1];
rz(-1.2073333) q[1];
sx q[1];
rz(-2.5055199) q[1];
rz(0.0089916747) q[2];
sx q[2];
rz(-1.7594382) q[2];
sx q[2];
rz(0.3450763) q[2];
rz(-1.5045139) q[3];
sx q[3];
rz(-0.87285973) q[3];
sx q[3];
rz(-1.9535337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
