OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.9935432) q[0];
sx q[0];
rz(-2.7033959) q[0];
sx q[0];
rz(-0.73102695) q[0];
rz(-2.464715) q[1];
sx q[1];
rz(-0.69898611) q[1];
sx q[1];
rz(-2.5875523) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9029516) q[0];
sx q[0];
rz(-1.1045505) q[0];
sx q[0];
rz(0.55517261) q[0];
x q[1];
rz(-0.24744769) q[2];
sx q[2];
rz(-1.1046865) q[2];
sx q[2];
rz(3.0296536) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.89008625) q[1];
sx q[1];
rz(-0.68822574) q[1];
sx q[1];
rz(-1.5176956) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1393806) q[3];
sx q[3];
rz(-2.3750288) q[3];
sx q[3];
rz(2.1994906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76250166) q[2];
sx q[2];
rz(-1.5753626) q[2];
sx q[2];
rz(-2.9762034) q[2];
rz(-2.9108544) q[3];
sx q[3];
rz(-0.24300353) q[3];
sx q[3];
rz(2.2802584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47925258) q[0];
sx q[0];
rz(-0.32821822) q[0];
sx q[0];
rz(1.3586556) q[0];
rz(-1.5232297) q[1];
sx q[1];
rz(-0.75444573) q[1];
sx q[1];
rz(2.3014136) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.019842783) q[0];
sx q[0];
rz(-2.5917834) q[0];
sx q[0];
rz(1.0417095) q[0];
x q[1];
rz(-0.39583037) q[2];
sx q[2];
rz(-2.2326222) q[2];
sx q[2];
rz(1.0368376) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3749373) q[1];
sx q[1];
rz(-1.8647653) q[1];
sx q[1];
rz(0.73237441) q[1];
x q[2];
rz(-0.76559098) q[3];
sx q[3];
rz(-1.4642707) q[3];
sx q[3];
rz(-1.2953143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75300616) q[2];
sx q[2];
rz(-1.9084088) q[2];
sx q[2];
rz(2.2071655) q[2];
rz(-1.8560483) q[3];
sx q[3];
rz(-2.3617187) q[3];
sx q[3];
rz(-0.88750315) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0282106) q[0];
sx q[0];
rz(-2.1206355) q[0];
sx q[0];
rz(1.5519979) q[0];
rz(-1.7734843) q[1];
sx q[1];
rz(-1.2721456) q[1];
sx q[1];
rz(1.1205463) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50110589) q[0];
sx q[0];
rz(-1.2418134) q[0];
sx q[0];
rz(1.0803243) q[0];
x q[1];
rz(2.0765614) q[2];
sx q[2];
rz(-2.6496844) q[2];
sx q[2];
rz(3.0658403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.016213266) q[1];
sx q[1];
rz(-2.1046791) q[1];
sx q[1];
rz(1.4631229) q[1];
rz(-1.2854101) q[3];
sx q[3];
rz(-0.42219463) q[3];
sx q[3];
rz(0.58870047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0418732) q[2];
sx q[2];
rz(-0.21328829) q[2];
sx q[2];
rz(0.29207692) q[2];
rz(-1.9000351) q[3];
sx q[3];
rz(-1.440666) q[3];
sx q[3];
rz(-2.6888473) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1555772) q[0];
sx q[0];
rz(-2.7283675) q[0];
sx q[0];
rz(2.7098932) q[0];
rz(2.9875634) q[1];
sx q[1];
rz(-1.5700211) q[1];
sx q[1];
rz(2.0808751) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4765324) q[0];
sx q[0];
rz(-1.2516343) q[0];
sx q[0];
rz(1.9299797) q[0];
rz(2.9633386) q[2];
sx q[2];
rz(-1.3298103) q[2];
sx q[2];
rz(-2.9566233) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4894197) q[1];
sx q[1];
rz(-1.8967581) q[1];
sx q[1];
rz(0.23613813) q[1];
rz(-0.19017724) q[3];
sx q[3];
rz(-1.4162283) q[3];
sx q[3];
rz(-3.0319253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.74956191) q[2];
sx q[2];
rz(-1.531484) q[2];
sx q[2];
rz(0.57382601) q[2];
rz(-3.0841893) q[3];
sx q[3];
rz(-1.6618238) q[3];
sx q[3];
rz(0.81006947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(2.5350128) q[0];
sx q[0];
rz(-2.8856394) q[0];
sx q[0];
rz(-0.25273299) q[0];
rz(0.17929721) q[1];
sx q[1];
rz(-1.7959692) q[1];
sx q[1];
rz(-1.4089233) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4985663) q[0];
sx q[0];
rz(-1.6448978) q[0];
sx q[0];
rz(1.8084433) q[0];
rz(2.7149212) q[2];
sx q[2];
rz(-1.2152078) q[2];
sx q[2];
rz(1.6030965) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.682413) q[1];
sx q[1];
rz(-2.3485567) q[1];
sx q[1];
rz(-0.17211087) q[1];
rz(-pi) q[2];
rz(-1.9305655) q[3];
sx q[3];
rz(-1.0749987) q[3];
sx q[3];
rz(-1.825765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0130284) q[2];
sx q[2];
rz(-1.7594254) q[2];
sx q[2];
rz(1.2010835) q[2];
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
rz(0.078868911) q[0];
sx q[0];
rz(-1.4987334) q[0];
sx q[0];
rz(2.3003182) q[0];
rz(1.9367283) q[1];
sx q[1];
rz(-1.8236022) q[1];
sx q[1];
rz(-2.1119609) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58517367) q[0];
sx q[0];
rz(-1.5812567) q[0];
sx q[0];
rz(2.0812698) q[0];
rz(-pi) q[1];
rz(1.835498) q[2];
sx q[2];
rz(-1.0695188) q[2];
sx q[2];
rz(-3.1304111) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8358546) q[1];
sx q[1];
rz(-1.5785913) q[1];
sx q[1];
rz(-2.9628997) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7577517) q[3];
sx q[3];
rz(-1.8690946) q[3];
sx q[3];
rz(-1.3307856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.131375) q[2];
sx q[2];
rz(-0.75675941) q[2];
sx q[2];
rz(-0.49794623) q[2];
rz(0.52305269) q[3];
sx q[3];
rz(-1.4191351) q[3];
sx q[3];
rz(-3.0150692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2987357) q[0];
sx q[0];
rz(-3.0099478) q[0];
sx q[0];
rz(2.329622) q[0];
rz(2.2506591) q[1];
sx q[1];
rz(-1.9782601) q[1];
sx q[1];
rz(-0.99937159) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10103664) q[0];
sx q[0];
rz(-0.66185364) q[0];
sx q[0];
rz(0.7579114) q[0];
rz(3.0452764) q[2];
sx q[2];
rz(-1.673398) q[2];
sx q[2];
rz(1.4136537) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4564086) q[1];
sx q[1];
rz(-2.051609) q[1];
sx q[1];
rz(-0.64688869) q[1];
rz(-pi) q[2];
rz(1.9346775) q[3];
sx q[3];
rz(-0.69857222) q[3];
sx q[3];
rz(1.9649803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.31819764) q[2];
sx q[2];
rz(-0.91788936) q[2];
sx q[2];
rz(-1.2981752) q[2];
rz(-0.038481742) q[3];
sx q[3];
rz(-1.1063856) q[3];
sx q[3];
rz(-1.5255671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-1.0549952) q[0];
sx q[0];
rz(-2.0107858) q[0];
sx q[0];
rz(2.8259592) q[0];
rz(1.9075958) q[1];
sx q[1];
rz(-0.71593586) q[1];
sx q[1];
rz(1.2589781) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18470284) q[0];
sx q[0];
rz(-1.8543921) q[0];
sx q[0];
rz(-0.75329874) q[0];
rz(-pi) q[1];
rz(-0.76541111) q[2];
sx q[2];
rz(-1.4452626) q[2];
sx q[2];
rz(1.5620934) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.96871829) q[1];
sx q[1];
rz(-2.5596788) q[1];
sx q[1];
rz(-0.7919148) q[1];
x q[2];
rz(-3.0837667) q[3];
sx q[3];
rz(-0.10198051) q[3];
sx q[3];
rz(2.0174055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9180318) q[2];
sx q[2];
rz(-0.74350244) q[2];
sx q[2];
rz(-0.027916748) q[2];
rz(0.74808407) q[3];
sx q[3];
rz(-1.7223765) q[3];
sx q[3];
rz(2.9885651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41392031) q[0];
sx q[0];
rz(-0.99069178) q[0];
sx q[0];
rz(-0.45898166) q[0];
rz(1.1363632) q[1];
sx q[1];
rz(-1.7106979) q[1];
sx q[1];
rz(2.9232025) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9680351) q[0];
sx q[0];
rz(-1.5941248) q[0];
sx q[0];
rz(-1.5541881) q[0];
rz(-pi) q[1];
x q[1];
rz(1.893546) q[2];
sx q[2];
rz(-0.40229978) q[2];
sx q[2];
rz(0.93076555) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0291044) q[1];
sx q[1];
rz(-2.6262754) q[1];
sx q[1];
rz(-2.8397933) q[1];
x q[2];
rz(0.13111643) q[3];
sx q[3];
rz(-1.1857486) q[3];
sx q[3];
rz(-0.39310405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2961262) q[2];
sx q[2];
rz(-0.25456905) q[2];
sx q[2];
rz(1.0477585) q[2];
rz(0.74538499) q[3];
sx q[3];
rz(-1.5630009) q[3];
sx q[3];
rz(1.6489702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1202241) q[0];
sx q[0];
rz(-0.57149082) q[0];
sx q[0];
rz(-2.7987203) q[0];
rz(2.951237) q[1];
sx q[1];
rz(-1.0089259) q[1];
sx q[1];
rz(-2.6284133) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8489055) q[0];
sx q[0];
rz(-2.2063037) q[0];
sx q[0];
rz(0.19750316) q[0];
x q[1];
rz(0.35731213) q[2];
sx q[2];
rz(-1.0917959) q[2];
sx q[2];
rz(1.8425531) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78402218) q[1];
sx q[1];
rz(-0.62767941) q[1];
sx q[1];
rz(2.4178972) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8796092) q[3];
sx q[3];
rz(-0.43918375) q[3];
sx q[3];
rz(-2.9787946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.43789431) q[2];
sx q[2];
rz(-0.63592211) q[2];
sx q[2];
rz(2.5305914) q[2];
rz(2.8827187) q[3];
sx q[3];
rz(-1.9301819) q[3];
sx q[3];
rz(-1.7620311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3872869) q[0];
sx q[0];
rz(-1.5504693) q[0];
sx q[0];
rz(-1.5691527) q[0];
rz(-0.98987956) q[1];
sx q[1];
rz(-1.9342593) q[1];
sx q[1];
rz(0.63607279) q[1];
rz(-1.382147) q[2];
sx q[2];
rz(-1.5619642) q[2];
sx q[2];
rz(-1.2240338) q[2];
rz(-3.0627904) q[3];
sx q[3];
rz(-2.4410421) q[3];
sx q[3];
rz(-2.0564612) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
