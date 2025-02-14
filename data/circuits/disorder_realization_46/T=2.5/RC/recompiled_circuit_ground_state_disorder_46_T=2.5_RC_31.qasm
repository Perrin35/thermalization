OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.0134861) q[0];
sx q[0];
rz(-0.82733893) q[0];
sx q[0];
rz(1.7537355) q[0];
rz(-2.2881621) q[1];
sx q[1];
rz(-2.7397459) q[1];
sx q[1];
rz(-2.0274577) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6809876) q[0];
sx q[0];
rz(-1.3189684) q[0];
sx q[0];
rz(-0.84427174) q[0];
x q[1];
rz(-1.763249) q[2];
sx q[2];
rz(-2.2666605) q[2];
sx q[2];
rz(0.025283289) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.6889894) q[1];
sx q[1];
rz(-1.5509948) q[1];
sx q[1];
rz(-0.0048586998) q[1];
x q[2];
rz(-0.57350343) q[3];
sx q[3];
rz(-0.57799229) q[3];
sx q[3];
rz(1.0740394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.16093291) q[2];
sx q[2];
rz(-0.42155835) q[2];
sx q[2];
rz(-1.4483615) q[2];
rz(0.24250044) q[3];
sx q[3];
rz(-2.226053) q[3];
sx q[3];
rz(-1.2601132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.563797) q[0];
sx q[0];
rz(-1.8255434) q[0];
sx q[0];
rz(-2.1412361) q[0];
rz(0.73161221) q[1];
sx q[1];
rz(-1.7121406) q[1];
sx q[1];
rz(-1.1801205) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7655752) q[0];
sx q[0];
rz(-1.4808324) q[0];
sx q[0];
rz(0.1040118) q[0];
x q[1];
rz(-1.8487471) q[2];
sx q[2];
rz(-1.6033844) q[2];
sx q[2];
rz(-1.1162835) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47467642) q[1];
sx q[1];
rz(-1.4366526) q[1];
sx q[1];
rz(-3.139702) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3445156) q[3];
sx q[3];
rz(-2.7555572) q[3];
sx q[3];
rz(0.71938709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0483027) q[2];
sx q[2];
rz(-0.90128171) q[2];
sx q[2];
rz(2.062659) q[2];
rz(-3.0537187) q[3];
sx q[3];
rz(-1.352997) q[3];
sx q[3];
rz(2.8972076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3946149) q[0];
sx q[0];
rz(-1.1629539) q[0];
sx q[0];
rz(1.278247) q[0];
rz(0.46319115) q[1];
sx q[1];
rz(-1.3563503) q[1];
sx q[1];
rz(-0.82122222) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44959497) q[0];
sx q[0];
rz(-2.2173081) q[0];
sx q[0];
rz(2.9086935) q[0];
rz(-1.6326007) q[2];
sx q[2];
rz(-1.9414177) q[2];
sx q[2];
rz(-0.48423094) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.48135172) q[1];
sx q[1];
rz(-1.7180182) q[1];
sx q[1];
rz(-2.6805356) q[1];
rz(-pi) q[2];
rz(2.2021095) q[3];
sx q[3];
rz(-2.2243613) q[3];
sx q[3];
rz(-0.77142772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.897573) q[2];
sx q[2];
rz(-1.579318) q[2];
sx q[2];
rz(2.139034) q[2];
rz(1.2282486) q[3];
sx q[3];
rz(-1.7398261) q[3];
sx q[3];
rz(-1.7206515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2380075) q[0];
sx q[0];
rz(-1.0500195) q[0];
sx q[0];
rz(2.8557657) q[0];
rz(-0.26126513) q[1];
sx q[1];
rz(-1.8087872) q[1];
sx q[1];
rz(3.1108943) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0945902) q[0];
sx q[0];
rz(-1.6436716) q[0];
sx q[0];
rz(0.26396172) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4723563) q[2];
sx q[2];
rz(-0.67924309) q[2];
sx q[2];
rz(-1.6427276) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.49467898) q[1];
sx q[1];
rz(-2.621197) q[1];
sx q[1];
rz(-1.0444276) q[1];
rz(-pi) q[2];
rz(3.0604393) q[3];
sx q[3];
rz(-1.2581244) q[3];
sx q[3];
rz(-2.816538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3280481) q[2];
sx q[2];
rz(-1.3406465) q[2];
sx q[2];
rz(-3.1385341) q[2];
rz(0.84117738) q[3];
sx q[3];
rz(-0.58924651) q[3];
sx q[3];
rz(0.35471788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4641814) q[0];
sx q[0];
rz(-2.2941636) q[0];
sx q[0];
rz(-1.4671951) q[0];
rz(-1.5122308) q[1];
sx q[1];
rz(-2.1758175) q[1];
sx q[1];
rz(2.1893952) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.061082) q[0];
sx q[0];
rz(-3.0709776) q[0];
sx q[0];
rz(0.10297601) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69779863) q[2];
sx q[2];
rz(-1.875617) q[2];
sx q[2];
rz(1.6483726) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.72022382) q[1];
sx q[1];
rz(-1.4671869) q[1];
sx q[1];
rz(-3.1006579) q[1];
x q[2];
rz(-2.5467993) q[3];
sx q[3];
rz(-0.99579558) q[3];
sx q[3];
rz(1.9146384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3516922) q[2];
sx q[2];
rz(-1.5396427) q[2];
sx q[2];
rz(2.6969625) q[2];
rz(-0.45186684) q[3];
sx q[3];
rz(-2.5022652) q[3];
sx q[3];
rz(0.062084559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(0.43279466) q[0];
sx q[0];
rz(-0.72467703) q[0];
sx q[0];
rz(0.92025796) q[0];
rz(-0.96011773) q[1];
sx q[1];
rz(-1.7476387) q[1];
sx q[1];
rz(-0.49930176) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1007784) q[0];
sx q[0];
rz(-1.4199323) q[0];
sx q[0];
rz(-2.8013381) q[0];
rz(0.22412207) q[2];
sx q[2];
rz(-2.4680228) q[2];
sx q[2];
rz(0.78896713) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0671406) q[1];
sx q[1];
rz(-0.66508355) q[1];
sx q[1];
rz(-2.6159899) q[1];
rz(-0.74120993) q[3];
sx q[3];
rz(-1.831372) q[3];
sx q[3];
rz(-1.2737361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.7509191) q[2];
sx q[2];
rz(-1.5385188) q[2];
sx q[2];
rz(-2.1092559) q[2];
rz(0.8647626) q[3];
sx q[3];
rz(-1.4356177) q[3];
sx q[3];
rz(0.56149703) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.432935) q[0];
sx q[0];
rz(-1.189804) q[0];
sx q[0];
rz(1.7565961) q[0];
rz(-2.2191091) q[1];
sx q[1];
rz(-1.9360767) q[1];
sx q[1];
rz(0.14499697) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7940878) q[0];
sx q[0];
rz(-1.4055168) q[0];
sx q[0];
rz(-1.3896304) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.029399662) q[2];
sx q[2];
rz(-0.85762944) q[2];
sx q[2];
rz(-2.5273163) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8783826) q[1];
sx q[1];
rz(-0.07941281) q[1];
sx q[1];
rz(-1.9694049) q[1];
rz(2.3649366) q[3];
sx q[3];
rz(-1.1478416) q[3];
sx q[3];
rz(-0.91931988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.21236913) q[2];
sx q[2];
rz(-2.9085458) q[2];
sx q[2];
rz(-1.93335) q[2];
rz(0.82331795) q[3];
sx q[3];
rz(-1.9602785) q[3];
sx q[3];
rz(-1.7821144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(3.0804444) q[0];
sx q[0];
rz(-0.7875945) q[0];
sx q[0];
rz(1.0795235) q[0];
rz(0.98006717) q[1];
sx q[1];
rz(-2.1213687) q[1];
sx q[1];
rz(-0.10990873) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85724505) q[0];
sx q[0];
rz(-1.2099414) q[0];
sx q[0];
rz(-1.2566823) q[0];
rz(-pi) q[1];
rz(-2.8632322) q[2];
sx q[2];
rz(-1.5621857) q[2];
sx q[2];
rz(1.6842764) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.67092451) q[1];
sx q[1];
rz(-1.9940071) q[1];
sx q[1];
rz(-1.898134) q[1];
x q[2];
rz(2.917116) q[3];
sx q[3];
rz(-1.0521476) q[3];
sx q[3];
rz(-2.8751837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.93426934) q[2];
sx q[2];
rz(-2.0685652) q[2];
sx q[2];
rz(0.70460021) q[2];
rz(-0.28471026) q[3];
sx q[3];
rz(-2.4094818) q[3];
sx q[3];
rz(-2.708191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1339742) q[0];
sx q[0];
rz(-1.3173137) q[0];
sx q[0];
rz(2.2229069) q[0];
rz(-0.48565117) q[1];
sx q[1];
rz(-0.44356569) q[1];
sx q[1];
rz(-2.0408911) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0328631) q[0];
sx q[0];
rz(-0.19685611) q[0];
sx q[0];
rz(-2.3545676) q[0];
x q[1];
rz(-2.9100239) q[2];
sx q[2];
rz(-1.2763888) q[2];
sx q[2];
rz(0.60588479) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6599258) q[1];
sx q[1];
rz(-1.1705913) q[1];
sx q[1];
rz(0.11379644) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0555567) q[3];
sx q[3];
rz(-2.5042296) q[3];
sx q[3];
rz(-1.5051382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47621581) q[2];
sx q[2];
rz(-2.2874338) q[2];
sx q[2];
rz(0.7594792) q[2];
rz(-0.94799834) q[3];
sx q[3];
rz(-1.6839323) q[3];
sx q[3];
rz(-1.6215526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.043902472) q[0];
sx q[0];
rz(-0.49711415) q[0];
sx q[0];
rz(2.8105766) q[0];
rz(-0.76200062) q[1];
sx q[1];
rz(-0.29251978) q[1];
sx q[1];
rz(2.5133572) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8549919) q[0];
sx q[0];
rz(-2.1341265) q[0];
sx q[0];
rz(-1.9024293) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30769257) q[2];
sx q[2];
rz(-1.3465704) q[2];
sx q[2];
rz(-0.7219519) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3207631) q[1];
sx q[1];
rz(-2.4270128) q[1];
sx q[1];
rz(-2.5933867) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2248433) q[3];
sx q[3];
rz(-1.1802434) q[3];
sx q[3];
rz(-1.2929163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.021775333) q[2];
sx q[2];
rz(-0.47406083) q[2];
sx q[2];
rz(0.21381703) q[2];
rz(0.98624936) q[3];
sx q[3];
rz(-1.8503559) q[3];
sx q[3];
rz(-3.0333062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66961359) q[0];
sx q[0];
rz(-1.6983953) q[0];
sx q[0];
rz(-1.4629913) q[0];
rz(-1.9646473) q[1];
sx q[1];
rz(-1.6564449) q[1];
sx q[1];
rz(0.0080531837) q[1];
rz(2.9672553) q[2];
sx q[2];
rz(-1.871289) q[2];
sx q[2];
rz(1.540779) q[2];
rz(-1.4900653) q[3];
sx q[3];
rz(-0.92083081) q[3];
sx q[3];
rz(3.010078) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
