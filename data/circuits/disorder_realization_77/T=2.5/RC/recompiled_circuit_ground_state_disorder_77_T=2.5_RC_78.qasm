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
rz(2.4426065) q[1];
sx q[1];
rz(12.01233) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0813538) q[0];
sx q[0];
rz(-1.0805655) q[0];
sx q[0];
rz(-2.1054563) q[0];
rz(2.0238766) q[2];
sx q[2];
rz(-2.6181698) q[2];
sx q[2];
rz(-0.40016178) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.89008625) q[1];
sx q[1];
rz(-0.68822574) q[1];
sx q[1];
rz(-1.623897) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0022121) q[3];
sx q[3];
rz(-0.7665638) q[3];
sx q[3];
rz(2.1994906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.379091) q[2];
sx q[2];
rz(-1.5662301) q[2];
sx q[2];
rz(-2.9762034) q[2];
rz(0.23073828) q[3];
sx q[3];
rz(-2.8985891) q[3];
sx q[3];
rz(0.86133426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6623401) q[0];
sx q[0];
rz(-0.32821822) q[0];
sx q[0];
rz(1.3586556) q[0];
rz(-1.6183629) q[1];
sx q[1];
rz(-2.3871469) q[1];
sx q[1];
rz(-0.84017909) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5206487) q[0];
sx q[0];
rz(-2.0387702) q[0];
sx q[0];
rz(0.2999938) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1112059) q[2];
sx q[2];
rz(-2.3860117) q[2];
sx q[2];
rz(2.7019175) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6908535) q[1];
sx q[1];
rz(-0.87635856) q[1];
sx q[1];
rz(-1.9574375) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4235557) q[3];
sx q[3];
rz(-0.81063945) q[3];
sx q[3];
rz(-0.37732201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.75300616) q[2];
sx q[2];
rz(-1.2331839) q[2];
sx q[2];
rz(0.93442717) q[2];
rz(1.8560483) q[3];
sx q[3];
rz(-0.779874) q[3];
sx q[3];
rz(-0.88750315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.0282106) q[0];
sx q[0];
rz(-2.1206355) q[0];
sx q[0];
rz(1.5895948) q[0];
rz(-1.3681083) q[1];
sx q[1];
rz(-1.2721456) q[1];
sx q[1];
rz(2.0210463) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6137973) q[0];
sx q[0];
rz(-0.58301914) q[0];
sx q[0];
rz(0.9436508) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1324544) q[2];
sx q[2];
rz(-1.8016644) q[2];
sx q[2];
rz(1.1924274) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6419598) q[1];
sx q[1];
rz(-1.4781535) q[1];
sx q[1];
rz(2.605162) q[1];
rz(-pi) q[2];
rz(3.0157959) q[3];
sx q[3];
rz(-1.9748944) q[3];
sx q[3];
rz(2.8640797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0997194) q[2];
sx q[2];
rz(-0.21328829) q[2];
sx q[2];
rz(-2.8495157) q[2];
rz(1.9000351) q[3];
sx q[3];
rz(-1.440666) q[3];
sx q[3];
rz(-0.45274538) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1555772) q[0];
sx q[0];
rz(-0.41322511) q[0];
sx q[0];
rz(2.7098932) q[0];
rz(-2.9875634) q[1];
sx q[1];
rz(-1.5715716) q[1];
sx q[1];
rz(2.0808751) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78846473) q[0];
sx q[0];
rz(-1.9110838) q[0];
sx q[0];
rz(-2.8022604) q[0];
x q[1];
rz(-0.94580146) q[2];
sx q[2];
rz(-0.29871395) q[2];
sx q[2];
rz(0.83160366) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1460675) q[1];
sx q[1];
rz(-1.3473099) q[1];
sx q[1];
rz(-1.9053915) q[1];
rz(-pi) q[2];
rz(0.68935518) q[3];
sx q[3];
rz(-0.24447799) q[3];
sx q[3];
rz(2.1357128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3920307) q[2];
sx q[2];
rz(-1.6101086) q[2];
sx q[2];
rz(-2.5677666) q[2];
rz(-0.057403322) q[3];
sx q[3];
rz(-1.4797689) q[3];
sx q[3];
rz(0.81006947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60657984) q[0];
sx q[0];
rz(-2.8856394) q[0];
sx q[0];
rz(0.25273299) q[0];
rz(2.9622954) q[1];
sx q[1];
rz(-1.3456234) q[1];
sx q[1];
rz(-1.4089233) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64302639) q[0];
sx q[0];
rz(-1.4966949) q[0];
sx q[0];
rz(1.8084433) q[0];
x q[1];
rz(-1.9581355) q[2];
sx q[2];
rz(-1.172386) q[2];
sx q[2];
rz(-0.18926316) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.21633575) q[1];
sx q[1];
rz(-2.3489423) q[1];
sx q[1];
rz(1.7429732) q[1];
x q[2];
rz(-1.9305655) q[3];
sx q[3];
rz(-1.0749987) q[3];
sx q[3];
rz(-1.825765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0130284) q[2];
sx q[2];
rz(-1.3821673) q[2];
sx q[2];
rz(1.9405091) q[2];
rz(-0.80738336) q[3];
sx q[3];
rz(-1.3599334) q[3];
sx q[3];
rz(2.0090296) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078868911) q[0];
sx q[0];
rz(-1.6428592) q[0];
sx q[0];
rz(0.8412745) q[0];
rz(1.9367283) q[1];
sx q[1];
rz(-1.8236022) q[1];
sx q[1];
rz(1.0296317) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1372927) q[0];
sx q[0];
rz(-0.51057112) q[0];
sx q[0];
rz(-1.5493896) q[0];
x q[1];
rz(-2.6961521) q[2];
sx q[2];
rz(-0.56156073) q[2];
sx q[2];
rz(-2.6391921) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8334044) q[1];
sx q[1];
rz(-0.17886111) q[1];
sx q[1];
rz(-0.043828242) q[1];
rz(-pi) q[2];
x q[2];
rz(1.383841) q[3];
sx q[3];
rz(-1.2724981) q[3];
sx q[3];
rz(1.3307856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.010217696) q[2];
sx q[2];
rz(-0.75675941) q[2];
sx q[2];
rz(2.6436464) q[2];
rz(2.61854) q[3];
sx q[3];
rz(-1.7224576) q[3];
sx q[3];
rz(0.12652346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.842857) q[0];
sx q[0];
rz(-3.0099478) q[0];
sx q[0];
rz(0.81197062) q[0];
rz(0.89093351) q[1];
sx q[1];
rz(-1.9782601) q[1];
sx q[1];
rz(-2.1422211) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82841372) q[0];
sx q[0];
rz(-2.0069608) q[0];
sx q[0];
rz(2.6266813) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4677202) q[2];
sx q[2];
rz(-1.4749881) q[2];
sx q[2];
rz(0.16703781) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2216894) q[1];
sx q[1];
rz(-1.0070325) q[1];
sx q[1];
rz(2.1497822) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2362367) q[3];
sx q[3];
rz(-1.3398583) q[3];
sx q[3];
rz(-3.0311716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.31819764) q[2];
sx q[2];
rz(-2.2237033) q[2];
sx q[2];
rz(-1.8434175) q[2];
rz(-0.038481742) q[3];
sx q[3];
rz(-1.1063856) q[3];
sx q[3];
rz(1.6160256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0549952) q[0];
sx q[0];
rz(-1.1308068) q[0];
sx q[0];
rz(-0.31563345) q[0];
rz(1.9075958) q[1];
sx q[1];
rz(-2.4256568) q[1];
sx q[1];
rz(-1.2589781) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0961279) q[0];
sx q[0];
rz(-2.3466131) q[0];
sx q[0];
rz(0.40277512) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76541111) q[2];
sx q[2];
rz(-1.4452626) q[2];
sx q[2];
rz(1.5620934) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.96871829) q[1];
sx q[1];
rz(-2.5596788) q[1];
sx q[1];
rz(2.3496778) q[1];
rz(-pi) q[2];
rz(1.564882) q[3];
sx q[3];
rz(-1.4689869) q[3];
sx q[3];
rz(1.9592782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.2235609) q[2];
sx q[2];
rz(-2.3980902) q[2];
sx q[2];
rz(-0.027916748) q[2];
rz(0.74808407) q[3];
sx q[3];
rz(-1.7223765) q[3];
sx q[3];
rz(-0.15302756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.1509009) q[0];
sx q[0];
rz(0.45898166) q[0];
rz(1.1363632) q[1];
sx q[1];
rz(-1.4308948) q[1];
sx q[1];
rz(-2.9232025) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9680351) q[0];
sx q[0];
rz(-1.5941248) q[0];
sx q[0];
rz(-1.5541881) q[0];
rz(0.13414975) q[2];
sx q[2];
rz(-1.9512259) q[2];
sx q[2];
rz(-2.5594437) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7228411) q[1];
sx q[1];
rz(-1.7178078) q[1];
sx q[1];
rz(0.49572368) q[1];
x q[2];
rz(1.9588542) q[3];
sx q[3];
rz(-1.6922631) q[3];
sx q[3];
rz(-2.0133919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8454664) q[2];
sx q[2];
rz(-2.8870236) q[2];
sx q[2];
rz(2.0938342) q[2];
rz(0.74538499) q[3];
sx q[3];
rz(-1.5785917) q[3];
sx q[3];
rz(-1.6489702) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1202241) q[0];
sx q[0];
rz(-0.57149082) q[0];
sx q[0];
rz(-0.34287232) q[0];
rz(-2.951237) q[1];
sx q[1];
rz(-1.0089259) q[1];
sx q[1];
rz(2.6284133) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8489055) q[0];
sx q[0];
rz(-2.2063037) q[0];
sx q[0];
rz(0.19750316) q[0];
rz(-pi) q[1];
rz(2.076976) q[2];
sx q[2];
rz(-1.2551885) q[2];
sx q[2];
rz(-0.10135717) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0963971) q[1];
sx q[1];
rz(-2.0264813) q[1];
sx q[1];
rz(-1.1229188) q[1];
rz(-pi) q[2];
rz(-1.8796092) q[3];
sx q[3];
rz(-2.7024089) q[3];
sx q[3];
rz(-0.16279804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7036983) q[2];
sx q[2];
rz(-2.5056705) q[2];
sx q[2];
rz(2.5305914) q[2];
rz(0.25887394) q[3];
sx q[3];
rz(-1.9301819) q[3];
sx q[3];
rz(-1.3795616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7543058) q[0];
sx q[0];
rz(-1.5504693) q[0];
sx q[0];
rz(-1.5691527) q[0];
rz(-0.98987956) q[1];
sx q[1];
rz(-1.9342593) q[1];
sx q[1];
rz(0.63607279) q[1];
rz(-1.7594457) q[2];
sx q[2];
rz(-1.5796285) q[2];
sx q[2];
rz(1.9175588) q[2];
rz(0.078802303) q[3];
sx q[3];
rz(-2.4410421) q[3];
sx q[3];
rz(-2.0564612) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
