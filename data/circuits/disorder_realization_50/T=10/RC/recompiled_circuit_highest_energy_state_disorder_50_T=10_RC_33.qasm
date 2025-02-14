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
rz(0.4515689) q[0];
sx q[0];
rz(1.4181925) q[0];
sx q[0];
rz(12.137492) q[0];
rz(-0.14313993) q[1];
sx q[1];
rz(-2.0506471) q[1];
sx q[1];
rz(-3.0615038) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6038136) q[0];
sx q[0];
rz(-1.5575557) q[0];
sx q[0];
rz(0.079188345) q[0];
rz(-1.8378788) q[2];
sx q[2];
rz(-2.7686229) q[2];
sx q[2];
rz(-2.3139062) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.45177191) q[1];
sx q[1];
rz(-2.2597888) q[1];
sx q[1];
rz(-0.34386872) q[1];
rz(-pi) q[2];
x q[2];
rz(2.56404) q[3];
sx q[3];
rz(-1.7087382) q[3];
sx q[3];
rz(2.4393314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.074778883) q[2];
sx q[2];
rz(-1.0513693) q[2];
sx q[2];
rz(-1.4220413) q[2];
rz(-1.7285796) q[3];
sx q[3];
rz(-1.541411) q[3];
sx q[3];
rz(-1.6933256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0195352) q[0];
sx q[0];
rz(-1.3966565) q[0];
sx q[0];
rz(2.8080217) q[0];
rz(0.16593274) q[1];
sx q[1];
rz(-2.797778) q[1];
sx q[1];
rz(0.75210345) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2488326) q[0];
sx q[0];
rz(-1.1247096) q[0];
sx q[0];
rz(-0.89527881) q[0];
rz(-pi) q[1];
rz(1.5971489) q[2];
sx q[2];
rz(-1.7894723) q[2];
sx q[2];
rz(1.6353351) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1811407) q[1];
sx q[1];
rz(-1.4616832) q[1];
sx q[1];
rz(1.0427598) q[1];
rz(-pi) q[2];
rz(-2.4346967) q[3];
sx q[3];
rz(-0.33875468) q[3];
sx q[3];
rz(-1.930463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.077945262) q[2];
sx q[2];
rz(-0.91219488) q[2];
sx q[2];
rz(2.4845541) q[2];
rz(0.71197236) q[3];
sx q[3];
rz(-0.65198055) q[3];
sx q[3];
rz(2.3967192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4350568) q[0];
sx q[0];
rz(-0.61921316) q[0];
sx q[0];
rz(-1.0414498) q[0];
rz(2.7667747) q[1];
sx q[1];
rz(-1.5033787) q[1];
sx q[1];
rz(2.5846438) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.108088) q[0];
sx q[0];
rz(-1.9226941) q[0];
sx q[0];
rz(-0.75893388) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1191423) q[2];
sx q[2];
rz(-1.6385534) q[2];
sx q[2];
rz(0.51167831) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5365689) q[1];
sx q[1];
rz(-0.67702952) q[1];
sx q[1];
rz(-0.1635464) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8076423) q[3];
sx q[3];
rz(-0.93333731) q[3];
sx q[3];
rz(-1.6528296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.53691429) q[2];
sx q[2];
rz(-2.2173209) q[2];
sx q[2];
rz(-2.7064145) q[2];
rz(-2.6168881) q[3];
sx q[3];
rz(-2.8080431) q[3];
sx q[3];
rz(3.0068126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0404469) q[0];
sx q[0];
rz(-2.0772159) q[0];
sx q[0];
rz(-1.6271628) q[0];
rz(-2.0928404) q[1];
sx q[1];
rz(-0.92275134) q[1];
sx q[1];
rz(1.4357766) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1483106) q[0];
sx q[0];
rz(-1.7459501) q[0];
sx q[0];
rz(1.6143285) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5982791) q[2];
sx q[2];
rz(-1.8738973) q[2];
sx q[2];
rz(1.4599279) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.32325778) q[1];
sx q[1];
rz(-2.0916953) q[1];
sx q[1];
rz(3.0677615) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1952329) q[3];
sx q[3];
rz(-2.588996) q[3];
sx q[3];
rz(-0.70878562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.40020308) q[2];
sx q[2];
rz(-1.3578537) q[2];
sx q[2];
rz(3.1164361) q[2];
rz(-1.4870421) q[3];
sx q[3];
rz(-2.7436723) q[3];
sx q[3];
rz(0.81620836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0461034) q[0];
sx q[0];
rz(-0.64952055) q[0];
sx q[0];
rz(-2.982614) q[0];
rz(-1.1051296) q[1];
sx q[1];
rz(-2.4159894) q[1];
sx q[1];
rz(2.9172858) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8858356) q[0];
sx q[0];
rz(-2.1191264) q[0];
sx q[0];
rz(-1.095039) q[0];
x q[1];
rz(2.131046) q[2];
sx q[2];
rz(-1.5429826) q[2];
sx q[2];
rz(-1.5459488) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6857924) q[1];
sx q[1];
rz(-0.53736254) q[1];
sx q[1];
rz(1.3440029) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57711011) q[3];
sx q[3];
rz(-1.464586) q[3];
sx q[3];
rz(2.4662227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9362681) q[2];
sx q[2];
rz(-1.3621796) q[2];
sx q[2];
rz(0.74860191) q[2];
rz(-0.11048206) q[3];
sx q[3];
rz(-1.2275077) q[3];
sx q[3];
rz(-2.7301679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1614138) q[0];
sx q[0];
rz(-2.7730589) q[0];
sx q[0];
rz(-2.4112716) q[0];
rz(-1.5018564) q[1];
sx q[1];
rz(-1.390099) q[1];
sx q[1];
rz(-0.23434848) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4880331) q[0];
sx q[0];
rz(-0.83439487) q[0];
sx q[0];
rz(-2.5279765) q[0];
x q[1];
rz(-2.3086771) q[2];
sx q[2];
rz(-2.3489522) q[2];
sx q[2];
rz(-1.329525) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1533656) q[1];
sx q[1];
rz(-1.7333789) q[1];
sx q[1];
rz(1.52262) q[1];
x q[2];
rz(-1.9803011) q[3];
sx q[3];
rz(-2.9886768) q[3];
sx q[3];
rz(0.59329882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.19330567) q[2];
sx q[2];
rz(-0.81393465) q[2];
sx q[2];
rz(-1.8602271) q[2];
rz(0.60503259) q[3];
sx q[3];
rz(-2.0421959) q[3];
sx q[3];
rz(1.9090778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.989711) q[0];
sx q[0];
rz(-1.879377) q[0];
sx q[0];
rz(-0.61304098) q[0];
rz(2.4609861) q[1];
sx q[1];
rz(-2.0304408) q[1];
sx q[1];
rz(1.8162762) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26215689) q[0];
sx q[0];
rz(-0.61816321) q[0];
sx q[0];
rz(2.6690041) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21317404) q[2];
sx q[2];
rz(-1.6521378) q[2];
sx q[2];
rz(-2.6590462) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5449574) q[1];
sx q[1];
rz(-1.4442634) q[1];
sx q[1];
rz(-2.5759349) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58206733) q[3];
sx q[3];
rz(-1.9110693) q[3];
sx q[3];
rz(0.64990625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.94720542) q[2];
sx q[2];
rz(-1.9401865) q[2];
sx q[2];
rz(-0.80805937) q[2];
rz(-2.4957073) q[3];
sx q[3];
rz(-2.8736727) q[3];
sx q[3];
rz(0.3033692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5658257) q[0];
sx q[0];
rz(-2.4149826) q[0];
sx q[0];
rz(-0.44276825) q[0];
rz(0.59204656) q[1];
sx q[1];
rz(-2.4202085) q[1];
sx q[1];
rz(0.17568849) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9239227) q[0];
sx q[0];
rz(-2.6852073) q[0];
sx q[0];
rz(0.66500591) q[0];
rz(-0.08205361) q[2];
sx q[2];
rz(-0.4528762) q[2];
sx q[2];
rz(2.1575515) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2776043) q[1];
sx q[1];
rz(-0.93011198) q[1];
sx q[1];
rz(-1.6246269) q[1];
rz(-0.91276987) q[3];
sx q[3];
rz(-0.72421779) q[3];
sx q[3];
rz(2.4778544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5438133) q[2];
sx q[2];
rz(-2.3976349) q[2];
sx q[2];
rz(0.90710863) q[2];
rz(-2.234327) q[3];
sx q[3];
rz(-1.7472569) q[3];
sx q[3];
rz(-3.0239014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9573145) q[0];
sx q[0];
rz(-2.2397581) q[0];
sx q[0];
rz(2.9515475) q[0];
rz(-1.5589145) q[1];
sx q[1];
rz(-1.546944) q[1];
sx q[1];
rz(1.83439) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35018539) q[0];
sx q[0];
rz(-2.7481347) q[0];
sx q[0];
rz(1.1883452) q[0];
rz(-pi) q[1];
rz(0.3385915) q[2];
sx q[2];
rz(-0.89950746) q[2];
sx q[2];
rz(3.1289575) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1550473) q[1];
sx q[1];
rz(-1.6035756) q[1];
sx q[1];
rz(1.4630586) q[1];
rz(-2.8966859) q[3];
sx q[3];
rz(-0.4462113) q[3];
sx q[3];
rz(3.0065766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0731395) q[2];
sx q[2];
rz(-2.0544724) q[2];
sx q[2];
rz(0.20115176) q[2];
rz(2.1189832) q[3];
sx q[3];
rz(-0.68251959) q[3];
sx q[3];
rz(1.539591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96711838) q[0];
sx q[0];
rz(-1.0459463) q[0];
sx q[0];
rz(-2.2651267) q[0];
rz(-2.5190952) q[1];
sx q[1];
rz(-0.21177706) q[1];
sx q[1];
rz(-0.34271398) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22607813) q[0];
sx q[0];
rz(-0.63185872) q[0];
sx q[0];
rz(1.4591239) q[0];
x q[1];
rz(-1.2340401) q[2];
sx q[2];
rz(-1.4554664) q[2];
sx q[2];
rz(-2.8591836) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5805768) q[1];
sx q[1];
rz(-2.3388712) q[1];
sx q[1];
rz(-0.64508857) q[1];
x q[2];
rz(1.8233772) q[3];
sx q[3];
rz(-1.1942721) q[3];
sx q[3];
rz(-1.2600217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9524625) q[2];
sx q[2];
rz(-1.6996982) q[2];
sx q[2];
rz(-2.6978037) q[2];
rz(-1.1187547) q[3];
sx q[3];
rz(-1.5464562) q[3];
sx q[3];
rz(-0.55383033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8583869) q[0];
sx q[0];
rz(-1.5571742) q[0];
sx q[0];
rz(-1.5207186) q[0];
rz(-3.0781147) q[1];
sx q[1];
rz(-1.3451481) q[1];
sx q[1];
rz(1.4048911) q[1];
rz(-1.7495161) q[2];
sx q[2];
rz(-2.8061493) q[2];
sx q[2];
rz(0.2827927) q[2];
rz(-0.24334454) q[3];
sx q[3];
rz(-1.018033) q[3];
sx q[3];
rz(1.3924934) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
