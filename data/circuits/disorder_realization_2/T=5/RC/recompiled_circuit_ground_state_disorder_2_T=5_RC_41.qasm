OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6096697) q[0];
sx q[0];
rz(-1.6165531) q[0];
sx q[0];
rz(-1.2582231) q[0];
rz(1.5594856) q[1];
sx q[1];
rz(-0.48696247) q[1];
sx q[1];
rz(0.43509405) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5814911) q[0];
sx q[0];
rz(-1.035443) q[0];
sx q[0];
rz(-2.125432) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2266394) q[2];
sx q[2];
rz(-1.0593057) q[2];
sx q[2];
rz(-2.34735) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9261998) q[1];
sx q[1];
rz(-1.493931) q[1];
sx q[1];
rz(1.8410626) q[1];
rz(2.052911) q[3];
sx q[3];
rz(-2.0687244) q[3];
sx q[3];
rz(2.9256647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2121928) q[2];
sx q[2];
rz(-1.2779028) q[2];
sx q[2];
rz(1.0431935) q[2];
rz(-2.7405401) q[3];
sx q[3];
rz(-1.4570823) q[3];
sx q[3];
rz(-2.5636087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9143518) q[0];
sx q[0];
rz(-2.9965017) q[0];
sx q[0];
rz(2.5703854) q[0];
rz(-2.1696443) q[1];
sx q[1];
rz(-0.74058878) q[1];
sx q[1];
rz(1.1245701) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1069431) q[0];
sx q[0];
rz(-2.1799583) q[0];
sx q[0];
rz(-2.8657495) q[0];
rz(2.6827181) q[2];
sx q[2];
rz(-0.70415184) q[2];
sx q[2];
rz(1.6353232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3414866) q[1];
sx q[1];
rz(-0.84034318) q[1];
sx q[1];
rz(3.0154702) q[1];
rz(-pi) q[2];
rz(1.589682) q[3];
sx q[3];
rz(-1.369993) q[3];
sx q[3];
rz(0.41057472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9198415) q[2];
sx q[2];
rz(-1.5849179) q[2];
sx q[2];
rz(2.705503) q[2];
rz(0.74887216) q[3];
sx q[3];
rz(-2.336453) q[3];
sx q[3];
rz(2.3124636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1083168) q[0];
sx q[0];
rz(-1.3949787) q[0];
sx q[0];
rz(-3.0535611) q[0];
rz(-2.2875359) q[1];
sx q[1];
rz(-0.66103926) q[1];
sx q[1];
rz(2.0565654) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4861761) q[0];
sx q[0];
rz(-1.2987505) q[0];
sx q[0];
rz(0.06975091) q[0];
x q[1];
rz(0.95698551) q[2];
sx q[2];
rz(-0.26273604) q[2];
sx q[2];
rz(1.764591) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26478085) q[1];
sx q[1];
rz(-1.1327883) q[1];
sx q[1];
rz(2.9609738) q[1];
rz(0.17739968) q[3];
sx q[3];
rz(-0.39493233) q[3];
sx q[3];
rz(-1.653914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6701086) q[2];
sx q[2];
rz(-1.7210311) q[2];
sx q[2];
rz(0.74618435) q[2];
rz(0.21607312) q[3];
sx q[3];
rz(-1.2553299) q[3];
sx q[3];
rz(-0.20572534) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3716607) q[0];
sx q[0];
rz(-3.0682204) q[0];
sx q[0];
rz(1.368847) q[0];
rz(-0.61028496) q[1];
sx q[1];
rz(-1.7758324) q[1];
sx q[1];
rz(2.761421) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26977045) q[0];
sx q[0];
rz(-2.5996723) q[0];
sx q[0];
rz(-0.63585441) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39003464) q[2];
sx q[2];
rz(-2.6594498) q[2];
sx q[2];
rz(-0.48393656) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.11503765) q[1];
sx q[1];
rz(-1.4241744) q[1];
sx q[1];
rz(-2.0762334) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1690456) q[3];
sx q[3];
rz(-1.6309079) q[3];
sx q[3];
rz(-1.7463804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2569106) q[2];
sx q[2];
rz(-2.1941049) q[2];
sx q[2];
rz(2.103629) q[2];
rz(0.3642309) q[3];
sx q[3];
rz(-2.3528152) q[3];
sx q[3];
rz(-0.68575931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45201856) q[0];
sx q[0];
rz(-1.4301393) q[0];
sx q[0];
rz(2.3110287) q[0];
rz(3.0914302) q[1];
sx q[1];
rz(-0.12718931) q[1];
sx q[1];
rz(-0.62279472) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4573097) q[0];
sx q[0];
rz(-1.3535392) q[0];
sx q[0];
rz(-1.6506249) q[0];
rz(-pi) q[1];
rz(0.41968996) q[2];
sx q[2];
rz(-2.6869535) q[2];
sx q[2];
rz(-1.1722133) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2166136) q[1];
sx q[1];
rz(-2.2053524) q[1];
sx q[1];
rz(0.44559176) q[1];
rz(-pi) q[2];
rz(-3.1386802) q[3];
sx q[3];
rz(-0.022449819) q[3];
sx q[3];
rz(-2.8685304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.043776) q[2];
sx q[2];
rz(-2.706683) q[2];
sx q[2];
rz(-2.3338976) q[2];
rz(-1.5378599) q[3];
sx q[3];
rz(-1.5065498) q[3];
sx q[3];
rz(-0.22411552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9648801) q[0];
sx q[0];
rz(-1.3310615) q[0];
sx q[0];
rz(-1.6656732) q[0];
rz(1.2337947) q[1];
sx q[1];
rz(-1.7314311) q[1];
sx q[1];
rz(2.9880611) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79685211) q[0];
sx q[0];
rz(-1.5670876) q[0];
sx q[0];
rz(3.0319503) q[0];
x q[1];
rz(-2.0664882) q[2];
sx q[2];
rz(-1.9736991) q[2];
sx q[2];
rz(0.35225454) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8922032) q[1];
sx q[1];
rz(-1.0703328) q[1];
sx q[1];
rz(-1.0618147) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0136387) q[3];
sx q[3];
rz(-2.4989486) q[3];
sx q[3];
rz(1.2817739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4393356) q[2];
sx q[2];
rz(-0.80436891) q[2];
sx q[2];
rz(0.55629936) q[2];
rz(0.016911658) q[3];
sx q[3];
rz(-2.1664679) q[3];
sx q[3];
rz(-2.4771966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21338129) q[0];
sx q[0];
rz(-3.0688372) q[0];
sx q[0];
rz(-0.67759204) q[0];
rz(1.821359) q[1];
sx q[1];
rz(-1.0878539) q[1];
sx q[1];
rz(2.5755612) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6460577) q[0];
sx q[0];
rz(-1.2355775) q[0];
sx q[0];
rz(2.6779102) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8458869) q[2];
sx q[2];
rz(-0.054271532) q[2];
sx q[2];
rz(-1.4473297) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.54864751) q[1];
sx q[1];
rz(-0.8742395) q[1];
sx q[1];
rz(-0.19754438) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5291847) q[3];
sx q[3];
rz(-1.0806314) q[3];
sx q[3];
rz(-3.0050638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.90938202) q[2];
sx q[2];
rz(-0.38572329) q[2];
sx q[2];
rz(0.56987008) q[2];
rz(-1.0424403) q[3];
sx q[3];
rz(-1.180155) q[3];
sx q[3];
rz(1.0038143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7637699) q[0];
sx q[0];
rz(-2.7144987) q[0];
sx q[0];
rz(-1.0231934) q[0];
rz(-2.2238253) q[1];
sx q[1];
rz(-2.387391) q[1];
sx q[1];
rz(-2.2775211) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8811778) q[0];
sx q[0];
rz(-2.8679608) q[0];
sx q[0];
rz(2.0405962) q[0];
x q[1];
rz(0.90166868) q[2];
sx q[2];
rz(-2.5732917) q[2];
sx q[2];
rz(1.2327884) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1601711) q[1];
sx q[1];
rz(-2.2653901) q[1];
sx q[1];
rz(-2.2194374) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72743639) q[3];
sx q[3];
rz(-1.8525278) q[3];
sx q[3];
rz(1.4294595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0920022) q[2];
sx q[2];
rz(-1.2668173) q[2];
sx q[2];
rz(-2.3769296) q[2];
rz(2.2406254) q[3];
sx q[3];
rz(-0.67470208) q[3];
sx q[3];
rz(-2.0425792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6242999) q[0];
sx q[0];
rz(-0.39532548) q[0];
sx q[0];
rz(-2.1584391) q[0];
rz(0.82955018) q[1];
sx q[1];
rz(-1.364565) q[1];
sx q[1];
rz(-0.57873908) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0522636) q[0];
sx q[0];
rz(-0.39396414) q[0];
sx q[0];
rz(-1.6039421) q[0];
rz(-pi) q[1];
x q[1];
rz(0.54429599) q[2];
sx q[2];
rz(-1.0981264) q[2];
sx q[2];
rz(-0.67677697) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9959401) q[1];
sx q[1];
rz(-2.4697587) q[1];
sx q[1];
rz(-2.2493659) q[1];
x q[2];
rz(1.5341358) q[3];
sx q[3];
rz(-1.9158773) q[3];
sx q[3];
rz(-2.7861129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4911554) q[2];
sx q[2];
rz(-1.8530242) q[2];
sx q[2];
rz(-0.11162652) q[2];
rz(-2.1857183) q[3];
sx q[3];
rz(-0.99146944) q[3];
sx q[3];
rz(1.2607695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2535506) q[0];
sx q[0];
rz(-2.084806) q[0];
sx q[0];
rz(-3.1274617) q[0];
rz(1.1125394) q[1];
sx q[1];
rz(-0.91347778) q[1];
sx q[1];
rz(-1.6003476) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6131564) q[0];
sx q[0];
rz(-0.19599685) q[0];
sx q[0];
rz(-0.71511786) q[0];
x q[1];
rz(-1.7377183) q[2];
sx q[2];
rz(-2.4903565) q[2];
sx q[2];
rz(1.2204352) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.75493357) q[1];
sx q[1];
rz(-1.6303829) q[1];
sx q[1];
rz(0.3646551) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3722911) q[3];
sx q[3];
rz(-0.23450867) q[3];
sx q[3];
rz(-2.086139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9666226) q[2];
sx q[2];
rz(-1.9113767) q[2];
sx q[2];
rz(-0.54565412) q[2];
rz(-0.086006554) q[3];
sx q[3];
rz(-1.6274446) q[3];
sx q[3];
rz(-0.29451323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1957112) q[0];
sx q[0];
rz(-1.042689) q[0];
sx q[0];
rz(-1.9324017) q[0];
rz(-2.53881) q[1];
sx q[1];
rz(-2.5310015) q[1];
sx q[1];
rz(-2.3878154) q[1];
rz(-2.2529765) q[2];
sx q[2];
rz(-1.4063514) q[2];
sx q[2];
rz(0.92264639) q[2];
rz(1.6794593) q[3];
sx q[3];
rz(-2.1908219) q[3];
sx q[3];
rz(0.11809668) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
