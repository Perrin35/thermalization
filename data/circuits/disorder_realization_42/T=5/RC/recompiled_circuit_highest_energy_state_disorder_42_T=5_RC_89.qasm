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
rz(0.49405721) q[0];
sx q[0];
rz(2.9394975) q[0];
sx q[0];
rz(9.1261368) q[0];
rz(2.5477297) q[1];
sx q[1];
rz(-1.9057823) q[1];
sx q[1];
rz(1.8699345) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8182871) q[0];
sx q[0];
rz(-0.092490993) q[0];
sx q[0];
rz(-2.9748067) q[0];
rz(-0.80940445) q[2];
sx q[2];
rz(-2.9706741) q[2];
sx q[2];
rz(-1.6910291) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3385816) q[1];
sx q[1];
rz(-2.525248) q[1];
sx q[1];
rz(-1.0907465) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1966977) q[3];
sx q[3];
rz(-0.44741073) q[3];
sx q[3];
rz(1.5844031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1539845) q[2];
sx q[2];
rz(-2.5488904) q[2];
sx q[2];
rz(0.1565557) q[2];
rz(2.9708059) q[3];
sx q[3];
rz(-0.6905061) q[3];
sx q[3];
rz(1.7295711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60468948) q[0];
sx q[0];
rz(-1.5645744) q[0];
sx q[0];
rz(0.94737303) q[0];
rz(2.0715711) q[1];
sx q[1];
rz(-2.1242296) q[1];
sx q[1];
rz(-1.0191466) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0120901) q[0];
sx q[0];
rz(-2.580205) q[0];
sx q[0];
rz(1.8719869) q[0];
x q[1];
rz(-1.972946) q[2];
sx q[2];
rz(-2.0216591) q[2];
sx q[2];
rz(2.8666945) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.7887914) q[1];
sx q[1];
rz(-1.4104801) q[1];
sx q[1];
rz(-0.43771863) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89313497) q[3];
sx q[3];
rz(-0.97468978) q[3];
sx q[3];
rz(0.79465468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3715839) q[2];
sx q[2];
rz(-1.2086478) q[2];
sx q[2];
rz(-0.51119512) q[2];
rz(-0.33358556) q[3];
sx q[3];
rz(-0.89190069) q[3];
sx q[3];
rz(0.79264486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.12673512) q[0];
sx q[0];
rz(-2.6124586) q[0];
sx q[0];
rz(-2.1281309) q[0];
rz(2.2765989) q[1];
sx q[1];
rz(-1.8987055) q[1];
sx q[1];
rz(-1.6046883) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8415926) q[0];
sx q[0];
rz(-2.9766797) q[0];
sx q[0];
rz(2.4191036) q[0];
rz(-2.2584707) q[2];
sx q[2];
rz(-0.97849023) q[2];
sx q[2];
rz(1.9229213) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5549965) q[1];
sx q[1];
rz(-1.0585016) q[1];
sx q[1];
rz(-1.1050455) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99076267) q[3];
sx q[3];
rz(-1.8988852) q[3];
sx q[3];
rz(-2.6716731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4889739) q[2];
sx q[2];
rz(-1.3451312) q[2];
sx q[2];
rz(-0.49749231) q[2];
rz(-0.8688212) q[3];
sx q[3];
rz(-0.38058623) q[3];
sx q[3];
rz(-3.0697921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3196165) q[0];
sx q[0];
rz(-0.0030586547) q[0];
sx q[0];
rz(2.7724566) q[0];
rz(2.4311851) q[1];
sx q[1];
rz(-1.0514759) q[1];
sx q[1];
rz(1.8416454) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6746691) q[0];
sx q[0];
rz(-0.82991582) q[0];
sx q[0];
rz(-3.039917) q[0];
rz(-0.31710386) q[2];
sx q[2];
rz(-2.6786227) q[2];
sx q[2];
rz(2.3635039) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.59349651) q[1];
sx q[1];
rz(-1.9106459) q[1];
sx q[1];
rz(2.7738394) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10792144) q[3];
sx q[3];
rz(-2.0156337) q[3];
sx q[3];
rz(-0.51714424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9909782) q[2];
sx q[2];
rz(-1.6852448) q[2];
sx q[2];
rz(-0.72188226) q[2];
rz(0.4374831) q[3];
sx q[3];
rz(-2.7719154) q[3];
sx q[3];
rz(-2.8371401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40780145) q[0];
sx q[0];
rz(-2.3423539) q[0];
sx q[0];
rz(0.21402982) q[0];
rz(0.80263823) q[1];
sx q[1];
rz(-1.1624348) q[1];
sx q[1];
rz(0.42073694) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20223742) q[0];
sx q[0];
rz(-1.6488243) q[0];
sx q[0];
rz(-1.3952888) q[0];
rz(2.7710469) q[2];
sx q[2];
rz(-0.59736138) q[2];
sx q[2];
rz(-1.0295402) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.85462697) q[1];
sx q[1];
rz(-1.4564774) q[1];
sx q[1];
rz(3.0752606) q[1];
rz(-pi) q[2];
rz(1.4368128) q[3];
sx q[3];
rz(-2.4089536) q[3];
sx q[3];
rz(2.028038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3558827) q[2];
sx q[2];
rz(-0.37675253) q[2];
sx q[2];
rz(-2.9316087) q[2];
rz(-2.0436132) q[3];
sx q[3];
rz(-2.4185601) q[3];
sx q[3];
rz(-0.22600225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8534773) q[0];
sx q[0];
rz(-2.2324039) q[0];
sx q[0];
rz(0.01532456) q[0];
rz(-0.76459908) q[1];
sx q[1];
rz(-1.2654283) q[1];
sx q[1];
rz(-2.9008124) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6103196) q[0];
sx q[0];
rz(-1.9833779) q[0];
sx q[0];
rz(0.32915769) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50159769) q[2];
sx q[2];
rz(-1.3469509) q[2];
sx q[2];
rz(2.1783592) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.937433) q[1];
sx q[1];
rz(-1.6361409) q[1];
sx q[1];
rz(-1.3974299) q[1];
rz(-pi) q[2];
rz(2.1606275) q[3];
sx q[3];
rz(-2.5348041) q[3];
sx q[3];
rz(-2.2555707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2601872) q[2];
sx q[2];
rz(-2.0892961) q[2];
sx q[2];
rz(-0.50797272) q[2];
rz(-2.9925665) q[3];
sx q[3];
rz(-0.38145724) q[3];
sx q[3];
rz(-1.4582567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7070865) q[0];
sx q[0];
rz(-2.2266946) q[0];
sx q[0];
rz(-1.9515422) q[0];
rz(0.03935628) q[1];
sx q[1];
rz(-1.0330361) q[1];
sx q[1];
rz(2.8796223) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9761115) q[0];
sx q[0];
rz(-1.2369177) q[0];
sx q[0];
rz(-1.6209442) q[0];
rz(-0.33032051) q[2];
sx q[2];
rz(-1.9322504) q[2];
sx q[2];
rz(-0.20955958) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8855454) q[1];
sx q[1];
rz(-2.6441473) q[1];
sx q[1];
rz(0.65620436) q[1];
x q[2];
rz(1.8197601) q[3];
sx q[3];
rz(-2.1002703) q[3];
sx q[3];
rz(-1.362669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1956341) q[2];
sx q[2];
rz(-2.8527263) q[2];
sx q[2];
rz(-1.757901) q[2];
rz(1.9687999) q[3];
sx q[3];
rz(-1.6371509) q[3];
sx q[3];
rz(2.8891532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70343542) q[0];
sx q[0];
rz(-2.1717838) q[0];
sx q[0];
rz(-0.54617149) q[0];
rz(-0.81900412) q[1];
sx q[1];
rz(-1.1823697) q[1];
sx q[1];
rz(-2.6738653) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0696164) q[0];
sx q[0];
rz(-1.3107976) q[0];
sx q[0];
rz(-0.051957038) q[0];
rz(0.39132705) q[2];
sx q[2];
rz(-1.9576009) q[2];
sx q[2];
rz(-2.4387173) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7744171) q[1];
sx q[1];
rz(-1.0760592) q[1];
sx q[1];
rz(-0.53337388) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1728014) q[3];
sx q[3];
rz(-2.8492619) q[3];
sx q[3];
rz(0.013127931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.24432448) q[2];
sx q[2];
rz(-1.6017598) q[2];
sx q[2];
rz(2.3648105) q[2];
rz(-1.6952093) q[3];
sx q[3];
rz(-1.2897688) q[3];
sx q[3];
rz(0.18073925) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1553797) q[0];
sx q[0];
rz(-1.7387094) q[0];
sx q[0];
rz(-3.0631995) q[0];
rz(-0.29131237) q[1];
sx q[1];
rz(-2.4261609) q[1];
sx q[1];
rz(-1.327347) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9531265) q[0];
sx q[0];
rz(-0.92102599) q[0];
sx q[0];
rz(2.3461525) q[0];
x q[1];
rz(-1.615165) q[2];
sx q[2];
rz(-1.67919) q[2];
sx q[2];
rz(2.8381009) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1256333) q[1];
sx q[1];
rz(-1.2984097) q[1];
sx q[1];
rz(-0.65954882) q[1];
rz(-2.6895608) q[3];
sx q[3];
rz(-0.31637438) q[3];
sx q[3];
rz(-1.3133501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.766481) q[2];
sx q[2];
rz(-1.3263308) q[2];
sx q[2];
rz(-0.43409902) q[2];
rz(0.33657524) q[3];
sx q[3];
rz(-1.7124636) q[3];
sx q[3];
rz(2.7460639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8074938) q[0];
sx q[0];
rz(-0.33768001) q[0];
sx q[0];
rz(2.7423774) q[0];
rz(-2.4576064) q[1];
sx q[1];
rz(-1.3147048) q[1];
sx q[1];
rz(1.7302053) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44784495) q[0];
sx q[0];
rz(-1.240792) q[0];
sx q[0];
rz(1.3838612) q[0];
x q[1];
rz(-1.3547276) q[2];
sx q[2];
rz(-2.2377295) q[2];
sx q[2];
rz(-1.9151646) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.017104024) q[1];
sx q[1];
rz(-0.75316562) q[1];
sx q[1];
rz(0.54211771) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78090043) q[3];
sx q[3];
rz(-1.077543) q[3];
sx q[3];
rz(-2.3515612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.6902249) q[2];
sx q[2];
rz(-1.6121612) q[2];
sx q[2];
rz(2.8239047) q[2];
rz(1.0090656) q[3];
sx q[3];
rz(-1.3721481) q[3];
sx q[3];
rz(-0.38869977) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4036564) q[0];
sx q[0];
rz(-1.4690348) q[0];
sx q[0];
rz(0.22672155) q[0];
rz(2.8059088) q[1];
sx q[1];
rz(-2.2102368) q[1];
sx q[1];
rz(-0.23063901) q[1];
rz(3.1138218) q[2];
sx q[2];
rz(-2.4774629) q[2];
sx q[2];
rz(-3.1240875) q[2];
rz(-0.54343358) q[3];
sx q[3];
rz(-2.4590878) q[3];
sx q[3];
rz(1.2868846) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
