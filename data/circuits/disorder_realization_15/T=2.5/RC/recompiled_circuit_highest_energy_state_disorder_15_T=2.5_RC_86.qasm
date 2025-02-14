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
rz(-2.5315142) q[0];
sx q[0];
rz(-0.22992034) q[0];
sx q[0];
rz(0.71969405) q[0];
rz(2.5201058) q[1];
sx q[1];
rz(-1.7401594) q[1];
sx q[1];
rz(0.12535867) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1272362) q[0];
sx q[0];
rz(-1.171021) q[0];
sx q[0];
rz(-2.6771604) q[0];
rz(-pi) q[1];
rz(-2.8530082) q[2];
sx q[2];
rz(-1.6615531) q[2];
sx q[2];
rz(-1.2707658) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92491313) q[1];
sx q[1];
rz(-1.5694111) q[1];
sx q[1];
rz(3.1411704) q[1];
rz(1.9699279) q[3];
sx q[3];
rz(-2.3775775) q[3];
sx q[3];
rz(-1.3543874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7980935) q[2];
sx q[2];
rz(-0.40830475) q[2];
sx q[2];
rz(0.84585345) q[2];
rz(0.79126233) q[3];
sx q[3];
rz(-0.013412272) q[3];
sx q[3];
rz(-3.0748034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4132408) q[0];
sx q[0];
rz(-2.6744196) q[0];
sx q[0];
rz(3.0979284) q[0];
rz(1.5751669) q[1];
sx q[1];
rz(-1.7685726) q[1];
sx q[1];
rz(1.6427737) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18066809) q[0];
sx q[0];
rz(-2.5378413) q[0];
sx q[0];
rz(-1.5636428) q[0];
x q[1];
rz(0.0129406) q[2];
sx q[2];
rz(-0.99943752) q[2];
sx q[2];
rz(0.017627942) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1825492) q[1];
sx q[1];
rz(-1.601853) q[1];
sx q[1];
rz(-3.0654383) q[1];
rz(-pi) q[2];
rz(1.6632027) q[3];
sx q[3];
rz(-0.3457717) q[3];
sx q[3];
rz(-1.913663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0955536) q[2];
sx q[2];
rz(-2.9914896) q[2];
sx q[2];
rz(-2.6138439) q[2];
rz(2.2857417) q[3];
sx q[3];
rz(-0.0015043613) q[3];
sx q[3];
rz(1.2214448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7432778) q[0];
sx q[0];
rz(-2.1687431) q[0];
sx q[0];
rz(1.995218) q[0];
rz(1.3829117) q[1];
sx q[1];
rz(-0.29255602) q[1];
sx q[1];
rz(-0.10398277) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87097528) q[0];
sx q[0];
rz(-1.7954602) q[0];
sx q[0];
rz(1.490834) q[0];
x q[1];
rz(3.1244794) q[2];
sx q[2];
rz(-1.8430954) q[2];
sx q[2];
rz(0.903331) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.96668916) q[1];
sx q[1];
rz(-2.3123868) q[1];
sx q[1];
rz(-1.8494383) q[1];
rz(0.29071112) q[3];
sx q[3];
rz(-1.4658648) q[3];
sx q[3];
rz(-0.19874979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9581703) q[2];
sx q[2];
rz(-3.1347745) q[2];
sx q[2];
rz(-2.5799694) q[2];
rz(-3.0567452) q[3];
sx q[3];
rz(-3.1361129) q[3];
sx q[3];
rz(0.00076278846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7247923) q[0];
sx q[0];
rz(-2.8775207) q[0];
sx q[0];
rz(0.52077878) q[0];
rz(2.9837823) q[1];
sx q[1];
rz(-0.66693711) q[1];
sx q[1];
rz(3.0657943) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6822081) q[0];
sx q[0];
rz(-2.5494116) q[0];
sx q[0];
rz(-2.9367052) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5694322) q[2];
sx q[2];
rz(-1.5714386) q[2];
sx q[2];
rz(0.13880348) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8677249) q[1];
sx q[1];
rz(-0.40479044) q[1];
sx q[1];
rz(1.2048436) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0337379) q[3];
sx q[3];
rz(-1.1198992) q[3];
sx q[3];
rz(-1.0552561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.51125222) q[2];
sx q[2];
rz(-3.1255836) q[2];
sx q[2];
rz(1.4207669) q[2];
rz(3.1347647) q[3];
sx q[3];
rz(-0.029191645) q[3];
sx q[3];
rz(-1.487287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1348949) q[0];
sx q[0];
rz(-0.20773523) q[0];
sx q[0];
rz(-2.8790706) q[0];
rz(0.94995704) q[1];
sx q[1];
rz(-0.078153178) q[1];
sx q[1];
rz(2.9238759) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9757468) q[0];
sx q[0];
rz(-1.6941771) q[0];
sx q[0];
rz(1.1749772) q[0];
x q[1];
rz(-0.08723549) q[2];
sx q[2];
rz(-1.6542742) q[2];
sx q[2];
rz(-1.491445) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7568183) q[1];
sx q[1];
rz(-0.17233241) q[1];
sx q[1];
rz(1.5583355) q[1];
rz(3.0721139) q[3];
sx q[3];
rz(-1.0437878) q[3];
sx q[3];
rz(-1.219092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.21506423) q[2];
sx q[2];
rz(-0.0096409163) q[2];
sx q[2];
rz(-2.8936774) q[2];
rz(1.7667814) q[3];
sx q[3];
rz(-3.091076) q[3];
sx q[3];
rz(-1.4352528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.0404102) q[0];
sx q[0];
rz(-1.269416) q[0];
sx q[0];
rz(-2.5293479) q[0];
rz(-2.9712037) q[1];
sx q[1];
rz(-3.0590765) q[1];
sx q[1];
rz(1.6498529) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1338324) q[0];
sx q[0];
rz(-0.88307021) q[0];
sx q[0];
rz(2.4439815) q[0];
x q[1];
rz(0.014530226) q[2];
sx q[2];
rz(-1.5749914) q[2];
sx q[2];
rz(-0.99630492) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.93894847) q[1];
sx q[1];
rz(-0.30057014) q[1];
sx q[1];
rz(-0.2267129) q[1];
rz(2.0759367) q[3];
sx q[3];
rz(-2.7759984) q[3];
sx q[3];
rz(-0.58581173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6461569) q[2];
sx q[2];
rz(-0.010224552) q[2];
sx q[2];
rz(-2.9025027) q[2];
rz(-0.80860364) q[3];
sx q[3];
rz(-3.1294398) q[3];
sx q[3];
rz(1.3330207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7072153) q[0];
sx q[0];
rz(-3.0476397) q[0];
sx q[0];
rz(-2.3648426) q[0];
rz(3.1306733) q[1];
sx q[1];
rz(-0.24591406) q[1];
sx q[1];
rz(-1.6687261) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1253788) q[0];
sx q[0];
rz(-0.74170602) q[0];
sx q[0];
rz(1.2522231) q[0];
rz(-1.5644349) q[2];
sx q[2];
rz(-1.5472876) q[2];
sx q[2];
rz(-2.9671217) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5653021) q[1];
sx q[1];
rz(-1.6716752) q[1];
sx q[1];
rz(-2.3957344) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3192572) q[3];
sx q[3];
rz(-2.8467379) q[3];
sx q[3];
rz(2.7406462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9933068) q[2];
sx q[2];
rz(-3.1243656) q[2];
sx q[2];
rz(-0.41603184) q[2];
rz(2.3038583) q[3];
sx q[3];
rz(-0.13844027) q[3];
sx q[3];
rz(1.4778888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1763022) q[0];
sx q[0];
rz(-3.1017113) q[0];
sx q[0];
rz(-0.95529977) q[0];
rz(1.2166066) q[1];
sx q[1];
rz(-0.33733264) q[1];
sx q[1];
rz(-0.40562707) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5914456) q[0];
sx q[0];
rz(-0.87021512) q[0];
sx q[0];
rz(-1.1484543) q[0];
x q[1];
rz(1.4816059) q[2];
sx q[2];
rz(-1.729166) q[2];
sx q[2];
rz(-1.0430481) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4279528) q[1];
sx q[1];
rz(-1.7312839) q[1];
sx q[1];
rz(-1.493371) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6576976) q[3];
sx q[3];
rz(-1.7830866) q[3];
sx q[3];
rz(-0.014537285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2847298) q[2];
sx q[2];
rz(-0.036981985) q[2];
sx q[2];
rz(0.82376087) q[2];
rz(0.13122261) q[3];
sx q[3];
rz(-0.28990144) q[3];
sx q[3];
rz(-0.32740617) q[3];
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
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4284978) q[0];
sx q[0];
rz(-2.9976124) q[0];
sx q[0];
rz(-2.4289828) q[0];
rz(0.54829848) q[1];
sx q[1];
rz(-2.8916841) q[1];
sx q[1];
rz(-2.9329494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0731733) q[0];
sx q[0];
rz(-1.1688587) q[0];
sx q[0];
rz(-2.7059024) q[0];
x q[1];
rz(-1.0717355) q[2];
sx q[2];
rz(-0.03943561) q[2];
sx q[2];
rz(-0.60434228) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1125579) q[1];
sx q[1];
rz(-1.0150954) q[1];
sx q[1];
rz(-3.1359926) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94722469) q[3];
sx q[3];
rz(-2.1125018) q[3];
sx q[3];
rz(0.052963363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6930406) q[2];
sx q[2];
rz(-0.081144944) q[2];
sx q[2];
rz(-0.21716675) q[2];
rz(-0.36755696) q[3];
sx q[3];
rz(-0.03438545) q[3];
sx q[3];
rz(-0.99678451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9656068) q[0];
sx q[0];
rz(-0.088986926) q[0];
sx q[0];
rz(2.9678645) q[0];
rz(1.732775) q[1];
sx q[1];
rz(-1.6830187) q[1];
sx q[1];
rz(1.4558314) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30688863) q[0];
sx q[0];
rz(-1.7346037) q[0];
sx q[0];
rz(-2.6192946) q[0];
rz(-pi) q[1];
rz(-0.98241229) q[2];
sx q[2];
rz(-0.2478226) q[2];
sx q[2];
rz(-0.10977015) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0152564) q[1];
sx q[1];
rz(-0.70955196) q[1];
sx q[1];
rz(0.1800329) q[1];
rz(-pi) q[2];
x q[2];
rz(0.10723857) q[3];
sx q[3];
rz(-1.3235561) q[3];
sx q[3];
rz(3.1013754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9659861) q[2];
sx q[2];
rz(-0.49314988) q[2];
sx q[2];
rz(1.7612339) q[2];
rz(-2.4557451) q[3];
sx q[3];
rz(-0.0018456056) q[3];
sx q[3];
rz(0.67870158) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3995517) q[0];
sx q[0];
rz(-2.4492332) q[0];
sx q[0];
rz(-1.5204182) q[0];
rz(1.5834658) q[1];
sx q[1];
rz(-1.635066) q[1];
sx q[1];
rz(0.20546694) q[1];
rz(0.016318446) q[2];
sx q[2];
rz(-1.6961799) q[2];
sx q[2];
rz(0.21519306) q[2];
rz(1.2027005) q[3];
sx q[3];
rz(-2.8218269) q[3];
sx q[3];
rz(0.091824986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
