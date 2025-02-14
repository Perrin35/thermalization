OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.6505721) q[0];
sx q[0];
rz(0.84492961) q[0];
sx q[0];
rz(6.9965811) q[0];
rz(2.9333935) q[1];
sx q[1];
rz(-0.18728501) q[1];
sx q[1];
rz(-2.0050144) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57713415) q[0];
sx q[0];
rz(-0.73720471) q[0];
sx q[0];
rz(1.2080824) q[0];
rz(3.0805262) q[2];
sx q[2];
rz(-0.89043987) q[2];
sx q[2];
rz(-0.55168569) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7311095) q[1];
sx q[1];
rz(-1.4315194) q[1];
sx q[1];
rz(2.9248505) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4376115) q[3];
sx q[3];
rz(-1.8719421) q[3];
sx q[3];
rz(2.4011321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.95181695) q[2];
sx q[2];
rz(-1.2729898) q[2];
sx q[2];
rz(0.58958685) q[2];
rz(-0.014852614) q[3];
sx q[3];
rz(-1.2735406) q[3];
sx q[3];
rz(-2.5428037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1349161) q[0];
sx q[0];
rz(-2.0966457) q[0];
sx q[0];
rz(2.6892804) q[0];
rz(2.2882838) q[1];
sx q[1];
rz(-2.6494458) q[1];
sx q[1];
rz(-0.36202994) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2483082) q[0];
sx q[0];
rz(-2.4442857) q[0];
sx q[0];
rz(-2.6261283) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8325927) q[2];
sx q[2];
rz(-2.6297792) q[2];
sx q[2];
rz(1.1771894) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7457742) q[1];
sx q[1];
rz(-0.89268301) q[1];
sx q[1];
rz(0.78435244) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70254247) q[3];
sx q[3];
rz(-1.6857393) q[3];
sx q[3];
rz(3.0963932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8028458) q[2];
sx q[2];
rz(-2.5669528) q[2];
sx q[2];
rz(-2.23526) q[2];
rz(1.4551506) q[3];
sx q[3];
rz(-0.95821277) q[3];
sx q[3];
rz(-1.5866535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.0872658) q[0];
sx q[0];
rz(-1.9540906) q[0];
sx q[0];
rz(-0.97386709) q[0];
rz(-2.5663238) q[1];
sx q[1];
rz(-1.560805) q[1];
sx q[1];
rz(1.5781933) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3546936) q[0];
sx q[0];
rz(-1.5074566) q[0];
sx q[0];
rz(3.1194434) q[0];
rz(0.44610666) q[2];
sx q[2];
rz(-2.2488454) q[2];
sx q[2];
rz(-2.7715832) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0525731) q[1];
sx q[1];
rz(-1.1088437) q[1];
sx q[1];
rz(-0.79098742) q[1];
rz(1.7922014) q[3];
sx q[3];
rz(-2.3292037) q[3];
sx q[3];
rz(-2.098907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8156208) q[2];
sx q[2];
rz(-1.1824563) q[2];
sx q[2];
rz(-1.6131442) q[2];
rz(-0.01550393) q[3];
sx q[3];
rz(-1.9663234) q[3];
sx q[3];
rz(-1.6541245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7148606) q[0];
sx q[0];
rz(-0.69789129) q[0];
sx q[0];
rz(-3.0791855) q[0];
rz(1.9852253) q[1];
sx q[1];
rz(-0.9461177) q[1];
sx q[1];
rz(0.83414331) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1301889) q[0];
sx q[0];
rz(-1.0727278) q[0];
sx q[0];
rz(2.8047436) q[0];
rz(0.049872204) q[2];
sx q[2];
rz(-1.4446752) q[2];
sx q[2];
rz(-1.0654861) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8319611) q[1];
sx q[1];
rz(-2.1681962) q[1];
sx q[1];
rz(-0.064927622) q[1];
x q[2];
rz(2.0888791) q[3];
sx q[3];
rz(-1.7527767) q[3];
sx q[3];
rz(-0.35634867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1674898) q[2];
sx q[2];
rz(-0.17593273) q[2];
sx q[2];
rz(-1.9195732) q[2];
rz(-2.7410298) q[3];
sx q[3];
rz(-1.0223072) q[3];
sx q[3];
rz(-2.9266973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0790613) q[0];
sx q[0];
rz(-2.489594) q[0];
sx q[0];
rz(0.25667152) q[0];
rz(1.7968862) q[1];
sx q[1];
rz(-1.2213629) q[1];
sx q[1];
rz(1.7324956) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042085083) q[0];
sx q[0];
rz(-1.2718023) q[0];
sx q[0];
rz(-0.87845699) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95375421) q[2];
sx q[2];
rz(-0.34025345) q[2];
sx q[2];
rz(-2.798852) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.46542922) q[1];
sx q[1];
rz(-1.5309285) q[1];
sx q[1];
rz(-2.7207147) q[1];
rz(-pi) q[2];
x q[2];
rz(0.6677169) q[3];
sx q[3];
rz(-1.8103292) q[3];
sx q[3];
rz(2.4706877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1008489) q[2];
sx q[2];
rz(-2.433653) q[2];
sx q[2];
rz(2.9948998) q[2];
rz(-1.3692859) q[3];
sx q[3];
rz(-0.36543235) q[3];
sx q[3];
rz(0.2379612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20928243) q[0];
sx q[0];
rz(-2.3754061) q[0];
sx q[0];
rz(-0.76195088) q[0];
rz(-2.2867098) q[1];
sx q[1];
rz(-1.8652752) q[1];
sx q[1];
rz(-0.68861047) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013278339) q[0];
sx q[0];
rz(-1.6546895) q[0];
sx q[0];
rz(1.1021139) q[0];
rz(-pi) q[1];
rz(-2.9607331) q[2];
sx q[2];
rz(-3.0263889) q[2];
sx q[2];
rz(1.2173139) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.411158) q[1];
sx q[1];
rz(-1.6741006) q[1];
sx q[1];
rz(1.5048601) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1375539) q[3];
sx q[3];
rz(-1.8747624) q[3];
sx q[3];
rz(-0.36813018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7410572) q[2];
sx q[2];
rz(-0.94509882) q[2];
sx q[2];
rz(-0.98089027) q[2];
rz(-0.26997057) q[3];
sx q[3];
rz(-1.910784) q[3];
sx q[3];
rz(3.1321757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
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
rz(0.16801676) q[0];
sx q[0];
rz(-1.4485757) q[0];
sx q[0];
rz(-2.5981405) q[0];
rz(1.5914397) q[1];
sx q[1];
rz(-2.2563069) q[1];
sx q[1];
rz(-2.1563704) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60632432) q[0];
sx q[0];
rz(-0.54161763) q[0];
sx q[0];
rz(-0.4569502) q[0];
rz(2.7706861) q[2];
sx q[2];
rz(-2.1044253) q[2];
sx q[2];
rz(0.4069538) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7581525) q[1];
sx q[1];
rz(-1.952127) q[1];
sx q[1];
rz(2.1995981) q[1];
x q[2];
rz(-1.9051294) q[3];
sx q[3];
rz(-2.0591934) q[3];
sx q[3];
rz(0.56570977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.33124179) q[2];
sx q[2];
rz(-2.6369429) q[2];
sx q[2];
rz(-1.8028367) q[2];
rz(-1.4618368) q[3];
sx q[3];
rz(-1.5906518) q[3];
sx q[3];
rz(-2.4707879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(0.93002334) q[0];
sx q[0];
rz(-2.0116563) q[0];
sx q[0];
rz(-1.2343963) q[0];
rz(1.3537539) q[1];
sx q[1];
rz(-0.45629558) q[1];
sx q[1];
rz(0.096045883) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9864842) q[0];
sx q[0];
rz(-2.2791515) q[0];
sx q[0];
rz(-2.1848795) q[0];
rz(-pi) q[1];
rz(1.646461) q[2];
sx q[2];
rz(-1.0319508) q[2];
sx q[2];
rz(-2.7494753) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9191285) q[1];
sx q[1];
rz(-1.1908731) q[1];
sx q[1];
rz(1.6187731) q[1];
rz(-pi) q[2];
rz(1.6173564) q[3];
sx q[3];
rz(-0.73107204) q[3];
sx q[3];
rz(-2.1995167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8058406) q[2];
sx q[2];
rz(-2.7554607) q[2];
sx q[2];
rz(-1.1248355) q[2];
rz(2.2624894) q[3];
sx q[3];
rz(-1.2259038) q[3];
sx q[3];
rz(0.44900289) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45824555) q[0];
sx q[0];
rz(-1.6793716) q[0];
sx q[0];
rz(0.71518389) q[0];
rz(-0.31582754) q[1];
sx q[1];
rz(-0.68604699) q[1];
sx q[1];
rz(-2.979523) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93955646) q[0];
sx q[0];
rz(-1.9959772) q[0];
sx q[0];
rz(2.3935938) q[0];
x q[1];
rz(0.58312441) q[2];
sx q[2];
rz(-1.4004502) q[2];
sx q[2];
rz(1.9161461) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9980364) q[1];
sx q[1];
rz(-0.85696917) q[1];
sx q[1];
rz(-1.803323) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81406419) q[3];
sx q[3];
rz(-2.6912675) q[3];
sx q[3];
rz(1.361477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38966236) q[2];
sx q[2];
rz(-1.1952362) q[2];
sx q[2];
rz(1.2644838) q[2];
rz(-1.8090931) q[3];
sx q[3];
rz(-1.8615362) q[3];
sx q[3];
rz(2.8247824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7952809) q[0];
sx q[0];
rz(-0.77487159) q[0];
sx q[0];
rz(2.0922022) q[0];
rz(-0.28114444) q[1];
sx q[1];
rz(-2.099359) q[1];
sx q[1];
rz(1.3220471) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78681386) q[0];
sx q[0];
rz(-0.6118868) q[0];
sx q[0];
rz(-2.6953681) q[0];
rz(-pi) q[1];
rz(-0.63788267) q[2];
sx q[2];
rz(-1.1670127) q[2];
sx q[2];
rz(2.11602) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.57105455) q[1];
sx q[1];
rz(-1.7295383) q[1];
sx q[1];
rz(-1.2069993) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4720143) q[3];
sx q[3];
rz(-1.9552468) q[3];
sx q[3];
rz(-0.41562322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37810024) q[2];
sx q[2];
rz(-1.8201479) q[2];
sx q[2];
rz(0.05149252) q[2];
rz(-1.817305) q[3];
sx q[3];
rz(-1.4824425) q[3];
sx q[3];
rz(1.9342559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0180203) q[0];
sx q[0];
rz(-1.3326895) q[0];
sx q[0];
rz(1.4154758) q[0];
rz(2.3208658) q[1];
sx q[1];
rz(-1.6976994) q[1];
sx q[1];
rz(1.2669947) q[1];
rz(2.5222798) q[2];
sx q[2];
rz(-1.7650585) q[2];
sx q[2];
rz(-1.4256918) q[2];
rz(-2.5883417) q[3];
sx q[3];
rz(-0.58589952) q[3];
sx q[3];
rz(0.70675969) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
