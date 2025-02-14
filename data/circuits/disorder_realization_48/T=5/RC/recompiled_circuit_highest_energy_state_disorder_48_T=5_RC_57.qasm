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
rz(0.91520619) q[0];
sx q[0];
rz(-0.45867607) q[0];
sx q[0];
rz(-2.8235478) q[0];
rz(-1.7755427) q[1];
sx q[1];
rz(-0.69861689) q[1];
sx q[1];
rz(3.0787992) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1695367) q[0];
sx q[0];
rz(-2.0436624) q[0];
sx q[0];
rz(2.1465143) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89102052) q[2];
sx q[2];
rz(-1.2212379) q[2];
sx q[2];
rz(-1.0206136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1789238) q[1];
sx q[1];
rz(-2.0589863) q[1];
sx q[1];
rz(-2.8250384) q[1];
rz(-pi) q[2];
rz(2.7460119) q[3];
sx q[3];
rz(-2.1115344) q[3];
sx q[3];
rz(0.21025019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.723145) q[2];
sx q[2];
rz(-1.8258839) q[2];
sx q[2];
rz(1.0203699) q[2];
rz(2.4127035) q[3];
sx q[3];
rz(-2.708669) q[3];
sx q[3];
rz(-1.5322022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6272524) q[0];
sx q[0];
rz(-1.9753375) q[0];
sx q[0];
rz(-0.038272055) q[0];
rz(-3.017784) q[1];
sx q[1];
rz(-0.47272155) q[1];
sx q[1];
rz(-1.5708057) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8077652) q[0];
sx q[0];
rz(-0.021919576) q[0];
sx q[0];
rz(1.0967361) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0395987) q[2];
sx q[2];
rz(-0.77384678) q[2];
sx q[2];
rz(0.71600658) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.96882183) q[1];
sx q[1];
rz(-1.570373) q[1];
sx q[1];
rz(-1.5718196) q[1];
rz(0.51608508) q[3];
sx q[3];
rz(-1.6760964) q[3];
sx q[3];
rz(2.1886231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6191285) q[2];
sx q[2];
rz(-1.8121441) q[2];
sx q[2];
rz(0.5788571) q[2];
rz(-1.045643) q[3];
sx q[3];
rz(-2.3561616) q[3];
sx q[3];
rz(-2.3845909) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4724562) q[0];
sx q[0];
rz(-2.0643015) q[0];
sx q[0];
rz(2.6876167) q[0];
rz(-0.16481915) q[1];
sx q[1];
rz(-0.77607981) q[1];
sx q[1];
rz(-1.6710612) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31904083) q[0];
sx q[0];
rz(-1.3994354) q[0];
sx q[0];
rz(1.9644072) q[0];
x q[1];
rz(-1.0376588) q[2];
sx q[2];
rz(-2.5488944) q[2];
sx q[2];
rz(-1.8225841) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9413621) q[1];
sx q[1];
rz(-0.75298568) q[1];
sx q[1];
rz(2.689207) q[1];
rz(-pi) q[2];
rz(-0.65916797) q[3];
sx q[3];
rz(-1.9999095) q[3];
sx q[3];
rz(1.3339213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7620324) q[2];
sx q[2];
rz(-1.6635868) q[2];
sx q[2];
rz(-1.7851768) q[2];
rz(-0.49946579) q[3];
sx q[3];
rz(-1.5288589) q[3];
sx q[3];
rz(-1.0511901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22719638) q[0];
sx q[0];
rz(-0.90407404) q[0];
sx q[0];
rz(-2.0830578) q[0];
rz(1.1610441) q[1];
sx q[1];
rz(-0.5608905) q[1];
sx q[1];
rz(3.0243691) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0441705) q[0];
sx q[0];
rz(-1.6099842) q[0];
sx q[0];
rz(-0.54255658) q[0];
rz(-pi) q[1];
rz(2.9449844) q[2];
sx q[2];
rz(-0.93485445) q[2];
sx q[2];
rz(-1.3351118) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8414054) q[1];
sx q[1];
rz(-2.3267728) q[1];
sx q[1];
rz(2.7542265) q[1];
rz(-pi) q[2];
rz(2.6113308) q[3];
sx q[3];
rz(-0.70928364) q[3];
sx q[3];
rz(-1.5269296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3566572) q[2];
sx q[2];
rz(-1.4484118) q[2];
sx q[2];
rz(1.3680722) q[2];
rz(1.28537) q[3];
sx q[3];
rz(-2.0338438) q[3];
sx q[3];
rz(2.4135597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042498978) q[0];
sx q[0];
rz(-0.92785257) q[0];
sx q[0];
rz(-1.4121144) q[0];
rz(2.1389351) q[1];
sx q[1];
rz(-2.899677) q[1];
sx q[1];
rz(1.5826506) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2501734) q[0];
sx q[0];
rz(-2.1262263) q[0];
sx q[0];
rz(2.069418) q[0];
rz(-pi) q[1];
rz(1.6539668) q[2];
sx q[2];
rz(-0.85415936) q[2];
sx q[2];
rz(-2.4019474) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4950558) q[1];
sx q[1];
rz(-0.55083129) q[1];
sx q[1];
rz(1.6158094) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2941957) q[3];
sx q[3];
rz(-2.8081354) q[3];
sx q[3];
rz(-1.6113623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5530508) q[2];
sx q[2];
rz(-1.43247) q[2];
sx q[2];
rz(0.39349619) q[2];
rz(-2.1816175) q[3];
sx q[3];
rz(-0.67592755) q[3];
sx q[3];
rz(-2.1737607) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8680962) q[0];
sx q[0];
rz(-0.12729004) q[0];
sx q[0];
rz(-2.9582276) q[0];
rz(-0.49194899) q[1];
sx q[1];
rz(-0.79194561) q[1];
sx q[1];
rz(1.9221745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3122092) q[0];
sx q[0];
rz(-1.6512733) q[0];
sx q[0];
rz(0.089437251) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6492789) q[2];
sx q[2];
rz(-1.6011136) q[2];
sx q[2];
rz(-0.91532367) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1508091) q[1];
sx q[1];
rz(-1.1649302) q[1];
sx q[1];
rz(1.957914) q[1];
x q[2];
rz(-1.0991092) q[3];
sx q[3];
rz(-1.6460643) q[3];
sx q[3];
rz(-1.690563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0101953) q[2];
sx q[2];
rz(-1.7051899) q[2];
sx q[2];
rz(-0.08237002) q[2];
rz(-2.2149337) q[3];
sx q[3];
rz(-0.72946531) q[3];
sx q[3];
rz(-1.5026708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6019186) q[0];
sx q[0];
rz(-2.813756) q[0];
sx q[0];
rz(-0.71640054) q[0];
rz(2.7377103) q[1];
sx q[1];
rz(-1.5733746) q[1];
sx q[1];
rz(1.7049888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8536896) q[0];
sx q[0];
rz(-2.7349964) q[0];
sx q[0];
rz(1.9415744) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9173849) q[2];
sx q[2];
rz(-0.53935236) q[2];
sx q[2];
rz(-1.8856848) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.75175367) q[1];
sx q[1];
rz(-1.1734278) q[1];
sx q[1];
rz(1.4629055) q[1];
rz(-pi) q[2];
rz(-0.72258805) q[3];
sx q[3];
rz(-2.1919985) q[3];
sx q[3];
rz(0.91739839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3639823) q[2];
sx q[2];
rz(-0.97475514) q[2];
sx q[2];
rz(3.0762365) q[2];
rz(-0.86723793) q[3];
sx q[3];
rz(-1.6403653) q[3];
sx q[3];
rz(-1.8170099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48677483) q[0];
sx q[0];
rz(-1.7205394) q[0];
sx q[0];
rz(0.73100334) q[0];
rz(0.77516088) q[1];
sx q[1];
rz(-2.8072) q[1];
sx q[1];
rz(-2.7269272) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1486007) q[0];
sx q[0];
rz(-2.1777657) q[0];
sx q[0];
rz(0.66120045) q[0];
rz(-0.35279775) q[2];
sx q[2];
rz(-0.55856201) q[2];
sx q[2];
rz(1.1272023) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9379256) q[1];
sx q[1];
rz(-1.3199184) q[1];
sx q[1];
rz(-0.30955065) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.88736262) q[3];
sx q[3];
rz(-0.21036869) q[3];
sx q[3];
rz(0.31926814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9097462) q[2];
sx q[2];
rz(-0.48603386) q[2];
sx q[2];
rz(-0.092378423) q[2];
rz(2.3653024) q[3];
sx q[3];
rz(-1.679136) q[3];
sx q[3];
rz(-2.9379454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6542776) q[0];
sx q[0];
rz(-1.4947816) q[0];
sx q[0];
rz(0.017729433) q[0];
rz(-2.0954633) q[1];
sx q[1];
rz(-1.4132063) q[1];
sx q[1];
rz(2.0808992) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3002787) q[0];
sx q[0];
rz(-1.5567686) q[0];
sx q[0];
rz(3.0728691) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1938042) q[2];
sx q[2];
rz(-0.62219884) q[2];
sx q[2];
rz(2.5925328) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2613147) q[1];
sx q[1];
rz(-2.3572817) q[1];
sx q[1];
rz(1.8159291) q[1];
rz(-pi) q[2];
rz(1.6298196) q[3];
sx q[3];
rz(-1.0477598) q[3];
sx q[3];
rz(-0.28429261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9866508) q[2];
sx q[2];
rz(-2.3255746) q[2];
sx q[2];
rz(2.4972534) q[2];
rz(-2.2996357) q[3];
sx q[3];
rz(-1.569845) q[3];
sx q[3];
rz(-2.2346066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3802721) q[0];
sx q[0];
rz(-0.43571061) q[0];
sx q[0];
rz(-0.45183387) q[0];
rz(-1.0093581) q[1];
sx q[1];
rz(-1.6296856) q[1];
sx q[1];
rz(1.5712646) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6752351) q[0];
sx q[0];
rz(-1.5769099) q[0];
sx q[0];
rz(-2.9714022) q[0];
rz(-pi) q[1];
rz(1.794593) q[2];
sx q[2];
rz(-0.48529709) q[2];
sx q[2];
rz(-1.8388302) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.45598509) q[1];
sx q[1];
rz(-2.7311014) q[1];
sx q[1];
rz(-1.2430771) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5630003) q[3];
sx q[3];
rz(-0.8484133) q[3];
sx q[3];
rz(-0.16297271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5181804) q[2];
sx q[2];
rz(-2.0698915) q[2];
sx q[2];
rz(-1.5706459) q[2];
rz(-3.0359641) q[3];
sx q[3];
rz(-0.94625866) q[3];
sx q[3];
rz(-0.51641974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.327772) q[0];
sx q[0];
rz(-2.1229424) q[0];
sx q[0];
rz(0.13737296) q[0];
rz(0.2188006) q[1];
sx q[1];
rz(-1.4004424) q[1];
sx q[1];
rz(-0.91011824) q[1];
rz(2.0135638) q[2];
sx q[2];
rz(-1.762286) q[2];
sx q[2];
rz(1.0049835) q[2];
rz(0.46047138) q[3];
sx q[3];
rz(-2.4302767) q[3];
sx q[3];
rz(-2.5828491) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
