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
rz(0.49993604) q[0];
sx q[0];
rz(-1.191782) q[0];
sx q[0];
rz(-1.3718104) q[0];
rz(-1.491188) q[1];
sx q[1];
rz(-2.2743382) q[1];
sx q[1];
rz(2.7452918) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9306878) q[0];
sx q[0];
rz(-0.22432029) q[0];
sx q[0];
rz(0.36422642) q[0];
rz(-1.9849586) q[2];
sx q[2];
rz(-2.7348619) q[2];
sx q[2];
rz(1.0309645) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.22452488) q[1];
sx q[1];
rz(-2.6049419) q[1];
sx q[1];
rz(1.3858252) q[1];
rz(-pi) q[2];
rz(1.8627152) q[3];
sx q[3];
rz(-1.9510428) q[3];
sx q[3];
rz(1.0060665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1402011) q[2];
sx q[2];
rz(-2.3346257) q[2];
sx q[2];
rz(-1.8001455) q[2];
rz(2.0724824) q[3];
sx q[3];
rz(-0.1461229) q[3];
sx q[3];
rz(1.7599546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5694201) q[0];
sx q[0];
rz(-2.2241346) q[0];
sx q[0];
rz(0.65895748) q[0];
rz(0.072703687) q[1];
sx q[1];
rz(-2.0854988) q[1];
sx q[1];
rz(-1.8812077) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9958268) q[0];
sx q[0];
rz(-0.46016177) q[0];
sx q[0];
rz(-1.3868679) q[0];
rz(-0.44620338) q[2];
sx q[2];
rz(-1.1006736) q[2];
sx q[2];
rz(-2.0044495) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.61140984) q[1];
sx q[1];
rz(-1.5059222) q[1];
sx q[1];
rz(-0.26574175) q[1];
rz(-pi) q[2];
rz(-2.1815592) q[3];
sx q[3];
rz(-1.6690147) q[3];
sx q[3];
rz(1.4032422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.8947123) q[2];
sx q[2];
rz(-2.3983045) q[2];
sx q[2];
rz(-1.7903719) q[2];
rz(0.61659914) q[3];
sx q[3];
rz(-2.1117881) q[3];
sx q[3];
rz(-0.29115796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63224822) q[0];
sx q[0];
rz(-2.6743439) q[0];
sx q[0];
rz(-2.9314991) q[0];
rz(3.0585152) q[1];
sx q[1];
rz(-1.2385085) q[1];
sx q[1];
rz(-3.074379) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7930896) q[0];
sx q[0];
rz(-1.5572892) q[0];
sx q[0];
rz(1.5828787) q[0];
rz(-2.1234197) q[2];
sx q[2];
rz(-0.45028307) q[2];
sx q[2];
rz(-1.0339111) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76308474) q[1];
sx q[1];
rz(-2.3649667) q[1];
sx q[1];
rz(1.08294) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70560734) q[3];
sx q[3];
rz(-1.5021245) q[3];
sx q[3];
rz(-2.9182383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82789603) q[2];
sx q[2];
rz(-1.6500429) q[2];
sx q[2];
rz(1.3649887) q[2];
rz(2.7374173) q[3];
sx q[3];
rz(-1.3173236) q[3];
sx q[3];
rz(-1.1990168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2760524) q[0];
sx q[0];
rz(-1.4006389) q[0];
sx q[0];
rz(-2.6633967) q[0];
rz(-2.1887691) q[1];
sx q[1];
rz(-2.9915504) q[1];
sx q[1];
rz(0.40963867) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38423019) q[0];
sx q[0];
rz(-0.91529796) q[0];
sx q[0];
rz(-1.537498) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6851321) q[2];
sx q[2];
rz(-1.6586411) q[2];
sx q[2];
rz(2.8473471) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8781902) q[1];
sx q[1];
rz(-0.95118427) q[1];
sx q[1];
rz(2.3248169) q[1];
x q[2];
rz(-1.2495991) q[3];
sx q[3];
rz(-1.7068591) q[3];
sx q[3];
rz(2.4709159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59565583) q[2];
sx q[2];
rz(-1.4384392) q[2];
sx q[2];
rz(2.0036009) q[2];
rz(-0.99500895) q[3];
sx q[3];
rz(-2.5840839) q[3];
sx q[3];
rz(-0.086624302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.397641) q[0];
sx q[0];
rz(-1.3084363) q[0];
sx q[0];
rz(-0.27807903) q[0];
rz(-1.8771578) q[1];
sx q[1];
rz(-1.8587298) q[1];
sx q[1];
rz(1.5637195) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0959543) q[0];
sx q[0];
rz(-1.6976893) q[0];
sx q[0];
rz(-1.9357827) q[0];
x q[1];
rz(0.78863849) q[2];
sx q[2];
rz(-1.83162) q[2];
sx q[2];
rz(2.3344085) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3914096) q[1];
sx q[1];
rz(-0.58779189) q[1];
sx q[1];
rz(-2.9255735) q[1];
rz(1.4117354) q[3];
sx q[3];
rz(-1.4791227) q[3];
sx q[3];
rz(-1.5605469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9331253) q[2];
sx q[2];
rz(-1.5310023) q[2];
sx q[2];
rz(-0.36766407) q[2];
rz(-1.095088) q[3];
sx q[3];
rz(-1.0235267) q[3];
sx q[3];
rz(0.17525214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.6943738) q[0];
sx q[0];
rz(-2.8419438) q[0];
sx q[0];
rz(-1.8023941) q[0];
rz(-0.12204349) q[1];
sx q[1];
rz(-1.782878) q[1];
sx q[1];
rz(-0.36062127) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80194762) q[0];
sx q[0];
rz(-0.80251283) q[0];
sx q[0];
rz(-2.5967477) q[0];
rz(1.6595938) q[2];
sx q[2];
rz(-0.42150233) q[2];
sx q[2];
rz(0.56383609) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.569446) q[1];
sx q[1];
rz(-2.0098369) q[1];
sx q[1];
rz(-2.8296058) q[1];
rz(-pi) q[2];
x q[2];
rz(3.120947) q[3];
sx q[3];
rz(-0.25293487) q[3];
sx q[3];
rz(-1.3731352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.7755166) q[2];
sx q[2];
rz(-1.4630432) q[2];
sx q[2];
rz(-2.0029081) q[2];
rz(2.8876997) q[3];
sx q[3];
rz(-2.4963899) q[3];
sx q[3];
rz(2.3937288) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89895407) q[0];
sx q[0];
rz(-0.88976088) q[0];
sx q[0];
rz(-0.97578543) q[0];
rz(-1.1294533) q[1];
sx q[1];
rz(-2.6934846) q[1];
sx q[1];
rz(1.7879558) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6393051) q[0];
sx q[0];
rz(-0.76984084) q[0];
sx q[0];
rz(-0.42686157) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3951621) q[2];
sx q[2];
rz(-0.43000107) q[2];
sx q[2];
rz(1.912009) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5013988) q[1];
sx q[1];
rz(-0.36194776) q[1];
sx q[1];
rz(0.37935235) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6725572) q[3];
sx q[3];
rz(-1.1877265) q[3];
sx q[3];
rz(-1.9214326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5740616) q[2];
sx q[2];
rz(-0.83063829) q[2];
sx q[2];
rz(2.2243824) q[2];
rz(-1.5926825) q[3];
sx q[3];
rz(-2.8062688) q[3];
sx q[3];
rz(-2.3622021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(0.2061283) q[0];
sx q[0];
rz(-1.6213106) q[0];
sx q[0];
rz(3.1167378) q[0];
rz(-1.8553597) q[1];
sx q[1];
rz(-1.356946) q[1];
sx q[1];
rz(1.53055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0162189) q[0];
sx q[0];
rz(-1.9589931) q[0];
sx q[0];
rz(2.6632705) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8288245) q[2];
sx q[2];
rz(-0.90094968) q[2];
sx q[2];
rz(-0.19425288) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3250585) q[1];
sx q[1];
rz(-1.6055425) q[1];
sx q[1];
rz(2.860022) q[1];
rz(-pi) q[2];
x q[2];
rz(0.99175158) q[3];
sx q[3];
rz(-2.8914521) q[3];
sx q[3];
rz(-2.5219805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1207235) q[2];
sx q[2];
rz(-2.1376231) q[2];
sx q[2];
rz(-1.2845767) q[2];
rz(-0.32127109) q[3];
sx q[3];
rz(-2.4925241) q[3];
sx q[3];
rz(0.2215213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.016841737) q[0];
sx q[0];
rz(-1.1469954) q[0];
sx q[0];
rz(2.2220213) q[0];
rz(-0.59182566) q[1];
sx q[1];
rz(-2.070919) q[1];
sx q[1];
rz(-0.33669534) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3172835) q[0];
sx q[0];
rz(-1.705184) q[0];
sx q[0];
rz(-1.8407673) q[0];
rz(-1.4496964) q[2];
sx q[2];
rz(-1.9316439) q[2];
sx q[2];
rz(-2.0845513) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7408161) q[1];
sx q[1];
rz(-2.0753521) q[1];
sx q[1];
rz(1.7187121) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3161427) q[3];
sx q[3];
rz(-1.4166204) q[3];
sx q[3];
rz(-1.7642989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1153974) q[2];
sx q[2];
rz(-2.7112609) q[2];
sx q[2];
rz(2.3587312) q[2];
rz(0.10910263) q[3];
sx q[3];
rz(-1.0166054) q[3];
sx q[3];
rz(1.8946064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03610177) q[0];
sx q[0];
rz(-1.4694659) q[0];
sx q[0];
rz(0.92392695) q[0];
rz(0.52664122) q[1];
sx q[1];
rz(-1.8662607) q[1];
sx q[1];
rz(0.99871666) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.625542) q[0];
sx q[0];
rz(-1.417451) q[0];
sx q[0];
rz(-1.6244066) q[0];
rz(-pi) q[1];
rz(1.8417712) q[2];
sx q[2];
rz(-0.30108115) q[2];
sx q[2];
rz(-1.4512514) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1371138) q[1];
sx q[1];
rz(-2.5527337) q[1];
sx q[1];
rz(-1.5084672) q[1];
rz(2.6941006) q[3];
sx q[3];
rz(-0.72967859) q[3];
sx q[3];
rz(3.0270769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1363137) q[2];
sx q[2];
rz(-0.88999358) q[2];
sx q[2];
rz(0.33373731) q[2];
rz(-2.2809095) q[3];
sx q[3];
rz(-1.4849097) q[3];
sx q[3];
rz(0.0066283289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-0.46398791) q[0];
sx q[0];
rz(-1.1900359) q[0];
sx q[0];
rz(0.66587454) q[0];
rz(2.5869276) q[1];
sx q[1];
rz(-1.278109) q[1];
sx q[1];
rz(0.14229933) q[1];
rz(-1.3210422) q[2];
sx q[2];
rz(-1.142923) q[2];
sx q[2];
rz(-0.10071071) q[2];
rz(-2.5139204) q[3];
sx q[3];
rz(-2.0439507) q[3];
sx q[3];
rz(2.4491257) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
