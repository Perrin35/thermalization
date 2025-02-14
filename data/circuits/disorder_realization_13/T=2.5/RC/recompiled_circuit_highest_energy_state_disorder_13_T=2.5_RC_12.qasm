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
rz(0.27999347) q[0];
sx q[0];
rz(-1.0323098) q[0];
sx q[0];
rz(-1.1601467) q[0];
rz(-0.87814826) q[1];
sx q[1];
rz(-2.703009) q[1];
sx q[1];
rz(2.261472) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11186799) q[0];
sx q[0];
rz(-0.97767413) q[0];
sx q[0];
rz(0.2082227) q[0];
rz(0.080332412) q[2];
sx q[2];
rz(-1.1690306) q[2];
sx q[2];
rz(0.062594819) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.82950006) q[1];
sx q[1];
rz(-2.1662427) q[1];
sx q[1];
rz(2.9445346) q[1];
x q[2];
rz(2.6875917) q[3];
sx q[3];
rz(-1.6780644) q[3];
sx q[3];
rz(1.1196746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1151513) q[2];
sx q[2];
rz(-0.5618962) q[2];
sx q[2];
rz(1.119841) q[2];
rz(-1.0037054) q[3];
sx q[3];
rz(-0.34990889) q[3];
sx q[3];
rz(2.2899535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72964662) q[0];
sx q[0];
rz(-0.8929407) q[0];
sx q[0];
rz(0.44806421) q[0];
rz(-2.0662775) q[1];
sx q[1];
rz(-1.0766462) q[1];
sx q[1];
rz(0.040464673) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0300022) q[0];
sx q[0];
rz(-0.69911041) q[0];
sx q[0];
rz(1.2707236) q[0];
rz(-pi) q[1];
rz(-1.984111) q[2];
sx q[2];
rz(-1.4885474) q[2];
sx q[2];
rz(0.27378191) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.81442388) q[1];
sx q[1];
rz(-1.8393687) q[1];
sx q[1];
rz(-1.5280523) q[1];
x q[2];
rz(1.5621885) q[3];
sx q[3];
rz(-1.2581902) q[3];
sx q[3];
rz(1.7695939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2761817) q[2];
sx q[2];
rz(-2.2409596) q[2];
sx q[2];
rz(-0.5033699) q[2];
rz(1.447575) q[3];
sx q[3];
rz(-1.9083128) q[3];
sx q[3];
rz(-0.90731049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68509787) q[0];
sx q[0];
rz(-0.95360294) q[0];
sx q[0];
rz(2.8431235) q[0];
rz(-1.586277) q[1];
sx q[1];
rz(-1.7319873) q[1];
sx q[1];
rz(1.5063937) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28373805) q[0];
sx q[0];
rz(-1.2728134) q[0];
sx q[0];
rz(1.4922754) q[0];
rz(1.9048018) q[2];
sx q[2];
rz(-2.4270505) q[2];
sx q[2];
rz(-1.3151907) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3909482) q[1];
sx q[1];
rz(-1.7034638) q[1];
sx q[1];
rz(-0.68718095) q[1];
x q[2];
rz(-1.7603586) q[3];
sx q[3];
rz(-1.624057) q[3];
sx q[3];
rz(1.9587751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7895268) q[2];
sx q[2];
rz(-2.3730998) q[2];
sx q[2];
rz(-2.3699769) q[2];
rz(1.4113034) q[3];
sx q[3];
rz(-1.227042) q[3];
sx q[3];
rz(-0.82967776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70374933) q[0];
sx q[0];
rz(-0.65173906) q[0];
sx q[0];
rz(-1.334345) q[0];
rz(1.5127381) q[1];
sx q[1];
rz(-1.6213497) q[1];
sx q[1];
rz(-3.0049862) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.893394) q[0];
sx q[0];
rz(-1.2635487) q[0];
sx q[0];
rz(0.91715468) q[0];
rz(-pi) q[1];
rz(1.2567472) q[2];
sx q[2];
rz(-0.32448353) q[2];
sx q[2];
rz(2.9195234) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.97164153) q[1];
sx q[1];
rz(-1.2593237) q[1];
sx q[1];
rz(-2.4258652) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9196687) q[3];
sx q[3];
rz(-1.2563946) q[3];
sx q[3];
rz(-1.507198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2887349) q[2];
sx q[2];
rz(-1.4790269) q[2];
sx q[2];
rz(0.23957254) q[2];
rz(2.1824956) q[3];
sx q[3];
rz(-1.9775016) q[3];
sx q[3];
rz(-1.0378999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.1921444) q[0];
sx q[0];
rz(-2.0862155) q[0];
sx q[0];
rz(1.5978093) q[0];
rz(2.0315309) q[1];
sx q[1];
rz(-1.9381356) q[1];
sx q[1];
rz(-1.69453) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.869782) q[0];
sx q[0];
rz(-0.29848924) q[0];
sx q[0];
rz(-3.117352) q[0];
rz(3.1008312) q[2];
sx q[2];
rz(-1.9136179) q[2];
sx q[2];
rz(1.1072323) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2867109) q[1];
sx q[1];
rz(-2.4457702) q[1];
sx q[1];
rz(2.8665101) q[1];
x q[2];
rz(2.5041833) q[3];
sx q[3];
rz(-2.4968409) q[3];
sx q[3];
rz(1.7774323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1453104) q[2];
sx q[2];
rz(-2.0029533) q[2];
sx q[2];
rz(-0.81926695) q[2];
rz(0.94775689) q[3];
sx q[3];
rz(-0.78815931) q[3];
sx q[3];
rz(1.88571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6325697) q[0];
sx q[0];
rz(-0.13795723) q[0];
sx q[0];
rz(0.55009681) q[0];
rz(-0.24712786) q[1];
sx q[1];
rz(-1.988966) q[1];
sx q[1];
rz(-2.1955042) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2193748) q[0];
sx q[0];
rz(-1.3117466) q[0];
sx q[0];
rz(-1.9402615) q[0];
rz(-pi) q[1];
rz(1.3576415) q[2];
sx q[2];
rz(-1.5098796) q[2];
sx q[2];
rz(-0.83535458) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.94080776) q[1];
sx q[1];
rz(-1.1012804) q[1];
sx q[1];
rz(1.5307309) q[1];
rz(1.4660353) q[3];
sx q[3];
rz(-1.8083866) q[3];
sx q[3];
rz(1.6692357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.75521022) q[2];
sx q[2];
rz(-2.5911665) q[2];
sx q[2];
rz(0.68106252) q[2];
rz(0.81124535) q[3];
sx q[3];
rz(-1.0439876) q[3];
sx q[3];
rz(-0.55027858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022245971) q[0];
sx q[0];
rz(-2.5211054) q[0];
sx q[0];
rz(2.3959809) q[0];
rz(1.9781808) q[1];
sx q[1];
rz(-2.5201576) q[1];
sx q[1];
rz(2.9643639) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7435869) q[0];
sx q[0];
rz(-1.2400125) q[0];
sx q[0];
rz(-0.28689204) q[0];
x q[1];
rz(0.96259768) q[2];
sx q[2];
rz(-0.11044914) q[2];
sx q[2];
rz(-0.77753528) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2043269) q[1];
sx q[1];
rz(-2.4961053) q[1];
sx q[1];
rz(-0.93475229) q[1];
x q[2];
rz(-0.69474649) q[3];
sx q[3];
rz(-1.5383895) q[3];
sx q[3];
rz(-0.8424527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.54997286) q[2];
sx q[2];
rz(-2.1592996) q[2];
sx q[2];
rz(1.4493235) q[2];
rz(2.5189404) q[3];
sx q[3];
rz(-0.52716523) q[3];
sx q[3];
rz(-1.3194293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1249238) q[0];
sx q[0];
rz(-1.387384) q[0];
sx q[0];
rz(1.5119875) q[0];
rz(-0.92388693) q[1];
sx q[1];
rz(-1.4199363) q[1];
sx q[1];
rz(0.61663827) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4409326) q[0];
sx q[0];
rz(-2.1671139) q[0];
sx q[0];
rz(0.18429549) q[0];
rz(-3.055279) q[2];
sx q[2];
rz(-2.6722135) q[2];
sx q[2];
rz(-2.1522107) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.60395998) q[1];
sx q[1];
rz(-2.5164971) q[1];
sx q[1];
rz(2.9949466) q[1];
rz(-pi) q[2];
x q[2];
rz(2.736896) q[3];
sx q[3];
rz(-0.62129489) q[3];
sx q[3];
rz(0.23970793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7743249) q[2];
sx q[2];
rz(-2.7820945) q[2];
sx q[2];
rz(-0.16505879) q[2];
rz(0.74602357) q[3];
sx q[3];
rz(-1.6227928) q[3];
sx q[3];
rz(0.042796854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54166722) q[0];
sx q[0];
rz(-0.2921108) q[0];
sx q[0];
rz(-2.7274729) q[0];
rz(0.43553022) q[1];
sx q[1];
rz(-0.51594096) q[1];
sx q[1];
rz(1.863265) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1762487) q[0];
sx q[0];
rz(-2.8581736) q[0];
sx q[0];
rz(1.0229646) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3971073) q[2];
sx q[2];
rz(-1.2635482) q[2];
sx q[2];
rz(-1.9510614) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.49474459) q[1];
sx q[1];
rz(-0.068915135) q[1];
sx q[1];
rz(-2.9685741) q[1];
rz(-pi) q[2];
rz(-0.06287136) q[3];
sx q[3];
rz(-2.0977763) q[3];
sx q[3];
rz(1.8450027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4182959) q[2];
sx q[2];
rz(-1.7593242) q[2];
sx q[2];
rz(0.74271512) q[2];
rz(-2.3501979) q[3];
sx q[3];
rz(-0.68358889) q[3];
sx q[3];
rz(-1.7323823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0940229) q[0];
sx q[0];
rz(-2.8122734) q[0];
sx q[0];
rz(-0.91162115) q[0];
rz(-2.1564663) q[1];
sx q[1];
rz(-0.42675012) q[1];
sx q[1];
rz(1.3381348) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1250004) q[0];
sx q[0];
rz(-2.5574698) q[0];
sx q[0];
rz(0.97726269) q[0];
rz(1.192044) q[2];
sx q[2];
rz(-1.2745924) q[2];
sx q[2];
rz(-0.51329324) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9836925) q[1];
sx q[1];
rz(-2.165598) q[1];
sx q[1];
rz(2.6524794) q[1];
x q[2];
rz(1.9444405) q[3];
sx q[3];
rz(-1.257953) q[3];
sx q[3];
rz(0.88176239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7868339) q[2];
sx q[2];
rz(-1.6666731) q[2];
sx q[2];
rz(1.1653398) q[2];
rz(1.3960086) q[3];
sx q[3];
rz(-1.952012) q[3];
sx q[3];
rz(-1.5497807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6373445) q[0];
sx q[0];
rz(-1.5550384) q[0];
sx q[0];
rz(0.73873781) q[0];
rz(2.5042116) q[1];
sx q[1];
rz(-1.5370054) q[1];
sx q[1];
rz(-0.76269033) q[1];
rz(1.2554854) q[2];
sx q[2];
rz(-0.44305319) q[2];
sx q[2];
rz(-2.6534266) q[2];
rz(0.55841726) q[3];
sx q[3];
rz(-1.9600086) q[3];
sx q[3];
rz(0.10127898) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
