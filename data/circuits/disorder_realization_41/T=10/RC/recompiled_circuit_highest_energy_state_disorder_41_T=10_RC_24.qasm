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
rz(0.32349411) q[0];
sx q[0];
rz(-0.20322023) q[0];
sx q[0];
rz(-2.8414371) q[0];
rz(2.7872941) q[1];
sx q[1];
rz(-1.2343255) q[1];
sx q[1];
rz(-0.77114463) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3572039) q[0];
sx q[0];
rz(-1.3805806) q[0];
sx q[0];
rz(1.3804382) q[0];
x q[1];
rz(-2.0448737) q[2];
sx q[2];
rz(-1.0869622) q[2];
sx q[2];
rz(-2.3729618) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.31509545) q[1];
sx q[1];
rz(-0.56964785) q[1];
sx q[1];
rz(0.82514067) q[1];
rz(-pi) q[2];
x q[2];
rz(0.067581108) q[3];
sx q[3];
rz(-1.1376808) q[3];
sx q[3];
rz(-1.4560631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99042088) q[2];
sx q[2];
rz(-0.46151084) q[2];
sx q[2];
rz(2.0578461) q[2];
rz(-0.57461965) q[3];
sx q[3];
rz(-1.5465522) q[3];
sx q[3];
rz(2.3641724) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3144749) q[0];
sx q[0];
rz(-2.5125393) q[0];
sx q[0];
rz(2.7431059) q[0];
rz(1.9972948) q[1];
sx q[1];
rz(-0.60595787) q[1];
sx q[1];
rz(0.95432895) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7643549) q[0];
sx q[0];
rz(-1.5713552) q[0];
sx q[0];
rz(3.1370381) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2284438) q[2];
sx q[2];
rz(-1.4962915) q[2];
sx q[2];
rz(-2.7721921) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.40007517) q[1];
sx q[1];
rz(-1.0055826) q[1];
sx q[1];
rz(2.8533622) q[1];
rz(-pi) q[2];
rz(-0.54150906) q[3];
sx q[3];
rz(-1.9810105) q[3];
sx q[3];
rz(2.4334986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2362471) q[2];
sx q[2];
rz(-0.50062847) q[2];
sx q[2];
rz(-0.11032571) q[2];
rz(-2.3162383) q[3];
sx q[3];
rz(-2.3767411) q[3];
sx q[3];
rz(0.74025214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9731307) q[0];
sx q[0];
rz(-0.2178807) q[0];
sx q[0];
rz(0.88200283) q[0];
rz(1.0285671) q[1];
sx q[1];
rz(-0.55362892) q[1];
sx q[1];
rz(2.2291768) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19174978) q[0];
sx q[0];
rz(-1.1198567) q[0];
sx q[0];
rz(-0.080206932) q[0];
x q[1];
rz(1.6912314) q[2];
sx q[2];
rz(-2.2856973) q[2];
sx q[2];
rz(3.1034865) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6710224) q[1];
sx q[1];
rz(-2.467944) q[1];
sx q[1];
rz(1.0224875) q[1];
rz(-1.0718143) q[3];
sx q[3];
rz(-2.6423965) q[3];
sx q[3];
rz(0.62455356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2622111) q[2];
sx q[2];
rz(-0.24803455) q[2];
sx q[2];
rz(-0.81825078) q[2];
rz(-2.7815871) q[3];
sx q[3];
rz(-1.8428948) q[3];
sx q[3];
rz(0.037809614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70571947) q[0];
sx q[0];
rz(-2.191045) q[0];
sx q[0];
rz(-1.1608423) q[0];
rz(2.1707161) q[1];
sx q[1];
rz(-2.2258046) q[1];
sx q[1];
rz(-0.11347778) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16100854) q[0];
sx q[0];
rz(-2.1652114) q[0];
sx q[0];
rz(-1.0779146) q[0];
rz(1.214049) q[2];
sx q[2];
rz(-0.037790701) q[2];
sx q[2];
rz(0.99810696) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.66333713) q[1];
sx q[1];
rz(-2.2436444) q[1];
sx q[1];
rz(0.70391432) q[1];
rz(-pi) q[2];
rz(2.5069019) q[3];
sx q[3];
rz(-1.029976) q[3];
sx q[3];
rz(1.6010923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.11161441) q[2];
sx q[2];
rz(-1.8054211) q[2];
sx q[2];
rz(0.90799904) q[2];
rz(-2.8152483) q[3];
sx q[3];
rz(-2.58674) q[3];
sx q[3];
rz(-2.7197796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57080764) q[0];
sx q[0];
rz(-0.76199216) q[0];
sx q[0];
rz(2.1245891) q[0];
rz(2.6550338) q[1];
sx q[1];
rz(-2.1162972) q[1];
sx q[1];
rz(3.0466381) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5240898) q[0];
sx q[0];
rz(-1.4727778) q[0];
sx q[0];
rz(-1.6785511) q[0];
rz(-pi) q[1];
rz(-0.20304154) q[2];
sx q[2];
rz(-0.83492324) q[2];
sx q[2];
rz(-1.0395713) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6330715) q[1];
sx q[1];
rz(-1.2969368) q[1];
sx q[1];
rz(-2.7809596) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8311555) q[3];
sx q[3];
rz(-1.2708029) q[3];
sx q[3];
rz(2.7315549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93593705) q[2];
sx q[2];
rz(-2.50694) q[2];
sx q[2];
rz(2.671799) q[2];
rz(-0.69822407) q[3];
sx q[3];
rz(-1.2263115) q[3];
sx q[3];
rz(-0.32018143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6093589) q[0];
sx q[0];
rz(-2.4025752) q[0];
sx q[0];
rz(-0.21482378) q[0];
rz(2.9048982) q[1];
sx q[1];
rz(-0.57699811) q[1];
sx q[1];
rz(-0.41928852) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.740991) q[0];
sx q[0];
rz(-3.0508578) q[0];
sx q[0];
rz(-2.0965133) q[0];
rz(1.0299267) q[2];
sx q[2];
rz(-1.3020143) q[2];
sx q[2];
rz(-0.93712419) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9671075) q[1];
sx q[1];
rz(-1.9000783) q[1];
sx q[1];
rz(2.9261173) q[1];
rz(-pi) q[2];
rz(-1.7691408) q[3];
sx q[3];
rz(-0.86776483) q[3];
sx q[3];
rz(-1.8198215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.99570167) q[2];
sx q[2];
rz(-3.006101) q[2];
sx q[2];
rz(-3.0502012) q[2];
rz(-0.35178301) q[3];
sx q[3];
rz(-0.7472977) q[3];
sx q[3];
rz(-0.465213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9212937) q[0];
sx q[0];
rz(-1.969368) q[0];
sx q[0];
rz(1.9660796) q[0];
rz(-2.562404) q[1];
sx q[1];
rz(-2.3349031) q[1];
sx q[1];
rz(0.18956345) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80081683) q[0];
sx q[0];
rz(-2.918266) q[0];
sx q[0];
rz(2.6912265) q[0];
x q[1];
rz(3.0547943) q[2];
sx q[2];
rz(-2.2703836) q[2];
sx q[2];
rz(-1.1530641) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.8622562) q[1];
sx q[1];
rz(-1.5631884) q[1];
sx q[1];
rz(2.3235282) q[1];
rz(-0.060154288) q[3];
sx q[3];
rz(-1.0511479) q[3];
sx q[3];
rz(2.1148256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.2438948) q[2];
sx q[2];
rz(-1.2995259) q[2];
sx q[2];
rz(0.50590903) q[2];
rz(-0.77887744) q[3];
sx q[3];
rz(-0.39611045) q[3];
sx q[3];
rz(-1.8469384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27095723) q[0];
sx q[0];
rz(-1.0519692) q[0];
sx q[0];
rz(1.9860995) q[0];
rz(0.0011860154) q[1];
sx q[1];
rz(-2.3540034) q[1];
sx q[1];
rz(-0.40518618) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7195994) q[0];
sx q[0];
rz(-2.036839) q[0];
sx q[0];
rz(2.5525722) q[0];
rz(-pi) q[1];
rz(1.1616906) q[2];
sx q[2];
rz(-1.4295345) q[2];
sx q[2];
rz(-2.2337332) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.59711528) q[1];
sx q[1];
rz(-1.6469203) q[1];
sx q[1];
rz(-0.80634574) q[1];
rz(-pi) q[2];
rz(-1.5022329) q[3];
sx q[3];
rz(-1.6985095) q[3];
sx q[3];
rz(2.8721916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.323641) q[2];
sx q[2];
rz(-1.9743071) q[2];
sx q[2];
rz(2.6200068) q[2];
rz(-0.096605435) q[3];
sx q[3];
rz(-2.1229027) q[3];
sx q[3];
rz(0.49083138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7535962) q[0];
sx q[0];
rz(-0.8803519) q[0];
sx q[0];
rz(-3.0269347) q[0];
rz(-1.8278587) q[1];
sx q[1];
rz(-1.6856245) q[1];
sx q[1];
rz(0.1350666) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11896597) q[0];
sx q[0];
rz(-2.8014136) q[0];
sx q[0];
rz(2.4426798) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58456771) q[2];
sx q[2];
rz(-0.46007878) q[2];
sx q[2];
rz(2.8321617) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5852063) q[1];
sx q[1];
rz(-1.7328246) q[1];
sx q[1];
rz(2.1955447) q[1];
rz(-pi) q[2];
x q[2];
rz(0.073487893) q[3];
sx q[3];
rz(-2.1284496) q[3];
sx q[3];
rz(2.776655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.31585285) q[2];
sx q[2];
rz(-1.2154546) q[2];
sx q[2];
rz(2.1960171) q[2];
rz(-2.8816176) q[3];
sx q[3];
rz(-2.1187449) q[3];
sx q[3];
rz(-2.0135349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.05397) q[0];
sx q[0];
rz(-2.051351) q[0];
sx q[0];
rz(1.7386275) q[0];
rz(-2.2258017) q[1];
sx q[1];
rz(-2.5251838) q[1];
sx q[1];
rz(0.91555196) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85957324) q[0];
sx q[0];
rz(-1.7012736) q[0];
sx q[0];
rz(3.0436712) q[0];
x q[1];
rz(-0.19277566) q[2];
sx q[2];
rz(-2.0844242) q[2];
sx q[2];
rz(2.7889268) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.97595313) q[1];
sx q[1];
rz(-2.0054066) q[1];
sx q[1];
rz(2.6967816) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7677018) q[3];
sx q[3];
rz(-0.71929661) q[3];
sx q[3];
rz(-0.48433892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0182858) q[2];
sx q[2];
rz(-2.4985963) q[2];
sx q[2];
rz(-2.7389738) q[2];
rz(-1.1766524) q[3];
sx q[3];
rz(-2.8549356) q[3];
sx q[3];
rz(-2.3736242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1422414) q[0];
sx q[0];
rz(-0.9724697) q[0];
sx q[0];
rz(-1.019626) q[0];
rz(2.453852) q[1];
sx q[1];
rz(-2.3383457) q[1];
sx q[1];
rz(1.8170423) q[1];
rz(0.68647142) q[2];
sx q[2];
rz(-2.1778637) q[2];
sx q[2];
rz(-1.7354497) q[2];
rz(-1.8263893) q[3];
sx q[3];
rz(-0.9556613) q[3];
sx q[3];
rz(-1.2355604) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
