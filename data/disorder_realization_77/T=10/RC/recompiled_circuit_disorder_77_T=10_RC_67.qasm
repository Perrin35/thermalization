OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(0.16790976) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(-0.056161031) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1355609) q[0];
sx q[0];
rz(-2.6132085) q[0];
sx q[0];
rz(-1.1392659) q[0];
x q[1];
rz(2.216823) q[2];
sx q[2];
rz(-1.7416818) q[2];
sx q[2];
rz(0.31121635) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9819298) q[1];
sx q[1];
rz(-0.61239457) q[1];
sx q[1];
rz(1.0931404) q[1];
rz(2.8817301) q[3];
sx q[3];
rz(-1.6136323) q[3];
sx q[3];
rz(0.44997893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7636259) q[2];
sx q[2];
rz(-2.8597735) q[2];
sx q[2];
rz(-0.4326694) q[2];
rz(1.9487322) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(-2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6137961) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(1.312785) q[0];
rz(0.20547543) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(1.9899433) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8173556) q[0];
sx q[0];
rz(-1.2756057) q[0];
sx q[0];
rz(2.2833707) q[0];
x q[1];
rz(-2.4685681) q[2];
sx q[2];
rz(-0.7012127) q[2];
sx q[2];
rz(0.53403026) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9559905) q[1];
sx q[1];
rz(-2.5163109) q[1];
sx q[1];
rz(-2.37466) q[1];
rz(-2.6395256) q[3];
sx q[3];
rz(-1.2289398) q[3];
sx q[3];
rz(-1.348192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1318704) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(0.22182626) q[2];
rz(2.7644073) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(-2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8310228) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(3.1047399) q[0];
rz(-0.82551461) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(-3.085014) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6037613) q[0];
sx q[0];
rz(-1.0147525) q[0];
sx q[0];
rz(-2.6105196) q[0];
rz(-pi) q[1];
rz(-0.25643202) q[2];
sx q[2];
rz(-1.4415381) q[2];
sx q[2];
rz(1.3916707) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6371582) q[1];
sx q[1];
rz(-1.27099) q[1];
sx q[1];
rz(-1.7791041) q[1];
rz(-pi) q[2];
rz(-2.4967381) q[3];
sx q[3];
rz(-1.2554902) q[3];
sx q[3];
rz(-0.23526084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.27292192) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(2.2154714) q[2];
rz(2.5849294) q[3];
sx q[3];
rz(-2.848048) q[3];
sx q[3];
rz(2.0986957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91519231) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(-0.74209374) q[0];
rz(2.0023951) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-2.6779968) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47953654) q[0];
sx q[0];
rz(-0.89131309) q[0];
sx q[0];
rz(2.7583073) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9651627) q[2];
sx q[2];
rz(-0.33761218) q[2];
sx q[2];
rz(-0.38844973) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5220118) q[1];
sx q[1];
rz(-1.3995692) q[1];
sx q[1];
rz(-2.3492858) q[1];
x q[2];
rz(-0.47834088) q[3];
sx q[3];
rz(-1.0676749) q[3];
sx q[3];
rz(-0.92418811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3670369) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(-3.0920933) q[2];
rz(-0.1285304) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(0.11894225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054935) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(2.8438925) q[0];
rz(2.659335) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(0.94435) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5807242) q[0];
sx q[0];
rz(-1.3618016) q[0];
sx q[0];
rz(-3.0955549) q[0];
x q[1];
rz(-1.1733426) q[2];
sx q[2];
rz(-2.3059418) q[2];
sx q[2];
rz(0.022692516) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7625092) q[1];
sx q[1];
rz(-0.06113872) q[1];
sx q[1];
rz(-1.6054543) q[1];
rz(-2.6258351) q[3];
sx q[3];
rz(-2.6532647) q[3];
sx q[3];
rz(-0.26667903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.258761) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(-3.0333701) q[2];
rz(-3.1392858) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(-0.32430696) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43679431) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(0.58445245) q[0];
rz(-2.2553518) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(-3.086673) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6829837) q[0];
sx q[0];
rz(-0.28124547) q[0];
sx q[0];
rz(-1.2693229) q[0];
rz(-pi) q[1];
x q[1];
rz(1.247585) q[2];
sx q[2];
rz(-1.5885457) q[2];
sx q[2];
rz(2.7780967) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1850486) q[1];
sx q[1];
rz(-1.946432) q[1];
sx q[1];
rz(0.59021414) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.92298569) q[3];
sx q[3];
rz(-2.4499948) q[3];
sx q[3];
rz(1.5036316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.75446689) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(-2.1248655) q[2];
rz(-2.5975442) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6761557) q[0];
sx q[0];
rz(-0.98452079) q[0];
sx q[0];
rz(0.28453919) q[0];
rz(-0.94447213) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(-0.91032666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7401687) q[0];
sx q[0];
rz(-1.6356042) q[0];
sx q[0];
rz(3.0868953) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.89332135) q[2];
sx q[2];
rz(-0.51572323) q[2];
sx q[2];
rz(0.90781462) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6701723) q[1];
sx q[1];
rz(-2.6176665) q[1];
sx q[1];
rz(2.7736204) q[1];
rz(0.72083731) q[3];
sx q[3];
rz(-2.0257054) q[3];
sx q[3];
rz(-0.37973675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7515144) q[2];
sx q[2];
rz(-0.063515924) q[2];
sx q[2];
rz(-0.92203036) q[2];
rz(-0.56728029) q[3];
sx q[3];
rz(-1.7276948) q[3];
sx q[3];
rz(2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(-3.0122053) q[0];
rz(-0.63240504) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(0.30050373) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26586543) q[0];
sx q[0];
rz(-0.89571307) q[0];
sx q[0];
rz(2.6595594) q[0];
rz(-pi) q[1];
rz(2.5007162) q[2];
sx q[2];
rz(-2.2705728) q[2];
sx q[2];
rz(1.7388294) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6943372) q[1];
sx q[1];
rz(-1.2303423) q[1];
sx q[1];
rz(-1.3854331) q[1];
rz(2.825533) q[3];
sx q[3];
rz(-2.3028767) q[3];
sx q[3];
rz(2.4162606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5552716) q[2];
sx q[2];
rz(-2.1466612) q[2];
sx q[2];
rz(2.3596181) q[2];
rz(2.590495) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5732116) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(-0.12776275) q[0];
rz(-2.5993775) q[1];
sx q[1];
rz(-0.95710373) q[1];
sx q[1];
rz(2.382747) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42620537) q[0];
sx q[0];
rz(-1.1466768) q[0];
sx q[0];
rz(0.018652648) q[0];
x q[1];
rz(2.5326469) q[2];
sx q[2];
rz(-1.7296089) q[2];
sx q[2];
rz(1.2283404) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2488333) q[1];
sx q[1];
rz(-1.4151238) q[1];
sx q[1];
rz(-0.70181904) q[1];
rz(-pi) q[2];
rz(-1.5893448) q[3];
sx q[3];
rz(-2.4759001) q[3];
sx q[3];
rz(1.9513643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0163429) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(-2.6515567) q[2];
rz(-1.4222493) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816417) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(-3.066257) q[0];
rz(-0.8967337) q[1];
sx q[1];
rz(-1.9672085) q[1];
sx q[1];
rz(2.5316701) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1760575) q[0];
sx q[0];
rz(-2.3241204) q[0];
sx q[0];
rz(-2.7479991) q[0];
rz(-pi) q[1];
rz(2.5209849) q[2];
sx q[2];
rz(-0.45244103) q[2];
sx q[2];
rz(-2.0394182) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.016644195) q[1];
sx q[1];
rz(-0.84608191) q[1];
sx q[1];
rz(1.2490586) q[1];
x q[2];
rz(-0.53819733) q[3];
sx q[3];
rz(-0.81848577) q[3];
sx q[3];
rz(1.8173816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.909409) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(-2.4278736) q[2];
rz(-0.37832007) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.338035) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(2.4907885) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(1.3748319) q[2];
sx q[2];
rz(-1.6008196) q[2];
sx q[2];
rz(0.65884789) q[2];
rz(2.0416904) q[3];
sx q[3];
rz(-1.7970016) q[3];
sx q[3];
rz(-2.5947528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
