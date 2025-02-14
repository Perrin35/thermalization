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
rz(0.723847) q[0];
sx q[0];
rz(-1.2011733) q[0];
sx q[0];
rz(0.49402753) q[0];
rz(1.5552893) q[1];
sx q[1];
rz(2.2145693) q[1];
sx q[1];
rz(8.5742843) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5828398) q[0];
sx q[0];
rz(-0.90769288) q[0];
sx q[0];
rz(2.9302674) q[0];
x q[1];
rz(-0.57588864) q[2];
sx q[2];
rz(-0.73969361) q[2];
sx q[2];
rz(2.0832555) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.32341084) q[1];
sx q[1];
rz(-1.2480624) q[1];
sx q[1];
rz(1.3363911) q[1];
rz(-pi) q[2];
rz(-2.9539724) q[3];
sx q[3];
rz(-2.7696262) q[3];
sx q[3];
rz(-0.99727977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6583307) q[2];
sx q[2];
rz(-0.61507812) q[2];
sx q[2];
rz(0.94493803) q[2];
rz(-3.1350709) q[3];
sx q[3];
rz(-2.3801453) q[3];
sx q[3];
rz(0.25750461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1061123) q[0];
sx q[0];
rz(-2.1529614) q[0];
sx q[0];
rz(0.60580564) q[0];
rz(2.2094191) q[1];
sx q[1];
rz(-1.4235539) q[1];
sx q[1];
rz(0.1056284) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8307997) q[0];
sx q[0];
rz(-0.89185878) q[0];
sx q[0];
rz(1.3722769) q[0];
rz(-pi) q[1];
rz(-1.0829686) q[2];
sx q[2];
rz(-1.8886856) q[2];
sx q[2];
rz(-2.6018104) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.82984207) q[1];
sx q[1];
rz(-2.0220322) q[1];
sx q[1];
rz(-1.2430906) q[1];
x q[2];
rz(1.9808045) q[3];
sx q[3];
rz(-1.2224397) q[3];
sx q[3];
rz(0.94545555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2848844) q[2];
sx q[2];
rz(-1.363089) q[2];
sx q[2];
rz(-1.523783) q[2];
rz(2.3273322) q[3];
sx q[3];
rz(-1.7819449) q[3];
sx q[3];
rz(0.8849357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098323671) q[0];
sx q[0];
rz(-0.79541484) q[0];
sx q[0];
rz(1.5650308) q[0];
rz(-0.99524975) q[1];
sx q[1];
rz(-2.1627656) q[1];
sx q[1];
rz(-0.78688041) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7835307) q[0];
sx q[0];
rz(-1.770505) q[0];
sx q[0];
rz(-3.003503) q[0];
rz(-pi) q[1];
rz(-2.9126105) q[2];
sx q[2];
rz(-0.35576421) q[2];
sx q[2];
rz(-2.7965429) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7461024) q[1];
sx q[1];
rz(-1.8999892) q[1];
sx q[1];
rz(1.9026865) q[1];
x q[2];
rz(-1.6305109) q[3];
sx q[3];
rz(-1.51126) q[3];
sx q[3];
rz(-1.2950031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3028822) q[2];
sx q[2];
rz(-1.6532712) q[2];
sx q[2];
rz(-2.5035456) q[2];
rz(-2.6436515) q[3];
sx q[3];
rz(-2.0697856) q[3];
sx q[3];
rz(1.0897442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8968673) q[0];
sx q[0];
rz(-1.3571955) q[0];
sx q[0];
rz(-0.063752739) q[0];
rz(-1.4490734) q[1];
sx q[1];
rz(-1.3056825) q[1];
sx q[1];
rz(-1.3099028) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5638014) q[0];
sx q[0];
rz(-0.088106958) q[0];
sx q[0];
rz(-2.7035575) q[0];
rz(1.7559577) q[2];
sx q[2];
rz(-1.3755535) q[2];
sx q[2];
rz(2.756292) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3514362) q[1];
sx q[1];
rz(-0.75439851) q[1];
sx q[1];
rz(-0.18138563) q[1];
x q[2];
rz(2.784666) q[3];
sx q[3];
rz(-0.91536575) q[3];
sx q[3];
rz(-0.41567685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0706851) q[2];
sx q[2];
rz(-2.0347774) q[2];
sx q[2];
rz(-3.0359388) q[2];
rz(-1.9581155) q[3];
sx q[3];
rz(-2.2453997) q[3];
sx q[3];
rz(0.67160523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3463109) q[0];
sx q[0];
rz(-0.91049894) q[0];
sx q[0];
rz(-1.3443391) q[0];
rz(0.67289871) q[1];
sx q[1];
rz(-1.8310603) q[1];
sx q[1];
rz(0.013462822) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1707662) q[0];
sx q[0];
rz(-1.5970236) q[0];
sx q[0];
rz(-0.88084014) q[0];
x q[1];
rz(-1.5384664) q[2];
sx q[2];
rz(-2.7271977) q[2];
sx q[2];
rz(-2.3737645) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8103761) q[1];
sx q[1];
rz(-1.5575641) q[1];
sx q[1];
rz(-0.41535901) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0125403) q[3];
sx q[3];
rz(-0.75296445) q[3];
sx q[3];
rz(-0.66679614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.54231918) q[2];
sx q[2];
rz(-2.4794674) q[2];
sx q[2];
rz(-1.2467965) q[2];
rz(-0.072619297) q[3];
sx q[3];
rz(-0.19425546) q[3];
sx q[3];
rz(-0.3092002) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4319864) q[0];
sx q[0];
rz(-1.1259587) q[0];
sx q[0];
rz(-3.0288938) q[0];
rz(0.20507774) q[1];
sx q[1];
rz(-1.3321184) q[1];
sx q[1];
rz(-2.233706) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4604707) q[0];
sx q[0];
rz(-2.2325071) q[0];
sx q[0];
rz(0.76028334) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.142318) q[2];
sx q[2];
rz(-2.1399763) q[2];
sx q[2];
rz(2.6250397) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1334539) q[1];
sx q[1];
rz(-0.63888237) q[1];
sx q[1];
rz(2.9879346) q[1];
x q[2];
rz(-1.3566293) q[3];
sx q[3];
rz(-2.2244144) q[3];
sx q[3];
rz(0.71398338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9997361) q[2];
sx q[2];
rz(-0.33679589) q[2];
sx q[2];
rz(3.0736308) q[2];
rz(-1.147602) q[3];
sx q[3];
rz(-1.2679029) q[3];
sx q[3];
rz(-1.2609153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9481908) q[0];
sx q[0];
rz(-1.8302487) q[0];
sx q[0];
rz(-0.46352682) q[0];
rz(2.9947128) q[1];
sx q[1];
rz(-1.2356707) q[1];
sx q[1];
rz(-2.1489977) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50336058) q[0];
sx q[0];
rz(-1.8192756) q[0];
sx q[0];
rz(-0.72297024) q[0];
rz(-pi) q[1];
rz(-0.067258283) q[2];
sx q[2];
rz(-1.2853299) q[2];
sx q[2];
rz(2.26254) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2414244) q[1];
sx q[1];
rz(-1.1196152) q[1];
sx q[1];
rz(-2.3996949) q[1];
x q[2];
rz(2.4086508) q[3];
sx q[3];
rz(-1.2228726) q[3];
sx q[3];
rz(0.80845736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85840449) q[2];
sx q[2];
rz(-1.4294727) q[2];
sx q[2];
rz(-0.44537133) q[2];
rz(-1.8177659) q[3];
sx q[3];
rz(-2.0711074) q[3];
sx q[3];
rz(2.0557192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0367947) q[0];
sx q[0];
rz(-0.0093655149) q[0];
sx q[0];
rz(2.985756) q[0];
rz(1.8644631) q[1];
sx q[1];
rz(-1.7946449) q[1];
sx q[1];
rz(0.015930463) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8317922) q[0];
sx q[0];
rz(-0.81912097) q[0];
sx q[0];
rz(-2.8099634) q[0];
x q[1];
rz(0.792298) q[2];
sx q[2];
rz(-1.3092666) q[2];
sx q[2];
rz(0.025394414) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5585352) q[1];
sx q[1];
rz(-2.0777521) q[1];
sx q[1];
rz(0.81113775) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63834493) q[3];
sx q[3];
rz(-1.9192291) q[3];
sx q[3];
rz(-1.3909457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5126123) q[2];
sx q[2];
rz(-3.0292558) q[2];
sx q[2];
rz(0.12574276) q[2];
rz(0.9564774) q[3];
sx q[3];
rz(-1.4230909) q[3];
sx q[3];
rz(2.8023348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15928957) q[0];
sx q[0];
rz(-0.24022261) q[0];
sx q[0];
rz(2.6851728) q[0];
rz(1.5727111) q[1];
sx q[1];
rz(-2.1379037) q[1];
sx q[1];
rz(0.68663866) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52138162) q[0];
sx q[0];
rz(-0.98582637) q[0];
sx q[0];
rz(-2.1725031) q[0];
rz(-1.3977658) q[2];
sx q[2];
rz(-1.695172) q[2];
sx q[2];
rz(-1.591452) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7788199) q[1];
sx q[1];
rz(-1.6906926) q[1];
sx q[1];
rz(-0.1480379) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37952642) q[3];
sx q[3];
rz(-2.9538493) q[3];
sx q[3];
rz(1.2953341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.75866428) q[2];
sx q[2];
rz(-2.0186581) q[2];
sx q[2];
rz(-2.0056966) q[2];
rz(-2.1618333) q[3];
sx q[3];
rz(-1.3330678) q[3];
sx q[3];
rz(-0.69825828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45282388) q[0];
sx q[0];
rz(-1.8166421) q[0];
sx q[0];
rz(0.45595566) q[0];
rz(-3.0988354) q[1];
sx q[1];
rz(-1.9330934) q[1];
sx q[1];
rz(2.4483689) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42825481) q[0];
sx q[0];
rz(-1.9549668) q[0];
sx q[0];
rz(-0.1066271) q[0];
rz(-pi) q[1];
rz(-2.736959) q[2];
sx q[2];
rz(-1.1961968) q[2];
sx q[2];
rz(2.7092779) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1566504) q[1];
sx q[1];
rz(-1.4727815) q[1];
sx q[1];
rz(-3.1344862) q[1];
x q[2];
rz(-0.59488036) q[3];
sx q[3];
rz(-2.6118546) q[3];
sx q[3];
rz(1.4929508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.91528714) q[2];
sx q[2];
rz(-1.8733571) q[2];
sx q[2];
rz(-2.477296) q[2];
rz(1.213446) q[3];
sx q[3];
rz(-0.72187859) q[3];
sx q[3];
rz(-0.87219316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6912457) q[0];
sx q[0];
rz(-0.84083122) q[0];
sx q[0];
rz(-1.3837411) q[0];
rz(1.4585523) q[1];
sx q[1];
rz(-0.28266193) q[1];
sx q[1];
rz(-2.2255486) q[1];
rz(3.0677879) q[2];
sx q[2];
rz(-2.7701785) q[2];
sx q[2];
rz(-2.2899173) q[2];
rz(-0.70590677) q[3];
sx q[3];
rz(-1.1242031) q[3];
sx q[3];
rz(-0.28225337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
