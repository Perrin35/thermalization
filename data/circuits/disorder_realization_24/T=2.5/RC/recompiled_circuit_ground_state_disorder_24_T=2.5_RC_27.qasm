OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3396575) q[0];
sx q[0];
rz(-0.8987838) q[0];
sx q[0];
rz(-1.838983) q[0];
rz(-3.0175735) q[1];
sx q[1];
rz(-1.7350585) q[1];
sx q[1];
rz(-1.5136493) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6684912) q[0];
sx q[0];
rz(-1.9051747) q[0];
sx q[0];
rz(0.3121434) q[0];
rz(-0.75998016) q[2];
sx q[2];
rz(-0.70356762) q[2];
sx q[2];
rz(-2.4576996) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3213734) q[1];
sx q[1];
rz(-2.7966683) q[1];
sx q[1];
rz(1.4089415) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.94515015) q[3];
sx q[3];
rz(-1.8813475) q[3];
sx q[3];
rz(-0.20598447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0836432) q[2];
sx q[2];
rz(-1.3336072) q[2];
sx q[2];
rz(2.7570214) q[2];
rz(-0.40258506) q[3];
sx q[3];
rz(-1.4339002) q[3];
sx q[3];
rz(-2.3334077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.758051) q[0];
sx q[0];
rz(-1.533968) q[0];
sx q[0];
rz(-2.5480399) q[0];
rz(-1.8202164) q[1];
sx q[1];
rz(-1.079419) q[1];
sx q[1];
rz(-3.1022601) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2101321) q[0];
sx q[0];
rz(-1.1173741) q[0];
sx q[0];
rz(2.7260145) q[0];
x q[1];
rz(1.6850542) q[2];
sx q[2];
rz(-2.2041177) q[2];
sx q[2];
rz(-2.3267724) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.63168282) q[1];
sx q[1];
rz(-1.0541774) q[1];
sx q[1];
rz(1.6277908) q[1];
rz(-pi) q[2];
rz(0.54941191) q[3];
sx q[3];
rz(-1.5447541) q[3];
sx q[3];
rz(2.5766635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3010657) q[2];
sx q[2];
rz(-1.4378005) q[2];
sx q[2];
rz(-2.5145516) q[2];
rz(-0.4380694) q[3];
sx q[3];
rz(-1.4984683) q[3];
sx q[3];
rz(2.0048678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(2.7196734) q[0];
sx q[0];
rz(-1.0872343) q[0];
sx q[0];
rz(1.7701953) q[0];
rz(2.5321391) q[1];
sx q[1];
rz(-2.2513159) q[1];
sx q[1];
rz(-1.2249464) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6951675) q[0];
sx q[0];
rz(-1.7621627) q[0];
sx q[0];
rz(-1.6686977) q[0];
rz(-pi) q[1];
rz(0.89307066) q[2];
sx q[2];
rz(-1.2532824) q[2];
sx q[2];
rz(1.7567321) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.7912253) q[1];
sx q[1];
rz(-1.8944728) q[1];
sx q[1];
rz(1.1926427) q[1];
rz(-pi) q[2];
rz(1.8277934) q[3];
sx q[3];
rz(-2.42958) q[3];
sx q[3];
rz(-2.3745927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8372535) q[2];
sx q[2];
rz(-0.86696583) q[2];
sx q[2];
rz(-0.85401094) q[2];
rz(-0.52418661) q[3];
sx q[3];
rz(-1.5093404) q[3];
sx q[3];
rz(-2.1715651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9147341) q[0];
sx q[0];
rz(-2.9840042) q[0];
sx q[0];
rz(0.7937113) q[0];
rz(3.1397676) q[1];
sx q[1];
rz(-0.55627126) q[1];
sx q[1];
rz(-0.8108286) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9017667) q[0];
sx q[0];
rz(-1.7525273) q[0];
sx q[0];
rz(2.316409) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0982207) q[2];
sx q[2];
rz(-2.3956991) q[2];
sx q[2];
rz(0.30442849) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.86321837) q[1];
sx q[1];
rz(-2.4469118) q[1];
sx q[1];
rz(1.0287697) q[1];
x q[2];
rz(0.14489095) q[3];
sx q[3];
rz(-1.9959183) q[3];
sx q[3];
rz(1.2910049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4578555) q[2];
sx q[2];
rz(-2.808414) q[2];
sx q[2];
rz(-1.4972756) q[2];
rz(-1.3582683) q[3];
sx q[3];
rz(-2.0446916) q[3];
sx q[3];
rz(1.8339405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7233647) q[0];
sx q[0];
rz(-1.4606553) q[0];
sx q[0];
rz(-0.25890589) q[0];
rz(3.018191) q[1];
sx q[1];
rz(-0.63107189) q[1];
sx q[1];
rz(0.48615989) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0927057) q[0];
sx q[0];
rz(-0.98711038) q[0];
sx q[0];
rz(0.58764579) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45275432) q[2];
sx q[2];
rz(-1.6932906) q[2];
sx q[2];
rz(-1.233135) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9100719) q[1];
sx q[1];
rz(-1.9329482) q[1];
sx q[1];
rz(-2.9400184) q[1];
rz(-0.93323243) q[3];
sx q[3];
rz(-0.8094111) q[3];
sx q[3];
rz(-2.2385067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2316124) q[2];
sx q[2];
rz(-2.6284802) q[2];
sx q[2];
rz(-1.7463589) q[2];
rz(-0.95476556) q[3];
sx q[3];
rz(-1.747811) q[3];
sx q[3];
rz(2.098293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2024277) q[0];
sx q[0];
rz(-1.2347777) q[0];
sx q[0];
rz(0.760461) q[0];
rz(0.83433926) q[1];
sx q[1];
rz(-1.1015247) q[1];
sx q[1];
rz(-1.1044501) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.051726159) q[0];
sx q[0];
rz(-1.7891004) q[0];
sx q[0];
rz(-2.1283243) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9707174) q[2];
sx q[2];
rz(-1.8900223) q[2];
sx q[2];
rz(0.17344061) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.901746) q[1];
sx q[1];
rz(-0.86565986) q[1];
sx q[1];
rz(-0.88259952) q[1];
rz(2.1451925) q[3];
sx q[3];
rz(-1.3076915) q[3];
sx q[3];
rz(-2.1679479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1793648) q[2];
sx q[2];
rz(-0.70415512) q[2];
sx q[2];
rz(3.1000225) q[2];
rz(-0.19503322) q[3];
sx q[3];
rz(-0.2427559) q[3];
sx q[3];
rz(-2.354505) q[3];
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
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70260173) q[0];
sx q[0];
rz(-0.85346237) q[0];
sx q[0];
rz(0.75337291) q[0];
rz(2.0808749) q[1];
sx q[1];
rz(-2.3371688) q[1];
sx q[1];
rz(-0.20739584) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8327861) q[0];
sx q[0];
rz(-2.1784003) q[0];
sx q[0];
rz(2.4101546) q[0];
x q[1];
rz(1.2815525) q[2];
sx q[2];
rz(-2.0194032) q[2];
sx q[2];
rz(1.9236652) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7809697) q[1];
sx q[1];
rz(-1.0252254) q[1];
sx q[1];
rz(-1.2956133) q[1];
rz(-pi) q[2];
rz(1.5630521) q[3];
sx q[3];
rz(-1.5317347) q[3];
sx q[3];
rz(-2.9970053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9994026) q[2];
sx q[2];
rz(-1.157607) q[2];
sx q[2];
rz(3.0863777) q[2];
rz(2.2986872) q[3];
sx q[3];
rz(-2.3838145) q[3];
sx q[3];
rz(2.260476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76999369) q[0];
sx q[0];
rz(-2.433625) q[0];
sx q[0];
rz(-1.9644894) q[0];
rz(-0.87989315) q[1];
sx q[1];
rz(-2.8113007) q[1];
sx q[1];
rz(3.0807307) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0165981) q[0];
sx q[0];
rz(-1.6654764) q[0];
sx q[0];
rz(1.4461317) q[0];
rz(1.251077) q[2];
sx q[2];
rz(-1.1100612) q[2];
sx q[2];
rz(-1.899903) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5338075) q[1];
sx q[1];
rz(-2.5609511) q[1];
sx q[1];
rz(0.21368475) q[1];
x q[2];
rz(1.8350321) q[3];
sx q[3];
rz(-1.4995534) q[3];
sx q[3];
rz(-0.83282214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.047478288) q[2];
sx q[2];
rz(-0.55314174) q[2];
sx q[2];
rz(-1.4609569) q[2];
rz(-0.85584062) q[3];
sx q[3];
rz(-0.97270054) q[3];
sx q[3];
rz(2.262825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0395373) q[0];
sx q[0];
rz(-0.85298959) q[0];
sx q[0];
rz(0.0041740388) q[0];
rz(-1.2511085) q[1];
sx q[1];
rz(-3.0079542) q[1];
sx q[1];
rz(-0.68152308) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7345372) q[0];
sx q[0];
rz(-0.87293599) q[0];
sx q[0];
rz(-0.38311335) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84284346) q[2];
sx q[2];
rz(-2.3327391) q[2];
sx q[2];
rz(0.12241546) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8261823) q[1];
sx q[1];
rz(-1.8528954) q[1];
sx q[1];
rz(1.764707) q[1];
rz(-2.8130346) q[3];
sx q[3];
rz(-2.3577981) q[3];
sx q[3];
rz(0.72600049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.556813) q[2];
sx q[2];
rz(-2.6656373) q[2];
sx q[2];
rz(-1.7468096) q[2];
rz(-2.7252588) q[3];
sx q[3];
rz(-1.2952015) q[3];
sx q[3];
rz(-2.7487315) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3861179) q[0];
sx q[0];
rz(-2.8135354) q[0];
sx q[0];
rz(-2.1395444) q[0];
rz(-2.877032) q[1];
sx q[1];
rz(-1.8160276) q[1];
sx q[1];
rz(1.4904259) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4571085) q[0];
sx q[0];
rz(-2.5837645) q[0];
sx q[0];
rz(-2.8929263) q[0];
rz(2.3333346) q[2];
sx q[2];
rz(-1.4604476) q[2];
sx q[2];
rz(2.7791952) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.38768735) q[1];
sx q[1];
rz(-1.6847982) q[1];
sx q[1];
rz(3.0579159) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0496554) q[3];
sx q[3];
rz(-1.5109533) q[3];
sx q[3];
rz(0.53093225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.010290535) q[2];
sx q[2];
rz(-1.4069858) q[2];
sx q[2];
rz(-1.4170125) q[2];
rz(-2.8085282) q[3];
sx q[3];
rz(-0.76996961) q[3];
sx q[3];
rz(1.3781579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604816) q[0];
sx q[0];
rz(-1.2662553) q[0];
sx q[0];
rz(1.8929831) q[0];
rz(-0.026451182) q[1];
sx q[1];
rz(-1.5298264) q[1];
sx q[1];
rz(-1.5419921) q[1];
rz(1.1556861) q[2];
sx q[2];
rz(-2.2397556) q[2];
sx q[2];
rz(-0.02469183) q[2];
rz(1.2534441) q[3];
sx q[3];
rz(-2.1614634) q[3];
sx q[3];
rz(0.95017016) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
