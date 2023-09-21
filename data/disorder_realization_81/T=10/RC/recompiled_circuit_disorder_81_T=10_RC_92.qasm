OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.73683357) q[0];
sx q[0];
rz(-1.3614549) q[0];
sx q[0];
rz(1.7629495) q[0];
rz(-0.8575851) q[1];
sx q[1];
rz(-1.4839988) q[1];
sx q[1];
rz(-2.690697) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3987797) q[0];
sx q[0];
rz(-3.0516041) q[0];
sx q[0];
rz(2.9773657) q[0];
rz(-1.3583463) q[2];
sx q[2];
rz(-0.48848402) q[2];
sx q[2];
rz(-2.8993895) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3528459) q[1];
sx q[1];
rz(-2.5874918) q[1];
sx q[1];
rz(-2.7098141) q[1];
rz(0.88702918) q[3];
sx q[3];
rz(-1.4966045) q[3];
sx q[3];
rz(-0.85899734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9840055) q[2];
sx q[2];
rz(-1.459534) q[2];
sx q[2];
rz(2.297304) q[2];
rz(-0.44101161) q[3];
sx q[3];
rz(-2.7859272) q[3];
sx q[3];
rz(-0.60602337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.59250295) q[0];
sx q[0];
rz(-1.9117768) q[0];
sx q[0];
rz(2.8785008) q[0];
rz(2.198055) q[1];
sx q[1];
rz(-2.5448006) q[1];
sx q[1];
rz(1.1862322) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037558) q[0];
sx q[0];
rz(-0.058996011) q[0];
sx q[0];
rz(-2.8179413) q[0];
rz(-pi) q[1];
rz(0.29859121) q[2];
sx q[2];
rz(-2.9332187) q[2];
sx q[2];
rz(1.5339472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.26145229) q[1];
sx q[1];
rz(-1.6958691) q[1];
sx q[1];
rz(-0.90420453) q[1];
x q[2];
rz(-2.4466189) q[3];
sx q[3];
rz(-1.4222099) q[3];
sx q[3];
rz(0.23651628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1295604) q[2];
sx q[2];
rz(-1.0027145) q[2];
sx q[2];
rz(-1.9821232) q[2];
rz(2.7705079) q[3];
sx q[3];
rz(-1.5044731) q[3];
sx q[3];
rz(-2.8306567) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7611258) q[0];
sx q[0];
rz(-2.0115871) q[0];
sx q[0];
rz(-2.3348715) q[0];
rz(-0.21356788) q[1];
sx q[1];
rz(-0.49626207) q[1];
sx q[1];
rz(0.82021964) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41528156) q[0];
sx q[0];
rz(-1.6831241) q[0];
sx q[0];
rz(2.2583654) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9224612) q[2];
sx q[2];
rz(-1.0872772) q[2];
sx q[2];
rz(-0.52106524) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.080938235) q[1];
sx q[1];
rz(-1.9229691) q[1];
sx q[1];
rz(2.0209795) q[1];
x q[2];
rz(0.1180325) q[3];
sx q[3];
rz(-2.0327838) q[3];
sx q[3];
rz(-0.0052009728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8308668) q[2];
sx q[2];
rz(-1.6409637) q[2];
sx q[2];
rz(-2.2107928) q[2];
rz(-2.9860949) q[3];
sx q[3];
rz(-1.6379387) q[3];
sx q[3];
rz(0.29155198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.528462) q[0];
sx q[0];
rz(-0.72137946) q[0];
sx q[0];
rz(0.91127515) q[0];
rz(0.43831929) q[1];
sx q[1];
rz(-1.3221909) q[1];
sx q[1];
rz(1.320425) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62990084) q[0];
sx q[0];
rz(-1.9803847) q[0];
sx q[0];
rz(-0.011261777) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7772061) q[2];
sx q[2];
rz(-1.1279391) q[2];
sx q[2];
rz(-2.4368311) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.98852324) q[1];
sx q[1];
rz(-1.8338025) q[1];
sx q[1];
rz(2.4371229) q[1];
rz(-pi) q[2];
x q[2];
rz(0.78762357) q[3];
sx q[3];
rz(-2.0260603) q[3];
sx q[3];
rz(-2.5256707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0358255) q[2];
sx q[2];
rz(-2.2129009) q[2];
sx q[2];
rz(-2.7992115) q[2];
rz(-2.9648182) q[3];
sx q[3];
rz(-2.7084559) q[3];
sx q[3];
rz(2.0006196) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3115561) q[0];
sx q[0];
rz(-2.4139068) q[0];
sx q[0];
rz(0.86529055) q[0];
rz(1.9150437) q[1];
sx q[1];
rz(-0.98926917) q[1];
sx q[1];
rz(-1.3006166) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0375992) q[0];
sx q[0];
rz(-0.35498699) q[0];
sx q[0];
rz(-0.07261891) q[0];
x q[1];
rz(2.775035) q[2];
sx q[2];
rz(-1.6379116) q[2];
sx q[2];
rz(1.7442489) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.76104858) q[1];
sx q[1];
rz(-0.64861464) q[1];
sx q[1];
rz(-2.5785239) q[1];
rz(2.6404068) q[3];
sx q[3];
rz(-2.8445344) q[3];
sx q[3];
rz(0.81851573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1317923) q[2];
sx q[2];
rz(-1.1494145) q[2];
sx q[2];
rz(2.664393) q[2];
rz(0.19208433) q[3];
sx q[3];
rz(-1.447907) q[3];
sx q[3];
rz(-2.208476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79648298) q[0];
sx q[0];
rz(-2.5273297) q[0];
sx q[0];
rz(3.1298424) q[0];
rz(0.55039644) q[1];
sx q[1];
rz(-1.7852716) q[1];
sx q[1];
rz(1.5531497) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4295411) q[0];
sx q[0];
rz(-2.6459604) q[0];
sx q[0];
rz(0.80722157) q[0];
x q[1];
rz(-3.1197238) q[2];
sx q[2];
rz(-1.2049897) q[2];
sx q[2];
rz(1.0794229) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8223871) q[1];
sx q[1];
rz(-0.86383312) q[1];
sx q[1];
rz(1.6824526) q[1];
rz(-pi) q[2];
rz(0.074514975) q[3];
sx q[3];
rz(-2.2722368) q[3];
sx q[3];
rz(2.6889192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6340296) q[2];
sx q[2];
rz(-2.4812249) q[2];
sx q[2];
rz(-1.2825512) q[2];
rz(-1.7717308) q[3];
sx q[3];
rz(-1.3953352) q[3];
sx q[3];
rz(1.1184568) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5605374) q[0];
sx q[0];
rz(-0.16796172) q[0];
sx q[0];
rz(2.4643331) q[0];
rz(-0.15180763) q[1];
sx q[1];
rz(-1.3744524) q[1];
sx q[1];
rz(0.97704926) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86977406) q[0];
sx q[0];
rz(-0.5973814) q[0];
sx q[0];
rz(3.0074044) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5706967) q[2];
sx q[2];
rz(-1.7028371) q[2];
sx q[2];
rz(-0.11167234) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.44927412) q[1];
sx q[1];
rz(-2.793503) q[1];
sx q[1];
rz(1.3432137) q[1];
rz(2.3903923) q[3];
sx q[3];
rz(-0.45414543) q[3];
sx q[3];
rz(-2.5129012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3892422) q[2];
sx q[2];
rz(-2.3198979) q[2];
sx q[2];
rz(-1.0127257) q[2];
rz(-1.9536473) q[3];
sx q[3];
rz(-2.0690737) q[3];
sx q[3];
rz(0.48721203) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1241207) q[0];
sx q[0];
rz(-3.108232) q[0];
sx q[0];
rz(2.4429328) q[0];
rz(-1.1220804) q[1];
sx q[1];
rz(-0.84609234) q[1];
sx q[1];
rz(-1.2493856) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65803618) q[0];
sx q[0];
rz(-1.5244966) q[0];
sx q[0];
rz(-0.21090837) q[0];
rz(-pi) q[1];
rz(-0.47301936) q[2];
sx q[2];
rz(-1.0968536) q[2];
sx q[2];
rz(2.3854286) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.57751209) q[1];
sx q[1];
rz(-3.1187594) q[1];
sx q[1];
rz(0.24502416) q[1];
rz(-0.28835339) q[3];
sx q[3];
rz(-0.95803146) q[3];
sx q[3];
rz(1.8307277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.32968783) q[2];
sx q[2];
rz(-0.78616443) q[2];
sx q[2];
rz(1.1784941) q[2];
rz(-1.684749) q[3];
sx q[3];
rz(-2.0791576) q[3];
sx q[3];
rz(-0.38213521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(-1.6417398) q[0];
sx q[0];
rz(-1.3795744) q[0];
sx q[0];
rz(-1.8485803) q[0];
rz(-1.4216084) q[1];
sx q[1];
rz(-1.0363818) q[1];
sx q[1];
rz(-0.59757772) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9874728) q[0];
sx q[0];
rz(-3.0010536) q[0];
sx q[0];
rz(-1.8494291) q[0];
x q[1];
rz(2.5848128) q[2];
sx q[2];
rz(-2.0282929) q[2];
sx q[2];
rz(0.58820398) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4328879) q[1];
sx q[1];
rz(-1.8678027) q[1];
sx q[1];
rz(1.2328641) q[1];
x q[2];
rz(2.5556373) q[3];
sx q[3];
rz(-1.5065985) q[3];
sx q[3];
rz(-0.5639329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9188345) q[2];
sx q[2];
rz(-1.4514048) q[2];
sx q[2];
rz(-1.9082327) q[2];
rz(0.90138609) q[3];
sx q[3];
rz(-0.12005761) q[3];
sx q[3];
rz(1.4982769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0937061) q[0];
sx q[0];
rz(-2.369635) q[0];
sx q[0];
rz(-3.1179324) q[0];
rz(-2.1854782) q[1];
sx q[1];
rz(-1.3095983) q[1];
sx q[1];
rz(0.67217174) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32165747) q[0];
sx q[0];
rz(-1.5655087) q[0];
sx q[0];
rz(1.5810285) q[0];
rz(-pi) q[1];
rz(0.51151885) q[2];
sx q[2];
rz(-1.5788955) q[2];
sx q[2];
rz(2.9715003) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.086239554) q[1];
sx q[1];
rz(-1.9085911) q[1];
sx q[1];
rz(2.4243381) q[1];
rz(-1.724733) q[3];
sx q[3];
rz(-1.142475) q[3];
sx q[3];
rz(1.1101013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6293634) q[2];
sx q[2];
rz(-1.2269292) q[2];
sx q[2];
rz(-2.771634) q[2];
rz(-1.6379179) q[3];
sx q[3];
rz(-2.2556997) q[3];
sx q[3];
rz(1.2009719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5794012) q[0];
sx q[0];
rz(-2.777522) q[0];
sx q[0];
rz(1.2072442) q[0];
rz(0.72369408) q[1];
sx q[1];
rz(-0.98725286) q[1];
sx q[1];
rz(-0.90686803) q[1];
rz(1.205668) q[2];
sx q[2];
rz(-0.42386133) q[2];
sx q[2];
rz(-0.78122666) q[2];
rz(-1.9772114) q[3];
sx q[3];
rz(-1.2953399) q[3];
sx q[3];
rz(-3.0084707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];