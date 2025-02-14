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
rz(-1.0903519) q[0];
sx q[0];
rz(-0.89972377) q[0];
sx q[0];
rz(-1.1296912) q[0];
rz(0.55454412) q[1];
sx q[1];
rz(-0.50538844) q[1];
sx q[1];
rz(2.5880421) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77614895) q[0];
sx q[0];
rz(-1.5834785) q[0];
sx q[0];
rz(-2.074341) q[0];
rz(3.1316544) q[2];
sx q[2];
rz(-1.9807837) q[2];
sx q[2];
rz(-0.25914295) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.62088242) q[1];
sx q[1];
rz(-1.1739879) q[1];
sx q[1];
rz(-0.75735332) q[1];
rz(-pi) q[2];
rz(-3.0985988) q[3];
sx q[3];
rz(-1.3519796) q[3];
sx q[3];
rz(0.18312632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0022137) q[2];
sx q[2];
rz(-1.4993818) q[2];
sx q[2];
rz(1.1688983) q[2];
rz(0.72541952) q[3];
sx q[3];
rz(-1.9099216) q[3];
sx q[3];
rz(-1.5554844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91189522) q[0];
sx q[0];
rz(-0.40922368) q[0];
sx q[0];
rz(2.0912066) q[0];
rz(1.7627675) q[1];
sx q[1];
rz(-1.4070815) q[1];
sx q[1];
rz(1.001766) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9501099) q[0];
sx q[0];
rz(-1.8284197) q[0];
sx q[0];
rz(-0.10665032) q[0];
rz(-1.1092735) q[2];
sx q[2];
rz(-1.2987198) q[2];
sx q[2];
rz(0.22914722) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.54047062) q[1];
sx q[1];
rz(-1.3408288) q[1];
sx q[1];
rz(0.44039393) q[1];
x q[2];
rz(-1.863594) q[3];
sx q[3];
rz(-1.1637886) q[3];
sx q[3];
rz(-1.3425105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0687781) q[2];
sx q[2];
rz(-0.38663703) q[2];
sx q[2];
rz(1.4978503) q[2];
rz(-2.1554598) q[3];
sx q[3];
rz(-0.96753263) q[3];
sx q[3];
rz(0.51924527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7943952) q[0];
sx q[0];
rz(-2.6674542) q[0];
sx q[0];
rz(3.0679605) q[0];
rz(2.4108048) q[1];
sx q[1];
rz(-1.7981139) q[1];
sx q[1];
rz(-2.1582019) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4016664) q[0];
sx q[0];
rz(-2.8186975) q[0];
sx q[0];
rz(-2.0773204) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86536256) q[2];
sx q[2];
rz(-1.9459122) q[2];
sx q[2];
rz(-0.6626216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.19304744) q[1];
sx q[1];
rz(-1.3009503) q[1];
sx q[1];
rz(-0.13040925) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7119104) q[3];
sx q[3];
rz(-0.50715461) q[3];
sx q[3];
rz(0.53559723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.44707766) q[2];
sx q[2];
rz(-2.647001) q[2];
sx q[2];
rz(0.05973235) q[2];
rz(-0.51074243) q[3];
sx q[3];
rz(-0.35664883) q[3];
sx q[3];
rz(1.7662778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8259976) q[0];
sx q[0];
rz(-0.70846486) q[0];
sx q[0];
rz(-0.84061709) q[0];
rz(-2.9261342) q[1];
sx q[1];
rz(-1.0865728) q[1];
sx q[1];
rz(-0.1159018) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2471231) q[0];
sx q[0];
rz(-0.43465675) q[0];
sx q[0];
rz(-2.2306395) q[0];
x q[1];
rz(-2.5253549) q[2];
sx q[2];
rz(-0.35348693) q[2];
sx q[2];
rz(-0.36187672) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.92479372) q[1];
sx q[1];
rz(-1.2844286) q[1];
sx q[1];
rz(-2.22424) q[1];
rz(-1.674599) q[3];
sx q[3];
rz(-1.1624059) q[3];
sx q[3];
rz(-0.30882747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0243715) q[2];
sx q[2];
rz(-1.2933967) q[2];
sx q[2];
rz(0.8563861) q[2];
rz(-0.74445009) q[3];
sx q[3];
rz(-0.88862935) q[3];
sx q[3];
rz(-0.3652679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15904467) q[0];
sx q[0];
rz(-0.58369517) q[0];
sx q[0];
rz(2.1481376) q[0];
rz(-2.0221201) q[1];
sx q[1];
rz(-1.9529724) q[1];
sx q[1];
rz(-0.022389222) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9291344) q[0];
sx q[0];
rz(-1.6957383) q[0];
sx q[0];
rz(1.4897299) q[0];
rz(2.1447626) q[2];
sx q[2];
rz(-1.0901684) q[2];
sx q[2];
rz(-0.89546452) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7931079) q[1];
sx q[1];
rz(-0.99103084) q[1];
sx q[1];
rz(-2.0391885) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56110142) q[3];
sx q[3];
rz(-0.77563876) q[3];
sx q[3];
rz(-2.5833094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8523954) q[2];
sx q[2];
rz(-0.22326938) q[2];
sx q[2];
rz(1.4168463) q[2];
rz(-1.019574) q[3];
sx q[3];
rz(-2.151078) q[3];
sx q[3];
rz(-2.8250601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1932909) q[0];
sx q[0];
rz(-0.010007771) q[0];
sx q[0];
rz(-2.365812) q[0];
rz(-1.4099482) q[1];
sx q[1];
rz(-1.8312788) q[1];
sx q[1];
rz(1.5595248) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55627953) q[0];
sx q[0];
rz(-1.5405498) q[0];
sx q[0];
rz(1.6325238) q[0];
x q[1];
rz(-2.1072793) q[2];
sx q[2];
rz(-2.5051077) q[2];
sx q[2];
rz(0.45490593) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3821545) q[1];
sx q[1];
rz(-1.702679) q[1];
sx q[1];
rz(-2.0410755) q[1];
rz(-1.9385064) q[3];
sx q[3];
rz(-1.0818521) q[3];
sx q[3];
rz(0.5209825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.55283028) q[2];
sx q[2];
rz(-1.4344119) q[2];
sx q[2];
rz(2.8130048) q[2];
rz(1.5594679) q[3];
sx q[3];
rz(-2.4170473) q[3];
sx q[3];
rz(0.72415486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5174823) q[0];
sx q[0];
rz(-2.3243853) q[0];
sx q[0];
rz(-0.66993129) q[0];
rz(-1.5176557) q[1];
sx q[1];
rz(-2.0567963) q[1];
sx q[1];
rz(-2.0673015) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5094604) q[0];
sx q[0];
rz(-1.355251) q[0];
sx q[0];
rz(-0.68595217) q[0];
rz(-1.7543484) q[2];
sx q[2];
rz(-1.2129686) q[2];
sx q[2];
rz(2.7821469) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.049557471) q[1];
sx q[1];
rz(-1.2566031) q[1];
sx q[1];
rz(-0.097196984) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89364918) q[3];
sx q[3];
rz(-1.7405677) q[3];
sx q[3];
rz(-0.61299619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39277789) q[2];
sx q[2];
rz(-2.1608976) q[2];
sx q[2];
rz(-1.7559715) q[2];
rz(2.4393926) q[3];
sx q[3];
rz(-0.32679138) q[3];
sx q[3];
rz(-2.0328111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4362519) q[0];
sx q[0];
rz(-1.5974644) q[0];
sx q[0];
rz(-0.77823773) q[0];
rz(-1.2938007) q[1];
sx q[1];
rz(-1.2575282) q[1];
sx q[1];
rz(3.0857118) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.125408) q[0];
sx q[0];
rz(-1.2795537) q[0];
sx q[0];
rz(1.4945369) q[0];
x q[1];
rz(1.3827902) q[2];
sx q[2];
rz(-1.3488608) q[2];
sx q[2];
rz(1.7772016) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5949616) q[1];
sx q[1];
rz(-1.43838) q[1];
sx q[1];
rz(-2.754171) q[1];
x q[2];
rz(-2.1382525) q[3];
sx q[3];
rz(-2.3262089) q[3];
sx q[3];
rz(1.5559097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0539661) q[2];
sx q[2];
rz(-2.1199333) q[2];
sx q[2];
rz(0.6915687) q[2];
rz(2.9223053) q[3];
sx q[3];
rz(-1.6577474) q[3];
sx q[3];
rz(1.1933491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9241065) q[0];
sx q[0];
rz(-0.036343887) q[0];
sx q[0];
rz(-2.0935667) q[0];
rz(-0.74608392) q[1];
sx q[1];
rz(-1.8522976) q[1];
sx q[1];
rz(2.7229436) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6857938) q[0];
sx q[0];
rz(-1.5869194) q[0];
sx q[0];
rz(-3.1324194) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2144958) q[2];
sx q[2];
rz(-1.3852556) q[2];
sx q[2];
rz(-2.5531921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.56142226) q[1];
sx q[1];
rz(-0.19374312) q[1];
sx q[1];
rz(-2.3592739) q[1];
x q[2];
rz(-0.04808406) q[3];
sx q[3];
rz(-0.86838956) q[3];
sx q[3];
rz(2.4772205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5082385) q[2];
sx q[2];
rz(-2.5373122) q[2];
sx q[2];
rz(2.5227127) q[2];
rz(-0.45281705) q[3];
sx q[3];
rz(-1.6524977) q[3];
sx q[3];
rz(0.49310163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.598269) q[0];
sx q[0];
rz(-3.0040574) q[0];
sx q[0];
rz(-1.2602873) q[0];
rz(0.29809412) q[1];
sx q[1];
rz(-2.1838078) q[1];
sx q[1];
rz(1.4990998) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84628203) q[0];
sx q[0];
rz(-2.3822226) q[0];
sx q[0];
rz(1.6037386) q[0];
rz(-pi) q[1];
rz(-0.46447931) q[2];
sx q[2];
rz(-0.3128007) q[2];
sx q[2];
rz(2.0020747) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3819653) q[1];
sx q[1];
rz(-1.5542728) q[1];
sx q[1];
rz(-1.7970867) q[1];
x q[2];
rz(-1.9630249) q[3];
sx q[3];
rz(-1.4334599) q[3];
sx q[3];
rz(-1.8399505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5160211) q[2];
sx q[2];
rz(-1.2837774) q[2];
sx q[2];
rz(0.18216356) q[2];
rz(-0.22570172) q[3];
sx q[3];
rz(-1.00939) q[3];
sx q[3];
rz(-1.3753447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43212776) q[0];
sx q[0];
rz(-1.3668677) q[0];
sx q[0];
rz(0.9912542) q[0];
rz(1.2200914) q[1];
sx q[1];
rz(-2.5972001) q[1];
sx q[1];
rz(-0.77144815) q[1];
rz(0.091339672) q[2];
sx q[2];
rz(-1.7484574) q[2];
sx q[2];
rz(-1.4469528) q[2];
rz(-0.78990084) q[3];
sx q[3];
rz(-0.98179437) q[3];
sx q[3];
rz(-2.475987) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
