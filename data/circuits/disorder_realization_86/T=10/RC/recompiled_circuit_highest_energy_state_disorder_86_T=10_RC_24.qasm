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
rz(-1.2404233) q[0];
sx q[0];
rz(1.9440396) q[0];
sx q[0];
rz(9.64111) q[0];
rz(0.95739111) q[1];
sx q[1];
rz(-2.4654145) q[1];
sx q[1];
rz(1.5114991) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6304361) q[0];
sx q[0];
rz(-1.5804218) q[0];
sx q[0];
rz(-0.55625347) q[0];
rz(-pi) q[1];
rz(2.6681603) q[2];
sx q[2];
rz(-1.8634999) q[2];
sx q[2];
rz(-2.6930489) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.87147698) q[1];
sx q[1];
rz(-1.2387215) q[1];
sx q[1];
rz(0.52712743) q[1];
x q[2];
rz(-2.0262296) q[3];
sx q[3];
rz(-2.3854227) q[3];
sx q[3];
rz(0.84634534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1067074) q[2];
sx q[2];
rz(-0.35197508) q[2];
sx q[2];
rz(2.0931639) q[2];
rz(-0.18167051) q[3];
sx q[3];
rz(-2.1766267) q[3];
sx q[3];
rz(0.65339965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.492391) q[0];
sx q[0];
rz(-0.94434706) q[0];
sx q[0];
rz(-2.6948068) q[0];
rz(0.88042879) q[1];
sx q[1];
rz(-1.3648405) q[1];
sx q[1];
rz(0.78278881) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.114417) q[0];
sx q[0];
rz(-1.8150629) q[0];
sx q[0];
rz(0.063112325) q[0];
x q[1];
rz(1.1433463) q[2];
sx q[2];
rz(-2.1437217) q[2];
sx q[2];
rz(-1.7597511) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35905169) q[1];
sx q[1];
rz(-2.5757669) q[1];
sx q[1];
rz(-2.5579909) q[1];
rz(-pi) q[2];
rz(0.47329013) q[3];
sx q[3];
rz(-0.9447228) q[3];
sx q[3];
rz(1.2155217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.11144) q[2];
sx q[2];
rz(-1.6609265) q[2];
sx q[2];
rz(0.97935575) q[2];
rz(2.7566946) q[3];
sx q[3];
rz(-1.9210457) q[3];
sx q[3];
rz(-0.16429193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.6556743) q[0];
sx q[0];
rz(-3.0874708) q[0];
sx q[0];
rz(-0.78980494) q[0];
rz(2.9549331) q[1];
sx q[1];
rz(-1.7402382) q[1];
sx q[1];
rz(-2.1479215) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8756722) q[0];
sx q[0];
rz(-1.8199931) q[0];
sx q[0];
rz(-0.56206352) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3445856) q[2];
sx q[2];
rz(-0.96983428) q[2];
sx q[2];
rz(0.41659875) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4504272) q[1];
sx q[1];
rz(-1.9948927) q[1];
sx q[1];
rz(1.0662088) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30672725) q[3];
sx q[3];
rz(-2.3969458) q[3];
sx q[3];
rz(-0.70017231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8197202) q[2];
sx q[2];
rz(-2.0237782) q[2];
sx q[2];
rz(0.37128386) q[2];
rz(-0.34879455) q[3];
sx q[3];
rz(-1.0943509) q[3];
sx q[3];
rz(2.546052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6868941) q[0];
sx q[0];
rz(-1.0201539) q[0];
sx q[0];
rz(-1.7373079) q[0];
rz(-0.67717254) q[1];
sx q[1];
rz(-1.155747) q[1];
sx q[1];
rz(-1.4926532) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.416788) q[0];
sx q[0];
rz(-1.3287523) q[0];
sx q[0];
rz(1.7927732) q[0];
rz(-pi) q[1];
rz(-2.1588232) q[2];
sx q[2];
rz(-2.6361536) q[2];
sx q[2];
rz(-2.7014521) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.75912913) q[1];
sx q[1];
rz(-1.3037523) q[1];
sx q[1];
rz(-1.3401396) q[1];
rz(-pi) q[2];
rz(1.221232) q[3];
sx q[3];
rz(-2.3858983) q[3];
sx q[3];
rz(-0.44103482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6877785) q[2];
sx q[2];
rz(-2.1731589) q[2];
sx q[2];
rz(-2.7933534) q[2];
rz(-1.4549152) q[3];
sx q[3];
rz(-1.429052) q[3];
sx q[3];
rz(-1.3454364) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9652047) q[0];
sx q[0];
rz(-2.5768953) q[0];
sx q[0];
rz(0.89214605) q[0];
rz(-0.46420321) q[1];
sx q[1];
rz(-1.8966388) q[1];
sx q[1];
rz(-1.2947327) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2703122) q[0];
sx q[0];
rz(-2.0593606) q[0];
sx q[0];
rz(-2.110092) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6991529) q[2];
sx q[2];
rz(-1.0367298) q[2];
sx q[2];
rz(0.72557025) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.66440873) q[1];
sx q[1];
rz(-0.55969612) q[1];
sx q[1];
rz(-2.2137292) q[1];
rz(-pi) q[2];
rz(-2.283461) q[3];
sx q[3];
rz(-2.3225124) q[3];
sx q[3];
rz(-1.8392966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3225473) q[2];
sx q[2];
rz(-0.61083856) q[2];
sx q[2];
rz(2.5757705) q[2];
rz(3.0602509) q[3];
sx q[3];
rz(-0.9809202) q[3];
sx q[3];
rz(2.5210023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6767947) q[0];
sx q[0];
rz(-3.0749574) q[0];
sx q[0];
rz(1.5555405) q[0];
rz(-2.0809035) q[1];
sx q[1];
rz(-1.5763177) q[1];
sx q[1];
rz(0.63180822) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6752967) q[0];
sx q[0];
rz(-1.6501556) q[0];
sx q[0];
rz(-0.23510374) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.86132974) q[2];
sx q[2];
rz(-1.1785186) q[2];
sx q[2];
rz(-2.5226468) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.374208) q[1];
sx q[1];
rz(-0.65113089) q[1];
sx q[1];
rz(-2.1506449) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4238213) q[3];
sx q[3];
rz(-2.2490361) q[3];
sx q[3];
rz(3.105046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1598728) q[2];
sx q[2];
rz(-1.1184511) q[2];
sx q[2];
rz(0.49883207) q[2];
rz(1.3129129) q[3];
sx q[3];
rz(-2.4098318) q[3];
sx q[3];
rz(-1.4341199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6071534) q[0];
sx q[0];
rz(-0.3648912) q[0];
sx q[0];
rz(-2.4420807) q[0];
rz(-0.46547678) q[1];
sx q[1];
rz(-2.270348) q[1];
sx q[1];
rz(2.2043998) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1904966) q[0];
sx q[0];
rz(-1.7825923) q[0];
sx q[0];
rz(1.7429211) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9019321) q[2];
sx q[2];
rz(-1.6722888) q[2];
sx q[2];
rz(3.0200151) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.921007) q[1];
sx q[1];
rz(-1.3594419) q[1];
sx q[1];
rz(1.2799954) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30234594) q[3];
sx q[3];
rz(-1.8519326) q[3];
sx q[3];
rz(0.1880364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.297544) q[2];
sx q[2];
rz(-1.301845) q[2];
sx q[2];
rz(2.76827) q[2];
rz(-1.9518055) q[3];
sx q[3];
rz(-0.50896421) q[3];
sx q[3];
rz(-2.6313307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74476403) q[0];
sx q[0];
rz(-1.0823534) q[0];
sx q[0];
rz(0.74791351) q[0];
rz(0.76639908) q[1];
sx q[1];
rz(-0.26873573) q[1];
sx q[1];
rz(3.1386197) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2073313) q[0];
sx q[0];
rz(-2.2483279) q[0];
sx q[0];
rz(2.8939061) q[0];
x q[1];
rz(1.6417129) q[2];
sx q[2];
rz(-1.5902963) q[2];
sx q[2];
rz(1.8766581) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6312478) q[1];
sx q[1];
rz(-1.6142134) q[1];
sx q[1];
rz(-1.4166635) q[1];
x q[2];
rz(1.3722754) q[3];
sx q[3];
rz(-1.7619507) q[3];
sx q[3];
rz(1.6894345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.030674) q[2];
sx q[2];
rz(-1.7577533) q[2];
sx q[2];
rz(2.0745011) q[2];
rz(-0.083960697) q[3];
sx q[3];
rz(-0.48560086) q[3];
sx q[3];
rz(-2.3626204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79754168) q[0];
sx q[0];
rz(-0.9386971) q[0];
sx q[0];
rz(-3.0352266) q[0];
rz(0.97995177) q[1];
sx q[1];
rz(-1.6500902) q[1];
sx q[1];
rz(0.76593691) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.04127114) q[0];
sx q[0];
rz(-2.092397) q[0];
sx q[0];
rz(1.0688416) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.19975234) q[2];
sx q[2];
rz(-1.482666) q[2];
sx q[2];
rz(-0.55921184) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.76141833) q[1];
sx q[1];
rz(-1.5360218) q[1];
sx q[1];
rz(0.67616391) q[1];
rz(-3.1274904) q[3];
sx q[3];
rz(-1.5670793) q[3];
sx q[3];
rz(0.93059082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.47436675) q[2];
sx q[2];
rz(-0.52046481) q[2];
sx q[2];
rz(-0.70029798) q[2];
rz(-0.66655603) q[3];
sx q[3];
rz(-1.8066112) q[3];
sx q[3];
rz(-1.0880281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67548442) q[0];
sx q[0];
rz(-2.1145144) q[0];
sx q[0];
rz(-1.745537) q[0];
rz(1.667977) q[1];
sx q[1];
rz(-1.8831848) q[1];
sx q[1];
rz(2.4748763) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7233492) q[0];
sx q[0];
rz(-0.94941345) q[0];
sx q[0];
rz(-2.3410812) q[0];
rz(-pi) q[1];
rz(2.780203) q[2];
sx q[2];
rz(-2.4494684) q[2];
sx q[2];
rz(-0.036594242) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.279505) q[1];
sx q[1];
rz(-0.34395978) q[1];
sx q[1];
rz(0.39744795) q[1];
x q[2];
rz(1.0239086) q[3];
sx q[3];
rz(-2.0738261) q[3];
sx q[3];
rz(1.5708814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.066976808) q[2];
sx q[2];
rz(-0.65698996) q[2];
sx q[2];
rz(0.89078772) q[2];
rz(2.7095419) q[3];
sx q[3];
rz(-1.0703577) q[3];
sx q[3];
rz(0.74705684) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73410949) q[0];
sx q[0];
rz(-1.643184) q[0];
sx q[0];
rz(-1.2947422) q[0];
rz(1.8941849) q[1];
sx q[1];
rz(-2.4220962) q[1];
sx q[1];
rz(-1.9069506) q[1];
rz(1.1449849) q[2];
sx q[2];
rz(-2.7663284) q[2];
sx q[2];
rz(-0.13079499) q[2];
rz(-2.6525146) q[3];
sx q[3];
rz(-1.5413956) q[3];
sx q[3];
rz(1.2914381) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
