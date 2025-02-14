OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.40795657) q[0];
sx q[0];
rz(3.3269296) q[0];
sx q[0];
rz(11.067631) q[0];
rz(-0.19417956) q[1];
sx q[1];
rz(-0.3436389) q[1];
sx q[1];
rz(2.763881) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.728017) q[0];
sx q[0];
rz(-2.8110162) q[0];
sx q[0];
rz(-2.3982993) q[0];
x q[1];
rz(0.80773662) q[2];
sx q[2];
rz(-0.57538549) q[2];
sx q[2];
rz(-0.45562109) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39241782) q[1];
sx q[1];
rz(-0.91411829) q[1];
sx q[1];
rz(-1.6139469) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1123766) q[3];
sx q[3];
rz(-1.3061285) q[3];
sx q[3];
rz(2.9801369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3898042) q[2];
sx q[2];
rz(-1.3499539) q[2];
sx q[2];
rz(2.8765615) q[2];
rz(-2.7122279) q[3];
sx q[3];
rz(-1.2484442) q[3];
sx q[3];
rz(-3.140669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4877743) q[0];
sx q[0];
rz(-2.9980897) q[0];
sx q[0];
rz(-1.300746) q[0];
rz(-0.067151345) q[1];
sx q[1];
rz(-0.90922272) q[1];
sx q[1];
rz(2.2815509) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22500998) q[0];
sx q[0];
rz(-0.89677484) q[0];
sx q[0];
rz(2.3300578) q[0];
rz(-pi) q[1];
rz(0.76171909) q[2];
sx q[2];
rz(-2.2713619) q[2];
sx q[2];
rz(-0.127244) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7262267) q[1];
sx q[1];
rz(-0.72390717) q[1];
sx q[1];
rz(2.2103146) q[1];
rz(-1.8710526) q[3];
sx q[3];
rz(-1.7466063) q[3];
sx q[3];
rz(0.45775698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.63090515) q[2];
sx q[2];
rz(-2.5248933) q[2];
sx q[2];
rz(-0.56646937) q[2];
rz(-0.72930068) q[3];
sx q[3];
rz(-2.2972378) q[3];
sx q[3];
rz(2.9384889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94839621) q[0];
sx q[0];
rz(-1.0862792) q[0];
sx q[0];
rz(-0.36619827) q[0];
rz(-1.5557479) q[1];
sx q[1];
rz(-1.0428753) q[1];
sx q[1];
rz(-3.0911456) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6696861) q[0];
sx q[0];
rz(-1.6499053) q[0];
sx q[0];
rz(0.0418274) q[0];
rz(-pi) q[1];
rz(-0.57951219) q[2];
sx q[2];
rz(-0.77634927) q[2];
sx q[2];
rz(0.001359847) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0512426) q[1];
sx q[1];
rz(-1.7900677) q[1];
sx q[1];
rz(-1.7419113) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9962323) q[3];
sx q[3];
rz(-1.1499377) q[3];
sx q[3];
rz(-0.03736729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0744276) q[2];
sx q[2];
rz(-2.0815492) q[2];
sx q[2];
rz(1.830843) q[2];
rz(1.7835167) q[3];
sx q[3];
rz(-2.570593) q[3];
sx q[3];
rz(1.3091492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27291372) q[0];
sx q[0];
rz(-2.9153115) q[0];
sx q[0];
rz(-1.3847466) q[0];
rz(1.9028496) q[1];
sx q[1];
rz(-1.9107995) q[1];
sx q[1];
rz(1.1741656) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096701972) q[0];
sx q[0];
rz(-0.47098038) q[0];
sx q[0];
rz(-1.1440008) q[0];
rz(-pi) q[1];
rz(-1.3210682) q[2];
sx q[2];
rz(-1.839723) q[2];
sx q[2];
rz(-0.099322546) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82249043) q[1];
sx q[1];
rz(-0.61885364) q[1];
sx q[1];
rz(-0.98554218) q[1];
rz(2.6888108) q[3];
sx q[3];
rz(-2.2119129) q[3];
sx q[3];
rz(2.6124817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.43526402) q[2];
sx q[2];
rz(-1.3202983) q[2];
sx q[2];
rz(-3.1026133) q[2];
rz(-2.4987761) q[3];
sx q[3];
rz(-2.1233852) q[3];
sx q[3];
rz(2.5207998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3103127) q[0];
sx q[0];
rz(-2.1273002) q[0];
sx q[0];
rz(-0.020462791) q[0];
rz(2.9488355) q[1];
sx q[1];
rz(-0.61734504) q[1];
sx q[1];
rz(-1.9761168) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.504731) q[0];
sx q[0];
rz(-0.76380542) q[0];
sx q[0];
rz(0.59815191) q[0];
x q[1];
rz(2.0203677) q[2];
sx q[2];
rz(-2.1940941) q[2];
sx q[2];
rz(-0.69823182) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.36042696) q[1];
sx q[1];
rz(-1.5910223) q[1];
sx q[1];
rz(-1.0670877) q[1];
x q[2];
rz(1.3806512) q[3];
sx q[3];
rz(-1.645429) q[3];
sx q[3];
rz(-2.7281705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8274902) q[2];
sx q[2];
rz(-0.9919439) q[2];
sx q[2];
rz(-2.6143383) q[2];
rz(-0.54316795) q[3];
sx q[3];
rz(-2.4542464) q[3];
sx q[3];
rz(0.34272042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0813893) q[0];
sx q[0];
rz(-2.9855766) q[0];
sx q[0];
rz(-0.40670893) q[0];
rz(-1.9806865) q[1];
sx q[1];
rz(-2.2229767) q[1];
sx q[1];
rz(-1.0103753) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5766616) q[0];
sx q[0];
rz(-1.0258249) q[0];
sx q[0];
rz(-3.0002563) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7356304) q[2];
sx q[2];
rz(-0.34892198) q[2];
sx q[2];
rz(-0.59970784) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1045253) q[1];
sx q[1];
rz(-0.43018331) q[1];
sx q[1];
rz(2.6135315) q[1];
x q[2];
rz(2.3361037) q[3];
sx q[3];
rz(-1.7542766) q[3];
sx q[3];
rz(-3.0940521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0668209) q[2];
sx q[2];
rz(-1.2819042) q[2];
sx q[2];
rz(1.034896) q[2];
rz(-1.3567989) q[3];
sx q[3];
rz(-0.59485888) q[3];
sx q[3];
rz(-0.6234197) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9676301) q[0];
sx q[0];
rz(-1.5246464) q[0];
sx q[0];
rz(0.3279283) q[0];
rz(-0.11218849) q[1];
sx q[1];
rz(-1.9393549) q[1];
sx q[1];
rz(0.97253886) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5060726) q[0];
sx q[0];
rz(-1.0131665) q[0];
sx q[0];
rz(-2.7235759) q[0];
rz(-pi) q[1];
rz(-2.8424524) q[2];
sx q[2];
rz(-1.1465596) q[2];
sx q[2];
rz(1.5787293) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.211208) q[1];
sx q[1];
rz(-2.5519785) q[1];
sx q[1];
rz(-3.1077887) q[1];
x q[2];
rz(0.48613207) q[3];
sx q[3];
rz(-0.68607578) q[3];
sx q[3];
rz(0.11763517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5358676) q[2];
sx q[2];
rz(-0.87656993) q[2];
sx q[2];
rz(-0.65315872) q[2];
rz(0.43290916) q[3];
sx q[3];
rz(-0.50986367) q[3];
sx q[3];
rz(2.9292817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2046278) q[0];
sx q[0];
rz(-0.84119868) q[0];
sx q[0];
rz(-2.9168108) q[0];
rz(-2.3460491) q[1];
sx q[1];
rz(-1.8276151) q[1];
sx q[1];
rz(1.068211) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58086181) q[0];
sx q[0];
rz(-0.95968548) q[0];
sx q[0];
rz(-2.1944502) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3722367) q[2];
sx q[2];
rz(-1.1072259) q[2];
sx q[2];
rz(-3.0887512) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.873899) q[1];
sx q[1];
rz(-1.5452478) q[1];
sx q[1];
rz(2.4539095) q[1];
rz(-pi) q[2];
rz(-1.3490476) q[3];
sx q[3];
rz(-1.6283247) q[3];
sx q[3];
rz(1.8309995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.75605679) q[2];
sx q[2];
rz(-1.3827366) q[2];
sx q[2];
rz(2.5059911) q[2];
rz(-0.33291891) q[3];
sx q[3];
rz(-1.0115441) q[3];
sx q[3];
rz(-0.46271589) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3391089) q[0];
sx q[0];
rz(-0.026263069) q[0];
sx q[0];
rz(3.0185757) q[0];
rz(1.7440375) q[1];
sx q[1];
rz(-1.3879644) q[1];
sx q[1];
rz(2.3618598) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62209475) q[0];
sx q[0];
rz(-2.4197398) q[0];
sx q[0];
rz(-1.9869542) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9626368) q[2];
sx q[2];
rz(-0.6935941) q[2];
sx q[2];
rz(0.041698448) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.18786622) q[1];
sx q[1];
rz(-2.7147016) q[1];
sx q[1];
rz(-2.0679451) q[1];
x q[2];
rz(1.6980096) q[3];
sx q[3];
rz(-1.4641342) q[3];
sx q[3];
rz(2.9677109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66703779) q[2];
sx q[2];
rz(-1.2660618) q[2];
sx q[2];
rz(3.0511268) q[2];
rz(-0.41680923) q[3];
sx q[3];
rz(-2.3495245) q[3];
sx q[3];
rz(0.99378234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4823293) q[0];
sx q[0];
rz(-1.8195131) q[0];
sx q[0];
rz(0.089476712) q[0];
rz(-2.3327475) q[1];
sx q[1];
rz(-2.6409812) q[1];
sx q[1];
rz(-2.083875) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9442975) q[0];
sx q[0];
rz(-2.5677201) q[0];
sx q[0];
rz(-0.7446592) q[0];
x q[1];
rz(-1.5616778) q[2];
sx q[2];
rz(-2.5895725) q[2];
sx q[2];
rz(-2.7087351) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48185086) q[1];
sx q[1];
rz(-2.5999477) q[1];
sx q[1];
rz(-2.2996344) q[1];
rz(-2.3848885) q[3];
sx q[3];
rz(-1.5528402) q[3];
sx q[3];
rz(-1.1566263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.40314254) q[2];
sx q[2];
rz(-2.6046627) q[2];
sx q[2];
rz(-2.7682847) q[2];
rz(-2.8849854) q[3];
sx q[3];
rz(-1.4213057) q[3];
sx q[3];
rz(3.0671425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8858717) q[0];
sx q[0];
rz(-1.5789565) q[0];
sx q[0];
rz(-1.4174905) q[0];
rz(-1.7120842) q[1];
sx q[1];
rz(-0.34965873) q[1];
sx q[1];
rz(0.8082334) q[1];
rz(1.5803554) q[2];
sx q[2];
rz(-2.0295967) q[2];
sx q[2];
rz(1.5461736) q[2];
rz(0.96961602) q[3];
sx q[3];
rz(-1.6027228) q[3];
sx q[3];
rz(3.0016196) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
