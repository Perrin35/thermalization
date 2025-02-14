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
rz(2.8494868) q[0];
sx q[0];
rz(-2.5295244) q[0];
sx q[0];
rz(3.0906313) q[0];
rz(1.2904957) q[1];
sx q[1];
rz(-1.2761071) q[1];
sx q[1];
rz(-0.84604231) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.030662) q[0];
sx q[0];
rz(-1.4101068) q[0];
sx q[0];
rz(2.7881505) q[0];
x q[1];
rz(2.8522367) q[2];
sx q[2];
rz(-1.932992) q[2];
sx q[2];
rz(-1.2928499) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.5979833) q[1];
sx q[1];
rz(-0.56386891) q[1];
sx q[1];
rz(2.6257674) q[1];
x q[2];
rz(3.0676431) q[3];
sx q[3];
rz(-3.0294501) q[3];
sx q[3];
rz(-2.6611947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.35752615) q[2];
sx q[2];
rz(-0.94782031) q[2];
sx q[2];
rz(2.2785462) q[2];
rz(2.0283902) q[3];
sx q[3];
rz(-0.92156711) q[3];
sx q[3];
rz(-0.79871261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6031826) q[0];
sx q[0];
rz(-0.69284678) q[0];
sx q[0];
rz(-2.8811654) q[0];
rz(2.8159091) q[1];
sx q[1];
rz(-2.1574056) q[1];
sx q[1];
rz(0.30776417) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1528252) q[0];
sx q[0];
rz(-0.77251311) q[0];
sx q[0];
rz(1.1940184) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.4701647) q[2];
sx q[2];
rz(-1.2855296) q[2];
sx q[2];
rz(2.4131052) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5883397) q[1];
sx q[1];
rz(-0.25811895) q[1];
sx q[1];
rz(1.939276) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.006853718) q[3];
sx q[3];
rz(-2.2061976) q[3];
sx q[3];
rz(1.8346661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8327568) q[2];
sx q[2];
rz(-2.2726161) q[2];
sx q[2];
rz(2.5858509) q[2];
rz(2.6293788) q[3];
sx q[3];
rz(-0.54849505) q[3];
sx q[3];
rz(0.51928025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.888805) q[0];
sx q[0];
rz(-0.35211173) q[0];
sx q[0];
rz(0.74932253) q[0];
rz(-1.0961756) q[1];
sx q[1];
rz(-1.7671894) q[1];
sx q[1];
rz(0.38415092) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6836943) q[0];
sx q[0];
rz(-1.3078469) q[0];
sx q[0];
rz(2.5431741) q[0];
x q[1];
rz(-1.3170502) q[2];
sx q[2];
rz(-0.98461005) q[2];
sx q[2];
rz(1.3519999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0047915) q[1];
sx q[1];
rz(-2.8344732) q[1];
sx q[1];
rz(-1.7476487) q[1];
x q[2];
rz(1.2189193) q[3];
sx q[3];
rz(-2.4259544) q[3];
sx q[3];
rz(1.1688237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.38203794) q[2];
sx q[2];
rz(-1.3052772) q[2];
sx q[2];
rz(1.6273512) q[2];
rz(-0.73005992) q[3];
sx q[3];
rz(-0.74159139) q[3];
sx q[3];
rz(-0.39456427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6877614) q[0];
sx q[0];
rz(-1.4241968) q[0];
sx q[0];
rz(-0.060039595) q[0];
rz(-2.2278348) q[1];
sx q[1];
rz(-0.91122952) q[1];
sx q[1];
rz(2.2027016) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37591463) q[0];
sx q[0];
rz(-1.5374827) q[0];
sx q[0];
rz(-3.0785962) q[0];
rz(-pi) q[1];
rz(-0.36146116) q[2];
sx q[2];
rz(-0.088080125) q[2];
sx q[2];
rz(2.9437906) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9524433) q[1];
sx q[1];
rz(-1.5464142) q[1];
sx q[1];
rz(-1.0827176) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2430211) q[3];
sx q[3];
rz(-2.13604) q[3];
sx q[3];
rz(-0.94796255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.011220304) q[2];
sx q[2];
rz(-2.8488686) q[2];
sx q[2];
rz(2.4893153) q[2];
rz(1.8782714) q[3];
sx q[3];
rz(-1.616856) q[3];
sx q[3];
rz(1.1711082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0703099) q[0];
sx q[0];
rz(-1.2391397) q[0];
sx q[0];
rz(0.5624482) q[0];
rz(1.3485472) q[1];
sx q[1];
rz(-2.8548073) q[1];
sx q[1];
rz(2.5602692) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7331478) q[0];
sx q[0];
rz(-0.83400351) q[0];
sx q[0];
rz(2.187378) q[0];
rz(-pi) q[1];
rz(-2.6810886) q[2];
sx q[2];
rz(-1.6031613) q[2];
sx q[2];
rz(0.89214395) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6186468) q[1];
sx q[1];
rz(-0.60877555) q[1];
sx q[1];
rz(2.2227915) q[1];
rz(-pi) q[2];
rz(-1.6147037) q[3];
sx q[3];
rz(-0.9287408) q[3];
sx q[3];
rz(1.2894443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.44744667) q[2];
sx q[2];
rz(-1.4299102) q[2];
sx q[2];
rz(-0.68667975) q[2];
rz(-0.40003362) q[3];
sx q[3];
rz(-1.8828078) q[3];
sx q[3];
rz(3.1299642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7069063) q[0];
sx q[0];
rz(-0.97395581) q[0];
sx q[0];
rz(-0.45027012) q[0];
rz(-0.88092583) q[1];
sx q[1];
rz(-1.2342781) q[1];
sx q[1];
rz(2.0251822) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8378431) q[0];
sx q[0];
rz(-2.3275796) q[0];
sx q[0];
rz(0.26153691) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7911602) q[2];
sx q[2];
rz(-1.3964126) q[2];
sx q[2];
rz(-0.36878219) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9444869) q[1];
sx q[1];
rz(-0.61054269) q[1];
sx q[1];
rz(-1.6050299) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93706705) q[3];
sx q[3];
rz(-2.6014199) q[3];
sx q[3];
rz(2.7989509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.66699666) q[2];
sx q[2];
rz(-0.23385364) q[2];
sx q[2];
rz(-0.65143603) q[2];
rz(0.78178072) q[3];
sx q[3];
rz(-1.4835446) q[3];
sx q[3];
rz(0.93530161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5485789) q[0];
sx q[0];
rz(-0.23603708) q[0];
sx q[0];
rz(-2.6759942) q[0];
rz(1.9904526) q[1];
sx q[1];
rz(-1.2460037) q[1];
sx q[1];
rz(1.3394855) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5440774) q[0];
sx q[0];
rz(-1.1664412) q[0];
sx q[0];
rz(-0.57283516) q[0];
rz(0.8546245) q[2];
sx q[2];
rz(-0.058283866) q[2];
sx q[2];
rz(-0.56266498) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3834584) q[1];
sx q[1];
rz(-1.3515359) q[1];
sx q[1];
rz(0.077139826) q[1];
rz(-pi) q[2];
rz(-2.7202456) q[3];
sx q[3];
rz(-2.192492) q[3];
sx q[3];
rz(2.5558215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1107948) q[2];
sx q[2];
rz(-1.75533) q[2];
sx q[2];
rz(0.01290713) q[2];
rz(3.0068093) q[3];
sx q[3];
rz(-1.2603935) q[3];
sx q[3];
rz(-0.24136647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6808692) q[0];
sx q[0];
rz(-0.091982059) q[0];
sx q[0];
rz(1.9665834) q[0];
rz(-2.3911632) q[1];
sx q[1];
rz(-1.3576327) q[1];
sx q[1];
rz(-2.7480385) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9703044) q[0];
sx q[0];
rz(-1.3960724) q[0];
sx q[0];
rz(0.98002429) q[0];
rz(-pi) q[1];
rz(1.6874763) q[2];
sx q[2];
rz(-0.91800729) q[2];
sx q[2];
rz(1.0972925) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8361349) q[1];
sx q[1];
rz(-1.4609481) q[1];
sx q[1];
rz(1.1396465) q[1];
rz(-pi) q[2];
rz(-2.2116488) q[3];
sx q[3];
rz(-2.5287147) q[3];
sx q[3];
rz(-1.5811435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.63341081) q[2];
sx q[2];
rz(-2.0076624) q[2];
sx q[2];
rz(-0.22937648) q[2];
rz(3.0502099) q[3];
sx q[3];
rz(-1.6978076) q[3];
sx q[3];
rz(-2.5858333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.9230187) q[0];
sx q[0];
rz(-2.3526683) q[0];
sx q[0];
rz(2.5564585) q[0];
rz(-1.8642037) q[1];
sx q[1];
rz(-1.2150512) q[1];
sx q[1];
rz(-0.1543943) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47772022) q[0];
sx q[0];
rz(-1.1363875) q[0];
sx q[0];
rz(1.4503195) q[0];
x q[1];
rz(-2.7150776) q[2];
sx q[2];
rz(-1.8898481) q[2];
sx q[2];
rz(-1.4501377) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8081863) q[1];
sx q[1];
rz(-1.8318614) q[1];
sx q[1];
rz(-0.33269791) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.102907) q[3];
sx q[3];
rz(-1.7280735) q[3];
sx q[3];
rz(0.09889557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.020891) q[2];
sx q[2];
rz(-2.2521844) q[2];
sx q[2];
rz(0.54023877) q[2];
rz(-2.7464416) q[3];
sx q[3];
rz(-2.4873147) q[3];
sx q[3];
rz(-1.7925709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2111135) q[0];
sx q[0];
rz(-1.2808639) q[0];
sx q[0];
rz(-1.6092009) q[0];
rz(-2.2156175) q[1];
sx q[1];
rz(-1.4264868) q[1];
sx q[1];
rz(1.6645974) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4173399) q[0];
sx q[0];
rz(-0.84463929) q[0];
sx q[0];
rz(-1.5717616) q[0];
rz(2.8468644) q[2];
sx q[2];
rz(-1.2090313) q[2];
sx q[2];
rz(0.28361646) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1077233) q[1];
sx q[1];
rz(-1.0230015) q[1];
sx q[1];
rz(-2.7366927) q[1];
x q[2];
rz(-0.98201237) q[3];
sx q[3];
rz(-1.1945621) q[3];
sx q[3];
rz(-0.8134977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.28107873) q[2];
sx q[2];
rz(-1.8033359) q[2];
sx q[2];
rz(0.64858428) q[2];
rz(0.72566882) q[3];
sx q[3];
rz(-2.9026493) q[3];
sx q[3];
rz(1.749595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7657179) q[0];
sx q[0];
rz(-0.97828843) q[0];
sx q[0];
rz(0.82213415) q[0];
rz(1.0646461) q[1];
sx q[1];
rz(-1.6864265) q[1];
sx q[1];
rz(2.0815157) q[1];
rz(-2.9072472) q[2];
sx q[2];
rz(-0.6498944) q[2];
sx q[2];
rz(1.9882974) q[2];
rz(-2.1172932) q[3];
sx q[3];
rz(-2.306675) q[3];
sx q[3];
rz(-2.6425895) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
