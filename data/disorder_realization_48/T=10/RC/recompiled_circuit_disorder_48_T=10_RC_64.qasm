OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.7735908) q[0];
sx q[0];
rz(3.9324023) q[0];
sx q[0];
rz(12.232236) q[0];
rz(2.6842527) q[1];
sx q[1];
rz(-2.1973124) q[1];
sx q[1];
rz(-1.2184719) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63872913) q[0];
sx q[0];
rz(-1.6067061) q[0];
sx q[0];
rz(-2.3907651) q[0];
rz(-pi) q[1];
rz(-1.5300418) q[2];
sx q[2];
rz(-1.1062804) q[2];
sx q[2];
rz(1.0711311) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7114746) q[1];
sx q[1];
rz(-0.700044) q[1];
sx q[1];
rz(-2.3858566) q[1];
x q[2];
rz(-0.18913194) q[3];
sx q[3];
rz(-1.7889708) q[3];
sx q[3];
rz(1.3355586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5228287) q[2];
sx q[2];
rz(-0.4814119) q[2];
sx q[2];
rz(2.5640326) q[2];
rz(1.9918359) q[3];
sx q[3];
rz(-1.7532319) q[3];
sx q[3];
rz(-2.4770588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98786551) q[0];
sx q[0];
rz(-0.58652121) q[0];
sx q[0];
rz(-2.7541449) q[0];
rz(-0.93915835) q[1];
sx q[1];
rz(-0.99717957) q[1];
sx q[1];
rz(-1.4025677) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1419066) q[0];
sx q[0];
rz(-1.5667331) q[0];
sx q[0];
rz(-0.02948972) q[0];
x q[1];
rz(1.7895133) q[2];
sx q[2];
rz(-2.1195076) q[2];
sx q[2];
rz(0.95900853) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1837511) q[1];
sx q[1];
rz(-1.1693923) q[1];
sx q[1];
rz(-0.5387696) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6658695) q[3];
sx q[3];
rz(-2.2007696) q[3];
sx q[3];
rz(0.41439393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42276057) q[2];
sx q[2];
rz(-1.7451124) q[2];
sx q[2];
rz(0.31769162) q[2];
rz(2.9348532) q[3];
sx q[3];
rz(-2.5419149) q[3];
sx q[3];
rz(0.81682214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3690255) q[0];
sx q[0];
rz(-1.7018397) q[0];
sx q[0];
rz(-1.7279708) q[0];
rz(0.47779045) q[1];
sx q[1];
rz(-1.7910035) q[1];
sx q[1];
rz(2.7405222) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0862335) q[0];
sx q[0];
rz(-1.2463352) q[0];
sx q[0];
rz(-0.1089036) q[0];
rz(1.7602073) q[2];
sx q[2];
rz(-2.2064798) q[2];
sx q[2];
rz(0.96178255) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8071825) q[1];
sx q[1];
rz(-1.5644329) q[1];
sx q[1];
rz(0.73015726) q[1];
x q[2];
rz(2.0173666) q[3];
sx q[3];
rz(-2.7105769) q[3];
sx q[3];
rz(-0.22180804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59427375) q[2];
sx q[2];
rz(-1.5218364) q[2];
sx q[2];
rz(-2.5857914) q[2];
rz(0.9764955) q[3];
sx q[3];
rz(-2.591811) q[3];
sx q[3];
rz(-2.3613789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34898409) q[0];
sx q[0];
rz(-1.6225092) q[0];
sx q[0];
rz(1.4439616) q[0];
rz(-1.5199039) q[1];
sx q[1];
rz(-0.65602055) q[1];
sx q[1];
rz(2.8881883) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93638203) q[0];
sx q[0];
rz(-2.10996) q[0];
sx q[0];
rz(-3.0217231) q[0];
x q[1];
rz(2.0565815) q[2];
sx q[2];
rz(-1.5554785) q[2];
sx q[2];
rz(2.3682396) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.036382347) q[1];
sx q[1];
rz(-1.6670334) q[1];
sx q[1];
rz(-2.5114245) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0074189) q[3];
sx q[3];
rz(-0.65376702) q[3];
sx q[3];
rz(1.3282446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.110934) q[2];
sx q[2];
rz(-1.7548283) q[2];
sx q[2];
rz(0.78732642) q[2];
rz(-0.91286719) q[3];
sx q[3];
rz(-2.392231) q[3];
sx q[3];
rz(-2.1319938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71516365) q[0];
sx q[0];
rz(-0.63868317) q[0];
sx q[0];
rz(0.062967904) q[0];
rz(0.12403034) q[1];
sx q[1];
rz(-2.3359559) q[1];
sx q[1];
rz(-2.6834992) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9034018) q[0];
sx q[0];
rz(-1.5684677) q[0];
sx q[0];
rz(0.44263126) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45256726) q[2];
sx q[2];
rz(-2.1308225) q[2];
sx q[2];
rz(2.3134311) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8372247) q[1];
sx q[1];
rz(-2.3886884) q[1];
sx q[1];
rz(1.7423082) q[1];
x q[2];
rz(-0.094893806) q[3];
sx q[3];
rz(-1.0720383) q[3];
sx q[3];
rz(0.14548485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9690341) q[2];
sx q[2];
rz(-0.92323747) q[2];
sx q[2];
rz(-2.9120973) q[2];
rz(-0.0028006639) q[3];
sx q[3];
rz(-0.87001785) q[3];
sx q[3];
rz(-1.8026479) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5795508) q[0];
sx q[0];
rz(-2.8631449) q[0];
sx q[0];
rz(0.91947412) q[0];
rz(0.062285034) q[1];
sx q[1];
rz(-1.0039763) q[1];
sx q[1];
rz(-1.2671635) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7125268) q[0];
sx q[0];
rz(-2.292284) q[0];
sx q[0];
rz(-2.6119786) q[0];
rz(-pi) q[1];
rz(-2.3976372) q[2];
sx q[2];
rz(-1.4641054) q[2];
sx q[2];
rz(-0.26847408) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8333203) q[1];
sx q[1];
rz(-0.38494021) q[1];
sx q[1];
rz(-1.0233378) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0165765) q[3];
sx q[3];
rz(-2.6187754) q[3];
sx q[3];
rz(-1.0558053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.9617812) q[2];
sx q[2];
rz(-0.87819019) q[2];
sx q[2];
rz(-0.58376694) q[2];
rz(-0.70872712) q[3];
sx q[3];
rz(-1.3112336) q[3];
sx q[3];
rz(-0.023199737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0631183) q[0];
sx q[0];
rz(-0.45409504) q[0];
sx q[0];
rz(-2.069058) q[0];
rz(-2.5947) q[1];
sx q[1];
rz(-1.2439589) q[1];
sx q[1];
rz(-2.0297208) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3355013) q[0];
sx q[0];
rz(-1.6643235) q[0];
sx q[0];
rz(-2.0557687) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56356168) q[2];
sx q[2];
rz(-0.60699082) q[2];
sx q[2];
rz(2.3129472) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6560981) q[1];
sx q[1];
rz(-1.6472677) q[1];
sx q[1];
rz(-0.2904201) q[1];
rz(-2.2124347) q[3];
sx q[3];
rz(-2.540179) q[3];
sx q[3];
rz(-1.8638368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.83773461) q[2];
sx q[2];
rz(-1.6656817) q[2];
sx q[2];
rz(-1.1676577) q[2];
rz(-1.6052823) q[3];
sx q[3];
rz(-1.6746018) q[3];
sx q[3];
rz(-1.7355828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7325571) q[0];
sx q[0];
rz(-1.5719825) q[0];
sx q[0];
rz(2.4107966) q[0];
rz(-0.90019512) q[1];
sx q[1];
rz(-2.3370445) q[1];
sx q[1];
rz(-2.3866167) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2728111) q[0];
sx q[0];
rz(-2.4065354) q[0];
sx q[0];
rz(1.3520157) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63379143) q[2];
sx q[2];
rz(-0.83665028) q[2];
sx q[2];
rz(-1.447669) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6983812) q[1];
sx q[1];
rz(-1.4581919) q[1];
sx q[1];
rz(2.5618782) q[1];
rz(2.6492277) q[3];
sx q[3];
rz(-1.0400606) q[3];
sx q[3];
rz(-2.9363971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76453152) q[2];
sx q[2];
rz(-1.7691282) q[2];
sx q[2];
rz(-0.63684741) q[2];
rz(2.8751255) q[3];
sx q[3];
rz(-2.0839432) q[3];
sx q[3];
rz(1.586097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72717845) q[0];
sx q[0];
rz(-1.1295015) q[0];
sx q[0];
rz(1.138858) q[0];
rz(0.75421929) q[1];
sx q[1];
rz(-2.8051839) q[1];
sx q[1];
rz(3.1220904) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50197983) q[0];
sx q[0];
rz(-1.4907955) q[0];
sx q[0];
rz(1.3638858) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3703913) q[2];
sx q[2];
rz(-1.7192171) q[2];
sx q[2];
rz(2.5401149) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1972678) q[1];
sx q[1];
rz(-1.7515344) q[1];
sx q[1];
rz(0.9428057) q[1];
x q[2];
rz(1.8832302) q[3];
sx q[3];
rz(-1.7282681) q[3];
sx q[3];
rz(1.0851932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1372244) q[2];
sx q[2];
rz(-1.7251816) q[2];
sx q[2];
rz(2.2926889) q[2];
rz(2.7539339) q[3];
sx q[3];
rz(-2.013423) q[3];
sx q[3];
rz(1.5415812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7983109) q[0];
sx q[0];
rz(-2.9738975) q[0];
sx q[0];
rz(0.48450255) q[0];
rz(-1.7548521) q[1];
sx q[1];
rz(-1.4258899) q[1];
sx q[1];
rz(1.9932995) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0120221) q[0];
sx q[0];
rz(-2.9992933) q[0];
sx q[0];
rz(-0.15108959) q[0];
x q[1];
rz(2.8414367) q[2];
sx q[2];
rz(-2.488392) q[2];
sx q[2];
rz(0.73210994) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9705829) q[1];
sx q[1];
rz(-2.2323771) q[1];
sx q[1];
rz(-0.2172825) q[1];
rz(-pi) q[2];
rz(-2.9989472) q[3];
sx q[3];
rz(-2.7890165) q[3];
sx q[3];
rz(0.25776097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.91223532) q[2];
sx q[2];
rz(-1.8476013) q[2];
sx q[2];
rz(-2.7764376) q[2];
rz(0.12864104) q[3];
sx q[3];
rz(-1.9059076) q[3];
sx q[3];
rz(-0.45583367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1098332) q[0];
sx q[0];
rz(-0.84072996) q[0];
sx q[0];
rz(1.6054556) q[0];
rz(-0.96314349) q[1];
sx q[1];
rz(-1.2711202) q[1];
sx q[1];
rz(-1.0585379) q[1];
rz(1.9671494) q[2];
sx q[2];
rz(-0.94559961) q[2];
sx q[2];
rz(-0.0018975817) q[2];
rz(2.1309489) q[3];
sx q[3];
rz(-1.602136) q[3];
sx q[3];
rz(-1.7020561) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
