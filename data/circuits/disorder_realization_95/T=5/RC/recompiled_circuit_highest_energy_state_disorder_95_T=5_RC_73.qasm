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
rz(1.3746102) q[0];
sx q[0];
rz(-2.1687825) q[0];
sx q[0];
rz(1.2306124) q[0];
rz(1.5988916) q[1];
sx q[1];
rz(-2.7886432) q[1];
sx q[1];
rz(-0.43500873) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15947882) q[0];
sx q[0];
rz(-1.6702939) q[0];
sx q[0];
rz(1.9051108) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6224766) q[2];
sx q[2];
rz(-0.26072219) q[2];
sx q[2];
rz(-0.97445011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.75499308) q[1];
sx q[1];
rz(-2.3191929) q[1];
sx q[1];
rz(1.1851285) q[1];
rz(-0.013641274) q[3];
sx q[3];
rz(-1.5392) q[3];
sx q[3];
rz(-2.7073366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.9445442) q[2];
sx q[2];
rz(-0.41434449) q[2];
sx q[2];
rz(2.22331) q[2];
rz(0.36327547) q[3];
sx q[3];
rz(-1.2493635) q[3];
sx q[3];
rz(1.8173119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0927703) q[0];
sx q[0];
rz(-0.73960441) q[0];
sx q[0];
rz(1.1356461) q[0];
rz(2.8593235) q[1];
sx q[1];
rz(-1.6259364) q[1];
sx q[1];
rz(0.20271066) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7736063) q[0];
sx q[0];
rz(-1.3895036) q[0];
sx q[0];
rz(-1.3884991) q[0];
rz(-pi) q[1];
rz(-1.545867) q[2];
sx q[2];
rz(-0.48183107) q[2];
sx q[2];
rz(0.78311759) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.91071002) q[1];
sx q[1];
rz(-2.6697914) q[1];
sx q[1];
rz(-2.8698026) q[1];
rz(1.1003157) q[3];
sx q[3];
rz(-3.1300448) q[3];
sx q[3];
rz(1.601416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.89381605) q[2];
sx q[2];
rz(-2.717369) q[2];
sx q[2];
rz(-0.068537863) q[2];
rz(-2.2872772) q[3];
sx q[3];
rz(-1.2645384) q[3];
sx q[3];
rz(-3.1356649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0900367) q[0];
sx q[0];
rz(-0.55110252) q[0];
sx q[0];
rz(-1.0968444) q[0];
rz(-0.54496533) q[1];
sx q[1];
rz(-2.4577591) q[1];
sx q[1];
rz(-0.1942689) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0476463) q[0];
sx q[0];
rz(-1.1585278) q[0];
sx q[0];
rz(0.38239371) q[0];
rz(-pi) q[1];
rz(-1.0610007) q[2];
sx q[2];
rz(-0.61697996) q[2];
sx q[2];
rz(-1.7327019) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.771505) q[1];
sx q[1];
rz(-1.2125101) q[1];
sx q[1];
rz(0.53046988) q[1];
rz(-pi) q[2];
rz(-1.9241289) q[3];
sx q[3];
rz(-1.4446745) q[3];
sx q[3];
rz(-0.68870163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37476173) q[2];
sx q[2];
rz(-1.7685879) q[2];
sx q[2];
rz(2.4135446) q[2];
rz(0.22045615) q[3];
sx q[3];
rz(-2.9900592) q[3];
sx q[3];
rz(-0.63173405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24882889) q[0];
sx q[0];
rz(-1.4735824) q[0];
sx q[0];
rz(-2.072075) q[0];
rz(-0.88691521) q[1];
sx q[1];
rz(-2.6353757) q[1];
sx q[1];
rz(-1.3217529) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76218433) q[0];
sx q[0];
rz(-0.43566049) q[0];
sx q[0];
rz(3.0489151) q[0];
x q[1];
rz(-0.23897929) q[2];
sx q[2];
rz(-1.5859446) q[2];
sx q[2];
rz(-2.3292484) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9001213) q[1];
sx q[1];
rz(-1.6423499) q[1];
sx q[1];
rz(2.009719) q[1];
rz(-pi) q[2];
rz(-0.87476829) q[3];
sx q[3];
rz(-2.416055) q[3];
sx q[3];
rz(1.0410795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6679823) q[2];
sx q[2];
rz(-1.0429635) q[2];
sx q[2];
rz(1.854151) q[2];
rz(2.1227664) q[3];
sx q[3];
rz(-0.5580709) q[3];
sx q[3];
rz(-2.4655925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92893112) q[0];
sx q[0];
rz(-2.9086845) q[0];
sx q[0];
rz(0.88053298) q[0];
rz(-0.13678837) q[1];
sx q[1];
rz(-0.22577481) q[1];
sx q[1];
rz(0.61019623) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5870415) q[0];
sx q[0];
rz(-2.4992538) q[0];
sx q[0];
rz(1.5917529) q[0];
x q[1];
rz(-2.8183646) q[2];
sx q[2];
rz(-0.45219496) q[2];
sx q[2];
rz(-1.5791073) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5194441) q[1];
sx q[1];
rz(-0.3560027) q[1];
sx q[1];
rz(0.58575969) q[1];
rz(2.7022252) q[3];
sx q[3];
rz(-1.4476848) q[3];
sx q[3];
rz(-1.1151552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5337164) q[2];
sx q[2];
rz(-0.72489649) q[2];
sx q[2];
rz(1.9113212) q[2];
rz(2.2599334) q[3];
sx q[3];
rz(-1.6917546) q[3];
sx q[3];
rz(-1.3151883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45516685) q[0];
sx q[0];
rz(-1.5436341) q[0];
sx q[0];
rz(2.0879188) q[0];
rz(1.3764489) q[1];
sx q[1];
rz(-2.0339298) q[1];
sx q[1];
rz(-2.460316) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2864425) q[0];
sx q[0];
rz(-1.5871829) q[0];
sx q[0];
rz(-0.70006242) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7679484) q[2];
sx q[2];
rz(-2.3609997) q[2];
sx q[2];
rz(-0.045750387) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44418844) q[1];
sx q[1];
rz(-2.5696074) q[1];
sx q[1];
rz(2.0818656) q[1];
rz(-3.104171) q[3];
sx q[3];
rz(-2.2624514) q[3];
sx q[3];
rz(-3.0841923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.072558746) q[2];
sx q[2];
rz(-0.37798887) q[2];
sx q[2];
rz(-2.6663713) q[2];
rz(1.2230988) q[3];
sx q[3];
rz(-2.1971072) q[3];
sx q[3];
rz(0.16330115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.908919) q[0];
sx q[0];
rz(-0.95218807) q[0];
sx q[0];
rz(2.9222144) q[0];
rz(2.1351922) q[1];
sx q[1];
rz(-1.718113) q[1];
sx q[1];
rz(-1.9361608) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4426098) q[0];
sx q[0];
rz(-0.97215334) q[0];
sx q[0];
rz(3.1151616) q[0];
rz(-2.5616165) q[2];
sx q[2];
rz(-2.6784228) q[2];
sx q[2];
rz(1.620174) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.58118966) q[1];
sx q[1];
rz(-2.9415574) q[1];
sx q[1];
rz(-0.086189857) q[1];
rz(-pi) q[2];
rz(3.0053776) q[3];
sx q[3];
rz(-2.0117617) q[3];
sx q[3];
rz(-1.3974401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.22126234) q[2];
sx q[2];
rz(-0.39322501) q[2];
sx q[2];
rz(0.66366759) q[2];
rz(-2.4237733) q[3];
sx q[3];
rz(-0.65687537) q[3];
sx q[3];
rz(2.9653911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32888907) q[0];
sx q[0];
rz(-2.1825574) q[0];
sx q[0];
rz(1.2411728) q[0];
rz(-0.42789704) q[1];
sx q[1];
rz(-1.5792081) q[1];
sx q[1];
rz(1.4235628) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54518979) q[0];
sx q[0];
rz(-1.3529142) q[0];
sx q[0];
rz(-1.6473739) q[0];
x q[1];
rz(-2.9828885) q[2];
sx q[2];
rz(-2.2825512) q[2];
sx q[2];
rz(-1.2325328) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7562508) q[1];
sx q[1];
rz(-1.5495117) q[1];
sx q[1];
rz(0.90275928) q[1];
rz(-pi) q[2];
rz(0.10382048) q[3];
sx q[3];
rz(-2.8049433) q[3];
sx q[3];
rz(2.3262084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.9972035) q[2];
sx q[2];
rz(-2.453697) q[2];
sx q[2];
rz(-0.29590657) q[2];
rz(-2.0603254) q[3];
sx q[3];
rz(-1.5238949) q[3];
sx q[3];
rz(-2.9937939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57701552) q[0];
sx q[0];
rz(-2.5595589) q[0];
sx q[0];
rz(-0.70593315) q[0];
rz(2.9544746) q[1];
sx q[1];
rz(-0.83693224) q[1];
sx q[1];
rz(-0.44609889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69639403) q[0];
sx q[0];
rz(-2.9893251) q[0];
sx q[0];
rz(3.1294786) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1438333) q[2];
sx q[2];
rz(-2.5964542) q[2];
sx q[2];
rz(1.0447431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46221387) q[1];
sx q[1];
rz(-1.0125547) q[1];
sx q[1];
rz(0.39887303) q[1];
x q[2];
rz(-0.25983747) q[3];
sx q[3];
rz(-2.5492955) q[3];
sx q[3];
rz(-3.1282476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9624761) q[2];
sx q[2];
rz(-2.1160782) q[2];
sx q[2];
rz(-1.7402488) q[2];
rz(2.5380747) q[3];
sx q[3];
rz(-2.9602435) q[3];
sx q[3];
rz(3.0028499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03183455) q[0];
sx q[0];
rz(-1.7367481) q[0];
sx q[0];
rz(3.0314714) q[0];
rz(1.8203863) q[1];
sx q[1];
rz(-0.92480129) q[1];
sx q[1];
rz(0.22470156) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8256841) q[0];
sx q[0];
rz(-1.0498739) q[0];
sx q[0];
rz(0.45895001) q[0];
x q[1];
rz(-0.16188936) q[2];
sx q[2];
rz(-1.2873374) q[2];
sx q[2];
rz(-3.0155011) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4555295) q[1];
sx q[1];
rz(-1.7357329) q[1];
sx q[1];
rz(0.73442028) q[1];
x q[2];
rz(-0.026894491) q[3];
sx q[3];
rz(-1.3984612) q[3];
sx q[3];
rz(0.89629345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3118185) q[2];
sx q[2];
rz(-2.8303787) q[2];
sx q[2];
rz(-2.6149926) q[2];
rz(0.38463587) q[3];
sx q[3];
rz(-1.8943818) q[3];
sx q[3];
rz(2.9068936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2100691) q[0];
sx q[0];
rz(-0.43545224) q[0];
sx q[0];
rz(-0.42238105) q[0];
rz(2.3474563) q[1];
sx q[1];
rz(-1.7105449) q[1];
sx q[1];
rz(1.7078043) q[1];
rz(0.83195291) q[2];
sx q[2];
rz(-0.83870287) q[2];
sx q[2];
rz(2.4626682) q[2];
rz(-0.41153367) q[3];
sx q[3];
rz(-1.2306662) q[3];
sx q[3];
rz(-1.6231404) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
