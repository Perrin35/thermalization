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
rz(1.6744094) q[0];
sx q[0];
rz(-1.5404584) q[0];
sx q[0];
rz(1.2920446) q[0];
rz(-1.0979103) q[1];
sx q[1];
rz(3.8976964) q[1];
sx q[1];
rz(10.1367) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5789509) q[0];
sx q[0];
rz(-1.6058086) q[0];
sx q[0];
rz(-2.9981311) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.759077) q[2];
sx q[2];
rz(-1.7393987) q[2];
sx q[2];
rz(1.027866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.10140534) q[1];
sx q[1];
rz(-0.95164548) q[1];
sx q[1];
rz(1.2907486) q[1];
rz(-pi) q[2];
rz(2.6308832) q[3];
sx q[3];
rz(-2.1132971) q[3];
sx q[3];
rz(-1.6245019) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0277675) q[2];
sx q[2];
rz(-1.0372459) q[2];
sx q[2];
rz(0.85565957) q[2];
rz(1.8850108) q[3];
sx q[3];
rz(-0.62658566) q[3];
sx q[3];
rz(1.1945061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7282495) q[0];
sx q[0];
rz(-3.0276868) q[0];
sx q[0];
rz(0.86694992) q[0];
rz(1.4504245) q[1];
sx q[1];
rz(-1.1879286) q[1];
sx q[1];
rz(0.11122045) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040505458) q[0];
sx q[0];
rz(-1.936628) q[0];
sx q[0];
rz(2.5554197) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84663518) q[2];
sx q[2];
rz(-0.66613448) q[2];
sx q[2];
rz(-0.29633477) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4091089) q[1];
sx q[1];
rz(-2.1883211) q[1];
sx q[1];
rz(-0.78583048) q[1];
x q[2];
rz(-0.86386776) q[3];
sx q[3];
rz(-1.4480531) q[3];
sx q[3];
rz(1.8254759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1078681) q[2];
sx q[2];
rz(-0.76216951) q[2];
sx q[2];
rz(2.6503358) q[2];
rz(1.2563502) q[3];
sx q[3];
rz(-1.4278853) q[3];
sx q[3];
rz(-0.45891941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97378174) q[0];
sx q[0];
rz(-1.7901006) q[0];
sx q[0];
rz(0.0034045086) q[0];
rz(0.34524521) q[1];
sx q[1];
rz(-2.7216941) q[1];
sx q[1];
rz(2.2664864) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9291139) q[0];
sx q[0];
rz(-2.1262332) q[0];
sx q[0];
rz(1.655913) q[0];
rz(0.11981182) q[2];
sx q[2];
rz(-1.7498328) q[2];
sx q[2];
rz(1.2546292) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9713953) q[1];
sx q[1];
rz(-1.1430234) q[1];
sx q[1];
rz(-1.0703767) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3703824) q[3];
sx q[3];
rz(-1.5516087) q[3];
sx q[3];
rz(0.33972188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.379999) q[2];
sx q[2];
rz(-2.9271409) q[2];
sx q[2];
rz(0.3271412) q[2];
rz(-1.3867779) q[3];
sx q[3];
rz(-1.1969457) q[3];
sx q[3];
rz(-3.0676214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.73579329) q[0];
sx q[0];
rz(-1.0013591) q[0];
sx q[0];
rz(-1.557198) q[0];
rz(-2.7994075) q[1];
sx q[1];
rz(-2.1125968) q[1];
sx q[1];
rz(-2.7142966) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9153684) q[0];
sx q[0];
rz(-1.5378204) q[0];
sx q[0];
rz(0.054126496) q[0];
x q[1];
rz(-2.7703366) q[2];
sx q[2];
rz(-1.0635738) q[2];
sx q[2];
rz(1.660342) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.72361833) q[1];
sx q[1];
rz(-2.2811541) q[1];
sx q[1];
rz(2.3478226) q[1];
x q[2];
rz(1.2702291) q[3];
sx q[3];
rz(-2.1168613) q[3];
sx q[3];
rz(0.48496839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7067318) q[2];
sx q[2];
rz(-2.2123983) q[2];
sx q[2];
rz(0.84563196) q[2];
rz(0.027806824) q[3];
sx q[3];
rz(-1.3549201) q[3];
sx q[3];
rz(-1.6644679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.528275) q[0];
sx q[0];
rz(-2.1149825) q[0];
sx q[0];
rz(3.0897621) q[0];
rz(-1.1593646) q[1];
sx q[1];
rz(-0.66957384) q[1];
sx q[1];
rz(-0.61028284) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0230816) q[0];
sx q[0];
rz(-1.7753965) q[0];
sx q[0];
rz(-1.3955297) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1716144) q[2];
sx q[2];
rz(-1.8464247) q[2];
sx q[2];
rz(2.3720019) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7597639) q[1];
sx q[1];
rz(-1.0914769) q[1];
sx q[1];
rz(-1.3717531) q[1];
rz(0.60812833) q[3];
sx q[3];
rz(-1.593489) q[3];
sx q[3];
rz(-0.80858999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3940238) q[2];
sx q[2];
rz(-1.3538066) q[2];
sx q[2];
rz(2.5547011) q[2];
rz(1.759257) q[3];
sx q[3];
rz(-2.0943677) q[3];
sx q[3];
rz(0.97572485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1047644) q[0];
sx q[0];
rz(-0.47127518) q[0];
sx q[0];
rz(2.9608744) q[0];
rz(1.7432927) q[1];
sx q[1];
rz(-1.2340052) q[1];
sx q[1];
rz(1.7416551) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0796868) q[0];
sx q[0];
rz(-0.8041412) q[0];
sx q[0];
rz(1.2595533) q[0];
rz(-pi) q[1];
rz(-1.9701875) q[2];
sx q[2];
rz(-1.9353494) q[2];
sx q[2];
rz(0.44427085) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8508262) q[1];
sx q[1];
rz(-1.6572448) q[1];
sx q[1];
rz(0.87169164) q[1];
x q[2];
rz(-0.65989699) q[3];
sx q[3];
rz(-2.0273034) q[3];
sx q[3];
rz(-2.5259301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7515298) q[2];
sx q[2];
rz(-2.2921102) q[2];
sx q[2];
rz(-2.0709822) q[2];
rz(-1.6603598) q[3];
sx q[3];
rz(-2.345572) q[3];
sx q[3];
rz(-1.7695561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2299131) q[0];
sx q[0];
rz(-2.3340618) q[0];
sx q[0];
rz(-2.5813778) q[0];
rz(2.920976) q[1];
sx q[1];
rz(-2.1897774) q[1];
sx q[1];
rz(-2.2645948) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91966832) q[0];
sx q[0];
rz(-3.1238365) q[0];
sx q[0];
rz(1.4229126) q[0];
rz(-pi) q[1];
rz(-1.8401519) q[2];
sx q[2];
rz(-1.1860086) q[2];
sx q[2];
rz(1.332217) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6798579) q[1];
sx q[1];
rz(-1.6731687) q[1];
sx q[1];
rz(1.3030679) q[1];
rz(-pi) q[2];
rz(-0.10538606) q[3];
sx q[3];
rz(-0.8871173) q[3];
sx q[3];
rz(2.2677641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.47473869) q[2];
sx q[2];
rz(-0.89484221) q[2];
sx q[2];
rz(-1.8043) q[2];
rz(2.9232591) q[3];
sx q[3];
rz(-2.1209013) q[3];
sx q[3];
rz(-1.7201503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.6364994) q[0];
sx q[0];
rz(-1.8890843) q[0];
sx q[0];
rz(0.1388347) q[0];
rz(-2.2598963) q[1];
sx q[1];
rz(-1.9813184) q[1];
sx q[1];
rz(-2.3169611) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6178095) q[0];
sx q[0];
rz(-1.6629761) q[0];
sx q[0];
rz(3.1193798) q[0];
rz(-pi) q[1];
rz(1.4662554) q[2];
sx q[2];
rz(-1.4741159) q[2];
sx q[2];
rz(-2.6226875) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4887374) q[1];
sx q[1];
rz(-2.0664381) q[1];
sx q[1];
rz(-1.7486217) q[1];
x q[2];
rz(0.10004754) q[3];
sx q[3];
rz(-2.2216019) q[3];
sx q[3];
rz(-1.5318961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3489939) q[2];
sx q[2];
rz(-0.94453347) q[2];
sx q[2];
rz(2.3080431) q[2];
rz(-0.27093568) q[3];
sx q[3];
rz(-1.1214333) q[3];
sx q[3];
rz(-2.2942395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18868294) q[0];
sx q[0];
rz(-1.2485349) q[0];
sx q[0];
rz(0.46440014) q[0];
rz(-0.69350997) q[1];
sx q[1];
rz(-0.92718569) q[1];
sx q[1];
rz(-0.69673353) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20737843) q[0];
sx q[0];
rz(-1.5153645) q[0];
sx q[0];
rz(1.4595281) q[0];
x q[1];
rz(0.040573434) q[2];
sx q[2];
rz(-1.3907205) q[2];
sx q[2];
rz(0.96873084) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.33524698) q[1];
sx q[1];
rz(-1.5434524) q[1];
sx q[1];
rz(-1.8198479) q[1];
rz(-1.5570627) q[3];
sx q[3];
rz(-0.93051118) q[3];
sx q[3];
rz(-0.69956799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36737475) q[2];
sx q[2];
rz(-3.0493224) q[2];
sx q[2];
rz(0.83887678) q[2];
rz(1.687626) q[3];
sx q[3];
rz(-1.5980915) q[3];
sx q[3];
rz(1.669223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0812747) q[0];
sx q[0];
rz(-1.769279) q[0];
sx q[0];
rz(-1.5244315) q[0];
rz(0.35628191) q[1];
sx q[1];
rz(-1.7226912) q[1];
sx q[1];
rz(1.3391395) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5863748) q[0];
sx q[0];
rz(-1.5065743) q[0];
sx q[0];
rz(-2.1703266) q[0];
x q[1];
rz(-1.6422604) q[2];
sx q[2];
rz(-1.2515837) q[2];
sx q[2];
rz(0.45305529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0799651) q[1];
sx q[1];
rz(-1.6929757) q[1];
sx q[1];
rz(-0.6811209) q[1];
x q[2];
rz(-1.7130835) q[3];
sx q[3];
rz(-1.5350966) q[3];
sx q[3];
rz(-1.3653099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39849207) q[2];
sx q[2];
rz(-2.3091381) q[2];
sx q[2];
rz(-1.8771646) q[2];
rz(0.33665952) q[3];
sx q[3];
rz(-0.94069898) q[3];
sx q[3];
rz(-1.3338859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7548512) q[0];
sx q[0];
rz(-1.516153) q[0];
sx q[0];
rz(-1.0990912) q[0];
rz(-2.5418538) q[1];
sx q[1];
rz(-0.44282423) q[1];
sx q[1];
rz(2.5437358) q[1];
rz(2.2221634) q[2];
sx q[2];
rz(-0.50242699) q[2];
sx q[2];
rz(-2.4519359) q[2];
rz(-2.609008) q[3];
sx q[3];
rz(-1.5859133) q[3];
sx q[3];
rz(1.5487158) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
