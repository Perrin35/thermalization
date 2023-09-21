OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(-0.17012574) q[0];
sx q[0];
rz(2.3556019) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(-2.4906467) q[1];
sx q[1];
rz(-0.63408607) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68361002) q[0];
sx q[0];
rz(-1.5994706) q[0];
sx q[0];
rz(-0.81575127) q[0];
rz(-pi) q[1];
rz(-0.91648957) q[2];
sx q[2];
rz(-1.620703) q[2];
sx q[2];
rz(-2.8061342) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0178535) q[1];
sx q[1];
rz(-0.98736073) q[1];
sx q[1];
rz(0.53637335) q[1];
rz(-pi) q[2];
x q[2];
rz(0.05539031) q[3];
sx q[3];
rz(-2.7408544) q[3];
sx q[3];
rz(-1.8000359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.2259851) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(1.8784286) q[2];
rz(-1.8566711) q[3];
sx q[3];
rz(-1.6804755) q[3];
sx q[3];
rz(2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15853515) q[0];
sx q[0];
rz(-1.2936658) q[0];
sx q[0];
rz(-2.7040226) q[0];
rz(0.63105398) q[1];
sx q[1];
rz(-2.8129306) q[1];
sx q[1];
rz(-0.22110573) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1572004) q[0];
sx q[0];
rz(-1.385014) q[0];
sx q[0];
rz(-1.1093344) q[0];
x q[1];
rz(1.2220076) q[2];
sx q[2];
rz(-1.4663327) q[2];
sx q[2];
rz(1.8609357) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4202538) q[1];
sx q[1];
rz(-1.6920648) q[1];
sx q[1];
rz(0.46055693) q[1];
rz(-pi) q[2];
rz(2.6392194) q[3];
sx q[3];
rz(-1.6564192) q[3];
sx q[3];
rz(1.8300213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.21800403) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(-2.5533) q[2];
rz(-0.44899392) q[3];
sx q[3];
rz(-0.42789999) q[3];
sx q[3];
rz(-1.0480405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56851971) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(1.0004689) q[0];
rz(-2.4160642) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(-0.75769889) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0962778) q[0];
sx q[0];
rz(-1.0490388) q[0];
sx q[0];
rz(1.9406712) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4865815) q[2];
sx q[2];
rz(-2.8821324) q[2];
sx q[2];
rz(-0.26817817) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0211027) q[1];
sx q[1];
rz(-2.0334525) q[1];
sx q[1];
rz(-2.2282269) q[1];
rz(-pi) q[2];
rz(-0.63852255) q[3];
sx q[3];
rz(-1.4811852) q[3];
sx q[3];
rz(1.5935957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2086601) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(-0.83646742) q[2];
rz(1.6992735) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(-0.20382717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02012415) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(-0.51112038) q[0];
rz(-3.064149) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(1.3285332) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9851345) q[0];
sx q[0];
rz(-1.4159604) q[0];
sx q[0];
rz(1.3550718) q[0];
x q[1];
rz(-2.6287597) q[2];
sx q[2];
rz(-2.3094258) q[2];
sx q[2];
rz(-1.8682478) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9424142) q[1];
sx q[1];
rz(-1.5702827) q[1];
sx q[1];
rz(1.8375988) q[1];
rz(-0.10947157) q[3];
sx q[3];
rz(-1.3502035) q[3];
sx q[3];
rz(2.4601439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.093734309) q[2];
sx q[2];
rz(-1.9231984) q[2];
sx q[2];
rz(-1.1070586) q[2];
rz(-0.47248653) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1012786) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(-2.8033946) q[0];
rz(-1.8473373) q[1];
sx q[1];
rz(-1.5599524) q[1];
sx q[1];
rz(-2.6370874) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11855928) q[0];
sx q[0];
rz(-2.1001864) q[0];
sx q[0];
rz(0.27079196) q[0];
rz(-pi) q[1];
rz(-2.4733876) q[2];
sx q[2];
rz(-2.5720111) q[2];
sx q[2];
rz(0.44797209) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.075715847) q[1];
sx q[1];
rz(-1.8738235) q[1];
sx q[1];
rz(-0.24995835) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2284331) q[3];
sx q[3];
rz(-1.0201766) q[3];
sx q[3];
rz(-3.1331568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0304886) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(1.4939235) q[2];
rz(1.4533639) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(1.3822314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(2.9313653) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(-2.0507623) q[0];
rz(2.6018654) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(-2.9528023) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.990373) q[0];
sx q[0];
rz(-1.9301206) q[0];
sx q[0];
rz(1.8381401) q[0];
rz(-pi) q[1];
rz(1.2994453) q[2];
sx q[2];
rz(-1.6711298) q[2];
sx q[2];
rz(3.0999822) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.97555679) q[1];
sx q[1];
rz(-1.3284725) q[1];
sx q[1];
rz(1.5029328) q[1];
rz(-pi) q[2];
rz(2.8998428) q[3];
sx q[3];
rz(-1.3931735) q[3];
sx q[3];
rz(1.0279946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9699036) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.8264654) q[2];
rz(1.5054437) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(-1.2833387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.090102) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(-0.32399696) q[0];
rz(-1.3011159) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(-1.3791929) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1349072) q[0];
sx q[0];
rz(-2.8838257) q[0];
sx q[0];
rz(1.2447312) q[0];
rz(-0.71521476) q[2];
sx q[2];
rz(-2.5589802) q[2];
sx q[2];
rz(1.4441393) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.52723253) q[1];
sx q[1];
rz(-0.34253201) q[1];
sx q[1];
rz(-0.50369461) q[1];
rz(-pi) q[2];
rz(-0.93382436) q[3];
sx q[3];
rz(-2.056042) q[3];
sx q[3];
rz(2.3683734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.032701187) q[2];
sx q[2];
rz(-1.6872493) q[2];
sx q[2];
rz(2.1984055) q[2];
rz(0.33106783) q[3];
sx q[3];
rz(-1.7739242) q[3];
sx q[3];
rz(-2.846472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-3.023495) q[0];
sx q[0];
rz(-2.5546615) q[0];
sx q[0];
rz(-0.76422894) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(-1.8483298) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18023597) q[0];
sx q[0];
rz(-1.1399674) q[0];
sx q[0];
rz(2.9234617) q[0];
rz(-0.13883491) q[2];
sx q[2];
rz(-0.60563696) q[2];
sx q[2];
rz(-2.5790737) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.075899374) q[1];
sx q[1];
rz(-1.432044) q[1];
sx q[1];
rz(-1.9779512) q[1];
x q[2];
rz(2.7420298) q[3];
sx q[3];
rz(-1.385653) q[3];
sx q[3];
rz(-3.1294587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.76688898) q[2];
sx q[2];
rz(-2.1024487) q[2];
sx q[2];
rz(0.58132201) q[2];
rz(2.2733722) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(1.3114312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-0.54810369) q[0];
sx q[0];
rz(-2.5618401) q[0];
sx q[0];
rz(-1.8909489) q[0];
rz(0.99682322) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.8539799) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8087013) q[0];
sx q[0];
rz(-1.5960072) q[0];
sx q[0];
rz(-2.0602134) q[0];
rz(-1.397923) q[2];
sx q[2];
rz(-2.0506952) q[2];
sx q[2];
rz(2.929504) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7848283) q[1];
sx q[1];
rz(-2.462466) q[1];
sx q[1];
rz(-1.2374415) q[1];
rz(2.3171114) q[3];
sx q[3];
rz(-1.609625) q[3];
sx q[3];
rz(-2.3076434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.320497) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(-1.256475) q[2];
rz(2.9750032) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2607516) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(0.32521954) q[0];
rz(1.1351769) q[1];
sx q[1];
rz(-0.67712855) q[1];
sx q[1];
rz(2.7744055) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4627535) q[0];
sx q[0];
rz(-3.084393) q[0];
sx q[0];
rz(2.5662533) q[0];
rz(2.514421) q[2];
sx q[2];
rz(-1.7998724) q[2];
sx q[2];
rz(2.0723745) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.50670934) q[1];
sx q[1];
rz(-0.92223972) q[1];
sx q[1];
rz(1.8222005) q[1];
rz(0.92948593) q[3];
sx q[3];
rz(-1.1747642) q[3];
sx q[3];
rz(0.091077387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8964768) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(2.0754576) q[2];
rz(0.079244763) q[3];
sx q[3];
rz(-2.5388122) q[3];
sx q[3];
rz(-1.4762896) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29522482) q[0];
sx q[0];
rz(-1.5548779) q[0];
sx q[0];
rz(1.5665733) q[0];
rz(0.29905839) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(-0.085061442) q[2];
sx q[2];
rz(-2.3903575) q[2];
sx q[2];
rz(0.059546197) q[2];
rz(-1.7442262) q[3];
sx q[3];
rz(-0.72931029) q[3];
sx q[3];
rz(1.9839877) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];