OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(0.064602764) q[0];
sx q[0];
rz(6.3048007) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(-1.3316766) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0776805) q[0];
sx q[0];
rz(-2.4644116) q[0];
sx q[0];
rz(-1.022524) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9900436) q[2];
sx q[2];
rz(-0.24818072) q[2];
sx q[2];
rz(-2.575945) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6380438) q[1];
sx q[1];
rz(-0.81144864) q[1];
sx q[1];
rz(1.391295) q[1];
rz(-0.94727912) q[3];
sx q[3];
rz(-1.8489269) q[3];
sx q[3];
rz(-2.5840685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8109479) q[2];
sx q[2];
rz(-1.4928879) q[2];
sx q[2];
rz(-0.60418207) q[2];
rz(-2.1172681) q[3];
sx q[3];
rz(-1.1573135) q[3];
sx q[3];
rz(2.0143051) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8169096) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(-1.0634134) q[0];
rz(-2.2564607) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(0.0016454776) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5434108) q[0];
sx q[0];
rz(-3.0953005) q[0];
sx q[0];
rz(1.8020736) q[0];
x q[1];
rz(1.8284945) q[2];
sx q[2];
rz(-0.94444599) q[2];
sx q[2];
rz(-1.4305654) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.590608) q[1];
sx q[1];
rz(-1.3270757) q[1];
sx q[1];
rz(0.79382146) q[1];
rz(-pi) q[2];
x q[2];
rz(0.9886338) q[3];
sx q[3];
rz(-2.8391264) q[3];
sx q[3];
rz(-0.31103381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1559747) q[2];
sx q[2];
rz(-1.6321471) q[2];
sx q[2];
rz(0.70811159) q[2];
rz(1.0937141) q[3];
sx q[3];
rz(-1.8023068) q[3];
sx q[3];
rz(2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22096069) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(2.3143342) q[0];
rz(-0.0050841252) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(2.0522096) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82747805) q[0];
sx q[0];
rz(-0.70686045) q[0];
sx q[0];
rz(2.2543738) q[0];
x q[1];
rz(2.543534) q[2];
sx q[2];
rz(-1.2130249) q[2];
sx q[2];
rz(1.9269112) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8100064) q[1];
sx q[1];
rz(-2.7467105) q[1];
sx q[1];
rz(-1.3473131) q[1];
x q[2];
rz(-2.1316707) q[3];
sx q[3];
rz(-2.2787333) q[3];
sx q[3];
rz(2.5885955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0638782) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(2.2568978) q[2];
rz(1.9583154) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(1.6515091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5692212) q[0];
sx q[0];
rz(-0.5643934) q[0];
sx q[0];
rz(-0.72682056) q[0];
rz(2.3379393) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(-0.35983905) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026304631) q[0];
sx q[0];
rz(-1.7548772) q[0];
sx q[0];
rz(-0.91362761) q[0];
x q[1];
rz(0.6511351) q[2];
sx q[2];
rz(-1.480181) q[2];
sx q[2];
rz(-3.0504984) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.29405669) q[1];
sx q[1];
rz(-1.9901853) q[1];
sx q[1];
rz(-1.9860752) q[1];
x q[2];
rz(-1.2218277) q[3];
sx q[3];
rz(-1.2811321) q[3];
sx q[3];
rz(-1.3834013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5741253) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(0.7652258) q[2];
rz(-2.3848173) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(-2.9005907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9969479) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(-1.7255406) q[0];
rz(-2.7744746) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(1.0353154) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2238335) q[0];
sx q[0];
rz(-2.1609554) q[0];
sx q[0];
rz(0.53442861) q[0];
rz(-pi) q[1];
rz(0.073280235) q[2];
sx q[2];
rz(-0.49256941) q[2];
sx q[2];
rz(-0.57893334) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1977735) q[1];
sx q[1];
rz(-0.10995956) q[1];
sx q[1];
rz(-2.6206559) q[1];
x q[2];
rz(1.7807998) q[3];
sx q[3];
rz(-1.8404507) q[3];
sx q[3];
rz(-0.3013914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3570024) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(-2.3727097) q[2];
rz(-2.8055577) q[3];
sx q[3];
rz(-2.3612645) q[3];
sx q[3];
rz(1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1795905) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(-1.1451716) q[0];
rz(2.0369453) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(-2.9343658) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.094039) q[0];
sx q[0];
rz(-2.5045966) q[0];
sx q[0];
rz(-0.23547049) q[0];
rz(-pi) q[1];
rz(2.0358762) q[2];
sx q[2];
rz(-2.1207223) q[2];
sx q[2];
rz(-1.20649) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.44240272) q[1];
sx q[1];
rz(-1.8132205) q[1];
sx q[1];
rz(1.3139903) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9793671) q[3];
sx q[3];
rz(-2.3174006) q[3];
sx q[3];
rz(-1.970286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8032288) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(-0.72171372) q[2];
rz(1.8668113) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(0.26003626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8969144) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(0.73202837) q[0];
rz(-3.1320944) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(-0.2917372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1222629) q[0];
sx q[0];
rz(-0.42207345) q[0];
sx q[0];
rz(0.12515573) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10642274) q[2];
sx q[2];
rz(-1.2692045) q[2];
sx q[2];
rz(1.0943204) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74624324) q[1];
sx q[1];
rz(-0.33848539) q[1];
sx q[1];
rz(-0.45792087) q[1];
rz(-pi) q[2];
rz(-1.2613867) q[3];
sx q[3];
rz(-2.072812) q[3];
sx q[3];
rz(-1.7878143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8774524) q[2];
sx q[2];
rz(-1.9984657) q[2];
sx q[2];
rz(2.3042802) q[2];
rz(1.9705747) q[3];
sx q[3];
rz(-1.6224909) q[3];
sx q[3];
rz(0.062019197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0432805) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(-2.4777381) q[0];
rz(-0.10617667) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(2.1829139) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52866919) q[0];
sx q[0];
rz(-1.0597502) q[0];
sx q[0];
rz(-2.5445166) q[0];
x q[1];
rz(-1.1218698) q[2];
sx q[2];
rz(-2.556483) q[2];
sx q[2];
rz(1.6758855) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2221335) q[1];
sx q[1];
rz(-0.53680116) q[1];
sx q[1];
rz(1.7807351) q[1];
rz(-0.58416768) q[3];
sx q[3];
rz(-1.468588) q[3];
sx q[3];
rz(0.16971961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0083996) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(0.91910249) q[2];
rz(-1.7637926) q[3];
sx q[3];
rz(-1.931124) q[3];
sx q[3];
rz(1.1834043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.33525) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(-0.43689716) q[0];
rz(0.70029744) q[1];
sx q[1];
rz(-1.7179787) q[1];
sx q[1];
rz(0.46554309) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8252194) q[0];
sx q[0];
rz(-1.6557949) q[0];
sx q[0];
rz(2.4337016) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8579673) q[2];
sx q[2];
rz(-2.7138777) q[2];
sx q[2];
rz(-2.5123793) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1411966) q[1];
sx q[1];
rz(-1.9138412) q[1];
sx q[1];
rz(0.47980206) q[1];
x q[2];
rz(-0.97165473) q[3];
sx q[3];
rz(-2.7145436) q[3];
sx q[3];
rz(0.077582434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(1.643606) q[2];
rz(0.20478976) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(2.8695316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64514226) q[0];
sx q[0];
rz(-0.56149879) q[0];
sx q[0];
rz(2.0196594) q[0];
rz(0.75138584) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(-0.5823935) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96888992) q[0];
sx q[0];
rz(-2.4318998) q[0];
sx q[0];
rz(-1.885528) q[0];
rz(-pi) q[1];
rz(-1.0803797) q[2];
sx q[2];
rz(-1.7047802) q[2];
sx q[2];
rz(-1.6756563) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6873684) q[1];
sx q[1];
rz(-2.9240989) q[1];
sx q[1];
rz(-2.1464159) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5428316) q[3];
sx q[3];
rz(-1.7879322) q[3];
sx q[3];
rz(-2.12487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0104684) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(-0.76254145) q[2];
rz(-1.7307581) q[3];
sx q[3];
rz(-2.2236731) q[3];
sx q[3];
rz(-1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0762155) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(1.8021884) q[1];
sx q[1];
rz(-1.6911472) q[1];
sx q[1];
rz(-1.4858248) q[1];
rz(-0.28152485) q[2];
sx q[2];
rz(-2.3073961) q[2];
sx q[2];
rz(-2.6632593) q[2];
rz(0.4776095) q[3];
sx q[3];
rz(-1.1160679) q[3];
sx q[3];
rz(3.06649) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];