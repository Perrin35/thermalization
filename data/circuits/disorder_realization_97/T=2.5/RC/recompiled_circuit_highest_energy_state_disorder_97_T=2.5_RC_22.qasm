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
rz(2.4801369) q[0];
sx q[0];
rz(-2.8602726) q[0];
sx q[0];
rz(1.9982279) q[0];
rz(0.65302628) q[1];
sx q[1];
rz(-1.3209359) q[1];
sx q[1];
rz(0.015425711) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3000473) q[0];
sx q[0];
rz(-0.76800387) q[0];
sx q[0];
rz(0.83243911) q[0];
rz(-pi) q[1];
rz(-0.77818971) q[2];
sx q[2];
rz(-2.3986849) q[2];
sx q[2];
rz(-1.6011432) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.676486) q[1];
sx q[1];
rz(-2.3576405) q[1];
sx q[1];
rz(0.98760651) q[1];
rz(0.32415819) q[3];
sx q[3];
rz(-0.56190842) q[3];
sx q[3];
rz(3.0747385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.36110863) q[2];
sx q[2];
rz(-0.0094272308) q[2];
sx q[2];
rz(-2.0375605) q[2];
rz(-0.55936724) q[3];
sx q[3];
rz(-0.95061022) q[3];
sx q[3];
rz(-2.2573788) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0246494) q[0];
sx q[0];
rz(-0.58498061) q[0];
sx q[0];
rz(0.90644932) q[0];
rz(0.34140423) q[1];
sx q[1];
rz(-0.87541348) q[1];
sx q[1];
rz(-2.7934449) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6432836) q[0];
sx q[0];
rz(-1.3237244) q[0];
sx q[0];
rz(2.8227795) q[0];
x q[1];
rz(1.6115336) q[2];
sx q[2];
rz(-1.5968644) q[2];
sx q[2];
rz(-3.0422826) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9056919) q[1];
sx q[1];
rz(-2.5254619) q[1];
sx q[1];
rz(1.7006763) q[1];
rz(-pi) q[2];
rz(-1.320789) q[3];
sx q[3];
rz(-0.75859447) q[3];
sx q[3];
rz(1.7649253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.97588313) q[2];
sx q[2];
rz(-0.3419064) q[2];
sx q[2];
rz(-3.0561786) q[2];
rz(3.0488221) q[3];
sx q[3];
rz(-2.2723891) q[3];
sx q[3];
rz(-2.2658277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77466011) q[0];
sx q[0];
rz(-0.36151883) q[0];
sx q[0];
rz(2.7969978) q[0];
rz(1.5498281) q[1];
sx q[1];
rz(-1.7235618) q[1];
sx q[1];
rz(1.6835015) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.944779) q[0];
sx q[0];
rz(-0.85374248) q[0];
sx q[0];
rz(1.6193006) q[0];
rz(1.4665274) q[2];
sx q[2];
rz(-2.6223256) q[2];
sx q[2];
rz(-2.7291275) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.56006735) q[1];
sx q[1];
rz(-1.7398283) q[1];
sx q[1];
rz(-1.8174816) q[1];
rz(-pi) q[2];
x q[2];
rz(2.273748) q[3];
sx q[3];
rz(-0.72753564) q[3];
sx q[3];
rz(-0.65505469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.496326) q[2];
sx q[2];
rz(-0.93924773) q[2];
sx q[2];
rz(-2.0110896) q[2];
rz(2.1988403) q[3];
sx q[3];
rz(-2.0944984) q[3];
sx q[3];
rz(-1.2070791) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7154253) q[0];
sx q[0];
rz(-1.796145) q[0];
sx q[0];
rz(1.6834393) q[0];
rz(-1.0484877) q[1];
sx q[1];
rz(-1.6494992) q[1];
sx q[1];
rz(-0.47007158) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1972407) q[0];
sx q[0];
rz(-1.8507804) q[0];
sx q[0];
rz(0.51843317) q[0];
x q[1];
rz(-2.8747005) q[2];
sx q[2];
rz(-2.0837657) q[2];
sx q[2];
rz(1.5266872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.25950228) q[1];
sx q[1];
rz(-0.90356245) q[1];
sx q[1];
rz(1.3253071) q[1];
rz(2.4660956) q[3];
sx q[3];
rz(-2.1515931) q[3];
sx q[3];
rz(1.8449933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5423535) q[2];
sx q[2];
rz(-2.3848644) q[2];
sx q[2];
rz(-0.92764628) q[2];
rz(-1.2841691) q[3];
sx q[3];
rz(-1.410306) q[3];
sx q[3];
rz(0.94902432) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.014378431) q[0];
sx q[0];
rz(-2.7365186) q[0];
sx q[0];
rz(2.6601484) q[0];
rz(2.1638347) q[1];
sx q[1];
rz(-0.6441741) q[1];
sx q[1];
rz(-2.0668623) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.510493) q[0];
sx q[0];
rz(-1.298615) q[0];
sx q[0];
rz(-2.0049606) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3911209) q[2];
sx q[2];
rz(-0.40130645) q[2];
sx q[2];
rz(1.6502146) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8751443) q[1];
sx q[1];
rz(-1.1992663) q[1];
sx q[1];
rz(-2.8846011) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67258622) q[3];
sx q[3];
rz(-2.6574316) q[3];
sx q[3];
rz(1.3369833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0398728) q[2];
sx q[2];
rz(-2.0247255) q[2];
sx q[2];
rz(0.57363415) q[2];
rz(1.4839858) q[3];
sx q[3];
rz(-2.0935757) q[3];
sx q[3];
rz(-2.535533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023733519) q[0];
sx q[0];
rz(-2.8890299) q[0];
sx q[0];
rz(-2.8299487) q[0];
rz(-2.812884) q[1];
sx q[1];
rz(-1.8054104) q[1];
sx q[1];
rz(0.74105826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8097654) q[0];
sx q[0];
rz(-0.96967319) q[0];
sx q[0];
rz(-2.580216) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7996721) q[2];
sx q[2];
rz(-1.6308349) q[2];
sx q[2];
rz(2.6031074) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0523129) q[1];
sx q[1];
rz(-1.6340173) q[1];
sx q[1];
rz(2.2526645) q[1];
rz(-1.5696008) q[3];
sx q[3];
rz(-2.5107267) q[3];
sx q[3];
rz(1.6945171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3662423) q[2];
sx q[2];
rz(-2.5516208) q[2];
sx q[2];
rz(0.74014202) q[2];
rz(-0.38672334) q[3];
sx q[3];
rz(-2.4155278) q[3];
sx q[3];
rz(2.6630785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9750403) q[0];
sx q[0];
rz(-2.1599025) q[0];
sx q[0];
rz(0.24328406) q[0];
rz(-2.2443306) q[1];
sx q[1];
rz(-1.5829395) q[1];
sx q[1];
rz(-0.74403393) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1420741) q[0];
sx q[0];
rz(-1.4878977) q[0];
sx q[0];
rz(-0.3905475) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7626708) q[2];
sx q[2];
rz(-2.1728467) q[2];
sx q[2];
rz(0.10720358) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0571088) q[1];
sx q[1];
rz(-1.9388102) q[1];
sx q[1];
rz(2.1294566) q[1];
rz(0.18455581) q[3];
sx q[3];
rz(-1.0519487) q[3];
sx q[3];
rz(2.542582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1870785) q[2];
sx q[2];
rz(-2.3176471) q[2];
sx q[2];
rz(-3.1336866) q[2];
rz(1.6451969) q[3];
sx q[3];
rz(-3.0683066) q[3];
sx q[3];
rz(0.50905281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.024648333) q[0];
sx q[0];
rz(-3.109303) q[0];
sx q[0];
rz(-0.13667983) q[0];
rz(-1.3487123) q[1];
sx q[1];
rz(-1.8486479) q[1];
sx q[1];
rz(-2.6332556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7901403) q[0];
sx q[0];
rz(-0.41787028) q[0];
sx q[0];
rz(-1.2463742) q[0];
x q[1];
rz(-2.2466564) q[2];
sx q[2];
rz(-1.3351262) q[2];
sx q[2];
rz(1.6431944) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1994902) q[1];
sx q[1];
rz(-0.9170031) q[1];
sx q[1];
rz(-0.4545757) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6898722) q[3];
sx q[3];
rz(-1.649646) q[3];
sx q[3];
rz(0.54311968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6792182) q[2];
sx q[2];
rz(-2.4269673) q[2];
sx q[2];
rz(3.083631) q[2];
rz(-0.95311779) q[3];
sx q[3];
rz(-2.0942196) q[3];
sx q[3];
rz(0.63547772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5304831) q[0];
sx q[0];
rz(-1.7246752) q[0];
sx q[0];
rz(2.2755151) q[0];
rz(1.3653612) q[1];
sx q[1];
rz(-1.5581286) q[1];
sx q[1];
rz(-2.3376215) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83528501) q[0];
sx q[0];
rz(-2.270335) q[0];
sx q[0];
rz(2.7650096) q[0];
x q[1];
rz(-2.3525904) q[2];
sx q[2];
rz(-1.8937292) q[2];
sx q[2];
rz(-0.6168405) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.70529363) q[1];
sx q[1];
rz(-1.1501405) q[1];
sx q[1];
rz(0.83479014) q[1];
rz(-pi) q[2];
rz(-1.1671806) q[3];
sx q[3];
rz(-1.5909373) q[3];
sx q[3];
rz(0.27475629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0552401) q[2];
sx q[2];
rz(-2.99282) q[2];
sx q[2];
rz(2.0894334) q[2];
rz(-2.2822028) q[3];
sx q[3];
rz(-0.79901564) q[3];
sx q[3];
rz(2.1990282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4683485) q[0];
sx q[0];
rz(-2.4788661) q[0];
sx q[0];
rz(0.41917875) q[0];
rz(-0.016050922) q[1];
sx q[1];
rz(-1.586986) q[1];
sx q[1];
rz(-3.0174461) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9036839) q[0];
sx q[0];
rz(-0.57460898) q[0];
sx q[0];
rz(-1.2875071) q[0];
x q[1];
rz(2.9356758) q[2];
sx q[2];
rz(-1.3214998) q[2];
sx q[2];
rz(0.60520335) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0090186) q[1];
sx q[1];
rz(-0.30840519) q[1];
sx q[1];
rz(-2.5243702) q[1];
rz(-1.8087555) q[3];
sx q[3];
rz(-0.66413022) q[3];
sx q[3];
rz(0.40011621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.19196709) q[2];
sx q[2];
rz(-2.4032335) q[2];
sx q[2];
rz(0.16853608) q[2];
rz(-1.4847697) q[3];
sx q[3];
rz(-0.46983457) q[3];
sx q[3];
rz(-0.43009871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1930502) q[0];
sx q[0];
rz(-2.6621303) q[0];
sx q[0];
rz(2.9051176) q[0];
rz(-2.3169658) q[1];
sx q[1];
rz(-1.6025447) q[1];
sx q[1];
rz(-1.2784169) q[1];
rz(0.58165089) q[2];
sx q[2];
rz(-0.86514513) q[2];
sx q[2];
rz(-1.9164597) q[2];
rz(-1.6705728) q[3];
sx q[3];
rz(-2.592516) q[3];
sx q[3];
rz(1.9883131) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
