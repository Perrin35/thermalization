OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6931273) q[0];
sx q[0];
rz(5.7603523) q[0];
sx q[0];
rz(8.8011959) q[0];
rz(-2.8514255) q[1];
sx q[1];
rz(-0.71915141) q[1];
sx q[1];
rz(-2.6410988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252285) q[0];
sx q[0];
rz(-2.541757) q[0];
sx q[0];
rz(-0.0068378011) q[0];
x q[1];
rz(3.0481553) q[2];
sx q[2];
rz(-1.7200025) q[2];
sx q[2];
rz(2.0051533) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9762293) q[1];
sx q[1];
rz(-0.86127087) q[1];
sx q[1];
rz(-2.707259) q[1];
rz(-pi) q[2];
rz(0.15668232) q[3];
sx q[3];
rz(-1.1374258) q[3];
sx q[3];
rz(-0.79790786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7913251) q[2];
sx q[2];
rz(-1.2056377) q[2];
sx q[2];
rz(-1.9187437) q[2];
rz(-1.6932999) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(-2.1712415) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7704849) q[0];
sx q[0];
rz(-1.6587695) q[0];
sx q[0];
rz(-2.0626542) q[0];
rz(-1.7547912) q[1];
sx q[1];
rz(-2.3290122) q[1];
sx q[1];
rz(2.4761377) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6404214) q[0];
sx q[0];
rz(-1.2873532) q[0];
sx q[0];
rz(-2.6675176) q[0];
x q[1];
rz(1.3999248) q[2];
sx q[2];
rz(-1.2534281) q[2];
sx q[2];
rz(-2.3651809) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8382149) q[1];
sx q[1];
rz(-1.8068552) q[1];
sx q[1];
rz(-0.49535507) q[1];
x q[2];
rz(-0.034025107) q[3];
sx q[3];
rz(-2.2799387) q[3];
sx q[3];
rz(2.2505086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5370496) q[2];
sx q[2];
rz(-1.4031354) q[2];
sx q[2];
rz(-0.075142168) q[2];
rz(1.4705307) q[3];
sx q[3];
rz(-1.1276779) q[3];
sx q[3];
rz(0.69141928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5866518) q[0];
sx q[0];
rz(-1.8692769) q[0];
sx q[0];
rz(2.3828322) q[0];
rz(1.2930019) q[1];
sx q[1];
rz(-1.2932152) q[1];
sx q[1];
rz(-1.4000777) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3103257) q[0];
sx q[0];
rz(-1.184549) q[0];
sx q[0];
rz(-1.8589742) q[0];
rz(-pi) q[1];
rz(-1.7064352) q[2];
sx q[2];
rz(-1.3856158) q[2];
sx q[2];
rz(0.26091012) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.61664596) q[1];
sx q[1];
rz(-2.4541306) q[1];
sx q[1];
rz(1.0487268) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16617822) q[3];
sx q[3];
rz(-1.8718534) q[3];
sx q[3];
rz(2.500246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.65163461) q[2];
sx q[2];
rz(-0.48214665) q[2];
sx q[2];
rz(2.4839694) q[2];
rz(-1.970132) q[3];
sx q[3];
rz(-1.6572584) q[3];
sx q[3];
rz(-1.7224147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.3477429) q[0];
sx q[0];
rz(-0.94809735) q[0];
sx q[0];
rz(-1.4452274) q[0];
rz(-1.4472648) q[1];
sx q[1];
rz(-1.6479965) q[1];
sx q[1];
rz(2.7935374) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845683) q[0];
sx q[0];
rz(-0.94385249) q[0];
sx q[0];
rz(1.9355965) q[0];
rz(1.1966755) q[2];
sx q[2];
rz(-1.7110363) q[2];
sx q[2];
rz(-1.3333048) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.12236774) q[1];
sx q[1];
rz(-2.5248563) q[1];
sx q[1];
rz(-1.7255746) q[1];
rz(-pi) q[2];
rz(1.7054889) q[3];
sx q[3];
rz(-1.7309942) q[3];
sx q[3];
rz(2.687541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41670123) q[2];
sx q[2];
rz(-1.3092224) q[2];
sx q[2];
rz(2.7187738) q[2];
rz(-2.4041798) q[3];
sx q[3];
rz(-0.80602065) q[3];
sx q[3];
rz(0.084658682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87930644) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(-2.6111531) q[0];
rz(-2.2166705) q[1];
sx q[1];
rz(-1.5735807) q[1];
sx q[1];
rz(-1.2984498) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7206551) q[0];
sx q[0];
rz(-1.7147831) q[0];
sx q[0];
rz(1.4074586) q[0];
x q[1];
rz(1.2735882) q[2];
sx q[2];
rz(-0.98515918) q[2];
sx q[2];
rz(-0.98779087) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4358208) q[1];
sx q[1];
rz(-2.768369) q[1];
sx q[1];
rz(-0.2758287) q[1];
x q[2];
rz(0.97814822) q[3];
sx q[3];
rz(-1.1453298) q[3];
sx q[3];
rz(0.076171906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9412781) q[2];
sx q[2];
rz(-1.986074) q[2];
sx q[2];
rz(-2.6679664) q[2];
rz(3.04223) q[3];
sx q[3];
rz(-1.8615581) q[3];
sx q[3];
rz(2.30106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6376003) q[0];
sx q[0];
rz(-1.7853328) q[0];
sx q[0];
rz(-0.9978869) q[0];
rz(-2.2672794) q[1];
sx q[1];
rz(-2.120178) q[1];
sx q[1];
rz(-2.6748437) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40690639) q[0];
sx q[0];
rz(-2.2285301) q[0];
sx q[0];
rz(-1.8440194) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6830707) q[2];
sx q[2];
rz(-0.47007559) q[2];
sx q[2];
rz(-2.2058861) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1280061) q[1];
sx q[1];
rz(-0.15833536) q[1];
sx q[1];
rz(-2.8360785) q[1];
rz(1.407133) q[3];
sx q[3];
rz(-1.6657077) q[3];
sx q[3];
rz(1.519219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2816887) q[2];
sx q[2];
rz(-0.47912654) q[2];
sx q[2];
rz(1.5768645) q[2];
rz(-2.5148897) q[3];
sx q[3];
rz(-1.7765216) q[3];
sx q[3];
rz(2.4600162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(1.3307813) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(-0.41982857) q[0];
rz(2.9201674) q[1];
sx q[1];
rz(-2.6629993) q[1];
sx q[1];
rz(1.5484757) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89824642) q[0];
sx q[0];
rz(-1.6582489) q[0];
sx q[0];
rz(-1.5792219) q[0];
x q[1];
rz(1.2862455) q[2];
sx q[2];
rz(-1.3853405) q[2];
sx q[2];
rz(-0.68765771) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.35510264) q[1];
sx q[1];
rz(-3.0349602) q[1];
sx q[1];
rz(0.14270466) q[1];
rz(-pi) q[2];
rz(2.0358622) q[3];
sx q[3];
rz(-1.2601868) q[3];
sx q[3];
rz(3.0748933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.55591136) q[2];
sx q[2];
rz(-2.6101117) q[2];
sx q[2];
rz(-1.7377724) q[2];
rz(0.79706556) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(1.2169303) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4306915) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(2.8523493) q[0];
rz(-2.5166683) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(1.3141059) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24459141) q[0];
sx q[0];
rz(-2.3750711) q[0];
sx q[0];
rz(-0.15890973) q[0];
x q[1];
rz(-1.9004702) q[2];
sx q[2];
rz(-2.2909819) q[2];
sx q[2];
rz(2.5230797) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9574979) q[1];
sx q[1];
rz(-3.0093319) q[1];
sx q[1];
rz(1.9930507) q[1];
rz(-pi) q[2];
rz(1.2460327) q[3];
sx q[3];
rz(-2.357558) q[3];
sx q[3];
rz(1.8079545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4079995) q[2];
sx q[2];
rz(-1.5635798) q[2];
sx q[2];
rz(-1.4286263) q[2];
rz(0.96380487) q[3];
sx q[3];
rz(-2.0644085) q[3];
sx q[3];
rz(2.111964) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52255094) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(1.2458941) q[0];
rz(-3.030581) q[1];
sx q[1];
rz(-1.1975892) q[1];
sx q[1];
rz(0.54661173) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0194191) q[0];
sx q[0];
rz(-1.1328567) q[0];
sx q[0];
rz(2.7287448) q[0];
rz(-0.39491744) q[2];
sx q[2];
rz(-1.2262218) q[2];
sx q[2];
rz(-3.092569) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6555772) q[1];
sx q[1];
rz(-1.7421107) q[1];
sx q[1];
rz(1.1737215) q[1];
rz(0.9162174) q[3];
sx q[3];
rz(-1.9700053) q[3];
sx q[3];
rz(2.6938714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0299915) q[2];
sx q[2];
rz(-0.84838715) q[2];
sx q[2];
rz(-1.6938422) q[2];
rz(1.1374121) q[3];
sx q[3];
rz(-1.8177989) q[3];
sx q[3];
rz(0.64731961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2820213) q[0];
sx q[0];
rz(-2.8347926) q[0];
sx q[0];
rz(-0.71722537) q[0];
rz(-1.9316797) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(-0.70770121) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4004138) q[0];
sx q[0];
rz(-0.78010633) q[0];
sx q[0];
rz(-1.0723423) q[0];
x q[1];
rz(1.0102347) q[2];
sx q[2];
rz(-0.26294225) q[2];
sx q[2];
rz(-1.5514785) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6245539) q[1];
sx q[1];
rz(-2.3561764) q[1];
sx q[1];
rz(-1.1167657) q[1];
x q[2];
rz(0.95146146) q[3];
sx q[3];
rz(-1.2776432) q[3];
sx q[3];
rz(2.4952863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6282965) q[2];
sx q[2];
rz(-2.4847023) q[2];
sx q[2];
rz(-2.2383402) q[2];
rz(1.603027) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(-0.85047754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.306504) q[0];
sx q[0];
rz(-2.7764414) q[0];
sx q[0];
rz(2.2055702) q[0];
rz(-2.3256336) q[1];
sx q[1];
rz(-2.7201256) q[1];
sx q[1];
rz(1.0526007) q[1];
rz(0.46335285) q[2];
sx q[2];
rz(-1.9909161) q[2];
sx q[2];
rz(-1.9009895) q[2];
rz(3.0872185) q[3];
sx q[3];
rz(-1.5297223) q[3];
sx q[3];
rz(-2.3755304) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
