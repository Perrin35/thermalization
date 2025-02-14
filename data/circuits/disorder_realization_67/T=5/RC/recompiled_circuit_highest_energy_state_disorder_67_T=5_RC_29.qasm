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
rz(-0.29210583) q[0];
sx q[0];
rz(-0.6120683) q[0];
sx q[0];
rz(0.050961343) q[0];
rz(-1.851097) q[1];
sx q[1];
rz(4.4176997) q[1];
sx q[1];
rz(7.1292276) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.660419) q[0];
sx q[0];
rz(-1.9194897) q[0];
sx q[0];
rz(1.3997188) q[0];
rz(-2.2162391) q[2];
sx q[2];
rz(-0.45956372) q[2];
sx q[2];
rz(1.1499576) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4197246) q[1];
sx q[1];
rz(-1.8375735) q[1];
sx q[1];
rz(0.50290033) q[1];
rz(1.5624763) q[3];
sx q[3];
rz(-1.4589615) q[3];
sx q[3];
rz(-2.73561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35752615) q[2];
sx q[2];
rz(-2.1937723) q[2];
sx q[2];
rz(2.2785462) q[2];
rz(2.0283902) q[3];
sx q[3];
rz(-0.92156711) q[3];
sx q[3];
rz(2.34288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.6031826) q[0];
sx q[0];
rz(-2.4487459) q[0];
sx q[0];
rz(-2.8811654) q[0];
rz(0.32568359) q[1];
sx q[1];
rz(-0.98418701) q[1];
sx q[1];
rz(0.30776417) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4474898) q[0];
sx q[0];
rz(-1.3110975) q[0];
sx q[0];
rz(-2.3070241) q[0];
rz(-pi) q[1];
rz(-1.8886052) q[2];
sx q[2];
rz(-1.121064) q[2];
sx q[2];
rz(-0.98435246) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.968281) q[1];
sx q[1];
rz(-1.8112343) q[1];
sx q[1];
rz(-3.0467826) q[1];
rz(0.93538385) q[3];
sx q[3];
rz(-1.5652802) q[3];
sx q[3];
rz(0.25980205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8327568) q[2];
sx q[2];
rz(-0.86897659) q[2];
sx q[2];
rz(0.55574179) q[2];
rz(-0.51221383) q[3];
sx q[3];
rz(-2.5930976) q[3];
sx q[3];
rz(-0.51928025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.888805) q[0];
sx q[0];
rz(-2.7894809) q[0];
sx q[0];
rz(0.74932253) q[0];
rz(-1.0961756) q[1];
sx q[1];
rz(-1.3744033) q[1];
sx q[1];
rz(2.7574417) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6836943) q[0];
sx q[0];
rz(-1.3078469) q[0];
sx q[0];
rz(-2.5431741) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5403156) q[2];
sx q[2];
rz(-1.3601175) q[2];
sx q[2];
rz(-2.7803068) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8763346) q[1];
sx q[1];
rz(-1.6240082) q[1];
sx q[1];
rz(-1.8734147) q[1];
rz(-2.8504653) q[3];
sx q[3];
rz(-2.2343221) q[3];
sx q[3];
rz(-2.4255668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7595547) q[2];
sx q[2];
rz(-1.3052772) q[2];
sx q[2];
rz(1.6273512) q[2];
rz(2.4115327) q[3];
sx q[3];
rz(-0.74159139) q[3];
sx q[3];
rz(-0.39456427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45383129) q[0];
sx q[0];
rz(-1.7173959) q[0];
sx q[0];
rz(0.060039595) q[0];
rz(-2.2278348) q[1];
sx q[1];
rz(-0.91122952) q[1];
sx q[1];
rz(2.2027016) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70907043) q[0];
sx q[0];
rz(-0.071252206) q[0];
sx q[0];
rz(2.6547316) q[0];
rz(-0.082415103) q[2];
sx q[2];
rz(-1.6019099) q[2];
sx q[2];
rz(-1.7331725) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.36870391) q[1];
sx q[1];
rz(-1.0828754) q[1];
sx q[1];
rz(0.027603961) q[1];
rz(2.4603953) q[3];
sx q[3];
rz(-1.0170611) q[3];
sx q[3];
rz(-2.9217229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.011220304) q[2];
sx q[2];
rz(-0.29272407) q[2];
sx q[2];
rz(0.65227738) q[2];
rz(1.8782714) q[3];
sx q[3];
rz(-1.616856) q[3];
sx q[3];
rz(-1.9704845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071282722) q[0];
sx q[0];
rz(-1.902453) q[0];
sx q[0];
rz(0.5624482) q[0];
rz(-1.7930454) q[1];
sx q[1];
rz(-2.8548073) q[1];
sx q[1];
rz(-0.58132344) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5347915) q[0];
sx q[0];
rz(-2.013371) q[0];
sx q[0];
rz(-2.3032196) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5346709) q[2];
sx q[2];
rz(-2.0310406) q[2];
sx q[2];
rz(-0.66260105) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5229458) q[1];
sx q[1];
rz(-0.60877555) q[1];
sx q[1];
rz(-2.2227915) q[1];
rz(-pi) q[2];
rz(-0.64251803) q[3];
sx q[3];
rz(-1.5356365) q[3];
sx q[3];
rz(-2.8865451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.44744667) q[2];
sx q[2];
rz(-1.4299102) q[2];
sx q[2];
rz(-2.4549129) q[2];
rz(0.40003362) q[3];
sx q[3];
rz(-1.2587849) q[3];
sx q[3];
rz(-0.011628477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43468633) q[0];
sx q[0];
rz(-2.1676368) q[0];
sx q[0];
rz(0.45027012) q[0];
rz(0.88092583) q[1];
sx q[1];
rz(-1.2342781) q[1];
sx q[1];
rz(1.1164104) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6754782) q[0];
sx q[0];
rz(-2.3494534) q[0];
sx q[0];
rz(-1.8380497) q[0];
x q[1];
rz(-2.7911602) q[2];
sx q[2];
rz(-1.74518) q[2];
sx q[2];
rz(-2.7728105) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3456381) q[1];
sx q[1];
rz(-1.5904203) q[1];
sx q[1];
rz(-0.96052891) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93706705) q[3];
sx q[3];
rz(-2.6014199) q[3];
sx q[3];
rz(-0.34264178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.474596) q[2];
sx q[2];
rz(-0.23385364) q[2];
sx q[2];
rz(-2.4901566) q[2];
rz(0.78178072) q[3];
sx q[3];
rz(-1.6580481) q[3];
sx q[3];
rz(2.206291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5485789) q[0];
sx q[0];
rz(-2.9055556) q[0];
sx q[0];
rz(0.4655984) q[0];
rz(-1.9904526) q[1];
sx q[1];
rz(-1.2460037) q[1];
sx q[1];
rz(1.8021072) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8663686) q[0];
sx q[0];
rz(-1.0491519) q[0];
sx q[0];
rz(-2.0417968) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6147827) q[2];
sx q[2];
rz(-1.5325452) q[2];
sx q[2];
rz(1.7234617) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7581343) q[1];
sx q[1];
rz(-1.3515359) q[1];
sx q[1];
rz(-3.0644528) q[1];
x q[2];
rz(2.2364113) q[3];
sx q[3];
rz(-1.9097121) q[3];
sx q[3];
rz(1.9012332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1107948) q[2];
sx q[2];
rz(-1.75533) q[2];
sx q[2];
rz(-3.1286855) q[2];
rz(-3.0068093) q[3];
sx q[3];
rz(-1.2603935) q[3];
sx q[3];
rz(0.24136647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6808692) q[0];
sx q[0];
rz(-0.091982059) q[0];
sx q[0];
rz(1.1750093) q[0];
rz(0.75042945) q[1];
sx q[1];
rz(-1.3576327) q[1];
sx q[1];
rz(0.39355412) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8581482) q[0];
sx q[0];
rz(-0.99020105) q[0];
sx q[0];
rz(-0.20943187) q[0];
rz(-pi) q[1];
rz(1.4541164) q[2];
sx q[2];
rz(-0.91800729) q[2];
sx q[2];
rz(-1.0972925) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8361349) q[1];
sx q[1];
rz(-1.4609481) q[1];
sx q[1];
rz(-1.1396465) q[1];
rz(-pi) q[2];
rz(2.0840852) q[3];
sx q[3];
rz(-1.9218787) q[3];
sx q[3];
rz(-2.5834192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5081818) q[2];
sx q[2];
rz(-1.1339302) q[2];
sx q[2];
rz(-2.9122162) q[2];
rz(0.091382787) q[3];
sx q[3];
rz(-1.4437851) q[3];
sx q[3];
rz(0.55575931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9230187) q[0];
sx q[0];
rz(-2.3526683) q[0];
sx q[0];
rz(-0.58513418) q[0];
rz(-1.2773889) q[1];
sx q[1];
rz(-1.2150512) q[1];
sx q[1];
rz(-2.9871984) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0994249) q[0];
sx q[0];
rz(-1.4615566) q[0];
sx q[0];
rz(-0.43719284) q[0];
rz(-pi) q[1];
rz(2.4677708) q[2];
sx q[2];
rz(-2.6148863) q[2];
sx q[2];
rz(-2.6580842) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8789492) q[1];
sx q[1];
rz(-2.7217138) q[1];
sx q[1];
rz(-2.4559569) q[1];
rz(-pi) q[2];
rz(-2.9657439) q[3];
sx q[3];
rz(-2.0324586) q[3];
sx q[3];
rz(1.3929183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.020891) q[2];
sx q[2];
rz(-2.2521844) q[2];
sx q[2];
rz(2.6013539) q[2];
rz(-0.39515105) q[3];
sx q[3];
rz(-0.65427798) q[3];
sx q[3];
rz(-1.7925709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2111135) q[0];
sx q[0];
rz(-1.2808639) q[0];
sx q[0];
rz(-1.5323918) q[0];
rz(-2.2156175) q[1];
sx q[1];
rz(-1.4264868) q[1];
sx q[1];
rz(-1.4769953) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4187936) q[0];
sx q[0];
rz(-0.72615756) q[0];
sx q[0];
rz(-0.0010869507) q[0];
rz(-pi) q[1];
rz(0.29472827) q[2];
sx q[2];
rz(-1.9325614) q[2];
sx q[2];
rz(0.28361646) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3850579) q[1];
sx q[1];
rz(-1.2278265) q[1];
sx q[1];
rz(-2.1567731) q[1];
rz(-2.1891045) q[3];
sx q[3];
rz(-0.68653169) q[3];
sx q[3];
rz(-1.2602947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.28107873) q[2];
sx q[2];
rz(-1.3382567) q[2];
sx q[2];
rz(2.4930084) q[2];
rz(2.4159238) q[3];
sx q[3];
rz(-2.9026493) q[3];
sx q[3];
rz(1.3919977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7657179) q[0];
sx q[0];
rz(-2.1633042) q[0];
sx q[0];
rz(-2.3194585) q[0];
rz(-2.0769465) q[1];
sx q[1];
rz(-1.6864265) q[1];
sx q[1];
rz(2.0815157) q[1];
rz(-0.63659414) q[2];
sx q[2];
rz(-1.7117715) q[2];
sx q[2];
rz(-2.911917) q[2];
rz(-1.0242994) q[3];
sx q[3];
rz(-0.83491769) q[3];
sx q[3];
rz(0.49900311) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
