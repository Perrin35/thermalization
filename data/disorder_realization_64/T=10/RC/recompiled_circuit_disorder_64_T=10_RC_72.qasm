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
rz(-3.0769899) q[0];
sx q[0];
rz(3.1199772) q[0];
rz(2.1463483) q[1];
sx q[1];
rz(-1.8145476) q[1];
sx q[1];
rz(-1.8099161) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0776805) q[0];
sx q[0];
rz(-2.4644116) q[0];
sx q[0];
rz(2.1190686) q[0];
rz(-pi) q[1];
rz(1.1515491) q[2];
sx q[2];
rz(-0.24818072) q[2];
sx q[2];
rz(2.575945) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0845619) q[1];
sx q[1];
rz(-1.7006526) q[1];
sx q[1];
rz(-2.3741541) q[1];
x q[2];
rz(2.8033923) q[3];
sx q[3];
rz(-0.97465289) q[3];
sx q[3];
rz(0.81830922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8109479) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(2.5374106) q[2];
rz(2.1172681) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(-1.1272875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3246831) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(1.0634134) q[0];
rz(2.2564607) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(-3.1399472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3666653) q[0];
sx q[0];
rz(-1.6158551) q[0];
sx q[0];
rz(0.010618322) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4992141) q[2];
sx q[2];
rz(-1.3628236) q[2];
sx q[2];
rz(3.1285398) q[2];
rz(-pi) q[3];
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
rz(1.8257636) q[3];
sx q[3];
rz(-1.7353188) q[3];
sx q[3];
rz(2.4428575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.985618) q[2];
sx q[2];
rz(-1.6321471) q[2];
sx q[2];
rz(-2.4334811) q[2];
rz(-1.0937141) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(-0.99350199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22096069) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(-0.0050841252) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(-1.089383) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0076865772) q[0];
sx q[0];
rz(-1.0431164) q[0];
sx q[0];
rz(-0.49467996) q[0];
rz(2.543534) q[2];
sx q[2];
rz(-1.9285678) q[2];
sx q[2];
rz(1.2146815) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.1091724) q[1];
sx q[1];
rz(-1.4854327) q[1];
sx q[1];
rz(-1.9568155) q[1];
x q[2];
rz(1.009922) q[3];
sx q[3];
rz(-2.2787333) q[3];
sx q[3];
rz(-0.55299711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0638782) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(0.88469488) q[2];
rz(1.9583154) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(-1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5692212) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(0.72682056) q[0];
rz(2.3379393) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(-2.7817536) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7774178) q[0];
sx q[0];
rz(-2.4628277) q[0];
sx q[0];
rz(-1.2749519) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14881046) q[2];
sx q[2];
rz(-2.4850922) q[2];
sx q[2];
rz(-1.3615001) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4543912) q[1];
sx q[1];
rz(-1.193421) q[1];
sx q[1];
rz(-2.6881933) q[1];
rz(0.30712819) q[3];
sx q[3];
rz(-1.2369452) q[3];
sx q[3];
rz(-0.083837282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56746733) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(0.7652258) q[2];
rz(2.3848173) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(2.9005907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9969479) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(-1.416052) q[0];
rz(0.36711806) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(-2.1062772) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028776289) q[0];
sx q[0];
rz(-2.007764) q[0];
sx q[0];
rz(-0.90941888) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49145047) q[2];
sx q[2];
rz(-1.5361668) q[2];
sx q[2];
rz(-2.2143242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.1086515) q[1];
sx q[1];
rz(-1.5161533) q[1];
sx q[1];
rz(-0.095468949) q[1];
rz(-1.3607929) q[3];
sx q[3];
rz(-1.8404507) q[3];
sx q[3];
rz(-0.3013914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3570024) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(-0.76888293) q[2];
rz(0.33603493) q[3];
sx q[3];
rz(-2.3612645) q[3];
sx q[3];
rz(-1.4679573) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96200213) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(-1.1451716) q[0];
rz(-2.0369453) q[1];
sx q[1];
rz(-1.2303753) q[1];
sx q[1];
rz(0.2072269) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.094039) q[0];
sx q[0];
rz(-2.5045966) q[0];
sx q[0];
rz(-0.23547049) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5404262) q[2];
sx q[2];
rz(-1.1784369) q[2];
sx q[2];
rz(-0.62077921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0761557) q[1];
sx q[1];
rz(-1.8199311) q[1];
sx q[1];
rz(-2.8912828) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3239922) q[3];
sx q[3];
rz(-1.6896276) q[3];
sx q[3];
rz(2.6314051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.8032288) q[2];
sx q[2];
rz(-0.88023606) q[2];
sx q[2];
rz(2.4198789) q[2];
rz(1.8668113) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(0.26003626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24467829) q[0];
sx q[0];
rz(-1.5667863) q[0];
sx q[0];
rz(2.4095643) q[0];
rz(-3.1320944) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(2.8498555) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1222629) q[0];
sx q[0];
rz(-0.42207345) q[0];
sx q[0];
rz(3.0164369) q[0];
rz(-pi) q[1];
rz(-1.2675915) q[2];
sx q[2];
rz(-1.469194) q[2];
sx q[2];
rz(0.50819699) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2597255) q[1];
sx q[1];
rz(-1.423466) q[1];
sx q[1];
rz(2.8357361) q[1];
x q[2];
rz(2.6187906) q[3];
sx q[3];
rz(-1.8409981) q[3];
sx q[3];
rz(-0.064388007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.26414028) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(0.83731246) q[2];
rz(1.9705747) q[3];
sx q[3];
rz(-1.6224909) q[3];
sx q[3];
rz(-3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0432805) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(-0.6638546) q[0];
rz(-0.10617667) q[1];
sx q[1];
rz(-2.5352434) q[1];
sx q[1];
rz(-2.1829139) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6129235) q[0];
sx q[0];
rz(-2.0818424) q[0];
sx q[0];
rz(2.5445166) q[0];
rz(-pi) q[1];
rz(1.0326951) q[2];
sx q[2];
rz(-1.8128464) q[2];
sx q[2];
rz(-0.48697105) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2221335) q[1];
sx q[1];
rz(-0.53680116) q[1];
sx q[1];
rz(1.3608576) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58416768) q[3];
sx q[3];
rz(-1.6730047) q[3];
sx q[3];
rz(-0.16971961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.133193) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(-0.91910249) q[2];
rz(1.7637926) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(-1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.33525) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(0.43689716) q[0];
rz(-0.70029744) q[1];
sx q[1];
rz(-1.7179787) q[1];
sx q[1];
rz(2.6760496) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7882999) q[0];
sx q[0];
rz(-0.71209891) q[0];
sx q[0];
rz(-3.0112991) q[0];
rz(-1.8579673) q[2];
sx q[2];
rz(-2.7138777) q[2];
sx q[2];
rz(-0.62921333) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1411966) q[1];
sx q[1];
rz(-1.2277514) q[1];
sx q[1];
rz(0.47980206) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2113308) q[3];
sx q[3];
rz(-1.8065479) q[3];
sx q[3];
rz(0.93709968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-2.6022544) q[2];
sx q[2];
rz(1.4979866) q[2];
rz(2.9368029) q[3];
sx q[3];
rz(-2.1361165) q[3];
sx q[3];
rz(2.8695316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4964504) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(2.0196594) q[0];
rz(2.3902068) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(0.5823935) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96888992) q[0];
sx q[0];
rz(-2.4318998) q[0];
sx q[0];
rz(-1.2560647) q[0];
rz(1.0803797) q[2];
sx q[2];
rz(-1.7047802) q[2];
sx q[2];
rz(1.6756563) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.45422428) q[1];
sx q[1];
rz(-2.9240989) q[1];
sx q[1];
rz(-2.1464159) q[1];
rz(0.21721812) q[3];
sx q[3];
rz(-1.5434885) q[3];
sx q[3];
rz(-2.581493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0104684) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(-2.3790512) q[2];
rz(1.7307581) q[3];
sx q[3];
rz(-2.2236731) q[3];
sx q[3];
rz(1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0762155) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(1.3394042) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(1.273524) q[2];
sx q[2];
rz(-0.77902972) q[2];
sx q[2];
rz(-3.0697889) q[2];
rz(2.0740261) q[3];
sx q[3];
rz(-1.9964841) q[3];
sx q[3];
rz(-1.869429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];