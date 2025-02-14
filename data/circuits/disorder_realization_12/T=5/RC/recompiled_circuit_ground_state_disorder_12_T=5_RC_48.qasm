OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0491068) q[0];
sx q[0];
rz(4.4409039) q[0];
sx q[0];
rz(9.470603) q[0];
rz(-1.8419645) q[1];
sx q[1];
rz(-1.3388495) q[1];
sx q[1];
rz(-1.2531228) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078431167) q[0];
sx q[0];
rz(-0.51126152) q[0];
sx q[0];
rz(-1.0529165) q[0];
rz(1.6145176) q[2];
sx q[2];
rz(-1.3695903) q[2];
sx q[2];
rz(1.9227366) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.040267) q[1];
sx q[1];
rz(-2.3803452) q[1];
sx q[1];
rz(2.9684307) q[1];
rz(-pi) q[2];
rz(1.4433483) q[3];
sx q[3];
rz(-1.5668673) q[3];
sx q[3];
rz(-1.1636101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8200298) q[2];
sx q[2];
rz(-1.4048615) q[2];
sx q[2];
rz(-0.31537867) q[2];
rz(0.086221181) q[3];
sx q[3];
rz(-3.0486221) q[3];
sx q[3];
rz(-2.2975547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.740199) q[0];
sx q[0];
rz(-1.711015) q[0];
sx q[0];
rz(-2.8061818) q[0];
rz(-2.7566578) q[1];
sx q[1];
rz(-2.3454869) q[1];
sx q[1];
rz(1.2892494) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7766472) q[0];
sx q[0];
rz(-1.9595265) q[0];
sx q[0];
rz(1.5035648) q[0];
x q[1];
rz(0.094893242) q[2];
sx q[2];
rz(-2.2101363) q[2];
sx q[2];
rz(2.4415093) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8131539) q[1];
sx q[1];
rz(-0.48373172) q[1];
sx q[1];
rz(2.4991922) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3277169) q[3];
sx q[3];
rz(-1.7020924) q[3];
sx q[3];
rz(-2.8158902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5424767) q[2];
sx q[2];
rz(-1.3834388) q[2];
sx q[2];
rz(1.6790338) q[2];
rz(0.8820495) q[3];
sx q[3];
rz(-0.65989152) q[3];
sx q[3];
rz(-1.5006458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7733234) q[0];
sx q[0];
rz(-0.48650807) q[0];
sx q[0];
rz(0.64273709) q[0];
rz(2.2679988) q[1];
sx q[1];
rz(-1.3106376) q[1];
sx q[1];
rz(0.98683039) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8514252) q[0];
sx q[0];
rz(-2.254524) q[0];
sx q[0];
rz(0.19398035) q[0];
rz(-pi) q[1];
rz(-1.4535232) q[2];
sx q[2];
rz(-2.6542695) q[2];
sx q[2];
rz(-2.8685958) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9171608) q[1];
sx q[1];
rz(-0.81714918) q[1];
sx q[1];
rz(2.3520873) q[1];
rz(-pi) q[2];
rz(1.6968459) q[3];
sx q[3];
rz(-0.37739649) q[3];
sx q[3];
rz(2.6871329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.70283908) q[2];
sx q[2];
rz(-3.1162016) q[2];
sx q[2];
rz(2.0849483) q[2];
rz(1.3422525) q[3];
sx q[3];
rz(-1.4547576) q[3];
sx q[3];
rz(-2.2630283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1037647) q[0];
sx q[0];
rz(-2.8570638) q[0];
sx q[0];
rz(-1.9985265) q[0];
rz(-2.4567538) q[1];
sx q[1];
rz(-1.3416483) q[1];
sx q[1];
rz(-0.24982223) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9651325) q[0];
sx q[0];
rz(-1.1145381) q[0];
sx q[0];
rz(0.81911032) q[0];
rz(-2.5673713) q[2];
sx q[2];
rz(-1.5674233) q[2];
sx q[2];
rz(-0.33004883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6982242) q[1];
sx q[1];
rz(-1.0920381) q[1];
sx q[1];
rz(0.67826346) q[1];
rz(-1.5644844) q[3];
sx q[3];
rz(-2.8996116) q[3];
sx q[3];
rz(0.93130563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.41986856) q[2];
sx q[2];
rz(-1.7896174) q[2];
sx q[2];
rz(-1.8458337) q[2];
rz(1.3755679) q[3];
sx q[3];
rz(-1.4294759) q[3];
sx q[3];
rz(-0.099055722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38882035) q[0];
sx q[0];
rz(-1.3714014) q[0];
sx q[0];
rz(-2.2950628) q[0];
rz(-1.1835774) q[1];
sx q[1];
rz(-1.1776244) q[1];
sx q[1];
rz(1.902098) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13279937) q[0];
sx q[0];
rz(-1.4504878) q[0];
sx q[0];
rz(-1.8254721) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7991512) q[2];
sx q[2];
rz(-1.7424118) q[2];
sx q[2];
rz(-1.9894853) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8628564) q[1];
sx q[1];
rz(-1.5148316) q[1];
sx q[1];
rz(3.1296594) q[1];
rz(-pi) q[2];
rz(0.059885177) q[3];
sx q[3];
rz(-0.97447534) q[3];
sx q[3];
rz(-2.8879762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5014629) q[2];
sx q[2];
rz(-1.5093466) q[2];
sx q[2];
rz(-2.787369) q[2];
rz(-0.7192449) q[3];
sx q[3];
rz(-0.57643276) q[3];
sx q[3];
rz(-2.332212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3767913) q[0];
sx q[0];
rz(-0.23623315) q[0];
sx q[0];
rz(-2.7948622) q[0];
rz(2.4193343) q[1];
sx q[1];
rz(-1.5303231) q[1];
sx q[1];
rz(2.9647656) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5520598) q[0];
sx q[0];
rz(-2.0314275) q[0];
sx q[0];
rz(0.15710196) q[0];
rz(0.093164154) q[2];
sx q[2];
rz(-2.0679065) q[2];
sx q[2];
rz(0.99592956) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4011098) q[1];
sx q[1];
rz(-1.675247) q[1];
sx q[1];
rz(0.20616504) q[1];
x q[2];
rz(-3.0400254) q[3];
sx q[3];
rz(-1.0103161) q[3];
sx q[3];
rz(-1.54502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.56883812) q[2];
sx q[2];
rz(-2.4080031) q[2];
sx q[2];
rz(1.0432165) q[2];
rz(2.9853232) q[3];
sx q[3];
rz(-2.0782491) q[3];
sx q[3];
rz(3.0472896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.25065502) q[0];
sx q[0];
rz(-1.0303048) q[0];
sx q[0];
rz(-1.913273) q[0];
rz(-0.43371513) q[1];
sx q[1];
rz(-2.3408196) q[1];
sx q[1];
rz(-1.0450276) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9244773) q[0];
sx q[0];
rz(-1.3682559) q[0];
sx q[0];
rz(3.0826352) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82257338) q[2];
sx q[2];
rz(-1.2116836) q[2];
sx q[2];
rz(2.1445779) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4853128) q[1];
sx q[1];
rz(-1.6846141) q[1];
sx q[1];
rz(-0.95358221) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6125836) q[3];
sx q[3];
rz(-0.85011357) q[3];
sx q[3];
rz(-1.0320272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.23994437) q[2];
sx q[2];
rz(-1.0239235) q[2];
sx q[2];
rz(2.955692) q[2];
rz(1.2402041) q[3];
sx q[3];
rz(-1.5408206) q[3];
sx q[3];
rz(-2.127229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1321201) q[0];
sx q[0];
rz(-1.2484231) q[0];
sx q[0];
rz(0.15952071) q[0];
rz(2.3347143) q[1];
sx q[1];
rz(-1.9069549) q[1];
sx q[1];
rz(3.1330915) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59183622) q[0];
sx q[0];
rz(-1.7248132) q[0];
sx q[0];
rz(1.9817673) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0372856) q[2];
sx q[2];
rz(-2.3416714) q[2];
sx q[2];
rz(1.4664354) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.3169049) q[1];
sx q[1];
rz(-1.7959692) q[1];
sx q[1];
rz(-1.4359739) q[1];
rz(-pi) q[2];
rz(-0.17592497) q[3];
sx q[3];
rz(-2.0485544) q[3];
sx q[3];
rz(1.1087518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.27935091) q[2];
sx q[2];
rz(-1.3631577) q[2];
sx q[2];
rz(-0.19374338) q[2];
rz(-1.6019609) q[3];
sx q[3];
rz(-1.1774457) q[3];
sx q[3];
rz(1.1295454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2368471) q[0];
sx q[0];
rz(-1.756825) q[0];
sx q[0];
rz(-0.73614502) q[0];
rz(0.88788095) q[1];
sx q[1];
rz(-0.45279756) q[1];
sx q[1];
rz(3.091541) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55287191) q[0];
sx q[0];
rz(-1.3952198) q[0];
sx q[0];
rz(-2.7953327) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4718945) q[2];
sx q[2];
rz(-1.2575694) q[2];
sx q[2];
rz(-3.1196496) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0298652) q[1];
sx q[1];
rz(-1.7741412) q[1];
sx q[1];
rz(-2.8257269) q[1];
rz(-pi) q[2];
rz(-2.5727651) q[3];
sx q[3];
rz(-1.0320559) q[3];
sx q[3];
rz(2.8728268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8686409) q[2];
sx q[2];
rz(-1.7014039) q[2];
sx q[2];
rz(-3.1325373) q[2];
rz(2.791259) q[3];
sx q[3];
rz(-2.83559) q[3];
sx q[3];
rz(0.30456021) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7250799) q[0];
sx q[0];
rz(-1.059499) q[0];
sx q[0];
rz(2.1685725) q[0];
rz(1.5618207) q[1];
sx q[1];
rz(-2.2353954) q[1];
sx q[1];
rz(-0.80950338) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3448787) q[0];
sx q[0];
rz(-2.2911706) q[0];
sx q[0];
rz(-2.0551089) q[0];
x q[1];
rz(-1.3415178) q[2];
sx q[2];
rz(-1.7588758) q[2];
sx q[2];
rz(0.089872472) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1228408) q[1];
sx q[1];
rz(-2.0301308) q[1];
sx q[1];
rz(-0.18045119) q[1];
rz(-0.9904434) q[3];
sx q[3];
rz(-1.1709612) q[3];
sx q[3];
rz(-0.70580182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1805264) q[2];
sx q[2];
rz(-2.3193391) q[2];
sx q[2];
rz(-0.52214885) q[2];
rz(-2.7850049) q[3];
sx q[3];
rz(-0.57727376) q[3];
sx q[3];
rz(0.063551158) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6575573) q[0];
sx q[0];
rz(-0.35150305) q[0];
sx q[0];
rz(1.2765314) q[0];
rz(-0.36318489) q[1];
sx q[1];
rz(-2.6138432) q[1];
sx q[1];
rz(1.6539727) q[1];
rz(2.5978418) q[2];
sx q[2];
rz(-0.85493543) q[2];
sx q[2];
rz(-2.2731352) q[2];
rz(2.5616931) q[3];
sx q[3];
rz(-1.2968466) q[3];
sx q[3];
rz(-2.235835) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
